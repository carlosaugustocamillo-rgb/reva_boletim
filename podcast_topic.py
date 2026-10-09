"""Podcast-only thematic entrypoint; never invokes campaign or calendar services."""
import json
import re
import uuid
import unicodedata
from pathlib import Path

import podcast_editorial


def validate_request(body):
    mode = body.get('modo_podcast', 'semanal')
    if mode not in ('semanal', 'tema'):
        raise ValueError('Modo de podcast inválido.')
    if mode == 'semanal':
        return {'modo_podcast': mode}
    theme = body.get('tema_podcast', '')
    instructions = body.get('orientacoes_podcast', '')
    if not isinstance(theme, str) or not 1 <= len(theme.strip()) <= 300:
        raise ValueError('Informe um tema com até 300 caracteres.')
    if not isinstance(instructions, str) or len(instructions) > 3000:
        raise ValueError('As orientações devem ter até 3000 caracteres.')
    selection = body.get('curadoria_tema_id')
    if not body.get('somente_curadoria') and not body.get('roteiro_aprovado_id'):
        if not isinstance(selection, str) or not re.fullmatch(r'[a-f0-9]{32}', selection):
            raise ValueError('Busque e selecione as referências do tema antes de gerar o roteiro.')
        ids = body.get('artigos_podcast_aprovados')
        if not isinstance(ids, list) or not 1 <= len(ids) <= 6 or any(not isinstance(x, str) for x in ids):
            raise ValueError('Selecione entre um e seis artigos para o podcast.')
    return {'modo_podcast': mode, 'tema_podcast': theme.strip(),
            'orientacoes_podcast': instructions.strip(), 'curadoria_tema_id': selection,
            'resumos': False, 'mailchimp': False, 'audio': False, 'firebase': False,
            'roteiro': True}


def _normalized_title(value):
    return ' '.join(re.findall(r'\w+', unicodedata.normalize('NFKC', value).casefold()))


def search_articles(theme):
    """Search only the supplied input, never weekly queries or publication dates.

    Resolve a pasted citation before general thematic expansion so the requested
    paper is not lost among the first 30 papers of a broader keyword search.
    """
    from Bio import Entrez
    from boletim_service import buscar_info_estruturada
    from revamais_service import gerar_query_pubmed_tema

    def search(query):
        with Entrez.esearch(db='pubmed', term=query, retmax=30, sort='relevance') as handle:
            ids = Entrez.read(handle)['IdList']
        return buscar_info_estruturada(ids) if ids else []

    value = theme.strip()
    pmid = re.fullmatch(r'(?:PMID:\s*|https?://pubmed\.ncbi\.nlm\.nih\.gov/)?(\d{1,9})/?', value, re.I)
    doi = re.fullmatch(r'(?:https?://doi\.org/|doi:\s*)?(10\.\d{4,9}/[^\s"]+)', value, re.I)
    if pmid or doi:
        query = f'{pmid.group(1)}[UID]' if pmid else f'"{doi.group(1)}"[AID]'
        articles = search(query)
    else:
        title = value.replace('"', '').strip().rstrip('.')
        # Full quoted titles may be absent from PubMed's phrase index. Match
        # title words first, then verify the complete normalized title locally.
        stopwords = {'a', 'an', 'the', 'of', 'to', 'in', 'and', 'with', 'for', 'on', 'by', 'from', 'at'}
        words = list(dict.fromkeys(w for w in _normalized_title(title).split() if w not in stopwords))
        if not words:
            raise ValueError('Informe um tema ou título com palavras específicas.')
        title_query = ' AND '.join(f'{word}[Title]' for word in words)
        exact = [a for a in search(title_query)
                 if _normalized_title(a.get('titulo', '')) == _normalized_title(title)]
        if exact:
            query, articles = title_query, exact
        else:
            query = gerar_query_pubmed_tema(value)
            articles = search(query)
    return query, [a for a in articles if a.get('resumo_original', '').strip()]


def run_topic(options, base_dir, client):
    options = {**options, **validate_request(options)}
    directory = Path(base_dir) / 'podcast_topics'
    directory.mkdir(parents=True, exist_ok=True)
    theme = options['tema_podcast']
    if options.get('somente_curadoria'):
        yield f'🔎 Nova busca exclusiva no PubMed: {theme} (sem limite de data; sem consultas do boletim semanal).'
        query, articles = search_articles(theme)
        yield f'Consulta utilizada: {query}'
        if not articles:
            yield 'Nenhum artigo com resumo encontrado. Ajuste o tema e busque novamente.'
        selection_id = uuid.uuid4().hex
        snapshot = {'id': selection_id, 'tema': theme, 'query': query, 'articles': articles}
        (directory / f'{selection_id}.json').write_text(json.dumps(snapshot, ensure_ascii=False), encoding='utf-8')
        yield {'tipo': 'curadoria_podcast', 'modo_podcast': 'tema', 'tema_podcast': theme,
               'orientacoes_podcast': options['orientacoes_podcast'],
               'curadoria_tema_id': selection_id, 'artigos_sugeridos': articles, 'limite': 6,
               'busca': {'origem': 'tema', 'consulta': query, 'sem_limite_data': True},
               'mensagem': 'Selecione até seis fontes para o episódio temático.' if articles else
                           'Nenhum artigo com resumo encontrado. Ajuste o tema e busque novamente.'}
        return
    path = directory / f"{options['curadoria_tema_id']}.json"
    if not path.is_file():
        raise ValueError('Curadoria não encontrada. Busque novamente as referências.')
    snapshot = json.loads(path.read_text(encoding='utf-8'))
    if snapshot['tema'] != theme:
        raise ValueError('O tema mudou. Busque novamente as referências.')
    ids = options['artigos_podcast_aprovados']
    by_id = {str(a['pmid']): a for a in snapshot['articles']}
    if len(set(ids)) != len(ids) or any(i not in by_id for i in ids):
        raise ValueError('A seleção não corresponde aos artigos desta curadoria.')
    articles = [by_id[i] for i in ids]
    context = options.get('contexto_pubmed_manual')
    if context:
        expected = set(ids)
        actual = {str(s.get('pmid_ancora', '')) for s in context.get('estudos', [])}
        if actual != expected:
            raise ValueError('Prepare novamente o contexto para os artigos selecionados.')
    else:
        from pubmed_related import enrich_episode, enabled_from_env
        yield '📚 Preparando contexto científico do podcast...'
        context = enrich_episode(articles, base_dir, client,
                                 enabled=enabled_from_env() and options.get('referencias_pubmed', True))
    draft = yield from podcast_editorial.generate_episode_steps(
        articles, context, client, base_dir=base_dir,
        episode_context={'modo': 'tema', 'tema': theme, 'orientacoes': options['orientacoes_podcast']})
    podcast_editorial.save_draft(base_dir, draft)
    review = podcast_editorial.review_payload(draft)
    yield {'modo_podcast': 'tema', 'tema_podcast': theme, 'roteiro_editorial': review,
           'roteiro_erro': None if review['can_approve'] else 'Confira as pendências no parecer científico.',
           'mailchimp': {'status': 'skipped', 'error': None}, 'referencias_pubmed': context}
