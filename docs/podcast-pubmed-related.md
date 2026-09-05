# Referências anteriores no podcast RevaCast Weekly

A alteração é exclusiva do **podcast**. As buscas, traduções, textos principal e
detalhado, HTML e agendamento do e-mail permanecem iguais. Não requer chave do
Semantic Scholar nem assinatura do Connected Papers.

## Uso

No painel, mantenha **Criar Roteiro do Podcast** e **Contextualizar o podcast com
estudos anteriores (PubMed)** marcados. O recurso acompanha a geração normal.
Na API, `referencias_pubmed=false` desativa a etapa por execução, tanto em
`POST /iniciar-boletim` como em `GET /rodar-boletim-stream`.

O limite inicial é de **dois estudos principais**, **12 candidatos por estudo**
e **até duas referências anteriores por estudo**. O total atual de estudos do
podcast não muda. A pesquisa semanal continua a mesma; referências históricas
não são adicionadas ao boletim por e-mail nem apresentadas como novidades.

1. Depois da seleção dos estudos do podcast, consulta `ELink` com
   `linkname=pubmed_pubmed` e `cmd=neighbor_score`.
2. Busca títulos, abstracts e metadados de candidatos em lote por `EFetch`.
3. Exclui duplicatas por PMID/DOI/título, artigos sem resumo, retratações/alertas,
   publicações inadequadas e estudos animais identificáveis nos descritores.
4. Só aceita como referência anterior quando a data mais tardia possível do
   candidato precede a data mais antiga possível do estudo principal. Datas
   parciais são tratadas conservadoramente. Não inventa dia/mês.
5. Usa uma chamada de triagem por estudo principal, com saída JSON estruturada:
   população, intervenção e desfecho devem ser explicitamente compatíveis;
   comparador, desenho e limitações são registrados. PMIDs e trechos literais
   que sustentam a decisão são conferidos contra as fontes recebidas.
6. Envia ao roteirista apenas as referências admitidas, com seus abstracts e
   atribuições. Solicita um bloco de 2–4 falas adicionais por estudo, mencionando
   ao menos uma referência admitida por autor/ano e explicando sua contribuição
   e uma limitação.

A triagem é automatizada, **não uma avaliação clínica definitiva**. O score de
similaridade não mede qualidade nem concordância. Achados nulos e discordantes
podem ser admitidos. O roteiro e as referências continuam sujeitos à revisão
editorial. Os dados se limitam aos abstracts; não pressupõem texto integral.

## Falhas e custos

- No máximo uma requisição NCBI por segundo neste cliente/processo; execuções
  paralelas em várias instâncias devem ser evitadas ou usar limitação compartilhada.
- Até três tentativas por chamada HTTP, timeout de conexão/leitura e backoff.
  `Retry-After` longo encerra a tentativa de enriquecimento daquela âncora.
- Cache por PMID por 24 horas (resposta vazia por uma hora), escrito atomicamente.
- Triagem: `gpt-4o`, temperatura zero, timeout de 35 segundos, sem retentativas
  automáticas; até duas chamadas adicionais por episódio na configuração padrão.
- As consultas NCBI não exigem pagamento; a triagem usa a conta OpenAI já
  configurada. Tokens reais ficam em `tokens_triagem` no relatório. A estimativa
  fixa de custo de texto do painel é legada e **não mede** esse gasto adicional.
- Falha, ausência de resultados ou triagem inválida: o roteiro continua baseado
  somente no estudo principal. Isso não significa ausência de literatura prévia.

## Configuração

Veja `.env.example`. Não é necessário alterar o `.env` para os padrões:

```dotenv
PODCAST_PUBMED_RELATED_ENABLED=true
PODCAST_PUBMED_RELATED_MAX_ANCHORS=2
PODCAST_PUBMED_RELATED_MAX_CANDIDATES=12
PODCAST_PUBMED_RELATED_MAX_REFERENCES=2
PODCAST_PUBMED_RELATED_CACHE_HOURS=24
PODCAST_PUBMED_RELATED_MODEL=gpt-4o
```

`ENTREZ_EMAIL` já existente identifica a aplicação. `NCBI_API_KEY` é opcional.
`PODCAST_PUBMED_RELATED_ENABLED=false` desativa globalmente, mesmo quando a opção
do painel estiver marcada. As alterações precisam estar no servidor em execução;
não há implantação automática por editar o repositório local.

## Rastreabilidade

- Relatório: `data/referencias/contexto_pubmed_YYYY-MM-DD.json` (última execução do
  dia com o recurso ativo), incluindo fontes recebidas, exclusões e triagem.
- Download: `/baixar-referencias-podcast/YYYY-MM-DD`, também acessível ao concluir
  pelo painel.
- Notas com links: anexadas ao arquivo de brief do Spotify. O RSS existente usa
  sua descrição padrão; não publica automaticamente o conteúdo do brief.
- Se o upload de roteiros estiver ativo, o contexto também acompanha o documento
  do roteiro no Firestore. Os arrays de diálogo mantêm o formato existente.

## Verificação

```sh
.venv/bin/python -m unittest test_pubmed_related test_podcast_pubmed_integration test_elevenlabs_utils -v
```

Os testes não enviam e-mail, não publicam, não consomem créditos e não precisam
de segredos. Cobrem HTTP/cache, deduplicação, cronologia, triagem conservadora,
falhas e integração. A comparação com o recurso ligado/desligado exige igualdade
byte a byte dos resumos e igualdade das chamadas de conteúdo/agendamento Mailchimp.

Referências técnicas: [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/),
[limites NCBI](https://www.ncbi.nlm.nih.gov/books/NBK25497/) e
[OpenAI Structured Outputs](https://developers.openai.com/api/docs/guides/structured-outputs).
