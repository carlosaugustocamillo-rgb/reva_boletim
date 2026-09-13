# Roteiro editorial do RevaCast Weekly

Esta mudança se aplica somente ao podcast. Consultas, traduções, conteúdo HTML/TXT
e envio do boletim por e-mail permanecem no fluxo anterior. O Reva+ também não muda.

## Modelo e etapas

`PODCAST_SCRIPT_MODEL=gpt-6-astra` é o padrão, confirmado como disponível na conta
durante a implementação. Pode ser substituído por `gpt-5.5` nessa variável sem
alterar os modelos de resumos, triagem PubMed ou brief Spotify. Não há fallback
silencioso para GPT-4o ou Gemini no roteiro. Falhas ficam visíveis no painel.

O cliente usa Chat Completions, saída JSON Schema estrita, sem `temperature` ou
`top_p`, sem alterar SDK/provedor, e esforço de raciocínio padrão do modelo.
São três chamadas por episódio (não três por artigo), até 12.000 tokens de saída
por chamada, timeout de 180 segundos e uma repetição de transporte pelo SDK.
O objetivo de palavras é 180 + 280 por estudo; o limite é 350 + 420 por estudo.
Isso é orçamento de texto, não promessa de duração exata de áudio.

1. `prompts/podcast/plan.txt`: pauta global, achados com trechos literais e avaliação
   crítica preliminar de cada estudo.
2. `prompts/podcast/write.txt`: conversa Ivo/Manu completa com abertura e fechamento,
   contexto complementar, transições e ressalvas integradas.
3. `prompts/podcast/audit.txt`: checagem contra as fontes. Pendências científicas
   bloqueiam aprovação; observações de estilo ficam para o revisor.

As etapas aparecem nos logs e oferecem pontos de cancelamento entre chamadas.
O modelo não navega: usa somente os artigos aprovados e as referências selecionadas.
O uso real de tokens por etapa fica no parecer. A estimativa antiga de custos
do boletim **não inclui** o novo roteiro; isso está explícito no retorno.

## Qualidade científica: o que é e o que não é

O material disponível é resumo original, tradução identificada ou só metadados.
BibTeX sem resumo não sustenta afirmações sobre resultados. A análise não se
apresenta como leitura do texto completo, RoB 2, ROBINS-I, AMSTAR ou GRADE.

O parecer distingue ponto favorável documentado, limitação documentada, informação
não relatada no resumo e limite inferencial do desenho. Não transforma ausência de
informação em falha metodológica. Não privilegia automaticamente ECR/RS nem descarta
qualitativos/observacionais. Similaridade não é concordância nem qualidade.

Validações locais checam PMIDs/ordem, fontes, existência literal dos trechos de apoio,
participação dos dois hosts, presença da ressalva nas falas e limite de palavras.
Isso não prova a validade científica da interpretação: a auditoria é automatizada
e a leitura/aprovação humana continua obrigatória.

## Fluxo no painel

1. Sugerir artigos → escolher até seis → selecionar referências por artigo.
2. Executar a geração do roteiro. Para um piloto, desmarcar Mailchimp e Firebase.
   Mesmo com áudio marcado, um **novo** roteiro aguarda aprovação.
3. Ler a conversa e abrir “Parecer científico e fontes por estudo”; baixar TXT se
   desejar. “Revisar último roteiro” recupera a versão persistida após recarregar.
4. Marcar a confirmação de leitura e clicar em “Aprovar e gerar prévia de áudio”.
   Essa ação não refaz busca/roteiro, não envia e-mail e não publica no RSS.
5. Ouvir o arquivo em “Baixar prévia de áudio”. O piloto de conteúdo com dois
   artigos reais deve ser revisado antes de publicação/deploy definitivo.

Não há editor manual de falas nesta versão. Havendo erro científico bloqueante,
conferir/ajustar as fontes e gerar um novo rascunho para revisão. Não foi acrescentada
publicação automática da prévia. Os controles de publicação existentes permanecem
separados; não iniciar outra síntese paga sem intenção explícita.

## Persistência e áudio

`data/editorial/<uuid>.json` contém fontes, pauta, falas, auditoria, modelo e tokens.
O TXT ao lado contém só as falas. A aprovação exige o SHA-256 da versão exibida,
que protege roteiro, pauta, fontes e auditoria. Alterações invalidam a versão.
Os JSON/TXT antigos de etapas 8/8.5 são aliases de compatibilidade; não são a fonte
de autorização para áudio. O áudio usa o JSON editorial aprovado.

Não existe mais introdução aleatória fora do roteiro nem substituição das
transições por frases genéricas. O TTS recebe o texto canônico; siglas/pronúncia
devem ser resolvidas na escrita, antes da revisão. A trilha musical fixa permanece.
Se faltar um bloco de áudio, o episódio parcial não é montado/publicado.

As vozes aprovadas permanecem: Ivo/Archer e Manu/Hope, Eleven Multilingual v2,
velocidade 1.1, estabilidade/similaridade 1.0, estilo 0.0. Ver `podcast-voices.md`.

## API e implantação

- `GET /ultimo-roteiro`: acrescenta `editorial`, preserva `conteudo` legado.
- `GET /podcast-roteiro/{id}`: texto, parecer, fontes e identificador de versão.
- `GET /podcast-roteiro/{id}/texto`: download legível.
- `POST /podcast-roteiro/{id}/aprovar`, corpo `{"sha256":"..."}`.
- `POST /iniciar-boletim`: áudio de versão aprovada exige `roteiro_aprovado_id` e
  `roteiro_aprovado_sha256`. Para prévia: áudio true; resumos, roteiro, Mailchimp,
  Firebase, referências e brief false.
- `GET /baixar-audio-podcast/{filename}`: somente nomes de episódios MP3 válidos.

Railway: manter `OPENAI_API_KEY`; acrescentar `PODCAST_SCRIPT_MODEL=gpt-6-astra`
se quiser configuração explícita (já é o padrão). Não alterar `OPENAI_TEXT_MODEL`
por causa do podcast. Manter as variáveis de voz aprovadas. Backend e frontend
precisam ser implantados juntos. A pasta `data` deve continuar persistida no volume.

Esta implementação não faz commit, push, deploy nem publicação automaticamente.

## Verificação

`python -m unittest test_podcast_editorial test_podcast_editorial_api test_podcast_pubmed_integration test_podcast_voice_settings test_elevenlabs_utils test_pubmed_related test_connected_papers test_pubmed_fetch_resilience`

Os testes usam fontes sintéticas e serviços simulados; não são um piloto clínico.
Há regressão do conteúdo do e-mail contra a versão commitada, falha editorial sem
impedir envio, aprovação/versionamento, texto exato no TTS e bloqueio de áudio parcial.

Documentação oficial consultada: https://developers.openai.com/api/docs/guides/latest-model
(modelo GPT-6 Astra, compatibilidade Chat Completions/Structured Outputs e parâmetros).
