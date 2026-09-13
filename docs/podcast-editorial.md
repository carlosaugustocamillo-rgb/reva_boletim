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
Se a ressalva registrada não corresponder às falas, há no máximo uma chamada
adicional de escrita com as mesmas fontes e o erro identificado. O resultado
passa novamente pela validação e auditoria; não há repetição ilimitada.
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
PDFs pesquisáveis podem ser anexados tanto ao artigo principal quanto às referências
semelhantes. O texto é extraído página a página, persistido no volume em
`data/podcast_pdfs` e identificado por hash. Quando há PDF, ele é a fonte científica
prioritária nas três etapas; sem PDF, permanece o resumo com sua limitação explícita.
O uso real de tokens por etapa fica no parecer. A estimativa antiga de custos
do boletim **não inclui** o novo roteiro; isso está explícito no retorno.

## Conversa e ritmo (`editorial-4-conversation`)

A pauta agora organiza dúvidas e limites que os dois colegas podem discutir.
A escrita pede que cada resposta desenvolva algo específico da fala anterior:
um resultado, uma ressalva ou uma interpretação. Ambos podem perguntar e explicar.
Não há número obrigatório de trocas, interjeições ou perguntas por estudo.

Para reduzir períodos longos, a referência de escrita é uma ideia por frase,
geralmente em 8–18 palavras. Acima de 25 palavras, o escritor reavalia onde separar
a explicação, sem cortar relações científicas. São referências de ritmo, não
limites impostos por código nem prova de naturalidade. Respostas breves convivem
com explicações maiores; o orçamento total de palavras continua igual.

`write.txt` inclui exemplos fictícios e uma revisão da interação dentro da mesma
chamada. O ajuste já existente de `spoken_caution` usa esse mesmo prompt, para
preservar o estilo ao corrigir uma divergência. A auditoria procura exposições
alternadas, perguntas de apresentador, frases densas e cautelas recitadas; indica
a localização como `note`, sem bloquear a aprovação por preferência de estilo.
Erros científicos continuam bloqueantes. A auditoria distingue PDF de resumo e
não interpreta informação ausente como procedimento que não foi realizado.

O fluxo ativo não usa os reescritores legados de transições nem uma etapa de
polimento depois da aprovação. `dialogues()` retira só os metadados e o TTS recebe
`text` exatamente como foi aprovado, inclusive respostas curtas e falas sucessivas
do mesmo locutor. Modelo, vozes, velocidade e pausas de mixagem não mudaram.

### Comparação controlada

`scripts/compare_podcast_conversation.py` executa as etapas editoriais antiga e nova
usando os mesmos dois estudos sintéticos e a mesma referência complementar de
`tests/fixtures/podcast_conversation.json`. Usa a mesma meta de palavras, modelo,
schema e parâmetros; só muda o diretório de prompts. É uma execução paga de texto
sob demanda, fora do pipeline de publicação, e não chama ElevenLabs/Firebase.

```
./.venv/bin/python scripts/compare_podcast_conversation.py \
  --before-prompts /caminho/snapshot-dos-prompts \
  --output /private/tmp/comparacao-nova
```

Os artefatos incluem fontes, prompts e hashes, roteiros, parecer, consumo real e
contagens descritivas. A avaliação humana deve conferir dependência entre falas,
contribuição dos dois, ritmo, transições e preservação dos fatos. Não usar apenas
a redução no tamanho das frases para decidir se o diálogo ficou natural.

Ver [o comparativo realizado em 13/09/2026](podcast-conversation-comparison.md),
com trechos reais das respostas da API sobre fontes sintéticas e limites da avaliação.

## Qualidade científica: o que é e o que não é

O material disponível é texto extraído de PDF, resumo original, tradução identificada
ou só metadados. BibTeX sem resumo/PDF não sustenta afirmações sobre resultados.
A análise identifica o material efetivamente fornecido; não atribui uma avaliação
formal de RoB 2, ROBINS-I, AMSTAR ou GRADE.

O parecer distingue ponto favorável documentado, limitação documentada, informação
não relatada no material disponível e limite inferencial do desenho. Não transforma ausência de
informação em falha metodológica. Não privilegia automaticamente ECR/RS nem descarta
qualitativos/observacionais. Similaridade não é concordância nem qualidade.

Validações locais checam PMIDs/ordem, fontes, existência literal dos trechos de apoio,
participação dos dois hosts, presença da ressalva nas falas e limite de palavras.
Isso não prova a validade científica da interpretação: a auditoria é automatizada
e a leitura/aprovação humana continua obrigatória.

A comparação da ressalva tolera maiúsculas, espaços e pontuação textual e pode
abranger falas consecutivas. Não usa similaridade aproximada nem ignora negações
ou diferenças numéricas. Divergências persistentes identificam o PMID e o texto
registrado e continuam bloqueando a aprovação.

## Fluxo no painel

1. Sugerir artigos → escolher até seis → selecionar referências por artigo.
2. Executar a geração do roteiro. Para um piloto, desmarcar Mailchimp e Firebase.
   Mesmo com áudio marcado, um **novo** roteiro aguarda aprovação.
3. Ler a conversa e abrir “Parecer científico e fontes por estudo”; baixar TXT se
   desejar. “Revisar último roteiro” recupera a versão persistida após recarregar.
4. Para revisar fora do painel, usar “Copiar roteiro” e depois “Editar / colar texto
   revisado”. Cada fala começa em nova linha com `Ivo:` ou `Manu:`; rótulos em
   negrito copiados do ChatGPT também são aceitos. Colar apenas as falas completas.
   Dois-pontos dentro das falas são preservados. A edição manual aceita até
   60.000 caracteres, sem aplicar o alvo de palavras por estudo da escrita automática.
5. “Salvar e conferir texto” cria outro UUID, preserva o original e inicia somente
   a auditoria das fontes em segundo plano. O painel acompanha essa conferência.
   A auditoria usa créditos de texto, mas não reescreve falas. PDF e resumos já
   persistidos continuam disponíveis; não é preciso reenviá-los.
6. Ler o novo parecer, confirmar a leitura e clicar em “Aprovar, gerar e publicar
   áudio”. Essa ação usa o texto salvo, sem nova escrita, e publica no Firebase/RSS.

A edição usa `script.edited_dialogue`, uma lista ordenada de falas. A associação
das afirmações aos estudos e a presença das ressalvas são verificadas pela auditoria
contra todas as fontes e a pauta original; a edição pode reorganizar os estudos.
O parecer da pauta é apresentado como planejado, e as observações da checagem
automática são refeitas para o texto editado. Essa checagem não garante ausência
de erro: permanece necessária a revisão humana.

Uma versão em `auditing` não pode ser aprovada. Se a checagem falhar, o texto
permanece salvo como `blocked` e pode ser editado/salvo novamente. Se o servidor
reiniciar durante a tarefa em segundo plano, a versão pode permanecer em `auditing`;
“Revisar último roteiro” recupera o texto, e salvar uma nova tentativa reinicia a
checagem. O SHA muda com a auditoria, exigindo aprovação da versão atual.

## Persistência e áudio

`data/editorial/<uuid>.json` contém fontes, pauta, falas, auditoria, modelo e tokens.
O TXT ao lado contém só as falas. A aprovação exige o SHA-256 da versão exibida,
que protege roteiro, pauta, fontes e auditoria. Alterações invalidam a versão.
Os JSON/TXT antigos de etapas 8/8.5 são aliases de compatibilidade; não são a fonte
de autorização para áudio. O áudio usa o JSON editorial aprovado.

Roteiros gerados que falhem na validação ou na auditoria são preservados como
`blocked`, com pendência visível no painel, sem substituir os aliases ativos ou
produzir um brief. Quando houver tentativa de ajuste, `original_script` guarda
a primeira versão no JSON. Sem uma resposta de escrita completa e estruturada,
não existe rascunho narrável a preservar.

Sem episódio completo aprovado desta execução, o pipeline não executa upload,
RSS, resgate de segmentos antigos ou rascunho WhatsApp. O envio de e-mail segue
as opções já selecionadas. O status terminal continua `completed` para encerrar
o polling, mas `outcome=partial` e a mensagem final explicitam pendências do
podcast; `outcome=pending_review` distingue a espera normal por aprovação.

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
- `POST /podcast-roteiro/{id}/editar`, corpo `{"sha256":"...", "text":"Ivo: ...\nManu: ..."}`:
  salva e retorna uma nova versão em `auditing`; consultar seu GET até encerrar.
- `POST /iniciar-boletim`: áudio de versão aprovada exige `roteiro_aprovado_id` e
  `roteiro_aprovado_sha256`. A execução do áudio aprovado publica automaticamente
  o MP3 no Firebase e atualiza o RSS; não há etapa manual adicional.
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
