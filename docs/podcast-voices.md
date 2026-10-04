# Vozes do RevaCast Weekly

Configuração solicitada em 11/09/2026, exclusiva do áudio do podcast:

| Apresentador | Voz | ID | Modelo | Velocidade | Estabilidade | Similaridade | Estilo |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Ivo | Voz brasileira selecionada | `NFmEzNOony1UsEJGXLth` | `eleven_v4` / Dialogue | 1.1 | 100% | 100% | 0% |
| Manu | Hope | `uYXf8XasLslADfZ2MB4u` | `eleven_v4` / Dialogue | 1.1 | 100% | 100% | 0% |

Após a primeira prévia, a velocidade de ambos foi aumentada para 1.1 a pedido do usuário.

O mesmo helper aplica os parâmetros à apresentação inicial e a todas as falas
dos estudos, incluindo a despedida. `use_speaker_boost=true` foi preservado.
Não há ajuste de velocidade posterior por edição do arquivo: `speed` é enviado
ao ElevenLabs em `voice_settings` para cada apresentador.

## Variáveis no Railway

Depois de publicar o código, atualizar as variáveis no serviço do backend:

```dotenv
ELEVEN_VOICE_ID_HOST=NFmEzNOony1UsEJGXLth
ELEVEN_VOICE_ID_COHOST=uYXf8XasLslADfZ2MB4u
ELEVEN_AUDIO_MODEL=eleven_v4
ELEVEN_AUDIO_DIALOGUE_ENABLED=true
ELEVEN_AUDIO_DIALOGUE_MODEL=eleven_v4
ELEVEN_AUDIO_FALLBACK_MODEL=eleven_multilingual_v2
ELEVEN_VOICE_SPEED_HOST=1.1
ELEVEN_VOICE_SPEED_COHOST=1.1
ELEVEN_VOICE_STABILITY=1.0
ELEVEN_VOICE_SIMILARITY_BOOST=1.0
ELEVEN_VOICE_STYLE=0.0
ELEVEN_AUDIO_LANGUAGE_CODE=pt
```

Esses valores são também os novos padrões do código. Variáveis já existentes
no Railway sobrescrevem os padrões: os IDs antigos precisam ser atualizados.
A edição do `.env` local não altera o Railway. Não é preciso mudar a chave de API.

`ELEVEN_AUDIO_MODEL=eleven_v4` e `ELEVEN_AUDIO_DIALOGUE_ENABLED=true` ativam o
Text to Dialogue para Ivo e Manu. O pipeline envia as duas vozes em cada bloco
para permitir pausas, alternância e maior variação emocional. Se a chamada v4
falhar, o estudo usa o `ELEVEN_AUDIO_FALLBACK_MODEL` para não perder a execução.
`ELEVEN_AUDIO_FALLBACK_MODEL` só seleciona o TTS alternativo no modo Dialogue;
não substitui o modelo principal.

Os percentuais são enviados na escala de **0 a 1**, e não de 0 a 100.
As velocidades aceitas pela configuração estão entre 0.7 e 1.2.

## Português e limitação da substituição de idioma

O idioma editorial continua sendo português. Porém, segundo a
[referência oficial do Text to Speech](https://elevenlabs.io/docs/api-reference/text-to-speech/convert),
`language_code` **não é suportado pelo Multilingual v2**. Por isso o pipeline
registra o idioma solicitado, avisa sobre a detecção automática pelo texto e
omite esse parâmetro na requisição v2. Não é uma garantia de sotaque brasileiro.
Não foi usado outro modelo nem texto adicional artificial para contornar a limitação.

## Verificação

```shell
.venv_test/bin/python -m unittest test_podcast_voice_settings test_elevenlabs_utils test_podcast_pubmed_integration
```

Os testes verificam payloads sem síntese paga. Configurações não regeneram
áudios já existentes; amostras e episódios antigos continuam com suas vozes originais.
