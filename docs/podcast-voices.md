# Vozes do RevaCast Weekly

Configuração solicitada em 11/09/2026, exclusiva do áudio do podcast:

| Apresentador | Voz | ID | Modelo | Velocidade | Estabilidade | Similaridade | Estilo |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Ivo | Archer | `L0Dsvb3SLTyegXwtm47J` | `eleven_multilingual_v2` | 1.1 | 100% | 100% | 0% |
| Manu | Hope | `uYXf8XasLslADfZ2MB4u` | `eleven_multilingual_v2` | 1.1 | 100% | 100% | 0% |

Após a primeira prévia, a velocidade de ambos foi aumentada para 1.1 a pedido do usuário.

O mesmo helper aplica os parâmetros à apresentação inicial e a todas as falas
dos estudos, incluindo a despedida. `use_speaker_boost=true` foi preservado.
Não há ajuste de velocidade posterior por edição do arquivo: `speed` é enviado
ao ElevenLabs em `voice_settings` para cada apresentador.

## Variáveis no Railway

Depois de publicar o código, atualizar as variáveis no serviço do backend:

```dotenv
ELEVEN_VOICE_ID_HOST=L0Dsvb3SLTyegXwtm47J
ELEVEN_VOICE_ID_COHOST=uYXf8XasLslADfZ2MB4u
ELEVEN_AUDIO_MODEL=eleven_multilingual_v2
ELEVEN_AUDIO_DIALOGUE_ENABLED=false
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

`ELEVEN_AUDIO_MODEL` seleciona o modelo principal. Uma antiga configuração
`ELEVEN_AUDIO_DIALOGUE_ENABLED=true` não ativa v3 enquanto o modelo selecionado
for Multilingual v2. Para habilitar v3 no futuro, é necessário selecionar
explicitamente `ELEVEN_AUDIO_MODEL=eleven_v3` e ativar a flag de diálogo.
`ELEVEN_AUDIO_FALLBACK_MODEL` só seleciona o TTS alternativo no modo v3;
não substitui o modelo principal v2.

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
