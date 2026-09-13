# Comparação do diálogo — 13/09/2026

Teste real da API com `gpt-6-astra`, usando exclusivamente os dois estudos e a
referência complementar **fictícios** de
[`podcast_conversation.json`](../tests/fixtures/podcast_conversation.json).
Os trechos abaixo não são evidência clínica nem um episódio para publicação.

As versões antiga e final receberam as mesmas fontes, meta de 740 palavras,
schema, modelo e parâmetros. Cada versão executou planejamento, escrita e
auditoria, sem precisar da tentativa adicional de ajuste da ressalva. Entre elas,
uma versão intermediária foi avaliada e refinada porque ainda continha respostas
longas com formato de exposição. São observações de um caso, não um benchmark.

| Medida descritiva | Prompt anterior | Prompt final |
| --- | ---: | ---: |
| Palavras do episódio | 669 | 664 |
| Falas | 15 | 25 |
| Palavras por frase, média aproximada | 14,2 | 10,5 |
| Maior frase, palavras | 31 | 21 |
| Frases acima de 25 palavras | 2 | 0 |

A redução no tamanho das frases não prova, sozinha, naturalidade. A leitura dos
trechos mostra também perguntas sobre algo recém-discutido, retomadas específicas
e conclusões construídas em mais de uma troca. A abertura e algumas sínteses
ainda podem ser refinadas pela avaliação humana do episódio real.

## Antes — trecho gerado com o prompt anterior

> Manu: Os relatos mostram aspectos diferentes da mesma experiência. Não precisar se deslocar facilitou o acesso. Ao mesmo tempo, conexões de internet instáveis interromperam a participação, e alguns participantes queriam mais contato com o terapeuta. A conveniência do atendimento em casa não elimina essas dificuldades relatadas.
>
> Ivo: Isso ajuda um serviço a formular perguntas sobre conectividade e contato durante o programa, sem prometer melhora clínica. O estudo não mediu eficácia comparativa, nem permite calcular com que frequência essas percepções aparecem entre todas as pessoas com DPOC. Seu resultado é uma descrição situada, não uma contagem populacional.
>
> Manu: Ouvir apenas concluintes de um centro não revela por que outros abandonaram nem demonstra benefício clínico. As razões de abandono não foram investigadas. Para aproveitar esses relatos em outro serviço, é preciso considerar o contexto local: internet e contato com o terapeuta são pontos a explorar, não uma receita universal.

## Depois — trecho gerado com o prompt final

> Ivo: Essa facilidade resolve uma barreira concreta. Mas realizar a reabilitação em casa eliminou as dificuldades de participação?
>
> Manu: Não. Conexões instáveis de internet interromperam a participação. Alguns também queriam mais contato com o terapeuta.
>
> Ivo: Então evitar deslocamento não basta para descrever a experiência como fácil. E esses problemas apareceram mesmo entre quem concluiu. Sabemos o que aconteceu com quem abandonou?
>
> Manu: A amostra incluiu apenas concluintes de um único centro. Os motivos de abandono não foram investigados.

Aqui, a observação sobre internet e contato leva Ivo a distinguir acesso de
experiência de participação. A informação de que todos concluíram o programa
abre a pergunta sobre abandonos, e Manu delimita o alcance dos dados. O desenho
da troca preserva o fato de que as razões de abandono não foram investigadas;
não o transforma em abandono causado pela internet.

## Conferência dos fatos e do contrato

Na versão final, o ensaio sintético durante quimioterapia manteve os 80
participantes, duas sessões semanais por oito semanas, diferença de −0,4 ponto
e IC95% de −1,1 a +0,3. Preservou as 12 avaliações finais incompletas, a ausência
de informação sobre seu tratamento e a distinção entre não demonstrar diferença
e demonstrar equivalência. A referência após quimioterapia manteve população,
dose, duração e intervalo distintos; a conversa não atribuiu o contraste a uma
causa testada. O estudo qualitativo continuou restrito a concluintes de um centro,
sem estimar frequência populacional ou eficácia clínica comparativa.

As duas auditorias automatizadas retornaram `issues: []`. Isso não substitui a
avaliação humana: os dados são sintéticos, e não foi realizado teste auditivo ou
geração de áudio neste comparativo. Não foi feito upload, aprovação de episódio,
envio de e-mail ou publicação no RSS.

Passaram 36 testes com:

```
./.venv/bin/python -m unittest test_podcast_editorial test_podcast_pubmed_integration test_podcast_voice_settings
```

A regressão de integração agora inclui respostas curtas, uma ressalva distribuída
entre duas falas e falas sucessivas do mesmo locutor. Verifica que o texto aprovado
chega exatamente ao TTS e que nenhuma nova chamada de redação é feita na síntese.
Os testes também cobrem o orçamento de chamadas editoriais, fontes, bloqueios
científicos, e-mail, geração parcial e configurações das vozes.

## Reprodução e artefatos

O comparativo pode ser repetido com
[`compare_podcast_conversation.py`](../scripts/compare_podcast_conversation.py).
O parâmetro `--only after` permite testar um refinamento sem repetir a geração
antiga. A execução é paga, apenas textual, e grava fontes, prompts, hashes,
respostas completas, auditoria e consumo real em uma pasta nova.

Artefatos locais desta execução:

- Anterior e intermediário: `/private/tmp/revacast-dialogue-comparison-20260913-live`.
- Final: `/private/tmp/revacast-dialogue-comparison-20260913-refined`.

Hashes SHA-256 do prompt de escrita:

- Anterior: `d855dfdf1da48bc34f5705f242f31fe1832d0130d64fe984bac0325feb0222dc`.
- Final: `3c2eb6e940fe8675f1310a982acdd92e5ee5bbf92f29ce1d5fa1e960a3f23006`.

A mudança de estilo segue a orientação oficial de especificar explicitamente
o estilo e revisar instruções concorrentes no
[guia do GPT-6 Astra](https://developers.openai.com/api/docs/guides/latest-model).
O modelo e as configurações de áudio foram mantidos. A aplicação da mudança
depende do deploy do backend e da geração de um novo roteiro; arquivos já
aprovados continuam vinculados ao texto que foi revisado.
