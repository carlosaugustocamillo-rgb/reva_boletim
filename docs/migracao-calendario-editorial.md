# Migração do calendário editorial Reva+

Data da auditoria: 2026-10-04

## Objetivo

Substituir o calendário mecânico de 299 linhas por um calendário ativo contendo somente:

1. conteúdos efetivamente enviados;
2. pautas novas aprovadas;
3. pautas planejadas para os próximos ciclos.

Temas antigos não executados deixam de ser agenda. Eles continuam disponíveis no histórico do Git e no banco de pautas para comparação editorial.

## Reconciliação das fontes

O Firestore registrava 46 índices concluídos em um calendário de 299 linhas, mas os índices foram acumulados durante diferentes versões do CSV. Como um mesmo índice apontava para títulos diferentes ao longo do histórico, o Firestore não foi usado isoladamente para reconstruir as publicações.

A fonte autoritativa adotada foi o Mailchimp:

- 49 campanhas enviadas para a audiência Reva+;
- 3 reenvios ou duplicações exatas;
- 46 campanhas únicas concluídas;
- cada campanha recebeu como ID estável o identificador `mc-<campaign_id>`;
- data histórica baseada no envio registrado pelo Mailchimp, convertida para `America/Sao_Paulo`.

O arquivo `calendario_editorial_150_semanas.csv` contém as 46 campanhas reconciliadas com status `Concluído`.

## Oito pautas originais do novo ciclo

| ID | Categoria | Decisão | Título de trabalho |
|---|---|---|---|
| 2026-10-03-D1 | Diabetes | Incluir | Uso caneta para diabetes e fico enjoado depois de comer: posso treinar mesmo assim? |
| 2026-10-03-D4 | Diabetes | Incluir | Apareceu uma bolha ou ferida no pé: posso continuar caminhando para controlar a glicemia? |
| 2026-10-03-D5 | Diabetes | Incluir reformulada | Dor abdominal intensa e vômitos usando caneta para diabetes: por que o treino deve esperar? |
| 2026-10-03-C1 | Cardiovascular | Incluir | Tenho insuficiência cardíaca e ganhei 2–3 kg em poucos dias: faço a caminhada ou aviso minha equipe? |
| 2026-10-03-C3 | Cardiovascular | Incluir | Na perimenopausa, palpitação e cansaço no treino são “só hormônios” ou merecem avaliação? |
| 2026-10-03-R1 | Respiratória | Incluir | Tenho Covid longa e pioro 12–48 horas depois do esforço: devo insistir no treino? |
| 2026-10-03-K1 | Câncer | Incluir | A dor articular do inibidor de aromatase está atrapalhando meu treino: como continuar com segurança? |
| 2026-10-03-K2 | Câncer | Incluir | Tenho port-a-cath: depois que cicatriza, posso movimentar o braço e fazer musculação? |

## Reavaliação das 12 pautas similares

O critério inicial comparou os candidatos com todas as 299 pautas antigas. Na migração, uma pauta antiga nunca publicada não bloqueia automaticamente uma pergunta atual e melhor formulada.

### Recuperadas para inclusão

| ID | Categoria | Título |
|---|---|---|
| 2026-10-03-D3 | Diabetes | Minha glicemia subiu depois da musculação: isso significa que o treino fez mal? |
| 2026-10-03-C2 | Cardiovascular | Meu relógio mede pressão arterial: posso usar esse número para decidir se treino hoje? |
| 2026-10-03-C4 | Cardiovascular | Tenho cardiopatia e o dia passou de 40 °C: devo reduzir ou adiar minha caminhada? |
| 2026-10-03-C5 | Cardiovascular | Uso anticoagulante e bati a cabeça numa queda: posso só observar e voltar a treinar? |
| 2026-10-03-R2 | Respiratória | A fumaça das queimadas está forte: quem tem asma ou DPOC deve trocar o treino ao ar livre por treino dentro de casa? |
| 2026-10-03-R4 | Respiratória | Depois de uma gripe ou piora da asma, quando é seguro voltar a treinar? |
| 2026-10-03-R5 | Respiratória | Tossi sangue depois do exercício: quando devo procurar atendimento imediatamente? |
| 2026-10-03-K3 | Câncer | Tenho metástase óssea: posso fazer musculação sem aumentar o risco de fratura? |

### Recuperadas com recorte mais específico

| ID | Categoria | Título revisado | Diferença editorial |
|---|---|---|---|
| 2026-10-03-R3 | Respiratória | O oxímetro caiu durante o treino, mas meu dedo estava frio: posso confiar no número? | O conteúdo enviado tratou queda de SpO₂; esta pauta trata confiabilidade da medida e artefatos. |
| 2026-10-03-K5 | Câncer | Já tenho linfedema: a manga compressiva é obrigatória durante a musculação? | O conteúdo enviado tratou força e risco de linfedema; esta pauta aborda a decisão individual sobre compressão. |

### Mantidas bloqueadas

| ID | Categoria | Título | Conteúdo enviado correspondente |
|---|---|---|---|
| 2026-10-03-D2 | Diabetes | Estou emagrecendo com caneta para diabetes: musculação ajuda a preservar massa muscular? | Perda de peso e ganho de massa magra no DM2: papel da musculação |
| 2026-10-03-K4 | Câncer | No dia da quimioterapia, é melhor treinar antes, depois ou não treinar? | Exercício durante quimioterapia; fadiga em dia de quimioterapia; náusea e exercício |

## Resultado

- 46 campanhas concluídas reconciliadas pelo Mailchimp.
- 8 pautas originais do novo ciclo.
- 10 das 12 pautas similares recuperadas.
- 2 pautas bloqueadas por repetição real.
- 18 pautas planejadas entre 2026-10-06 e 2026-12-04.
- 64 linhas no calendário ativo.

## Arquivos persistentes

- `banco_de_pautas.csv`: preserva os 20 candidatos e a decisão original da rodada; as colunas de revisão registram a decisão final, o título efetivamente adotado, a data e o formato planejados.
- `pautas_pendentes.csv`: contém somente as 18 pautas aprovadas, na mesma ordem cronológica do calendário.
- `calendario_editorial_150_semanas.csv`: contém 46 publicações concluídas e as mesmas 18 pautas planejadas.

Os arquivos originais em `Downloads` não foram alterados.

## Migração do estado

O código passou a aceitar `completed_ids`, mantendo compatibilidade temporária com `completed_indices`. No primeiro carregamento após a publicação conjunta do código e do novo CSV, as linhas marcadas como `Concluído` serão a fonte autoritativa e o documento `revamais_state/editorial_calendar` será recalculado automaticamente com:

- `completed_ids`: os 46 IDs iniciados por `mc-`;
- `completed_indices`: índices de 0 a 45;
- `next_index`: 46;
- `total_rows`: 64;
- `next_title`: primeira pauta planejada;
- `next_format`: formato da primeira pauta planejada.

O Firestore não foi alterado durante a preparação local da migração.
