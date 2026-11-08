# Frente 2a — validação externa com Dias Pinto & da Rocha (2011) (08/11/2026)

**Natureza:** registro documental. Nenhum código, teste, CI ou arquivo em `inputs/` muda. **Base:** `develop` em
`6433b932988db9323b005f91f085bf74efe57d37`. Implementa D-014 ([replanejamento](2026-11-07-replanning-decisions.md)).

## 1. Objetivo e método

Comparar o LECTK com os diagramas médios por fase (Figs. 6, 10 e 14) de Dias Pinto & da Rocha (2011),
*J. Geophys. Res.*, 116, D14112, doi:10.1029/2011JD016217, para os ciclones Reg1, Reg2 e Reg3.

- Dados: NCEP/NCAR R1, 2,5°, 6 h.
- Caixas fixas das legendas da Fig. 2; fases das legendas das Figs. 6, 10 e 14
  ([arquivos](../frente-2a/README.md)).
- Topo efetivo de 100 hPa: o ômega do R1 tem só 12 níveis (1000–100 hPa).
- Sensibilidade de topo: 17 níveis com ômega = 0 de 70 a 10 hPa; envelope E = |principal − sensibilidade|.
- Transcrição dos diagramas do artigo independente das saídas do LECTK.

## 2. Pré-registro (aprovado pelo responsável em 08/11/2026, texto integral)

1. 16 termos por diagrama (∂AZ/∂t, ∂AE/∂t, ∂KZ/∂t, ∂KE/∂t; CZ, CA, CE, CK; BAZ, BAE, BKZ, BKE; RGZ, RGE, RKZ, RKE),
   média por fase, 3 casos = 144 valores, em W m⁻².
2. Convenções: CZ AZ→KZ; CA AZ→AE; CE AE→KE; CK KE→KZ. Transcrição: positivo no sentido do painel (a);
   B e RG positivos se a seta entra na caixa; RK positivo se entra em K.
3. Transcrição dupla e independente; em cada caixa, |tendência − soma das setas| ≤ 0,1. Resolver as
   discrepâncias ANTES de abrir as saídas do LECTK.
4. Envelope de topo E por termo = |principal − sensibilidade (omega = 0 de 70 a 10 hPa)|.
5. Concordância: mesmo sinal e |LECTK − artigo| ≤ max(0,25·|artigo|; 0,5) + E.
   Discordância explicada: fora disso, mas coberta por E ou pela troca da derivada centrada por diferença avançada.
   Falha: sinal oposto com |artigo| ≥ 1, ou fora da faixa sem explicação. |artigo| < 0,5: só o sinal, informativo.
6. PASS numérico: ≥ 80% dos termos com |artigo| ≥ 1 em concordância ou discordância explicada E mecanismos
   dominantes iguais (Reg1 formação: CK < 0 e relevante; Reg2 maturidade: CE dominante; Reg3: CA e CE dominantes).
   FAIL: < 60% ou mecanismo dominante invertido. Caso contrário, inconclusivo.
7. Controles (Reg1 com a caixa +15°E; Reg3 com dados de 2007) devem ser majoritariamente falha; se não forem,
   o veredito global é inconclusivo.
8. PORTÃO FINAL: o veredito numérico é relatado como calculado; a aprovação final depende da inspeção visual
   do responsável. Divergências entre os dois são registradas, nunca corrigidas.

**Interpretações operacionais aceitas pelo responsável** ([regras congeladas](../frente-2a/REGRAS_OPERACIONAIS_P5.md),
escritas antes de abrir as saídas):

1. Classificação na ordem literal do item 5: sinal oposto com |A| ≥ 1 é falha sem possibilidade de explicação.
2. “Relevante” (Reg1, CK) significa |CK| ≥ 1 W m⁻².
3. O mecanismo do Reg3 (sem fase indicada) é exigido nas três fases.
4. Os controles usam o E do caso correspondente; a regra do item 7 vale também sobre FAIL.

## 3. Desvio de ordem

As execuções do LECTK (casos, sensibilidade e controles) foram feitas **antes** da aprovação do pré-registro.
O agente declarou por escrito que não calculou médias por fase nem comparou com as Figs. 6, 10 e 14 antes da
aprovação. Antes disso, só leu nomes de colunas e contagens de NaN das saídas, e não usou digitalizações
anteriores do Reg1. As regras operacionais e seus hashes foram congelados antes da primeira leitura de valores.

## 4. Resultado numérico: **FAIL global**

| Escopo | Aceitáveis (\|A\| ≥ 1) | % | Mecanismo dominante | Veredito |
|---|---|---|---|---|
| Reg1 | 11/27 | 40,7 | igual | FAIL |
| Reg2 | 27/34 | 79,4 | igual | INCONCLUSIVO |
| Reg3 | 15/28 | 53,6 | igual | FAIL |
| **Global** | **53/89** | **59,6** (limiar 60) | 3/3 iguais | **FAIL** |

Controles: falha em 92,6 % (Reg1 +15°E) e 64,3 % (Reg3 com 2007), ambos majoritariamente falha.

Por grupo (aceitáveis/termos com |A| ≥ 1): conversões 22/24; tendências 18/24; BA 4/9; RG 5/12; BK 1/6; RK 3/14.

## 5. Interpretação

- **Fronteira × resíduo:** os erros dos termos de fronteira e dos resíduos se compensam (correlação dos erros
  LECTK − artigo nas 9 fases: r ≈ −0,98 no ramo zonal, BAZ×RGZ e BKZ×RKZ; r ≈ −0,91 em KE, BKE×RKE). As
  tendências e as conversões fecham. Portanto, a **partição entre fronteira e resíduo diverge**. A causa não é
  determinável: o artigo não especifica os níveis, o tratamento do ômega acima de 100 hPa nem a formulação
  discreta dos termos de fronteira.
- **Conversões:** mecanismos reproduzidos nos três casos; magnitudes próximas em Reg1 e Reg2. No Reg3, CA e CE
  são sistematicamente 19–25 % menores no LECTK, em todas as fases, no limite da tolerância.
  **Hipótese não testada:** níveis abaixo do relevo andino na caixa do Reg3. Tema excluído do escopo
  atual (D-015); registrado como evidência caso o responsável decida reabri-lo.

## 6. Limitações

- Resultado no limiar: um termo a mais aceitável daria 60,7 % (INCONCLUSIVO).
- O controle do Reg3 discrimina pouco (35,7 % aceitável, contra 53,6 % no caso real).
- Comparar com outra implementação não é validação física independente; concordância ou divergência não diz
  qual das duas está certa.
- Portão visual (item 8): inspeção do responsável em 08/11 concorda que as divergências principais estão nas
  fronteiras; registra desconforto com diferenças nas conversões, em especial o viés sistemático de CA e CE no Reg3.

## 7. Verificação independente

A conversa-filha informou: reclassificação dos 144 termos 144/144 igual; fechamento da transcrição ≤ 0,07 W m⁻²
nas 36 caixas. As duas transcrições (sessão principal e subagente isolado) coincidiram em 144/144 valores.

## 8. Evidências externas (fora do Git, retenção sob GOV-RET-01)

Manifesto consolidado `MANIFESTO_SHA256.tsv` (438 arquivos), SHA-256
`1b268bffb12e5fba8342ecd0bbcbff5305a697c497841664c03f0fdc647c6916`. Principais itens:

| Arquivo | SHA-256 |
|---|---|
| `STATUS_2026-11-08.md` | `e749769632941e6aef4f73bcddbea71ab74d84779ff8732e4986e66c4c9a5e57` |
| `REGRAS_OPERACIONAIS_P5.md` | `248d57ad1fbbd1b4f30b969d9fe3e163d50a55800ce611af9f4631a945ed66d9` |
| `logs/p5_congelamento.txt` | `9c215273f3cc6fb1170f7ed88149f2cd547fefce2483078e99e5ee5c7f8d37fd` |
| `resultados/p5_tabela_144.csv` | `a055ded1ee4dd819a46a00da42a508d6d8d034755bce9d91c2686f86c7895271` |
| `resultados/p5_vereditos.md` | `2ed383e6482ca114b7124065438f637e885b4118a11773f77ab554f467f717a5` |
| `resultados/p5_controles.csv` | `b245ca85e16e07078d8bea9f12674efe2d157af8a38c8e713b1fcef1ac518dc2` |
| `resultados/p5_robustez.txt` | `f6dbf72478c593da1a023f27ea3084b80b4bdbed72adf40efca56bd1c1c30575` |
| `transcricoes/transcricao_A.csv` (= `transcricao_final.csv`) | `1a8689e48809da687d8e46625ab69230efcb034f57a691e7eff0d510b64a7b63` |
| `transcricoes/transcricao_B.csv` | `cc926535a684d7e469b922a66952f99a2a7f779198a63ad37993a7b8e751693c` |
| `transcricoes/teste_fechamento.csv` | `dc70a47dda93bd342f82aafce79f38082da3b306a33a4caaa3063fccf3b878aa` |
| `scripts/p5b_analise.py` | `e2300a42871686f1e46aa19e1ca9f3187af4f1eb7bebd83b356c88d8eb9aaacf` |
| `scripts/p5c_figuras.py` | `e7d3891fb9f9d6871dde29a261085d42a173fd5159b3116d44b1782caf43cd68` |

Figuras, recortes do artigo, dados R1, logs e saídas brutas não são versionados.

## 9. Próximos passos propostos (não abertos)

- Verificação interna dos termos de fronteira do LECTK.
- Veiga (2008) e Revegoca (2026) apenas como checagem de coerência; o código de Revegoca não está disponível.
