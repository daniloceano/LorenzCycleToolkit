# Etapa 3 — Níveis verticais: resultados da medição simples (27/11/2026)

**Natureza:** registro documental. Torna permanentes os resultados da etapa 3 (D-016), executada em 27/11/2026
conforme o [pré-registro](2026-11-27-stage-3-vertical-levels-preregistration.md). **Nada executável muda no
repositório:** nenhum código, teste, dependência, padrão de níveis ou CI. Dados, scripts, tabelas e selos ficam fora
do Git, na Swell. `main`, tags e releases continuam inalterados.

**Base:** `develop` em `754ec2e` (PR #92).

**Estado:** etapa 3 CONCLUÍDA. Amostra reduzida: 10 casos, Atlântico Sul, 1979–1980, caixa fixa. Não é a hipótese K.

## Resumo em linguagem simples

**Pergunta:** quanto muda o ciclo de energia de Lorenz no LECTK quando ele usa só os 13 níveis da grade ARCO em vez
dos 32 níveis ERA5 que aproveita por padrão (1000–10 hPa)?

**Método:** um download WEkEO por caso com os 37 níveis, nas caixas seladas da 2b; os conjuntos de níveis são
recortes do mesmo arquivo, escritos pelo mesmo script; o LECTK rodou em caixa fixa com resíduos (`-f -r`).

**Resultado:** na mediana dos 10 casos, a troca de 32 para 13 níveis muda os termos em 0,2 % a 5,7 % da sua escala.
Nenhum reservatório ou conversão inverte de sinal, e todos têm correlação temporal ≥ 0,99 com a referência. Dos 26
termos do CSV, 13 ficam na classe b (pequena), 12 na c (moderada) e 1 na d (relevante): o termo M, cuja magnitude é
próxima de zero. Os 13 níveis tendem a dar reservatórios um pouco menores (Az −3,1 %, Kz −1,6 %, Ae −1,5 %) e Ca
um pouco menor (−2,0 %). Em casos isolados, termos de fronteira e resíduos mudam bem mais (até 37–54 %).

## 1. Execução

- **Pasta (fora do Git):** `/p1-swell/danilocs/lectk_stage3_levels_20261127/`.
- **Código do LECTK:** clone de `develop` em `754ec2e4a924b8a7d0a6f9ec5ac18227e338a7eb`, árvore limpa antes e
  depois (D-040); SHA gravado em cada `run_manifest.json`. Ambiente `lorenz-4d` (D-043); para o download, venv
  `--system-site-packages` sobre ele com `hda==2.40`.
- **Scripts da etapa (fora do Git):** `scripts/stage3.py` (sha256 `2ab66775…`; download, recortes, LEC, C1) e
  `scripts/compare3.py` (`dbc0655e…`; C2, C3 e tabelas). A conversão para a grade comum e a consulta WEkEO são cópias
  das funções do harness da 2b (`bench2b.py`, `cc1538a0…`), lido sem alterar.
- **Download:** 13 consultas WEkEO (`EO:ECMWF:DAT:REANALYSIS_ERA5_PRESSURE_LEVELS`), 37 níveis, `t z u v w`,
  00/06/12/18 UTC, caixas seladas (`boxes.json`, selo conferido). 11:43–13:26 UTC, mais uma repetição às 16:23–16:25
  UTC (§4). Tráfego: 3,99 GB, contando a tentativa perdida (teto: 10 GB). Disco novo: 18,9 GB (teto: 20 GB).
- **Recortes:** 31 arquivos (L32, L29 e L13 dos 10 casos; L37 só do caso de controle). Todos completos, sem valores
  não finitos e com faixas plausíveis; o arquivo bruto tem 37 níveis em todos os casos.
- **LEC:** 43 execuções, todas com código de retorno 0: 30 dos conjuntos, 2 do C1, 1 do C2 e 10 do C3. Níveis efetivos
  pelo log: L32 → 32 (topo 10 hPa), L29 → 29 e L13 → 13 (topo 50 hPa); L37 → 37 na entrada, 5 descartados, 32 usados.
- **Selos:** `report/SHA256SUMS` com 1.493 arquivos (brutos, recortes, execuções, tabelas, scripts, estado e logs);
  sha256 do próprio arquivo `c5c3aad92d5ebb2407c8a984400ca2137f19ab66ca1a0438070a79dcd2b43b23`.
- **Itens retidos da 2b:** só lidos. Os 10 arquivos de 13 níveis do grupo {CDS, WEkEO} conferem com os sha256
  selados nos manifestos da fase 2.

## 2. Controles

| Controle | Resultado |
|---|---|
| C1, repetição (caso 19790205, L32 e L13 duas vezes) | PASS: CSVs idênticos byte a byte (L32 `7212f199…`, L13 `935ed785…`) |
| C2, recorte (L37 reescrito pelo script × L32) | PASS: CSV idêntico byte a byte (`7212f199…`) |
| C3, fonte (L13 novo × arquivo de 13 níveis da 2b, grupo {CDS, WEkEO}) | 10 de 10 casos: arquivos de entrada idênticos byte a byte; LEC com o mesmo código idêntico byte a byte; LEC idêntico byte a byte também ao CSV de referência da fase 2 da 2b (código `068c3a9`, ambiente `lorenz`) |

O C1 passou antes de qualquer outra execução do LEC. O C3 mostra que a nova entrega WEkEO, em 37 níveis, é a mesma
da 2b nos 13 níveis comuns.

## 3. Resultados

Escala `S = max(|média_ref|, desvio-padrão temporal_ref)`; `d = |Δ médio| / S`; classe pela mediana dos 10 casos
(pré-registro, §5). "Sentido" é o sinal mediano de Δ médio / |média_ref| (L13 menor ou maior que a referência).

### 3.1 Comparação principal: L32 → L13

| Termo | d mediano | Faixa de d | Sentido (mediana) | r mínimo | Inversões | r < 0,9 | Classe |
|---|---|---|---|---|---|---|---|
| Az | 3,08 % | 1,91–6,26 % | −3,1 % (10/10 menores) | 0,998 | 0 | 0 | c |
| Ae | 1,45 % | 0,30–2,73 % | −1,5 % (10/10) | 0,988 | 0 | 0 | b |
| Kz | 1,62 % | 0,16–16,60 % | −1,6 % (10/10) | 0,995 | 0 | 0 | b |
| Ke | 0,24 % | 0,04–3,67 % | −0,2 % (9/10) | 0,998 | 0 | 0 | b |
| Cz | 2,04 % | 0,48–3,67 % | −0,3 % (5/10) | 0,997 | 0 | 0 | c |
| Ca | 2,01 % | 0,28–7,22 % | −2,0 % (9/10) | 0,990 | 0 | 0 | c |
| Ck | 0,57 % | 0,22–5,19 % | +0,6 % (4/10 menores) | 0,991 | 0 | 0 | b |
| Ce | 0,66 % | 0,19–3,49 % | −0,5 % (9/10) | 0,998 | 0 | 0 | b |
| C_overturning | 1,74 % | 0,31–7,22 % | — | 0,996 | 1 | 0 | c |
| BAz | 1,84 % | 0,07–4,16 % | — | 0,996 | 0 | 0 | b |
| BAe | 1,97 % | 0,20–3,81 % | — | 0,997 | 0 | 0 | b |
| BKz | 2,46 % | 0,46–36,89 % | — | 0,989 | 0 | 0 | c |
| BKe | 1,35 % | 0,03–9,79 % | — | 0,997 | 0 | 0 | b |
| BΦZ | 2,96 % | 0,26–47,08 % | — | 0,856 | 0 | 1 | c |
| BΦE | 3,00 % | 0,73–29,61 % | — | 0,810 | 0 | 1 | c |
| M | 14,88 % | 0,97–66,07 % | — | 0,103 | 1 | 9 | **d** |
| Gz | 4,43 % | 0,02–16,29 % | −2,0 % (8/10) | 0,991 | 0 | 0 | c |
| Ge | 5,06 % | 0,19–11,99 % | +2,8 % (3/10) | 0,988 | 0 | 0 | c |
| ∂Az/∂t | 0,89 % | 0,37–2,61 % | — | 0,994 | 0 | 0 | b |
| ∂Ae/∂t | 0,89 % | 0,09–2,81 % | — | 0,993 | 0 | 0 | b |
| ∂Kz/∂t | 1,65 % | 0,10–9,82 % | — | 0,988 | 0 | 0 | b |
| ∂Ke/∂t | 0,40 % | 0,12–2,22 % | — | 0,997 | 0 | 0 | b |
| RGz | 5,74 % | 1,94–25,25 % | — | 0,969 | 0 | 0 | c |
| RKz | 2,35 % | 0,02–53,80 % | — | 0,919 | 0 | 0 | c |
| RGe | 5,61 % | 0,54–16,86 % | — | 0,984 | 1 | 0 | c |
| RKe | 1,66 % | 0,04–12,47 % | — | 0,980 | 0 | 0 | b |

Classes: a = 0, b = 13, c = 12, d = 1.

### 3.2 Decomposição: topo (L32 → L29) e discretização (L29 → L13)

| Termo | d mediano, topo | d mediano, discretização | d mediano, total |
|---|---|---|---|
| Az | 0,41 % | 2,03 % | 3,08 % |
| Ae | 0,32 % | 0,87 % | 1,45 % |
| Kz | 1,77 % | 0,44 % | 1,62 % |
| Ke | 0,66 % | 0,36 % | 0,24 % |
| Ca | 0,09 % | 1,89 % | 2,01 % |
| Ce | 0,17 % | 0,50 % | 0,66 % |
| Ck | 0,42 % | 0,83 % | 0,57 % |
| Cz | 0,63 % | 0,89 % | 2,04 % |
| Gz | 0,35 % | 4,93 % | 4,43 % |
| Ge | 0,40 % | 4,34 % | 5,06 % |
| M | 1,75 % | 13,59 % | 14,88 % |

Classes: topo — b = 21, c = 4 (C_overturning, BΦZ, RKz, RGe), d = 1 (M); discretização — b = 17, c = 8, d = 1 (M).
As medianas não se somam (são medianas de |Δ| sobre casos diferentes e os efeitos podem ter sinais opostos).

### 3.3 Interpretação **[PRELIMINAR — amostra reduzida]**

- **A discretização domina** Az, Ca, Gz, Ge e os resíduos de geração: na mesma camada 1000–50 hPa, os 13 níveis dão
  Az menor em 10 de 10 casos (−2,0 %) e Ca menor em 9 de 10 (−1,9 %).
- **O topo (10–50 hPa) domina Kz:** cortar os níveis acima de 50 hPa reduz Kz em 10 de 10 casos (−1,8 %) e Ke em
  10 de 10 (−0,7 %); a energia cinética dessa camada deixa de entrar no volume de controle.
- **Ck e Ce** mudam menos de 1 % na mediana; são os termos mais estáveis, como na LZV.
- **M** (d ≈ 15 %, r mediano 0,50) é um termo de magnitude ~10⁻⁵; a normalização é dominada por ruído e a classe d
  diz pouco sobre o ciclo. Fica como ressalva, não como achado físico.
- **Casos isolados** com fronteiras e resíduos de 30–54 % (BKz, BΦZ, BΦE, RKz) mostram que a mediana esconde casos
  sensíveis; os termos de fronteira de geopotencial são os mais afetados.
- **Antecedente da LZV:** no caso 19791150, aqui em caixa fixa, Ca muda −0,3 % (L32 → L13), contra −5,1 % na LZV
  (modo móvel, 29 → 13 níveis). Configurações diferentes; não se comparam diretamente.
- **Sem generalização** para outras regiões, estações ou o modo móvel. Mudar o padrão de níveis do LECTK está fora
  do escopo; se for considerado, vira recomendação para decisão separada.

## 4. Desvios e interpretações (registrados antes do cálculo, salvo o 2)

1. **13 consultas, não 10 pedidos.** Três casos atravessam a virada do mês (19790644, 19791150, 19800256) e a
   consulta é por mês, como na 2b. Continua um pedido por caso e por mês.
2. **Arquivo truncado, sem erro do cliente.** A primeira entrega do caso 19790166 tinha 300.292.096 de 521.610.410
   bytes (HDF5 truncado); `hda` 2.40 não acusou erro. Detectado ao abrir o arquivo antes dos recortes; a consulta foi
   repetida (arquivo completo) e a tentativa perdida ficou em `state/download.json` (`previous_attempts`).
3. **Faixa de plausibilidade de z:** limite superior de 2,2 × 10⁵ (2b, topo 50 hPa) para 5,0 × 10⁵ m² s⁻² (1 hPa).
4. **L37 só no caso de controle** (19790205, escolhido antes de qualquer cálculo), por causa do teto de disco.
5. **Termos:** o CSV do LECTK não tem Dz/De; tem RGz, RKz, RGe, RKe, C_overturning e M. Foram usadas todas as 26
   colunas.
6. **Classe b com uma inversão:** "sem inversão de sinal" foi lido como nenhuma inversão nos 10 casos; com d ≤ 2 % e
   uma inversão, o termo vai para c (C_overturning e RGe).
7. **C3 ampliado:** além da comparação dos campos, o LEC rodou com o código desta etapa sobre os 10 arquivos retidos
   da 2b, e os CSVs foram comparados com os da fase 2.

## 5. Fechamento (D-044)

Aprovado pelo responsável, PR #93.

## 6. Nota de correção (27/11/2026)

O texto acima fica como foi aprovado; esta nota corrige duas frases do resumo e acrescenta uma ressalva, à luz da tabela do §3.1:

- "Nenhum reservatório ou conversão inverte de sinal": vale para Az, Ae, Kz, Ke, Ca, Ce, Ck e Cz. **C_overturning**,
  listado com as conversões na tabela, inverte o sinal da média em **1** dos 10 casos.
- "Todos têm correlação temporal ≥ 0,99" (reservatórios e conversões): o r mínimo de **Ae** é **0,988** (e o de Ck, 0,991, arredonda para 0,99).
  O correto é "r ≥ 0,988".
- A faixa "0,2 % a 5,7 %" da mediana **exclui o termo M** (mediana de 14,9 %, classe d); a mesma ressalva vale para a
  linha de 27/11/2026 em [overview.md](../../overview.md).

Encerramento e retenção: [registro de retenção da etapa 3](2026-11-27-stage-3-retention.md).
