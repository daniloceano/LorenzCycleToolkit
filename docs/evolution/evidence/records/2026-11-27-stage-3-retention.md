# Etapa 3 — retenção dos dados na Swell (GOV-RET-01) e encerramento (27/11/2026)

**Estado operacional: ENCERRADA — RETENÇÃO VALIDADA.**

**Natureza:** só documentação no repositório; nenhum código, teste, dependência, padrão do LECTK ou CI mudou. Na Swell,
as exclusões foram feitas exclusivamente pelos manifestos aprovados. `main`, tags e releases ficam inalterados.

**Base:** `develop` em `b9f8ea4` (PR #93). Branch única do encerramento (D-042): `docs/stage-3-closure-2026-11-27`.

## Resumo em linguagem simples

A medição da etapa 3 tinha deixado **19,11 GB** na pasta `/p1-swell/danilocs/lectk_stage3_levels_20261127/`. Foram
**apagados 19,09 GB**: os downloads brutos da WEkEO, os arquivos recortados por conjunto de níveis, a venv do cliente
`hda` e o clone do código. Ficam **15,2 MB**: os resultados do LEC (CSVs, logs, manifestos de execução), os scripts,
as tabelas, o estado do download e os selos. Os números da etapa podem ser conferidos com o que ficou; para refazer o
cálculo a partir dos dados é preciso baixar de novo da WEkEO (§6).

## 1. Escopo e autorização

- **Política:** [GOV-RET-01](../../retention-policy.md); D-034 (Swell com consulta prévia), D-042 (branch única) e
  D-044 (fechamento).
- **Pasta:** só `/p1-swell/danilocs/lectk_stage3_levels_20261127/`. Os retidos da 2b, da etapa 4 e da LZV não foram
  tocados; dois arquivos selados da 2b foram apenas lidos para conferir entradas (§2).
- **Parte 1:** inventário só de leitura, sem seguir links. Manifesto apresentado ao responsável.
- **Parte 2:** **aprovação expressa do responsável em 27/11/2026**, incluindo o descarte de `code/`, que estava em
  INVESTIGAR.

## 2. Manifesto

`MANIFESTO_EXCLUSAO_stage3.tsv`, na raiz da pasta (fora da área apagada), somente leitura;
sha256 `2e9ab8c1335439b16db794844976eb37013630a1d1ce366d9bf34e771ef9d602`. Gerado por
`scripts/build_manifest_stage3.py` (`cd12ab40…`). Todo item a descartar tem linha própria com sha256; os retidos
menores que 1 MB aparecem agrupados por contagem e bytes.

| Classe | Entradas | Bytes | Conteúdo |
|---|---|---|---|
| DESCARTAR | 1.005 arquivos, 4 links, 172 pastas | 18.933.680.589 | 13 brutos WEkEO (3.690.860.574 B), 31 recortes L37/L32/L29/L13 (15.231.529.857 B), venv do `hda` (875 arquivos e 4 links), `runs/*/inputs/` (86 arquivos) |
| RETER até o Marco B | 1.365 arquivos | 14.728.931 | `runs/` sem `inputs/` (43 CSVs de resultado, 989 CSVs por nível, `run_manifest.json`, `lec_metadata.json`, `run_status.json`, `domain_selection.csv`, logs), `scripts/`, `report/` (inclui `SHA256SUMS`), `state/`, `logs/` |
| INVESTIGAR | 489 arquivos | 158.368.905 | `code/lectk_develop`: clone de `develop` em `754ec2e`, com `.git` |

- **`runs/*/inputs/`:** recriáveis. Os 43 `namelist` têm sha256 idêntico ao do original da 2b
  (`lectk_bench2b_phase2_20261119/runs/ref/cds_19790166/inputs/namelist`); os 43 `box_limits` têm conteúdo idêntico
  ao derivado das caixas seladas (`lectk_bench2b_20261109/boxes/boxes.json`). Nenhuma divergência.
- **`boxes.json`:** a pasta não tinha cópia; os scripts leram o original selado da 2b, que continua retido.
- **`code/` (INVESTIGAR → DESCARTAR por decisão do responsável):** como o manifesto principal só tinha a contagem dos
  arquivos pequenos de `code/`, foi gerado um adendo com uma linha e um sha256 por arquivo, depois de conferir que o
  clone estava limpo em `754ec2e`: `MANIFESTO_EXCLUSAO_stage3_code.tsv`, sha256
  `dd2cc76aee8bfda1c84098142fdbc4f183d77d96ac4cdba996b24f0d1cc97599` (489 arquivos, 70 pastas, 158.368.905 B).

## 3. Execução (27/11/2026, 20:55–20:56 UTC)

Script `scripts/execute_exclusion_stage3.py` (`24ebd394…`):

1. conferiu os selos dos dois manifestos;
2. para cada arquivo: tipo regular, tamanho e sha256 calculado imediatamente antes de apagar; para cada link: tipo e
   alvo (nunca seguido);
3. pastas por último, das mais profundas para as mais rasas, só se vazias.

**Resultado:** 1.498 arquivos e links e 242 pastas removidos, **19.092.049.494 bytes** (soma exata dos dois
manifestos), **sem nenhuma divergência**. `du` da pasta: 19.107.052.617 B antes, 15.154.889 B depois; espaço livre de
`/p1-swell` aumentou 19.096.477.696 B. Log em `logs/exclusion_stage3.log`; resumo em `state/exclusion_stage3.json`.

## 4. O que ficou (até o Marco B, D-013)

| O que | Onde (`lectk_stage3_levels_20261127/`) | Arquivos | Bytes |
|---|---|---|---|
| Execuções do LEC (sem `inputs/`) | `runs/` | 1.333 | 13.786.388 |
| Scripts da etapa e do encerramento | `scripts/` (`stage3.py`, `compare3.py`, `launch_lectk.py`, `build_manifest_stage3.py`, `execute_exclusion_stage3.py`) | 5 | 39.472 |
| Tabelas, controles e selos | `report/` (`summary*`, `table_*`, `C1–C3.json`, `levels_effective.json`, `SHA256SUMS`) | 12 | 337.203 |
| Estado do download, recortes e exclusão | `state/` | 4 | 74.876 |
| Logs | `logs/` | 14 | 497.260 |
| Manifestos e selos | raiz | 4 | ≈ 0,42 MB |

O `report/SHA256SUMS` continua listando os arquivos apagados; serve para conferir um novo download (§6).

## 5. Fora deste encerramento

Retidos da 2b, da etapa 4 e da LZV (inclusive as caixas seladas e os arquivos de 13 níveis usados no C3): inalterados.

## 6. Ressalva de reprodução

Os dados de entrada **não ficam guardados**. Refazer o cálculo exige um novo download WEkEO com as consultas
registradas em `state/download.json` (13 consultas, 37 níveis, caixas seladas da 2b) e o código em `754ec2e`. Na
etapa 3, a entrega WEkEO foi idêntica, byte a byte nos 13 níveis comuns, à da 2b (C3), mas a identidade de uma entrega
futura não é garantida; um novo download deve ser conferido contra os sha256 de `report/SHA256SUMS` e
`state/download.json`.

## 7. Fechamento (D-044)

_(preenchido no commit de encerramento: aprovado pelo responsável, PR #N)_
