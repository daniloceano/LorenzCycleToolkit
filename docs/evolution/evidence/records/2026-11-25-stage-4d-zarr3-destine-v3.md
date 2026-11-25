# Etapa 4d — ambiente com zarr 3 e DestinE no store Zarr v3 (25/11/2026)

**Natureza:** código, dependências, testes, documentação de usuário e registros. Cálculo do LEC inalterado
(`src/analysis`, `src/frameworks` sem mudança). `main`, tags e releases inalterados.

**Base:** `develop` em `91dfbed` (PR #87); `main` em `48fe8e2`. Branch única da etapa: `feat/stage-4d-zarr3` (D-042).

**Estado deste registro:** Partes 1 a 3 concluídas (regressão sem download e testes reais com o commit `f923d80`,
D-040). Encerrado em 25/11/2026: PR #88, merge `436279b`, CI verde; ver §9.

## Resumo em linguagem simples

O ambiente do LECTK passa a usar zarr 3 com um xarray novo, porque o store da DestinE que continua sendo atualizado
(Zarr v3) só é lido com zarr 3, e a operacionalização futura do LEC precisará desse ambiente de qualquer modo. Só a
camada de leitura de dados foi atualizada; o núcleo numérico (numpy, pandas, scipy, netCDF4) ficou igual, para que o
LEC dê resultados idênticos bit a bit sobre os mesmos arquivos. `--destine` passa a ler o store v3 por padrão; o
store v2 da 2b continua disponível com `--destine-store v2` enquanto for servido. No teste real, os dois stores
deram **valores idênticos bit a bit** e o mesmo LEC: os dados não mudaram, só os blocos. Os metadados do v3
declaram arredondamento para 13 bits de mantissa, mas os valores têm os mesmos 10 bits do v2 (§3, §7).

## 1. Decisão D-043 (responsável, 25/11/2026)

Opção A: atualizar as dependências para zarr 3 e ler a ERA5 da DestinE no store Zarr v3 atual. O store v2 está
congelado desde 06/07/2026 e pode ser retirado sem aviso; a operacionalização do LEC precisará do ambiente
atualizado. Escolhas aprovadas na Etapa 0 (25/11/2026, "Sim, as 4 recomendações"):

1. conjunto conservador de versões (só a camada de E/S muda);
2. destino do v2: mantido como alternativa fora do padrão (`--destine-store v2`), v3 padrão; sai quando a Earth Data
   Hub retirar o v2;
3. zarr, numcodecs e aiohttp no `requirements.txt`, para o CI rodar os testes Zarr (revoga a escolha 6 da 4b·2 e da
   4b·3);
4. sondagem dos metadados do v3 na Parte 2 (≈ 10 requisições, < 1 MB).

## 2. Versões

| Pacote | Mac, CI e Swell antes (`lorenz`) | Depois (`requirements.txt`, Swell `lorenz-4d`) | Motivo |
|---|---|---|---|
| Python | 3.12.12 (Mac, Swell); 3.12.16 (CI) | igual | — |
| zarr / numcodecs | — | 3.4.2 / 0.17.1 | D-043 |
| xarray | 2024.2.0 | 2026.9.0 | primeiro que lê com zarr 3 nas versões disponíveis |
| fsspec | 2024.6.1 | 2026.9.0 | acesso HTTP |
| dask | 2024.7.0 | 2026.8.0 | mínimo do xarray 2026.9 é 2025.2; 2026.11.0 (23/11) evitado por ser recente demais |
| packaging / typing_extensions | 24.1 / 4.12.2 | 26.3 / 4.16.0 | mínimos do xarray e do zarr |
| aiohttp e dependências | só no extra | 3.14.5; aiohappyeyeballs 2.7.1, aiosignal 1.4.0, frozenlist 1.8.0, multidict 7.2.1, propcache 0.5.4, yarl 1.25.1 | CI roda os testes Zarr |
| donfig, google-crc32c, msgspec | — | 0.8.1.post1, 1.9.0, 0.22.0 | dependências do zarr 3 |
| numpy, pandas, scipy, netCDF4, h5py, h5netcdf, MetPy, Cartopy, matplotlib | 2.0.0, 2.2.2, 1.14.0, 1.7.1.post1, 3.11.0, 1.3.0, 1.6.2, 0.23.0, 3.9.1 | **iguais** | identidade bit a bit |

Ambiente `lorenz-4d` na Swell: clone de `lorenz` + `pip install -r requirements.txt` da branch; `pip check` sem
problemas; diferença de `pip freeze` = exatamente as linhas acima. `lorenz` não foi alterado.

`pyproject.toml`: extras `[arco]` e `[destine]` = `zarr>=3.4.2,<4`, `numcodecs>=0.17.1`, `fsspec>=2026.9.0`,
`aiohttp>=3.9`; limites inferiores de xarray (2026.9.0) e dask (2026.8.0) = versões testadas. Os workflows de CI não
mudaram: instalam `requirements.txt`, que agora inclui o zarr.

## 3. Sondagem do store v3 (Swell, 25/11/2026, token da 2b via `DESTINE_EDH_TOKEN`)

6 requisições, ≈ 94 KB, todas HTTP 200: `zarr.json` raiz do v3 (49 KB), `.zmetadata` do v2 (42 KB), coordenadas
`isobaricInhPa`, `latitude`, `longitude` do v3 e um HEAD de `valid_time`.

| Item | v3 `era5-pressure-levels-v0.zarr` | v2 `reanalysis-era5-pressure-levels-v0.zarr` |
|---|---|---|
| Formato / metadados | Zarr 3, consolidados no `zarr.json` raiz, sem sharding | Zarr 2, `.zmetadata` |
| Dimensões | `valid_time`, `isobaricInhPa` (int32, hPa, 19 níveis 1000→1), `latitude` (90→−90), `longitude` (0→359,75) | mesmas |
| Bloco dos campos | 1440 × 1 × 60 × 60 (float32) | 120 × 1 × 256 × 256 |
| Codecs | `numcodecs.bitround` keepbits **13**, bytes, blosc zstd 5 shuffle | bitround keepbits **10**, blosc zstd 5 shuffle |
| Tempo | 761 544 horas (1940-01-01 a 2026-11-15 23 UTC), um bloco (≈ 0,8 MB comprimido) | congelado em meados de 2026 |

**Leitura na Parte 2 e correção na Parte 3:** os metadados do v2 confirmam a causa dos degraus da DestinE na 2b
("não verificada" na 4c): arredondamento para 10 bits de mantissa (T: meio degrau 0,0625 K abaixo de 256 K e 0,125 K
acima; z: 64 m² s⁻² até 262 144). Na Parte 2 li o `keepbits 13` do v3 como precisão 8 × maior e escrevi isso no
aviso, na ajuda e nas docs. **O teste real (§7) desmentiu essa leitura:** os valores do v3 são idênticos aos do v2 e
os 13 bits menos significativos da mantissa são zero em 100 % dos valores, ou seja, os dados foram copiados do v2 já
arredondados para 10 bits e o codec de 13 bits não os altera. Os textos foram corrigidos depois do commit `f923d80`
(§7.3). Dados acrescentados ao v3 depois de 07/2026 não foram verificados.

## 4. Implementação

- `src/utils/arco_zarr.py`: `load_packages` exige zarr ≥ 3; `_is_metadata` reconhece `zarr.json`; novo `InFlight`
  (no máximo 4 leituras ao mesmo tempo, pico registrado); `CountingStore` inalterado na contagem, limites, tentativas
  e recusa de bloco ausente; novo adaptador `MappingStore` (store do zarr 3, só leitura, cada chave lida do mapeamento
  contado numa thread) e `open_store_dataset(store, zarr_format)` (formato explícito: nenhuma chave do outro formato é
  pedida); store com sharding recusado na estimativa; proveniência com `max_parallel_reads`.
- `src/utils/destine_edh.py`: `STORES` (v3 padrão, v2), `--destine-store`; aviso de precisão por store; mensagem de
  horário fora do store por store; proveniência com `store_version` e `zarr_format`; mesmas proteções (tetos de bytes
  e requisições antes de ler, contagem a cada tentativa, 429/401/403 sem nova tentativa, cada bloco lido uma vez, sem
  dask).
- `lorenzcycletoolkit.py`: `--cdsapi` sem "(experimental)"; `--time-resolution` vale para todas as rotas, com as
  restrições de cada uma; `--destine-store {v3,v2}` (só com `--destine`); ajuda de `--destine` atualizada.
- Rotas kerchunk, Get a Subset, NCSS, WEkEO e CDS: nenhum código mudou.

## 5. Testes offline

- Stores sintéticos gravados com zarr 3 (formato 2 para ARCO e DestinE v2, formato 3 para DestinE v3) e convertidos
  em dicionários de bytes; `lengthen_time_axis` portado.
- Novos: conteúdo e contrato do arquivo com v3 e v2; leitura única de cada bloco com store curto e longo, v3 e v2
  (teste do grafo); store com sharding recusado antes de ler blocos; opção v2 e store desconhecido; no máximo 4
  leituras simultâneas (ARCO e DestinE); chaves de metadados dos dois formatos; recusa de zarr 2; ajuda da linha de
  comando (`tests/test_cli_help.py`).
- Swell, `lorenz-4d`, árvore de trabalho: testes Zarr e de ajuda 100 passaram. Suíte completa: ver §5.1.

### 5.1 Suíte completa

Swell, 25/11/2026, `CARTOPY_DATA_DIR` = cache Natural Earth da Swell (sem ele, os testes de mapas dão erro de
ambiente, como nas suítes da 4b·3):

| Código / ambiente | Resultado |
|---|---|
| `develop` `91dfbed` / `lorenz` (antes) | 745 passaram, 34 ignorados (testes Zarr: sem zarr), 0 falhas |
| branch, árvore de trabalho / `lorenz-4d` | 792 passaram, 1 falha, 0 ignorados |
| idem, depois da correção abaixo | `test_packaging.py`: 4 passaram |

A falha era o teste de empacotamento: o adaptador passou a ter um `from zarr.abc.store import Store` explícito, e o
teste exige que todo import de terceiros esteja declarado. O zarr foi incluído em `OPTIONAL_IMPORTS` (é extra
`[arco]`/`[destine]`, como `hda` para `[wekeo]`). Aviso novo só em testes: o xarray 2026.9 avisa sobre dimensão
ilimitada no `to_netcdf` de um auxiliar de `test_input_cleanup.py` (5 testes); nenhum efeito no código do toolkit.
flake8 com os erros que bloqueiam o CI (E9, F63, F7, F82): 0.

## 6. Regressão sem download (Parte 3.1)

**Código:** `git archive` do commit `f923d80` da branch (D-040), em `lectk_stage4d_zarr3_20261125/code_f923d80`;
ambiente `lorenz-4d`. **Execução:** `regress_4d.py` (na pasta da etapa), 8 rodadas em paralelo; cada uma numa pasta
nova `regress_f923d80/<id>/` com cópia de `inputs/` e da entrada (sha256 da cópia = da origem), os mesmos argumentos
da referência, `OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=2`, `PYTHONHASHSEED=0` (como o harness da fase 2), credenciais
removidas do ambiente e lançador `tests/run_cli_offline.py` (rede bloqueada). Nenhuma pasta de referência foi
alterada; as cópias das entradas foram apagadas depois de cada rodada (sobram 105 MB de saídas).

**Resultado: 78 de 78 CSVs idênticos bit a bit às referências; código de saída 0 em todas; nenhuma tentativa de
rede; entradas inalteradas; 68 de 68 entradas da fase 2 iguais ao sha256 selado.**

| Referência (caso 19790205) | Entrada | CSV referência = novo |
|---|---|---|
| 4a.1, comando principal `--ncar` | `6881a228` | `ffad993f` |
| 4a, grade comum kerchunk | `567ccdce` | `ac5c7eec` |
| 4b·1, comando principal `--wekeo` | `cf2ea3a2` | `8c919a68` |
| 4b·1, grade comum WEkEO nova | `5e497c48` | `0fe175d2` |
| 4b·1, grade comum CDS 2b (34 h) | `7d124089` | `d46d82ba` |
| 4b·2, comando principal `--arco` (execução 1) | `99aebc35` | `810d45e4` |
| 4b·2, comando principal `--arco` (execução 2) | `7c0fe5b8` | `810d45e4` |
| 4b·2, grade comum ARCO nova | `6fc9c51c` | `810d45e4` |
| 4b·3, comando principal `--destine` | `ee5e8463` | `7ba757d2` |
| 4b·3, grade comum DestinE nova | `dd7d837f` | `ced56f3e` |

Fase 2 da 2b: 68 rodadas (10 `ref`, 10 `c1`, 1 piloto, 47 por rota: ARCO, DestinE, kerchunk e WEkEO 10 cada, Get a
Subset 7), todas idênticas ao `outputs_sha256` dos manifestos. Tabela completa: `regress_f923d80/regression.tsv`.
Tempo somado 1 816 s.

**Observação:** o xarray 2026.9 emite `FutureWarning` em `src/utils/preprocessing.py:428` e `:430` (`drop` →
`drop_vars`). Sem efeito nos resultados; fica como pendência para uma etapa futura (não é ajuste necessário agora).

## 7. Testes reais (Parte 3.2 e 3.3)

Limites aprovados pelo responsável em 25/11/2026. Código: export de `f923d80` (sem `.git`, por isso
`lectk_commit` é nulo na proveniência; o commit está em `code_f923d80/COMMIT` e nos scripts `run_destine_v3.sh` e
`run_arco.sh`), ambiente `lorenz-4d`, mesmo track e caixa da 4b·2/4b·3 (`inputs/` copiados), pastas novas
`real_destine_v3/` e `real_arco/`. As duas corridas foram seguidas, sem sobreposição.

### 7.1 DestinE, store v3 (padrão)

**Estimativa antes de ler blocos** (`--max-gib 0.01`): 5 700 blocos, ≈ 5 707 requisições, teto 111,82 GiB
(120 067 829 703 B), 1 754,9 MiB de valores úteis; recusada sem ler nenhum bloco. Uma primeira tentativa com
`--max-gib 0.001` parou já nos metadados (a coordenada `valid_time` do v3 tem ≈ 0,8 MB), também sem ler bloco.

**Comando** (`run_destine_v3.sh`): `lorenzcycletoolkit.py case_19790205.nc -f -r --destine --trackfile inputs/track
--time-resolution 6 --max-requests 6000 --max-gib 112 -o destine_v3_19790205_main`.

| Item | v3 (4d) | v2 (4b·3) |
|---|---|---|
| Requisições | 5 707 (7 de metadados + 5 700 blocos), todas HTTP 200, 0 exceções, 0 novas tentativas | 1 528 |
| Recebido | 26 485 904 273 B = 24,67 GiB (2,48 MB de metadados) | 8,79 GiB |
| Teto conferido antes | 111,82 GiB | 45,23 GiB |
| Leituras em paralelo (pico) | 4 | — |
| Blocos lidos mais de uma vez | 0 | 0 |
| Transferência / total | 3 050 s / 51,6 min (18:15:25–19:07:00 UTC) | 1 010 s / 17,8 min |
| Arquivo | 1 840 155 497 B, sha256 `4310c8e4…` | `ee5e8463…` |
| LEC | rodou até o fim; CSV `7ba757d2…` | `7ba757d2…` |

**Dados v3 × v2 (arquivo inteiro: 34 horários × 19 níveis × 225 × 633, 92 006 550 valores por campo):**
coordenadas iguais e **valores idênticos bit a bit em t, z, u, v e w**. O sha256 do arquivo difere só porque o atributo
`history` traz a hora e `source` diz "v3". Em 100 % dos valores os 13 bits menos significativos da mantissa são zero:
os valores têm 10 bits, embora os metadados do v3 declarem `keepbits 13` (dados copiados do v2 já arredondados).

**Grade comum da 2b** (`to_common_grid` do harness da rodada 2, 34 horários × 13 níveis): arquivo idêntico ao da 4b·3
(sha256 `dd7d837f…`); LEC com o código `f923d80` (rede bloqueada): CSV `ced56f3e…`, **idêntico** ao da 4b·3. Contra o
CDS da 2b: diferença máxima 0,125 K em T e 64 m² s⁻² em z, a mesma da 2b, ou seja, o mesmo meio-degrau de 10 bits.
Comparação com o CDS de 34 horários portanto igual à da 4b·3 (só BΦE acima de 1 %, 1,21 %).

**Conclusão:** "os dados não mudam, os blocos sim" confirmado para este caso, inclusive no LEC. O v3 custa 3,7 ×
mais requisições e 2,8 × mais bytes neste caso curto e é o único atualizado. Dados acrescentados ao v3 depois de
07/2026 não foram verificados (podem ter de fato 13 bits).

### 7.2 ARCO

**Comando** (`run_arco.sh`, limites da 4b·2): `lorenzcycletoolkit.py case_19790205.nc -f -r --arco --trackfile
inputs/track --time-resolution 6 --max-gib 20 -o arco_19790205`.

| Item | 4d (zarr 3) | 4b·2 execução 2 (zarr 2) |
|---|---|---|
| Teto conferido antes | 13 053 214 720 B (12,16 GiB) | igual |
| Recebido | 4 944 004 372 B (metadados 59 535 B) | 4 944 004 336 B |
| Leituras | 4 433, 0 novas tentativas, pico de 4 em paralelo | 4 431 |
| Transferência / total | 978 s / 16,9 min | 793 s / 13,8 min |
| Valores do arquivo | **idênticos** (u, v, t, w, z e coordenadas; só `history` difere) | — |
| CSV do LEC | `810d45e4…` | `810d45e4…` (**idêntico**) |

As 2 leituras e 36 bytes a mais são metadados lidos pelo zarr 3 ao abrir o store v2; os blocos são os mesmos.

### 7.3 Correção depois do commit `f923d80`

O aviso de precisão do v3, a ajuda de `--destine`/`--destine-store`, as docs (§2.7, §3.1, limitações), o CHANGELOG, a
D-043 e este registro diziam que o v3 guarda 13 bits e é 8 × mais preciso. Corrigidos para o que foi medido (10 bits nos
dois stores, mesmos valores). O aviso de precisão passou a ser o mesmo para v3 e v2 (`PRECISION_WARNING`, o texto da
4b·3). Testes afetados (DestinE, ajuda, ARCO, empacotamento): 104 passaram na Swell.

## 8. Pendências para a retenção

Nada foi apagado na Swell. Pasta da etapa `/p1-swell/danilocs/lectk_stage4d_zarr3_20261125` (4,6 GB):

- `real_destine_v3/` (arquivo v3 1,84 GB, proveniência, LEC), `real_arco/` (arquivo 1,26 GB, proveniência, LEC),
  `cmp/` (grade comum v3, LEC, `cmp_v3.out`), `regress_f923d80/` (78 rodadas, 105 MB), `code_f923d80/`, `dev/`,
  `dev2/`, `base/` (cópias de trabalho para as suítes), `env/` (freeze e sondagem), `logs/`, `logs_real/`, scripts;
- ambiente conda `lorenz-4d` (clone de `lorenz` + pins da 4d);
- as pastas das etapas 4a–4b·3 e da 2b só foram lidas.

Para a retenção (D-028): os arquivos v3 e v2 do caso são idênticos em valores, então o v3 não é evidência nova além da
proveniência e do LEC; `dev/`, `dev2/`, `base/` e os `.tgz` são descartáveis; decidir o destino de `lorenz` × `lorenz-4d`
depois do merge.

## 9. Nota de fechamento (25/11/2026): aprovação, PR, CI, merge e releitura

- **Aprovação:**
  - o responsável inspecionou os diffs no VS Code e autorizou o commit `f923d80`, que é o código dos testes da Parte 3;
  - depois autorizou o commit `dd1bdf0` (correção dos 13 bits, §7.3) e o PR.
- **Commits e PR:**
  - `f923d80` e `dd1bdf0` na branch única da etapa (`feat/stage-4d-zarr3`, D-042);
  - PR único para a `develop`: [#88](https://github.com/daniloceano/LorenzCycleToolkit/pull/88).
- **CI do PR:** verde. Os dois `build` do workflow "Python application" passaram (15 min 13 s e 11 min 22 s), e o
  `ci/circleci: build_test` também. O workflow instala o `requirements.txt` com zarr 3, então rodou a suíte inteira com
  os testes Zarr depois da correção.
- **Merge:** feito pelo responsável em 25/11/2026 às 20:54 UTC; o commit de merge é `436279b`
  (`436279b70b8f18adb2bb76b4fb201bd771647c4f`).
- **Releitura no SHA do merge:**
  - a árvore do merge é idêntica à do commit `dd1bdf0` (diff vazio).
  - Entre `f923d80` (o código dos testes reais) e o merge, o código só muda em texto:
    - a constante do aviso de precisão (`PRECISION_WARNING_V3` removida; o v3 usa o mesmo `PRECISION_WARNING` do v2);
    - a ajuda de `--destine` e `--destine-store`;
    - o docstring de `destine_edh.py`.

    Nenhuma leitura, seleção ou cálculo mudou, então os testes reais e a regressão valem para o código do merge.
  - O build do Sphinx (`-n -W`, mesmo ambiente isolado) a partir de um `git archive` de `436279b` terminou sem avisos.
  - A página `remote_data_acquisition` diz "10 bits of mantissa". "13 bits" só aparece onde se explica que é o valor
    declarado nos metadados do v3.
  - Em `docs/source` não há caminho da Swell nem credencial.
- **Pendências da etapa 4:** retenção dos dados na Swell (D-028; ver §8) e resumo de transferência. Fica registrado
  também o `FutureWarning` de `drop` → `drop_vars` em `src/utils/preprocessing.py` (linhas 428 e 430) visto na
  regressão (§6); ele não afeta o resultado.
