# Etapa 4d — ambiente com zarr 3 e DestinE no store Zarr v3 (25/11/2026)

**Natureza:** código, dependências, testes, documentação de usuário e registros. Cálculo do LEC inalterado
(`src/analysis`, `src/frameworks` sem mudança). `main`, tags e releases inalterados.

**Base:** `develop` em `91dfbed` (PR #87); `main` em `48fe8e2`. Branch única da etapa: `feat/stage-4d-zarr3` (D-042).

**Estado deste registro:** Partes 1 e 2 concluídas; **Parte 3 (Swell) pendente** — as seções 6 e 7 serão
preenchidas depois da regressão e dos testes reais, com o commit usado (D-040).

## Resumo em linguagem simples

O ambiente do LECTK passa a usar zarr 3 com um xarray novo, porque o store da DestinE que continua sendo atualizado
(Zarr v3) só é lido com zarr 3, e a operacionalização futura do LEC precisará desse ambiente de qualquer modo. Só a
camada de leitura de dados foi atualizada; o núcleo numérico (numpy, pandas, scipy, netCDF4) ficou igual, para que o
LEC dê resultados idênticos bit a bit sobre os mesmos arquivos. `--destine` passa a ler o store v3 por padrão; o
store v2 da 2b continua disponível com `--destine-store v2` enquanto for servido. Os dois stores **não** guardam os
mesmos valores: o v2 arredonda a mantissa para 10 bits, o v3 para 13 (confirmado nos metadados).

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

**Consequência:** a causa dos degraus da DestinE na 2b ("não verificada" na 4c) está confirmada: arredondamento para
10 bits de mantissa no v2 (T: meio degrau 0,0625 K abaixo de 256 K e 0,125 K acima; z: 64 m² s⁻² até 262 144). No
v3, com 13 bits, os erros máximos caem 8 ×: T ≤ 0,0156 K, z ≤ 16 m² s⁻². Portanto o critério "dados iguais ao v2 na
região comum" não pode valer bit a bit; a diferença esperada v3 − v2 está limitada pela soma dos meios degraus
(§7).

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

*(pendente: commit, plano aprovado, tabela por sha256.)*

## 7. Testes reais (Parte 3.2 e 3.3)

*(pendente: estimativa aprovada, DestinE v3 no caso 19790205, comparação com v2 e com o LEC da 4b·3; ARCO curto com
os limites da 4b·2.)*

## 8. Pendências para a retenção

*(preencher no fim.)*
