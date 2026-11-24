# Etapa 4b·3 — cliente DestinE nativo (`--destine`), com a Etapa 0 (24/11/2026)

**Natureza:** registro de implementação. Acrescenta a rota `--destine` ao programa principal (D-041). Nenhum valor de
dado de outra rota muda, e o cálculo do LEC também não muda. `main`, tags e releases inalterados.

**Base:** `develop` em `eca3864` (merge do encerramento da 4b·2, PR #83, conferido antes da branch
`feat/stage-4b3-destine`); `main` em `48fe8e2`.

**Estado:** implementação, testes offline e teste real concluídos na Swell (§4, §5). Os dois critérios foram atendidos:
dados idênticos à DestinE da 2b e, contra o CDS, só BΦE acima de 1 % (1,21 %), como na 2b para este caso.

## Resumo em linguagem simples

O usuário passa a poder baixar o ERA5 da DestinE (Earth Data Hub) direto pelo LECTK, com `--destine`. O arquivo entregue
ao LEC é igual ao do `--cdsapi` em nomes, unidades e ordem; o pedido usa a mesma caixa e os mesmos horários. O programa
avisa, no início e na proveniência, as duas limitações do produto: **19 níveis próprios** (não os 37) e **valores
guardados com precisão reduzida** ("degraus"), que na 2b mudaram 18 de 240 termos do ciclo entre 1 % e 2,7 % em relação
ao CDS. Cada bloco lido é **uma requisição HTTP** contra a **cota mensal de 500 000** da conta: antes de ler, o programa
mostra quantas requisições vai fazer e recusa se passar de `--max-requests` (padrão 10 000). Pela D-015, o LECTK informa e
não proíbe.

## 1. Etapa 0

### 1.1 Harness da 2b na Swell (só leitura, avisado ao responsável em 24/11, D-034)

**Encontrado** em `/p1-swell/danilocs/lectk_bench2b_20261109/harness/bench2b.py` (o registro da 4b·2 disse "não
encontrado", porque procurou em outras pastas). Credenciais não lidas nem copiadas.

| Item | Harness da 2b (`route_zarr`, `_zarr_open`, `_zarr_select`) |
|---|---|
| Store | `https://data.earthdatahub.destine.eu/era5/reanalysis-era5-pressure-levels-v0.zarr` (Zarr v2, metadados consolidados) |
| Autenticação | HTTP *basic*, usuário `edh` e o token como senha (`aiohttp.BasicAuth`); token em `~/.config/lectk_bench/destine_token` |
| Bibliotecas | zarr 2.18.4, numcodecs 0.15.1, fsspec 2024.6.1, aiohttp 3.14.5, xarray 2024.2.0, dask 2024.7.0 |
| Paralelismo | 4 (`gather_batch_size` e `num_workers`) |
| Recorte | `sel` por horários, níveis e fatias de latitude/longitude; longitude 0..360 dividida em dois arcos quando cruza Greenwich |
| Contagem | *trace* do aiohttp: requisições (`on_request_end`), corpos (`on_response_chunk_received`) e códigos de estado |
| Conversão | `_zarr_select` + `load()` com dask (`chunks={}` sobre o store inteiro) |
| Regra própria | teto de 500 GB e parada em 150 000 requisições no mês (cota 500 000) |

Manifesto do caso 19790205 (rep. 2, caixa de 5°, 13 níveis, 35 horários de 6 h): nível `isobaricInhPa`, tempo
`valid_time`, variáveis `u v t w z`; blocos de **120 h × 1 nível × 256 × 256** (31,5 MB sem compressão); **1 560 blocos,
1 569 requisições** (todas HTTP 200), **10 030 009 482 B recebidos** para um limite superior de 49 073 356 800 B (20 %),
861 s de aquisição. As três repetições receberam os mesmos bytes e requisições. Arquivo de grade comum
`runs/campaign/rep2/destine/19790205_case/common_destine_19790205_case.nc` (sha256 `c3672d3c…`).

### 1.2 Documentação atual da Earth Data Hub (consultada em 24/11/2026)

Fontes: `earthdatahub.destine.eu/docs/data-access.md`, `…/help/troubleshooting.md`, `…/whats-new.md`, catálogo STAC
público. Nada exigiu login; nenhuma conta foi criada.

| Item | Resultado |
|---|---|
| Token | "API key" pessoal, gratuita, na página `earthdatahub.destine.eu/quota-api-keys` (conta da plataforma DestinE) |
| Onde fica | `~/.netrc` recomendado: `machine data.earthdatahub.destine.eu` + `password <key>` (sem linha `login`), e o mesmo para `api.earthdatahub.destine.eu`; alternativa "para teste rápido": a chave na URL (`https://edh:<key>@…`). **Nenhuma variável de ambiente documentada** |
| Cota | **500 000 requisições por mês** por usuário, renovadas à meia-noite do dia 1º; cada bloco e cada arquivo de metadados é uma requisição; acima disso, **HTTP 429**; uso visível na página da conta. Sem limite documentado por minuto nem por bytes |
| **Store** | **Em 06/07/2026 o ERA5 passou para stores Zarr v3** (`era5-pressure-levels-v0.zarr`, blocos de 1440 h × 60 × 60, zarr-python ≥ 3); "os dados não mudam, os blocos sim"; **o store v2 da 2b "não é mais atualizado"**. Os dois respondem HTTP 401 sem token (24/11) |
| Níveis | 19: 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10, 5, 1 hPa |
| Frequência | 1 h; v3 até 2026-11-15 23 UTC (STAC de 22/11); atualização mensal |
| Grade | 0,25°; longitude 0..359,75; latitude 90 → −90 |
| Armazenamento | **a documentação não diz** o tipo, os codecs nem se há compressão com perda; os metadados por variável exigem token |
| Licença | CC-BY 4.0 (ERA5 via CDS) |

**Sobre os degraus:** a causa continua **não verificada**. A documentação não descreve quantização; o tamanho dos
degraus da 2b (0,0625/0,125 K em T, até 64 m² s⁻² em z, constante por faixa de níveis) é compatível com valores guardados
com menos bits de mantissa. A leitura dos `.zarray` com token, no teste real, mostra o tipo e os filtros.

### 1.3 Estimador de bytes e de requisições, conferido contra a 2b

O limite superior é o da ARCO (blocos inteiros tocados, sem compressão) **mais uma margem de codec** por bloco
(1/64 do bloco + 4 KiB, acima do que Blosc ou zstd acrescentam a um bloco incompressível) **mais os bytes já recebidos**
dos metadados e coordenadas. Requisições previstas: as já feitas mais uma por bloco tocado (sem novas tentativas).

Conferência com a 2b (caso 19790205): o estimador da 2b deu 1 560 blocos e 49,07 GB; foram recebidos 10,03 GB em 1 569
requisições (1 560 blocos + 9 de metadados). A regra acima dá o mesmo número de blocos e um limite ≥ 49,07 GB: **limite
superior real** (recebido = 20 % dele). Os "18–20×" e os "3,1 GB por caso" da tarefa são da fase 1 inteira (mediana
entre casos); para este caso o excesso foi 14× e o volume, 10 GB.

**Para o pedido do teste real** (track da 4a.1/4b·1/4b·2, caixa de 15° N−24 W−82 S−80 E76, 34 horários de 6 h, 19
níveis), contado supondo a grade da 2b (latitude 90 → −90, longitude 0..360, blocos de 120 h a partir de 1940-01-01): **16 blocos por campo e nível**, 1 520 blocos, **≈ 1 530
requisições**, limite superior **≈ 45,2 GiB** (44,5 GiB de blocos + margem). Pela razão da 2b (20 %), o recebido deve
ficar perto de 9 GiB. Valores úteis: 1,71 GiB.

### 1.4 Reaproveitamento

| Peça | Origem | Uso |
|---|---|---|
| Caixa, horários, teto padrão, commit, recusa de 180° | `ncar_gdex`, `track_alignment.planned_request_times` | reutilizados sem cópia |
| Seleção por posições, nomes de dimensões (inclui `isobaricInhPa` e `valid_time`), duas convenções de longitude, estimativa por blocos inteiros, sha256, JSON | `arco_zarr` (`select`, `estimate_bytes`, `_dimension`, `_sha256`, `_write_json`) | importados; erros da ARCO convertidos em erros da DestinE com a mesma mensagem |
| Abrir **sem dask** e recortar antes de ler | correção da 4b·2 | mantido, e levado adiante: **não há dask nenhum**. Cada bloco da DestinE guarda 120 horas; a escrita da ARCO (uma tarefa por horário) leria o mesmo bloco até 34 vezes. Aqui cada leitura é uma seleção ortogonal por campo e nível sobre todos os horários, e o zarr lê cada bloco tocado uma só vez |
| Store contado, tentativas, bloco ausente, 401/403 | `arco_zarr.CountingStore` | reescrito na DestinE para contar também requisições, parar no limite de requisições e tratar 429 sem nova tentativa |
| Escrita do NetCDF na thread principal, checagem final, proveniência | `arco_zarr` | mesmo desenho |

## 2. Decisões

**Consulta ao responsável (24/11/2026), porque a Etapa 0 contrariou o padrão:** o store v2 da 2b está congelado desde
06/07/2026; o v3 atual exige zarr ≥ 3, que não funciona com o xarray 2024.2.0 do projeto e conflita com o `[arco]`. Para o
caso de teste, o v3 custaria ≈ 5 700 requisições e 110 GiB de limite superior, contra ≈ 1 530 e 45 GiB no v2. **Respostas
do responsável:** usar o **store v2 da 2b**; padrão de **10 000 requisições** por execução (2 % da cota mensal).

**D-041 (aprovada pelo responsável em 24/11/2026, com o diff; commit `cc9729f`):** rota DestinE (etapa 4b·3): `--destine`, mutuamente exclusivo com
`--cdsapi`, `--wekeo`, `--ncar` e `--arco`; store **Zarr v2 da 2b** (`reanalysis-era5-pressure-levels-v0.zarr`), congelado
desde 06/07/2026, com recusa antes de qualquer leitura de horários fora dele; níveis do store (19), sem interpolar;
qualquer `--time-resolution` em horas inteiras; mesma caixa, mesmos horários e mesmo contrato de arquivo do CDS
(`inputs/namelist_ERA5-cdsapi`); avisos (a) dos níveis e (b) da precisão reduzida no início, no atributo do arquivo e na
proveniência; token da Earth Data Hub em `~/.netrc` (como a documentação recomenda) ou `DESTINE_EDH_TOKEN` (variável do
LECTK, não da EDH), autenticação *basic* com usuário `edh`, nunca registrado; extra `[destine]` (os pacotes do `[arco]`);
**`--max-requests`** (padrão 10 000), conferido antes de ler blocos sobre as requisições previstas e durante a leitura
sobre cada tentativa; `--max-gib` (padrão 20 GiB) sobre os blocos inteiros + margem de codec e depois sobre os bytes
recebidos; 3 tentativas por bloco, 4 em paralelo; 401, 403 e 429 sem nova tentativa; requisições HTTP reais contadas pelo
*trace* do aiohttp e registradas.

Escolhas de implementação, para a revisão:

1. **Sem dask.** Ver §1.4. Um teste garante que cada bloco é lido uma só vez e que o número de blocos lidos não muda
   com o tamanho do store (96 horas ou 400 dias).
2. **Os 19 níveis são baixados**, embora o LEC descarte 5 e 1 hPa (corte abaixo de 10 hPa, como para todas as fontes):
   cerca de 10 % de bytes e requisições a mais, pela regra "os níveis do produto".
3. **Contagem dupla das requisições:** `requests_store` (cada tentativa de cada chave, a que o limite usa) e
   `requests_http` (o que o aiohttp realmente enviou). Diferenças entre os dois apareceriam no teste real.
4. **Margem de codec** no limite superior, para que ele continue limite mesmo para um bloco incompressível.
5. **Não há registro local do uso mensal** da cota: o LECTK não sabe o que a conta já gastou em outras execuções ou
   ferramentas; o uso está na página da conta.
6. Como na 4b·2, `zarr` e `aiohttp` **não entram no `requirements.txt`**: na CI, os testes que escrevem ou leem Zarr
   ficam ignorados.

## 3. Arquivos

Novos: `src/utils/destine_edh.py`, `tests/test_destine_edh.py`, este registro. Alterados: `lorenzcycletoolkit.py`
(`--destine`, `--max-requests`, exclusões, `--max-gib` e `--delete-downloaded-input`), `src/utils/preprocessing.py` (ramo
`--destine`, `namelist_for_run`), `src/utils/input_cleanup.py`, `src/utils/run_timing.py` (comentários), `pyproject.toml`
(extra `[destine]`), `requirements.txt` (comentário), `tests/test_run_namelist.py` (`--destine` nos casos),
`docs/source/remote_data_acquisition.rst` (nota no topo e aviso na §2.6), `docs/source/configuration.rst`,
`CHANGELOG.md`, `docs/evolution/queue.md`, `docs/evolution/decisions.md`. As mensagens de `--max-gib` e
`--delete-downloaded-input` ganharam "--destine" mantendo os trechos que os testes anteriores procuram.

## 4. Testes

Na **Swell** (D-034), pasta `/p1-swell/danilocs/lectk_stage4b3_destine_20261124/`: `base/` (export de `eca3864`) e `work/`
(árvore de trabalho da branch, copiada por `tar`/`scp`); ambiente `lorenz` (Python 3.12.12, numpy 2.0.0, xarray 2024.2.0,
sem zarr, como a CI) e um `venv` sobre ele com zarr 2.18.4, numcodecs 0.15.1, aiohttp 3.14.5 e fsspec 2024.6.1. Suíte
completa, sem rede:

| | passou | ignorado | falhou | erro |
|---|---|---|---|---|
| Antes (`eca3864`), sem zarr | 699 | 16 | 0 | 16 |
| Antes (`eca3864`), com zarr | 715 | 0 | 0 | 16 |
| Depois, sem zarr (como a CI) | **729** | **34** | 0 | 16 |
| Depois, com zarr | **763** (+46 da DestinE, +2 de `--destine` no namelist) | 0 | 0 | 16 |

Os 16 erros são o cache Natural Earth ausente (`CARTOPY_DATA_DIR`), **a mesma lista antes e depois** (md5 da lista
`6db6a0d6…` nos dois casos). Sem zarr, 18 dos 46 testes novos ficam ignorados e 28 rodam. Nenhum teste de outra rota
mudou de resultado.

**Dois incidentes, registrados:**

1. A primeira rodada foi tentada **neste Mac** (contra a D-034, sem consulta) e o travou: o teste "o grafo não cresce"
   montava um store sintético de 400 dias com a grade inteira em memória (~99 GB). O teste foi reescrito: o store longo é
   alongado **só nos metadados** (`resize` do zarr), e os stores dos demais testes foram reduzidos a 13 horas e à caixa;
   o arquivo inteiro roda com no máximo 0,7 GB. Depois disso, todo o trabalho de teste passou para a Swell, por
   orientação do responsável.
2. Na primeira rodada da Swell, `test_packaging::test_every_imported_distribution_is_declared` falhou por causa de
   32 arquivos `._*.py` (atributos estendidos do macOS) criados pelo `tar` da cópia, não pelo código; apagados das cópias,
   as duas suítes "depois" foram refeitas (números acima).

Cobertura de `tests/test_destine_edh.py`, com store Zarr v2 sintético no formato do real (`valid_time`, `isobaricInhPa`,
19 níveis, 1 h, latitude decrescente, longitude 0..360, blocos de 120 h × 1 nível × uma placa): plano igual ao do CDS;
Greenwich com o store em 0..360; recusa de 180°; horas inteiras (1, 3 e 6 h servidas; 0, −3 e 1,5 recusadas); contrato do
arquivo e valores idênticos aos do store nas duas convenções de longitude; níveis em Pa e em outro número; namelist
`inputs/namelist_ERA5-cdsapi`; avisos (a) e (b) no log, no atributo do arquivo e na proveniência; **cada bloco lido uma
só vez e o mesmo número de blocos para um store de 13 horas e de 400 dias**; nenhum array dask; estimativa ≥ recebido e
requisições previstas = feitas; recusa acima de `--max-requests` e de `--max-gib` sem ler bloco nenhum, e o valor dito na
recusa basta; disco; contagem real que para no limite de bytes e no de requisições; tentativas limitadas e contadas
como requisições; 401, 403 e 429 sem nova tentativa; bloco ausente; horário fora do store congelado; token pelo `.netrc`
(formato da EDH, sem `login`), por `NETRC` e por variável; token ausente antes de qualquer conexão; pacote ausente e
zarr 3; autenticação *basic* e *trace* só no *mapper*; contador HTTP; arquivo existente reaproveitado; CLI (exclusões,
`--max-requests`, `--delete-downloaded-input`, opções aceitas, texto do `--help`); o token não aparece em log nem na
proveniência.

## 5. Teste real na Swell (24/11/2026)

**Autorizações do responsável (24/11):** diff aprovado e commit `cc9729f` antes do teste; credencial: **o token da 2b**
(`~/.config/lectk_bench/destine_token`), passado ao processo por `DESTINE_EDH_TOKEN` dentro do comando, sem ser lido,
copiado ou impresso; limites **`--max-requests 2000`** e **`--max-gib 46`**, aprovados sobre a estimativa da §1.3.

**Ambiente (pasta nova, fora do repositório):** `/p1-swell/danilocs/lectk_stage4b3_destine_20261124/` (2,8 GB). **Código:
clone git da branch no commit `cc9729f2ae76aa4f0da2976906c765d7d844675a`**, sem alterações locais (`git status` vazio);
a proveniência grava esse SHA em `lectk_commit` (D-040 atendida sem comparação de hashes). `venv-zarr` sobre o ambiente
`lorenz` (Python 3.12.12, numpy 2.0.0, xarray 2024.2.0) mais zarr 2.18.4, numcodecs 0.15.1, aiohttp 3.14.5 e fsspec
2024.6.1. Entradas copiadas de `lectk_stage4b2_arco_20261124/run2/inputs` (track `9da9d809…`, o selado da 2b; caixa e
namelist, sha256 iguais). Nenhuma pasta existente foi alterada (2b, fase 2, 4b·1 e 4b·2 só lidas).

**Comando** (`run_main.sh 46 2000`): `lorenzcycletoolkit.py case_19790205.nc -f -r --destine --trackfile inputs/track
--time-resolution 6 --max-gib 46 --max-requests 2000 -o destine_19790205`. `--time-resolution 6` como a 2b; o track é de
3 h e foi alinhado aos 34 horários de 6 h (mensagem do programa).

### 5.1 Execução

| Item | Resultado |
|---|---|
| Início / fim (UTC) | 20:32:22 / 20:50:11; código de saída 0; LEC `success` |
| Tempo total | **1 066 s (17 min 46 s)**: aquisição 1 010,4 s, preparação 0,09 s, LEC 55,8 s |
| Previsão antes de ler (log) | **1 520 blocos, 1 528 requisições, até 45,23 GiB** (48 568 570 372 B); igual à §1.3 |
| Requisições | **1 528** contadas pelo store e **1 528** enviadas pelo aiohttp, todas **HTTP 200**; 8 de metadados + 1 520 blocos; **0 novas tentativas**; nenhum bloco lido duas vezes |
| Bytes recebidos | **9 435 028 228 B (8,79 GiB)**, dos quais 168 452 B de metadados; **19,4 % do limite superior** (2b: 20 %); 5,4 × os valores úteis (1 754,9 MiB) |
| Níveis | 19 confirmados no store (1000 → 1 hPa), avisos (a) e (b) no log, no atributo do arquivo e na proveniência |
| Arquivo | 1 840 155 398 B; sha256 `ee5e846388d892573e34bfaf9faa92b8681a706268f7c1529c46e5b1ca258add`; 34 horários (1979-03-18 12 UTC a 03-26 18 UTC) × 19 níveis × 225 × 633; N−24 W−82 S−80 E76 |
| CSV do LEC (34 × 26) | sha256 `7ba757d2a4def7b830e492205a6f88ade932959e40082d09da2df5138da8b18c` |
| Token | não aparece em nenhum arquivo de `run/` e `logs_real/` (busca feita na Swell sem imprimir o valor) |

O uso desta execução na cota de novembro: 1 528 das 500 000 requisições (0,3 %).

### 5.2 Critério 1: dados na região comum

Grade comum da 2b (34 horários × 13 níveis × 142 × 550, 34 520 200 valores por campo), com `to_common_grid` do harness
da rodada 2 (verificações `ok`); `cmp_data.py`, só leitura das pastas da 2b:

| Comparado com | `t` | `z` | `u`, `v`, `w` |
|---|---|---|---|
| **DestinE da 2b, rep. 1, 2 e 3** | **idêntico** | **idêntico** | **idênticos** |
| CDS da 2b | máx. 0,125 K, 99,99 % dos pontos | máx. 64 m² s⁻², 99,9 % | máx. 0,031 / 0,031 m s⁻¹ e 0,006 Pa s⁻¹ |
| WEkEO e ARCO da 2b | iguais ao CDS | iguais ao CDS | iguais ao CDS |

**Critério 1: atendido.** Mesmo produto que a DestinE da 2b, bit a bit; contra o CDS, a diferença é a da fase 2 (máx.
0,125 K e 64 m² s⁻², classe "relevante" nos campos de entrada; u, v, ω com diferença pequena). Os degraus são do produto,
não do cliente.

### 5.3 Critério 2: LEC

O LEC rodou até o fim no comando principal. Para comparar com o CDS, a grade comum do arquivo novo foi processada com o
**mesmo export de código da 4b·1 e da 4b·2** (`068c3a9`, `launch_lectk.py`, `-f -r`, caixa e namelist da 4b·2, sha256
iguais) e comparada com o LEC do CDS da 2b recortado nos mesmos 34 horários (`d46d82ba…`, o da 4b·1, só lido),
**com a classificação da fase 2** (`classify_term` do harness: e_max = máx.|Δ| / máx(média|ref|, F), "relevante" acima
de 1 %); `cmp_lec.py`:

| | Resultado |
|---|---|
| CSV da grade comum | `ced56f3e…` (difere do CDS `d46d82ba…`) |
| Termos (24 classificados) | **23 pequena, 1 relevante** (+2 auxiliares sem classe) |
| Acima de 1 % | **só BΦE: 1,21 %** |
| Seguintes | Ge 0,91 %, RGe 0,75 %, BΦZ 0,71 %, RGz 0,56 %, Gz 0,55 %, BAe 0,35 % |
| Reservatórios e conversões | todos "pequena": Az 0,010 %, Ae 0,026 %, Kz 0,005 %, Ke 0,001 %; Cz 0,020 %, Ca 0,029 %, Ck 0,027 %, Ce 0,030 % |
| Na 2b, mesmo caso (35 h) | **só BΦE relevante, 1,24 %** |

**Critério 2: atendido.** Os mesmos termos passam de 1 % que na 2b para este caso: **só BΦE**, com o mesmo tamanho
(1,21 % × 1,24 %). Os 18 de 240 da 2b somam os 10 casos; neste, era um. Reservatórios e conversões não mudam (≤ 0,03 %).
**Não se declara equivalência com o CDS.**

**Observação [PRELIMINAR], fora dos critérios:** o CSV do comando principal (19 níveis, até 10 hPa no LEC) difere do da
grade comum (13 níveis, até 50 hPa), ambos DestinE, na mesma caixa fixa: médias de Kz +2,1 %, Ke +0,9 %, Cz +1,0 % e
**Ck +12,5 %**; Ca e Ce < 0,2 %. É o efeito dos níveis 70, 30, 20 e 10 hPa neste caso, **não** o efeito de 19 contra 37
níveis, que não foi medido aqui (a LZV mediu outro caso por subamostragem).

### 5.4 Onde estão os resultados

Pasta `/p1-swell/danilocs/lectk_stage4b3_destine_20261124/`: `code/` (clone em `cc9729f`), `run/` e `logs_real/` (teste
real), `cmp/` (grade comum, LEC da comparação, `lec_by_term_destine_new_vs_cds34.csv`, saídas), `cmp_data.py`,
`cmp_lec.py`, `run_main.sh`; `base/`, `base_z/`, `work/`, `work_z/` e `logs/` (suítes da §4). O destino da pasta fica
para o fim da etapa 4 (D-028).

## 6. Limitações (leia antes de usar `--destine`)

- **Store congelado:** o v2 não é atualizado desde 06/07/2026 e pode ser retirado sem aviso; datas posteriores são
  recusadas. O v3 atual não é suportado (exigiria zarr 3 e outro xarray).
- **19 níveis** e **precisão reduzida**: ver os avisos; o efeito só foi medido na 2b (13 níveis comuns, 10 casos).
- **Cota:** 500 000 requisições por mês na conta; o LECTK limita cada execução, não o mês.
- **Volume:** blocos de 120 h × 256 × 256; no teste real, 8,79 GiB recebidos (5,4 × o útil, 19 % do limite superior);
  17 min 46 s no total (uma medição).
- **Um só caso** (19790205) e uma execução real.
- **Causa dos degraus:** não verificada (§1.2).

## 7. Acréscimos de 24/11/2026 (encerramento documental; o texto acima não foi reescrito)

### 7.1 Integração e CI

PR #84 integrado em `develop` pelo responsável em 24/11/2026 (23:42 UTC); merge
`83a10be7d182b758da702a9356dab04cd1b8378a`, com pais `eca3864` (develop) e `07bc624` (cabeça da branch). Os três checks
passaram (`build` do GitHub Actions, duas vezes, e `ci/circleci: build_test`). `main` ficou em `48fe8e2`, sem tag nem
release; nenhuma tag contém `cc9729f`.

### 7.2 Identidade do código testado na Swell

**O que foi comparado.** O teste real (§5) rodou num clone git da branch no commit `cc9729f`, sem alterações locais
(§5, `git status` vazio), e `lectk_commit` da proveniência registra esse SHA. Na releitura (24/11/2026), comparou-se esse
commit com o merge pelo Git:

- as árvores de `07bc624` e de `83a10be` são idênticas;
- entre `cc9729f` e `83a10be` mudam só três documentos (`docs/evolution/decisions.md`, `docs/evolution/queue.md` e este
  registro), todos do commit `07bc624`; nenhum código, teste, dependência ou página do Read the Docs muda;
- o sha256 de `src/utils/destine_edh.py` no merge é `601f57eb…`, o mesmo do commit testado.

**Limite da afirmação.** A identidade vale para os arquivos rastreados; o ambiente da Swell (`lorenz` mais `venv-zarr`,
versões na §5) não é parte da comparação. Não houve nova leitura na Swell para este encerramento.

### 7.3 Pendências para a 4c

- Reescrever o §2.6 de `remote_data_acquisition.rst` (hoje com uma nota de atualização; o formato de bloco ali é o do v3).
- Decidir se o LECTK passa a ler os stores Zarr v3 da EDH (exige zarr 3 e outro xarray).
- Retenção da pasta da Swell (2,8 GB) no fim da etapa 4 (D-028).

