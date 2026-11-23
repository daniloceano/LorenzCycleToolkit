# Etapa 4b·1 — cliente WEkEO nativo (`--wekeo`), com a Etapa 0 (4b·0) (23/11/2026)

**Natureza:** registro de implementação. Acrescenta a rota `--wekeo` ao programa principal (D-036). Nenhum valor de
dado de outra rota muda, e o cálculo do LEC também não muda. `main`, tags e releases continuam inalterados.

**Base:** `develop` em `1a714dd` (conferido antes da branch); `main` em `48fe8e2`.

**Estado:** implementação, testes offline e teste real na Swell concluídos (§5). O teste real **não demonstrou
identidade bit a bit** com a rota CDS da 2b: `u`, `v` e `w` são idênticos; `t` e `z` diferem por arredondamento; o LEC
difere em no máximo 0,007 %.

## Resumo em linguagem simples

O usuário passa a poder baixar o ERA5 da WEkEO direto pelo LECTK, com `--wekeo`, no lugar de `--cdsapi`. O pedido é o
mesmo do CDS: mesma caixa, mesmos níveis e mesmos horários. O arquivo final também é igual, nos nomes, nas unidades, na
ordem das dimensões e na convenção de longitude. A biblioteca `hda` é opcional e só precisa ser instalada por quem usa a
WEkEO.

## 1. Etapa 0 (4b·0): investigação

### 1.1 Harness da 2b na Swell (só leitura, D-034)

O responsável foi avisado no início da leitura. Nenhuma credencial foi aberta, copiada ou citada; só se verificou que
`~/.hdarc` existe.

- **Código:** `lectk_bench2b_20261109/harness/bench2b.py`, função `route_wekeo`; o mesmo desenho nas versões da rodada
  2.
- **Biblioteca:** `hda` 2.40, no `venv` da campanha.
- **Dataset:** `EO:ECMWF:DAT:REANALYSIS_ERA5_PRESSURE_LEVELS`.
- **Pedido:** um por mês do caso, com `product_type`, `variable`, `year`, `month`, `day`, `time` e `pressure_level` no
  vocabulário do CDS, mais `data_format: netcdf` e `download_format: unarchived`. A caixa vai em
  `bbox` = [W, S, E, N] (ordem confirmada no pré-voo de 10/11).
- **Combinação de dias e horários:** o pedido mensal cobre todos os dias com todos os horários. No caso 19790205
  chegaram 40 horários para os 35 do caso, e a grade comum da 2b descartou os excedentes.
- **Arquivo recebido:** um NetCDF por mês, no mesmo formato do CDS: `valid_time`, `pressure_level`, `latitude`,
  `longitude`, com as coordenadas `number` e `expver`, float32 e longitudes negativas preservadas (−71,75 a 65,5 no caso
  19790205). Na 2b, CDS e WEkEO entregaram arquivos idênticos (fase 2: 240/240 termos idênticos).
- **Download:** o harness usava `SearchResults.download()`.

### 1.2 Documentação atual da HDA

- **Versão estável:** `hda` 2.40 (PyPI, 07/09/2026; também no conda-forge), Python ≥ 3.10. Depende de `requests` e
  `tqdm`.
- **Credencial:** `~/.hdarc` com as linhas `user:` e `password:`, ou as variáveis `HDA_USER` e `HDA_PASSWORD`. O
  código 2.40 também lê `HDA_RC` (outro arquivo) e `HDA_URL`. Desde a 2.40, contas EUMETSAT e CMEMS também são aceitas.
- **Limites:** 500 itens baixáveis por hora, sem cota mensal. Repetir um download já feito não conta de novo. A
  documentação não cita limite de tamanho nem de pedidos por minuto, nem descreve o HTTP 429. Pelo código, o 429 gera
  `QuotaReachedError`.
- **Fila:** pelo código, o pedido (`POST dataaccess/download`) devolve um `download_id`, e o cliente consulta o estado
  com `HEAD` (202 = em preparo, 200 = pronto). **Não há operação de cancelar ou apagar pedidos no cliente, e a
  documentação não diz se um pedido expira.**
- **Formato:** com `data_format=netcdf` e `download_format=unarchived`, chega um NetCDF; a documentação não trata do
  caso zip.

**Problemas da biblioteca que o LECTK evita:**

- `SearchResults.download()` registra a falha no log e **não lança erro**: um download falho passaria despercebido.
- O padrão é de **500 tentativas** HTTP.
- A espera pelo pedido **não tem prazo**.
- A biblioteca registra o token de acesso no nível DEBUG.

### 1.3 O que foi reaproveitado

| Parte | Origem | Uso na WEkEO |
|---|---|---|
| Caixa (15° de margem, arredondada para graus inteiros) | `ncar_gdex.area_from_track`, mesma regra de `tools.get_cdsapi_data` | igual |
| Horários de cada dia | `ncar_gdex.daily_hours`, mesma regra do CDS | igual; os excedentes do pedido mensal são descartados na montagem |
| Greenwich | o CDS repassa as longitudes como vêm do track | igual: o `bbox` usa os mesmos W e E |
| Checagem prévia dos horários do track | `preprocessing.check_planned_times` | igual |
| Montagem | concatenação em `valid_time` e `to_netcdf`, como no CDS | igual, mais a seleção dos horários exatos |
| Namelist | `inputs/namelist_ERA5-cdsapi` | igual |
| Teto padrão de 20 GiB e `--max-gib` | NCAR (D-035) | igual |
| Commit do LECTK na proveniência | `ncar_gdex._lectk_commit` | igual |

O `ByteBudget` do NCAR não foi reaproveitado, porque ele conta bytes no socket do `urllib`, e a `hda` usa `requests`.
A WEkEO conta os bytes do corpo de cada resposta.

## 2. Decisões de projeto (D-036)

Todas seguem os padrões da tarefa; a Etapa 0 não contrariou nenhum.

1. **Interface:** `--wekeo`, mutuamente exclusivo com `--cdsapi` e `--ncar`. Valem `--trackfile`, `--time-resolution`,
   `--max-gib` (agora com `--ncar` ou `--wekeo`) e `--delete-downloaded-input`.
2. **Credencial:** quem lê é só a `hda` (`Configuration()`). O LECTK nunca lê, registra nem mostra os valores. Sem
   credencial, a mensagem explica o `~/.hdarc` e as variáveis de ambiente e aponta o cadastro
   (`https://www.wekeo.eu/register`), sem repetir o texto da exceção. O logger `hda` fica em WARNING para não registrar
   o token.
3. **Dependência:** extra `[wekeo]` no `pyproject.toml` (`hda>=2.40` e `requests>=2.5.0`, este importado diretamente pelo
   módulo). `requirements.txt` e `environment.yml` só trazem um comentário, porque a dependência não é obrigatória. A
   `hda` é importada só quando `--wekeo` é usado; sem ela, a mensagem dá o comando de instalação.
4. **Pedido:** um por mês; variáveis e níveis iguais aos do CDS; mesma caixa e mesmo `bbox` [W, S, E, N].
5. **Saída:** mesmo contrato do CDS; proveniência em `<arquivo>.provenance.json` (rota, dataset, versão da `hda`,
   commit do LECTK, caixa, horários, teto, bytes, pedidos com id, bytes e sha256 de cada resposta, sha256 do arquivo).
6. **Proteções:**
   - **Volume:** antes de qualquer pedido, o limite superior (float32 sem compressão, mais 1 MiB por resposta) é
     comparado com o teto e com o disco livre (o dobro da estimativa). Na recusa, a mensagem dá o `--max-gib`
     necessário. Durante o download, uma resposta maior que o saldo não é lida, e a leitura para quando o saldo
     acaba. Bytes de tentativas com falha também contam.
   - **Tentativas:** 3 tentativas HTTP no cliente (em vez de 500) e no máximo 2 de download por pedido pronto, sempre
     no mesmo pedido. A espera na fila tem prazo de 1 h.
   - **Pedidos órfãos:** proteção proporcional à D-033. O registro `<arquivo>.wekeo-orders.json` é gravado **antes**
     de cada envio e recebe o id quando a WEkEO responde. Uma nova execução no mesmo destino retoma os pedidos
     registrados sem reenviar. Um envio registrado sem id (estado remoto desconhecido) bloqueia a rota até o usuário
     conferir a conta e apagar o registro. No fim, o registro e as respostas mensais são apagados, e a proveniência
     guarda os ids.
7. **Reuso:** um destino que já existe é reaproveitado sem conferência, como no `--cdsapi`.

## 3. Arquivos alterados

- `src/utils/wekeo_hda.py` (novo): plano, cliente, registro de pedidos, download, montagem e proveniência.
- `lorenzcycletoolkit.py`: `--wekeo`, exclusão mútua, `--max-gib` e `--delete-downloaded-input` com `--wekeo`.
- `src/utils/preprocessing.py`: ramo `--wekeo` em `get_data` e namelist ERA5 com `--wekeo`.
- `src/utils/input_cleanup.py`, `src/utils/run_timing.py`: `--wekeo` reconhecido (exclusão do arquivo e origem
  `wekeo`).
- `pyproject.toml`, `requirements.txt`, `environment.yml`: extra opcional `[wekeo]`.
- `tests/test_wekeo_hda.py` (novo, 20 testes) e `tests/test_packaging.py` (extras opcionais aceitos; teste novo de que a
  `hda` não é obrigatória).
- `docs/source/remote_data_acquisition.rst`: nota curta de disponibilidade (a seção completa fica para a 4c).
- `CHANGELOG.md`, `docs/evolution/decisions.md` (D-036), `docs/evolution/queue.md` e este registro.

## 4. Testes (Mac Mini, ambiente `lorenz`, Python 3.12)

| | Antes (`1a714dd`) | Depois |
|---|---|---|
| Passaram | 644 | 665 (644 + 20 novos da WEkEO + 1 novo de empacotamento) |
| Falharam | 1 (`test_installed_metadata_matches_pyproject`) | a mesma 1 |
| Erros | 16 (fixture do Natural Earth, arquivos ausentes no Mac; ficam para o CI) | os mesmos 16 |

A falha de `test_installed_metadata_matches_pyproject` vem de uma instalação local desatualizada do pacote (faltam
`cmocean`, `dask`, `shapely`, `cyclophaser` e `lorenz-phase-space` nos metadados instalados). Ela já existia antes desta
etapa e não depende dela.

Os testes novos são todos offline, com a `hda` e o `cdsapi` simulados:

- **pedido:** caixa e horários iguais aos do CDS, divisão por mês com horários excedentes, Greenwich, campos e níveis
  conferidos contra os pedidos que o próprio `get_cdsapi_data` faz;
- **equivalência:** o mesmo ERA5 sintético passa pelas duas rotas, e os arquivos resultantes são idênticos
  (`xr.testing.assert_identical`: valores, nomes, unidades, atributos, dimensões, ordem e longitudes);
- **pacote ou credencial ausente:** mensagens de instalação e de cadastro; 3 tentativas; logger da `hda` em WARNING;
- **volume:** estimativa acima do teto (nenhuma chamada ao cliente), resposta acima do saldo com e sem
  `Content-Length`, falta de disco;
- **pedidos:** registro gravado antes do envio; pedido registrado retomado sem reenvio; envio sem id bloqueia sem
  nenhuma chamada; pedido com falha; prazo da fila; resposta que não é NetCDF;
- **linha de comando:** combinações inválidas saem com código 2.

## 5. Teste real na Swell (23/11/2026)

**Autorização:** o responsável aprovou o diff (commit `cb54b69` na branch `feat/stage-4b1-wekeo`) e a comparação na
região comum. O teto de 1 GiB foi revisto após consulta: o limite superior calculado era de 3,54 GiB, e o responsável
avaliou o risco como baixo e liberou mais volume. O teste rodou com `--max-gib 4`.

**Montagem:**

- pasta nova `/p1-swell/danilocs/lectk_stage4b1_wekeo_20261123/`, fora do repositório (≈ 3 GB);
- código: `git archive` de `cb54b69`;
- `venv` com `--system-site-packages` sobre o ambiente `lorenz` e `hda` 2.40;
- track `9da9d809…` (o selado da 2b) e caixa e namelist copiados da pasta da 4a.1.

Nenhum dado existente foi alterado.

**Comando:**
`lorenzcycletoolkit.py case_19790205.nc -f -r --wekeo --trackfile inputs/track --time-resolution 6 --max-gib 4 -o wekeo_19790205`.

### 5.1 Tentativas

| Tentativa | Resultado |
|---|---|
| 1 (19:38 UTC) | Parou em 3 s, sem rede: faltava `inputs/namelist_ERA5-cdsapi` na pasta de execução (erro de montagem meu; o `--cdsapi` se comporta igual). |
| 2 (19:39–20:03 UTC) | Download e montagem concluídos. O LEC falhou com `KeyError: 'time'`, porque `run_lec_analysis` relê `inputs/namelist`, que era o namelist da fase 2 (`time`/`level`). É um defeito anterior a esta etapa, que também afeta o `--cdsapi` (§7). |
| 3 (20:04 UTC) | `inputs/namelist_ERA5-cdsapi` copiado para `inputs/namelist`, como manda `configuration.rst`. O arquivo foi reaproveitado (sem rede) e o LEC rodou até o fim: código 0, 148 s. |

**Vigia de tráfego:** não funcionou nas tentativas 1 e 2. O `pgrep` devolveu um PID que não era o do programa, e o vigia
terminou no primeiro ciclo. O tráfego da interface de rede **não foi medido**; o único contador é o do LECTK.

### 5.2 Medidas (tentativa 2)

| Grandeza | Valor |
|---|---|
| Pedidos | 1 (1979-03; 9 dias × 4 horários = 36, dos quais ficam 34); pedido `6b04965ff1874199ed1aa81f` |
| Espera na fila da WEkEO | 1 017 s |
| Transferência | 273 s |
| Montagem | 135 s |
| Aquisição total | 19:39:05 a 20:02:57 UTC (≈ 24 min) |
| Bytes recebidos (contador do LECTK) | 1 719 116 096 B (1,60 GiB), abaixo dos 4 GiB; limite superior previsto: 3,54 GiB |
| Resposta mensal | sha256 `4390fcab…` |
| Arquivo montado | 1 497 346 700 B; sha256 `cf2ea3a20f528640d094336ffaa0aef29b5aef16bbd33da5c2a47e00b80c5008`; 34 horários × 37 níveis × 225 × 633; N−24 W−82 S−80 E76 |
| CSV do LEC do comando principal | 34 × 27; sha256 `8c919a68…` (sem referência anterior: a caixa de dados de 15° é outra) |

A proveniência registra `lectk_commit: null`, porque o export do `git archive` não tem `.git`. O código é o do commit
`cb54b69`.

### 5.3 Critério 1: dados na região comum

O sha256 não é comparável por construção: a caixa de 15° e os 37 níveis do comando principal diferem da caixa de 5° e
dos 13 níveis da 2b. A comparação foi feita valor a valor na região comum (34 horários × 13 níveis × 142 × 550, 5
variáveis).

| Comparação | u, v, w | t | z |
|---|---|---|---|
| Novo × WEkEO da 2b (arquivo bruto) | idênticos | máx. \|Δ\| 1,2·10⁻⁴ K | máx. \|Δ\| 0,0029 m² s⁻² |
| Novo × CDS da 2b (grade comum, a referência da fase 2) | idênticos | diferem | diferem |
| Novo × kerchunk da 4a.1 (mesma caixa de 15°, 37 níveis) | idênticos | máx. 2,4·10⁻⁴ K, em 2,2 % dos pontos | máx. 0,0039 m² s⁻², em 0,7 % dos pontos |

**Critério 1 não atendido no sentido estrito:** não há identidade bit a bit com o arquivo CDS da 2b. Os arquivos CDS e
WEkEO da 2b, porém, continuam idênticos entre si.

**[HIPÓTESE, não verificada]** A diferença viria da compactação GRIB do ERA5, refeita pelo serviço para cada caixa
pedida, e não da rota. Ela ficaria só em `t` e `z`, abaixo de 2,5·10⁻⁴ K e 0,004 m² s⁻². Para separar o efeito da rota
do efeito da caixa, seria preciso um pedido CDS com a mesma caixa de 15°, que não foi feito.

### 5.4 Critério 2: LEC

O CSV da fase 2 (`935ed785…`, 35 horários) **não é reproduzível**: o comando principal pede 34 horários (até 26/03
18 UTC), e o pedido selado da 2b tinha 35. Em vez disso, o LEC da fase 2 (export de `develop` em `068c3a9`,
`launch_lectk.py`, `-f -r`, caixa selada e namelist da fase 2) rodou sobre dois arquivos de grade comum de 34 horários:

- o arquivo CDS da 2b, recortado nos 34 horários;
- o arquivo WEkEO novo, convertido com `to_common_grid` do harness da rodada 2.

| | CDS da 2b (34 h) | WEkEO novo (34 h) |
|---|---|---|
| Código de saída | 0 | 0 |
| CSV (34 × 26) | `d46d82ba…` | `0fe175d2…` |

- **Diferença:** os CSVs não são idênticos. A maior diferença relativa à escala do termo, como o e_max da fase 2, é de
  **0,007 %** (RGz); depois vêm ∂Az/∂t com 0,006 % e RGe com 0,005 %.
- **Sinal:** nenhum termo trocou de sinal.

**Critério 2: atendido com a ressalva de que a igualdade não é exata.** O LEC roda até o fim. A diferença em relação à
rota CDS fica na classe "pequena" da fase 2, como ARCO e kerchunk, muito abaixo do limiar de 1 %. Não se declara
equivalência bit a bit.

Scripts e saídas: `compare.py` e `lec34/` na pasta do teste. O destino da pasta (≈ 3 GB) fica para o fim da etapa 4,
junto com as demais (D-028).

## 6. Limitações

- **Cache na medida da 2b:** os pedidos das três repetições da 2b eram idênticos, e a WEkEO e o CDS podem ter servido
  resultados guardados. O tempo de ~1 min por caso pode ser só o tempo de entrega de um resultado já pronto, e não
  deve ser lido como tempo de fila de um pedido novo.
- **Estado remoto:** a `hda` não oferece cancelamento. Um pedido interrompido fica na WEkEO, mas o registro guarda o
  seu id, e a nova execução o retoma.
- **Reenvio HTTP:** o `POST` do pedido passa pelas tentativas internas da `hda` (até 3, em erros 5xx). Uma resposta
  perdida pode criar um segundo pedido no servidor sem id registrado. O risco é baixo e não foi testado.
- **Greenwich:** o tratamento é o mesmo do CDS. Um track em 0–360° que cruza 0° geraria W > E nas duas rotas, e isso
  não muda aqui.
- **Concorrência:** duas execuções no mesmo destino ao mesmo tempo não são impedidas. O registro só protege execuções
  sucessivas.
- **Estimativa de volume:** é um limite superior sem compressão, e o volume real costuma ser menor.
- **Destinos existentes:** são reaproveitados sem conferência, como no CDS.

## 7. Achado fora do escopo: namelist do LEC com `--cdsapi`/`--wekeo`

`prepare_data` troca para `inputs/namelist_ERA5-cdsapi` com `--cdsapi` (e agora com `--wekeo`). Já `run_lec_analysis`
(`lorenzcycletoolkit.py`) relê sempre `inputs/namelist`. Quando os dois namelists diferem, o LEC falha depois do
download (`KeyError: 'time'`). `configuration.rst` diz que o usuário não precisa configurar o namelist com
`--cdsapi`, mas o exemplo manda copiar. **Não foi corrigido nesta etapa;** fica como pendência para a 4c ou para uma
correção própria.

## Fora deste registro

ARCO e DestinE, robustez do kerchunk, conferência dos 34 × 35 horários e reescrita da documentação (4c). Merge, CI e
releitura entram depois, por nota.
