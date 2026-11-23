# Etapa 4b·1 — cliente WEkEO nativo (`--wekeo`), com a Etapa 0 (4b·0) (23/11/2026)

**Natureza:** registro de implementação. Acrescenta a rota `--wekeo` ao programa principal (D-036). Nenhum valor de
dado de outra rota muda, e o cálculo do LEC também não muda. `main`, tags e releases continuam inalterados.

**Base:** `develop` em `1a714dd` (conferido antes da branch); `main` em `48fe8e2`.

**Estado:** implementação e testes offline concluídos; **teste real na Swell pendente de consulta ao responsável**
(§5).

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

## 5. Teste real na Swell

**Não executado.** Aguarda consulta ao responsável (tarefa, Parte 3).

**Comparabilidade, a esclarecer antes do teste:** o comando principal pede a caixa do track com 15° de margem e os 37
níveis. O arquivo WEkEO/CDS da 2b usou a caixa selada com 5° de margem e 13 níveis, num pedido do harness. Por
construção, o sha256 não pode coincidir com o da 2b. O que se pode verificar é a igualdade valor a valor na região
comum (13 níveis, caixa da 2b, horários do caso), como foi feito na 4a.1 (§7.3 do registro da 4a). Com isso, o CSV do LEC
da fase 2 (`-f -r` na caixa selada, 13 níveis) só é reproduzível rodando o LEC sobre o recorte da região comum. A
proposta fica para a consulta.

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

## Fora deste registro

ARCO e DestinE, robustez do kerchunk, conferência dos 34 × 35 horários e reescrita da documentação (4c). Merge, CI e
releitura entram depois, por nota.
