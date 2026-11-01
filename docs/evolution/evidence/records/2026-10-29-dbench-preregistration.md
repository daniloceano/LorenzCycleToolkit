# D-Bench — pré-registro do benchmark CDS × NCAR (29/10/2026)

**Estado:** escrito e selado **antes de qualquer aquisição real da campanha**. Critérios, tolerâncias, ordem, orçamento e regras de parada não mudam depois de observados os resultados; um desvio necessário é registrado como desvio, com motivo, e os resultados afetados são reportados separadamente. O selo é o commit de `develop` que integrar este arquivo junto com o código que o executa (`benchmarks/d_bench/dbench.py`, `src/utils/cds_guard.py`). Plano de referência: [front-D-bench-cds-ncar.md](../../fronts/front-D-bench-cds-ncar.md).

## 1. População e configuração equivalente

- **Casos:** as dez trajetórias reais de 1979–1980 fornecidas pelo responsável, com 438 posições a cada 3 h. Os SHA-256 registrados no plano são verificados pelo runner antes de cada execução, e um arquivo divergente é recusado. Nenhuma trajetória sintética.
- **Fontes:** CDS `reanalysis-era5-pressure-levels` × NCAR GDEX `d633000`.
- **Campos:** temperatura, geopotencial, omega e vento zonal e meridional, nos 37 níveis de 1 a 1000 hPa, grade de 0,25°, passo de 3 h (`--time-resolution 3`).
- **Área:** a mesma regra nos dois caminhos: extremos do track inteiro, margem de 15°, arredondamento para inteiro.
- **LECTK:** `-t -r` (framework móvel com resíduos) e o mesmo track nos dois provedores.
- **Volume estimado:** 64 pedidos diários CDS e 25,67 GiB de valores float32 por passagem NCAR (piloto: 0,50 GiB).
- **Namelist:** cada execução roda num diretório privado novo, com `inputs/namelist` e `inputs/namelist_ERA5-cdsapi` iguais ao preset versionado `inputs/namelist_ERA5-cdsapi` (SHA-256 `bc24c2048c99d2dc779c7645fadaf4c77cd517d01002d05986e3a7339f69d0c7`, idêntico a `inputs/namelist.default`).
  - **Namelist de trabalho do servidor:** o `inputs/namelist` local (SHA-256 `57faf1bc…b13f`) segue a convenção WRF (`ght`, `omega`, `lon`, `lat`, `levels`) e **não é usado**; o arquivo global também não é alterado. Com ele, a preparação `--cdsapi` e o cálculo LEC leriam mapeamentos diferentes, e o `--ncar` nomearia a saída de outro modo. Por isso a correção é só de configuração, sem código.
  - **Nomes nas saídas:** com o mesmo namelist, as saídas CDS e NCAR usam os mesmos nomes (`t, z, w, u, v`; `valid_time`, `pressure_level`, `latitude`, `longitude`). Diferenças de convenção de longitude (−180…180 × 0…360) e de ordem dos eixos são normalizadas na comparação; não são diferenças físicas.
  - **Verificação offline:** o arquivo de amostra ERA5 renomeado para a convenção CDS e calculado com o preset dá CSV LEC **byte a byte idêntico** ao cálculo com o preset `namelist_ERA5`.

## 2. Aquisição fria e reuso

- **Aquisição fria:** diretório de execução novo e exclusivo (`ROOT/runs/<caso>/<provedor>/r<n>`, recusado se existir), destino ausente e nenhuma resposta parcial anterior.
  - CDS em modo protegido: recusa um destino existente.
  - NCAR: as respostas parciais ficam no próprio diretório novo.
- **Reuso de cache:** é um estrato separado, não planejado nesta campanha, e nunca é contado como aquisição.
- **Preservação dos NetCDF:** todos são mantidos (sem `--delete-downloaded-input`) para auditoria e reprocessamento.

## 3. Medidas

- **Primária:** `timing.phases.acquisition.seconds` do `run_status.json` (relógio monotônico, da chamada de aquisição até o NetCDF aceito pelo LECTK).
- **Secundárias:**
  - preparação, LEC e total interno;
  - relógio externo do runner (ponta a ponta, com a inicialização do Python);
  - pico de RSS e CPU (`/usr/bin/time -v`) e load average do servidor antes e depois;
  - bytes e SHA-256 do NetCDF;
  - totais NCAR da proveniência (respostas, bytes, tentativas);
  - por pedido CDS: horários de submissão e download no ledger privado.
- Tempo por pedido é diagnóstico, nunca substitui o tempo do caso completo.

## 4. Ordem, repetições e orçamento

- **Ordem:** piloto `19800296` primeiro, depois os outros nove em ordem pré-aleatorizada (semente `20261029`, `dbench.py plan`), com a ordem dos provedores sorteada por caso. Execução **sequencial**: um provedor e um caso por vez, sem concorrência CDS.
- **Repetições:**
  - **r1** dos dez casos nos dois provedores é a amostra primária.
  - **r2** (opcional) só depois que os **dez pares r1 terminarem em `success`** (verificado pelo runner), na mesma ordem e dentro do mesmo teto de 96 h.
  - **r3** não está planejada.
  - A análise primária usa só r1; r2 é descritiva (variabilidade).
- **Limites por caso e provedor:**
  - prazo global de aquisição de **12 h**, monotônico;
  - CDS: 3 tentativas HTTP, timeout de rede de 120 s, um pedido ativo por ledger **e por conta CDS neste servidor** (trava por credencial, nomeada por hash truncado; execuções em outros hosts não ficam cobertas), ledger antes de cada submissão e delete best-effort do próprio job após o prazo (ver §7);
  - NCAR: limites nativos (8 GiB por execução, 4 tentativas, timeout de 300 s);
  - parada externa se a execução passar do prazo em 30 min.
- **Campanha:** teto de segurança de **96 h**, contado desde o primeiro início e persistido em `ROOT` (uma retomada não reinicia o relógio). Uma execução só começa se o prazo dela mais a margem de parada (12 h 30 min) ainda couberem no teto. Os limites nunca são aumentados automaticamente. Se as 96 h não bastarem, a campanha é registrada como **parcial**, e os pares faltantes não são apresentados como concluídos.
- **Recursos verificados em 29/10:** 127 TiB livres no volume de trabalho, 239 GiB livres no disco raiz (não usado: `TMPDIR` fica no diretório da execução), 1 TiB de RAM e 112 núcleos com load em torno de 110 (servidor compartilhado: a carga é registrada, e preparação e LEC são sensíveis a ela).
- **Armazenamento previsto:** até cerca de 52 GiB por passagem (NCAR mais CDS, sem contar a compressão) e cerca de 105 GiB com r2. Antes de cada execução, o espaço livre precisa ser de pelo menos 3 × a estimativa + 5 GiB.

## 5. Integridade e tolerâncias (por par CDS × NCAR)

Os dois arquivos vêm do mesmo ERA5 final de 1979–1980, empacotado em GRIB de 16 bits por campo e nível. A expectativa é igualdade, ou diferença de no máximo um passo de quantização. Tolerâncias absolutas, nas unidades do arquivo:

| Campo | Tolerância |
|---|---|
| `t` (K) | 0,005 |
| `z` (m² s⁻²) | 2,0 |
| `w` (Pa s⁻¹) | 0,001 |
| `u`, `v` (m s⁻¹) | 0,01 |

O par é **não equivalente** se qualquer um destes ocorrer:

- conjuntos de horários diferentes, ou algum horário do track ausente;
- níveis diferentes dos 37;
- latitudes ou longitudes diferentes (após normalizar a longitude);
- campo ausente;
- unidades diferentes;
- valor não finito;
- `max |Δ|` acima da tolerância;
- CSV LEC com tempos ou colunas diferentes, ou alguma coluna fora da regra abaixo.

**Regra das colunas do CSV LEC:**

- **Escala:** `escala = max(max|CDS|, max|NCAR|)`, simétrica e ≥ 0.
- **Tolerância:** `tol = max(10⁻⁴ × escala, 10⁻⁹)`, nas unidades da coluna. O piso absoluto de 10⁻⁹ é a regra explícita para colunas nulas ou quase nulas.
- **Aprovação:** a coluna passa se `max|Δ| ≤ tol`.
- **Valores ausentes ou infinitos:** NaN nas mesmas posições é aceito; NaN em posições diferentes ou qualquer valor infinito torna a coluna não equivalente.

Também são reportadas a fração de valores idênticos e a média de |Δ|. Igualdade bit a bit não é exigida, e nenhuma normalização oculta é feita.

## 6. Falhas e pares válidos

- **Resultado de cada execução:**
  - `success`: status `success` e saída 0;
  - `censored`: prazo excedido;
  - `failed`: qualquer outro caso.
  - Censuras e falhas nunca viram zero nem são imputadas.
- **Par válido:** os dois provedores com `success` e o par equivalente.
- **CDS:**
  - nada é submetido se a listagem de pedidos pendentes da conta falhar ou mostrar algum pendente;
  - na primeira falha, censura ou interrupção, o CDS para pelo resto da campanha (o ledger também bloqueia), até verificação manual e `resolve`.
- **NCAR:** para após duas falhas consecutivas.
- **Medições unilaterais:** se só um provedor puder continuar, os resultados válidos dele são preservados como diagnóstico unilateral, nunca como comparação pareada.

## 7. Portão do piloto para a campanha

O piloto `19800296` só libera a campanha se todas as condições valerem:

1. os dois provedores terminam com `success`;
2. o par é equivalente segundo §5;
3. não houve evento de segurança CDS (bloqueio, pedido pendente sobrando, erro de credencial, cota ou política);
4. o pico de RSS ficou abaixo de 25% da RAM e o espaço livre segue dentro da regra de §4;
5. aquisição, preparação e LEC têm tempo medido nas duas execuções.

Se o CDS falhar por limite da conta ou do serviço, a parte NCAR pode continuar só como diagnóstico. O delete best-effort (`--cds-delete-on-stop`) entra no runner **somente com a aprovação do responsável**. Ele não é garantia de cancelamento, e o efeito registrado é o que o próprio CDS informar.

## 8. Análise pré-especificada (`dbench.py summary`)

- **Razão por caso:** `r = aquisição NCAR / aquisição CDS` em r1, nos pares válidos.
- **Estatística:** mediana de `r` com intervalo livre de distribuição, por estatísticas de ordem binomial (o mesmo fundamento do teste do sinal), com cobertura de **pelo menos** 95%. Para n = 10, o intervalo é (x₍₂₎, x₍₉₎), com cobertura exata de 97,9%. A cobertura exata é reportada. O intervalo descreve a mediana da razão para essas trajetórias, nas condições observadas. **Não é teste de significância**, e com dez pares é largo.
- **Regra de decisão:**
  - **vantagem prática NCAR:** n ≥ 6, mediana ≤ 0,80 e limite superior < 1;
  - **vantagem prática CDS:** n ≥ 6, mediana ≥ 1,25 e limite inferior > 1;
  - caso contrário: **nenhuma superioridade demonstrada**, o que não prova velocidade igual;
  - com n < 6, nenhuma comparação de velocidade.
- O tempo total ponta a ponta e as fases são reportados como secundários. O resultado vale para essas dez trajetórias e para as condições de serviço observadas no período, não como lei geral.

## 9. Riscos residuais conhecidos

- Interromper o cliente não cancela o pedido remoto. O delete best-effort atua só no próprio job, por ID; seu efeito sobre um job já em execução não foi demonstrado (verificação offline em §10).
- O `resolve` do ledger é uma **atestação humana**: exige `--confirm-no-pending` e uma nota, é recusado enquanto uma aquisição detém a trava e não consulta o CDS. A campanha ainda lista os pedidos pendentes da conta antes de cada caso CDS.
- Limite de pedidos na fila por dataset e por conta, e as "Data Stores Service protection measures" anunciadas pelo CDS em 06/10/2026, podem rejeitar pedidos.
- Servidor compartilhado e carregado.
- População de um período e de uma região restritos.

## 10. Verificação offline de `Remote.delete()` (ecmwf-datastores-client 0.4.2)

- **O que o código faz:** `Remote.delete()` envia **um** `DELETE` para a URL do próprio job (`…/retrieve/v1/jobs/<id>`), conforme simulado com uma sessão falsa e sem rede. A própria biblioteca usa esse método na limpeza automática (`cleanup`) e trata `dismissed`/`deleted` como estados terminais.
- **O padrão:** a OGC API – Processes, à qual o CDS diz ter alinhado seus estados, define `DELETE /jobs/{id}` como cancelamento ou dispensa do job.
- **Conclusão:**
  - **seletivo por construção;**
  - efeito provável sobre pedidos `accepted` (na fila);
  - **efeito não demonstrado** sobre um pedido já `running`.
  - O guard só envia o delete com credencial isolada, apenas para o objeto devolvido pela submissão do próprio processo (nunca para um ID tirado do ledger ou de listagem) e apenas enquanto o status lido for `accepted` ou `running`. Um job concluído, falho ou de estado desconhecido não recebe delete (`delete_skipped`), e seus resultados nunca são removidos.
  - O guard registra a resposta e o status lido depois, e mantém o bloqueio em qualquer caso.
- Nenhum cancelamento real nem teste autenticado foi feito.

## 11. Incidente de 29/10/2026 (antes deste pré-registro)

Numa verificação offline do agente, um arquivo de entrada ausente levou o modo `--cdsapi` convencional a submeter **um** pedido real com a credencial global do servidor, de outra conta, não a da D-Bench.

- **Resultado:** o pedido foi rejeitado pelo CDS por fila cheia daquela conta. Nada foi baixado, nenhum job ficou pendente e o processo terminou.
- **Causas:**
  - comandos encadeados sem parada em erro;
  - a exclusão automática do NetCDF pelo framework móvel, já corrigida: os arquivos agora são preservados e há uma flag explícita de exclusão;
  - a leitura implícita de `~/.cdsapirc` no modo convencional.
- **Medidas:**
  - verificações manuais da CLI com `--cdsapi`/`--ncar` só pelo lançador offline;
  - o modo protegido exige credencial explícita;
  - teste de regressão: entrada ausente com `--cdsapi` sem rede e sem `~/.cdsapirc` alcançável.
