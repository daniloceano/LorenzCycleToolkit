# Frente C — avaliação documental de fontes ERA5 e requisitos do LECTK

**Data da investigação:** 2026-10-25. **Base do código inspecionado:** `develop` `e1ce2bdaefa9bff86e6443bbcc03f66f0e12f9d1`. **Natureza:** investigação documental, sem novas execuções numéricas nem integração de provedores. **Estado deste documento:** **HOMOLOGADO — PASS DOCUMENTAL, Frente C encerrada em 25/10/2026**, após revisão aprovada, PR #36 integrado em `develop` (`abaab7add04c21c694a4aa90da38c7d1465e5add`) e releitura de três arquivos com hashes coincidentes. Não representa benchmark executado nem aprovação de arquitetura.

## 1. Objetivo, universo e critério de aceitação

Esta frente caracteriza o contrato atual de entrada do LorenzCycleToolkit (LECTK), separa provedor, produto, método de acesso e formato, e classifica a elegibilidade documental de fontes ERA5 como insumo para a futura Frente D. Respeita o `infile` NetCDF local como fluxo principal e a CDS API já existente como mecanismo opcional.

**Classes predefinidas:** (1) **compatível documentalmente** — satisfaz documentalmente o fluxo sem novo adaptador; (2) **compatível sob condições explícitas** — pode cumprir o contrato mediante seleção, conversão, configuração e/ou verificação adicional especificadas; (3) **incompatível com requisitos essenciais** — o produto isolado analisado não disponibiliza requisito imprescindível; (4) **evidência insuficiente** — a documentação não sustenta classificação conclusiva. As classes referem-se ao *produto específico analisado*, não a uma proibição permanente do provedor.

**PASS DOCUMENTAL técnico** exige: contrato rastreável; caracterização de NetCDF e CDS; avaliação do universo delimitado com fontes e lacunas; separação entre ciência, operação e casos históricos; ausência de alegação sem prova de equivalência/rapidez/custo; encaminhamento fundamentado para D. A integração documental e o encerramento administrativo exigem, adicionalmente, revisão e portões do fluxo Git. A ausência de um provedor elegível não seria, por si, FAIL da frente.

## 2. Contrato de entrada atualmente observado

| Dimensão | Requisito ou comportamento do código | Natureza / limite |
|---|---|---|
| Arquivo | Argumento posicional `infile` contendo campos meteorológicos; `xr.open_dataset(infile)`; `namelist` identifica os nomes das variáveis/coordenadas. | Interface operacional atual; novos formatos remotos demandam preparação para `infile` ou projeto de integração futuro. |
| Variáveis | Temperatura `T` (K), geopotencial `z` (m² s⁻²; há conversão de altura geopotencial em m), componentes do vento `u` e `v` (m s⁻¹), velocidade vertical de pressão `omega/w` (Pa s⁻¹). | Requisito meteorológico essencial; nomes exatos parametrizáveis no `namelist`, definições/unidades não intercambiáveis indiscriminadamente. |
| Eixos | Tempo, pressão isobárica, latitude e longitude em grade regular; níveis convertidos internamente para Pa; longitude e ordenação de coordenadas normalizadas. | Operacional e científico; dados em níveis de modelo não são entrada direta equivalente. Unidades verticais ausentes podem disparar suposição de hPa. |
| Pressão | Dois ou mais níveis de pressão utilizáveis; harmonização entre tempos/domínios requer níveis completos comuns, sem preencher lacunas internas de dados ausentes. O pré-processamento limita seleção superior nominal a 10 hPa quando disponível (`slice(1000, lowest_level)` em Pa). | Um subconjunto vertical deliberadamente esparso não é, por si só, lacuna de dados; limites efetivos de integração devem ser registrados caso a caso. Não há exigência universal de 37 níveis. |
| Tempo | Instantes ordenados e coerentes; tendências dos termos requerem séries com instantes suficientes; modo móvel calcula derivada temporal da temperatura. | A frequência de 3 h é referência histórica de uso, não requisito universal. Não presumir comparabilidade de 3 h e 6 h. |
| Espaço | Campo regular com cobertura do domínio fixo ou dos domínios móveis ao longo da trajetória, incluindo suporte necessário para gradientes e termos de borda. | Resolução e tamanho de caixa dependem do caso; não há limiar universal empiricamente aprovado de resolução. |
| Modos | `--fixed` seleciona caixa definida; `--track` utiliza trajetória; `--choose` depende de seleção de domínios pelo usuário. Nos modos móveis, diagnósticos adicionais usam **850 hPa**. | Compartilham módulos de termos energéticos, mas os fluxos não são idênticos (incluindo tratamento de tendência/aquecimento e definição de domínios). `--fixed` isolado não valida integralmente os outros dois. |
| Artefatos auxiliares | Trajetória com `time`, `Lat`, `Lon` no formato esperado pelo validador; configuração de `namelist` e limites de caixa pertinentes. | Conferir nomes, formatos e datas pela implementação, pois exemplos documentais podem divergir. |

**Rastreabilidade primária (base inspecionada):** [`lorenzcycletoolkit.py`](../../../lorenzcycletoolkit.py), [`src/utils/preprocessing.py`](../../../src/utils/preprocessing.py), [`src/utils/box_data.py`](../../../src/utils/box_data.py), [`src/utils/vertical_harmonization.py`](../../../src/utils/vertical_harmonization.py), [`src/frameworks/lec_fixed_framework.py`](../../../src/frameworks/lec_fixed_framework.py), [`src/frameworks/lec_moving_framework.py`](../../../src/frameworks/lec_moving_framework.py), [`src/utils/validation.py`](../../../src/utils/validation.py) e [`src/utils/tools.py`](../../../src/utils/tools.py). Todos se referem ao SHA-base acima, e não a uma suposta versão numericamente validada.

**Observação a conferir em etapa executiva:** `prepare_data` escolhe `inputs/namelist_ERA5-cdsapi` com `--cdsapi`, enquanto `run_lec_analysis` volta a ler `inputs/namelist`; isso indica possível divergência de configuração, **não uma falha demonstrada por execução**. Não corrigir nesta frente.

## 3. Matriz de fontes, produtos e elegibilidade

**Verificação de catálogos em 2026-10-25.** Não foram realizadas autenticação, chamadas API, download, medições de volume, execução científica, comparação de campos ou auditoria independente das coleções remotas. Coberturas finais e termos comerciais podem mudar e requerem consulta atualizada na D.

| Produto / método | Variáveis e níveis | Tempo, espaço, formato / recorte | Acesso e classificação documental |
|---|---|---|---|
| **NetCDF local** (`infile`) | Cinco campos nas definições/unidades exigidas; quantidade de níveis dependente do arquivo. | Grade regular, horários e região controlados pelo usuário; NetCDF aberto por `xarray`. | **Compatível documentalmente**, somente se o arquivo individual cumprir o contrato e o `namelist`; baseline operacional. |
| **CDS API — `reanalysis-era5-pressure-levels`** (aquisição existente) | Cinco campos; seleção de até 37 níveis isobáricos. | ERA5 horário, recorte por variáveis/área/instantes/níveis; opção NetCDF; código atual solicita dias separadamente e concatena. | **Compatível sob condições explícitas**: credenciais, disponibilidade da API, seleção correta de níveis e coerência de `namelist`, arquivo e horários; não pressupor fluxo validado nesta etapa. |
| **ECMWF CDS ARCO Zarr** (subconjunto do ERA5 pressão) | Cinco campos presentes na seleção de variáveis; **13 níveis**: 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 50 hPa. | **Seis-horário** segundo o guia do produto ARCO; 0,25°; Zarr com organização espacial e temporal em chunks, recortes via leitor apropriado. | **Compatível sob condições explícitas**: chave CDS, seleção Zarr, preparação de NetCDF ou futuro leitor, horários/níveis aceitos; não equivale a 37 níveis/3 h. |
| **CDS — `reanalysis-era5-pressure-levels-timeseries`** (serviço associado, distinto do store ARCO) | Cinco campos de interesse em **13 níveis**. | 0,25°, **6 h** (00,06,12,18 UTC), NetCDF/CSV; serviço orientado a séries de pontos/áreas pequenas, não presumir geração ideal de caixas 2D. | **Compatível sob condições explícitas** apenas caso a saída realmente forme uma grade/domínio suficiente; **não** tomar download de ponto como entrada 3D pronta do LECTK. Não é a mesma interface de acesso Zarr. |
| **WEkEO — `EO:ECMWF:DAT:REANALYSIS_ERA5_PRESSURE_LEVELS`** | Catálogo lista geopotencial, temperatura, `u`, `v` e `omega`, correspondendo ao ERA5 em níveis de pressão. | Produto horário, recorte/download via HDA API; confirmar níveis efetivamente selecionáveis, contrato do arquivo exportado, disponibilidade e dimensões. | **Compatível sob condições explícitas**: autenticação, requisição de subdomínio, NetCDF compatível/normalização, verificação de limites e cotas; não foi demonstrada velocidade. |
| **DestinE Earth Data Hub — `era5-pressure-levels-v0.zarr`** | **19 níveis**: 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10, 5, 1 hPa; variáveis `t`, `z`, `u`, `v`, `w`. | Horário, 0,25° global, Zarr v3; recorte espacial e temporal por Xarray, organização de chunks documentada. | **Compatível sob condições explícitas**: chave EDH, biblioteca Zarr v3, conversão/seleção para `infile`, gestão de chunks, cotas e credenciais. Cota documentada de **500.000 requisições por mês**, não taxa de download verificada. |
| **Google ARCO ERA5 — `full_37-1h-0p25deg-chunk-1.zarr-v3`** | 37 níveis, `temperature`, `u_component_of_wind`, `v_component_of_wind`, `vertical_velocity`; o inventário examinado lista `geopotential_at_surface`, **não geopotencial em níveis de pressão**. | Horário, 0,25°, Zarr público na infraestrutura GCS; regradeado. | **Incompatível com requisitos essenciais como fonte isolada diretamente utilizável**, pois falta comprovação do campo `z` 3D no produto específico; derivação/complementação constitui outra proposta, não avaliada. |
| **Google Earth Engine — `ECMWF/ERA5/HOURLY_PRESSURE_LEVELS`** | 37 níveis, bandas por variável e nível; o inventário consultado **não contém geopotencial de pressão**. | Horário, grade em bandas por nível, exportação e reestruturação seriam necessárias; projeto Earth Engine e cotas aplicáveis. | **Incompatível como fonte isolada diretamente utilizável** no produto consultado; eventual combinação de coleções ou derivação é trabalho novo, não autorizado. |

**Observação transversal:** 13 ou 19 níveis não tornam automaticamente a fonte inelegível. Configurações verticais reduzidas podem ser cientificamente suficientes para aplicações escolhidas pelo usuário, mas quantificar viés/convergência não faz parte da Frente C. Fontes que já incluam 850 hPa conservam pelo menos esse pré-requisito dos diagnósticos móveis atuais. A existência de cinco variáveis no catálogo não demonstra que todos os passos de adaptação funcionem.

### 3.1. Referências públicas oficiais

- **CDS ERA5 níveis de pressão (produto completo e ARCO):** https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels?tab=analysis_ready_data
- **ECMWF, guia do ARCO em pressão (seis-horário, chunking):** https://confluence.ecmwf.int/spaces/CKB/pages/639215289/ERA5+Analysis+Ready+Cloud+Optimised+ARCO+data+on+pressure+levels+from+1940+to+present+Product+User+Guide+PUG
- **CDS séries temporais de pressão:** https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels-timeseries?tab=overview
- **WEkEO, catálogo do produto:** https://data.wekeo.copernicus.eu/data?dataset=EO%3AECMWF%3ADAT%3AREANALYSIS_ERA5_PRESSURE_LEVELS&view=layers
- **WEkEO, acesso harmonizado (HDA):** https://wekeo.copernicus.eu/scaling-up-data-access
- **DestinE EDH, produto ERA5 de pressão:** https://earthdatahub.destine.eu/collections/era5/datasets/era5-pressure-levels
- **DestinE EDH, autenticação e cota:** https://earthdatahub.destine.eu/getting-started
- **Google ARCO ERA5, inventário e armazenamento:** https://github.com/google-research/arco-era5/blob/main/README.md
- **Google Earth Engine, coleção de pressão:** https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_HOURLY_PRESSURE_LEVELS

**Questões não resolvidas/documentação insuficiente:** desempenho real e custos por unidade de dado; limites atualizados de requisição, paralelismo e faturamento, especialmente WEkEO/Google; integridade efetiva de exportações; compatibilidade de cabeçalhos/unidades dos arquivos retornados; ausências pontuais de dados; equivalência com arquivos CDS originais. Não adotar números de desempenho/ganho como fatos.

## 4. Escopo deliberadamente separado: níveis verticais e histórico publicado

O responsável informou manter **resultados de milhares de ciclones**, obtidos originalmente do CDS/ERA5 e processados pelo LECTK, em servidor remoto, incluindo resultados publicados. **Esses arquivos não foram acessados ou inventariados nesta frente.** Podem servir de referências científicas históricas presumidamente confiáveis para reprodução de cenários semelhantes, sem presumir validação física independente universal, e sem exigir revalidação indiscriminada de resultados publicados.

A análise de possíveis vieses médios, dispersão das diferenças, sensibilidade e convergência por número/distribuição de níveis é uma **investigação científica separada** (backlog K). Uma eventual documentação futura desses resultados terá valor próprio, mas **não condiciona a classificação documental da Frente C** e não constitui benchmark já executado.

As **49 referências numéricas ERA5** em `tests/data/reference/` preservam regressões dos seus cenários documentados, não equivalência entre provedores. A ausência de referência numérica NCEP-R2 aceita permanece registrada pela Frente B; a necessidade disso será reavaliada antes do desenho final de D, sem reabertura automática de B.

## 5. Encaminhamento para benchmarks controlados (Frente D, ainda não aberta)

**Prioridade de investigação, sem compromisso de implementação:** (i) baseline NetCDF local e CDS atual; (ii) DestinE EDH 19 níveis, CDS ARCO 13 níveis, WEkEO; (iii) alternativas condicionais somente se novo produto demonstrar todos os campos, inclusive geopotencial 3D.

**Desenho mínimo proposto e não executado:**

1. **Pré-verificação sem download volumoso:** autenticação e metadados; cinco variáveis corretas; coordenadas, unidades, níveis, tempo, cobertura e formato; definir prioriamente medidas e tolerâncias relevantes.
2. **Aquisição mínima e `--fixed`:** período/caixa reduzidos mas suficientes para gradientes e tendências; preservar os mesmos níveis/horários/domínio quando a pergunta for equivalência de acesso, e registrar quando isso não for possível.
3. **Cobertura distinta dos modos:** breve caso `--track` que exercite 850 hPa, seleção móvel e tratamento temporal; `--choose` apenas na medida dos comportamentos não cobertos, sem assumir identidade integral com `--fixed`.
4. **Métricas em blocos:** tempo de fila/transferência, bytes transferidos, tempo de conversão/preparação, memória, cálculo, escrita, custo/cotas e integridade; reportar reprodutibilidade numérica apenas com controles e tolerâncias pré-declarados.
5. **Acervo histórico:** após autorização específica, inventariar em modo somente leitura o servidor e selecionar poucos casos representativos; não copiar dezenas de milhares de ciclones. Comparações contra resultados publicados devem preservar produto, configuração, domínio e modo, sem igualar `--track` histórico a um novo `--fixed` arbitrário.

**Fora de escopo:** nova arquitetura E, adaptadores, downloads/benchmarks agora, mutação de dados científicos, estudo de convergência K, alteração das equações, credenciais/serviços cobrados, `main`, tags, releases e publicação.

## 6. Estado, riscos e portões para homologação

**Avaliação técnica documental:** **PASS DOCUMENTAL** nos critérios pré-declarados, com limites explícitos. **Homologação e integração verificadas em 25/10/2026:** o responsável aprovou o conteúdo e a integração documental; o PR [#36](https://github.com/daniloceano/LorenzCycleToolkit/pull/36) alterou somente este relatório e dois índices; GitHub Actions (push e PR) e CircleCI terminaram com sucesso; `develop` passou ao merge commit `abaab7add04c21c694a4aa90da38c7d1465e5add`; os três arquivos integrados foram relidos e seus blobs coincidiram com os textos aprovados. **Conclusão administrativa: Frente C encerrada por PASS DOCUMENTAL**, não por execução de comparações de provedores. A Frente D não foi aberta automaticamente.

**Riscos remanescentes:** documentação de produto pode mudar; versões de bibliotecas e acesso/quotas devem ser verificadas antes de chamadas reais; a classificação de produto isolado não cobre arquiteturas combinadas; diferenças de níveis/tempo/grade afetam interpretação científica e não foram mensuradas aqui. No fechamento revisar retenção mínima de artefatos operacionais, sem exclusões automáticas.

## 7. Retenção e transferência

Não houve coleta remota de ERA5, acesso ao servidor de ciclones ou geração de dados científicos nesta frente. Este relatório, as fontes citadas, a fila, o histórico do PR #36 e os hashes do repositório são as evidências permanentes; não há artefato operacional exclusivo identificado para retenção nesta atuação documental. Nenhuma exclusão ou limpeza de arquivos locais/remotos foi autorizada ou realizada. Aplicar a revisão de retenção aos eventuais materiais de agentes ou ambientes externos quando inventariados, sem presumir que os arquivos do servidor foram preservados pelo Git.

**Transferência à conversa-mãe:** apresentar o handoff automático com resultado, SHA de `develop`, PR #36, ressalvas, backlog K e fila; a mãe deve reler o registro integrado e não abrir D sem decisão contextual expressa.
