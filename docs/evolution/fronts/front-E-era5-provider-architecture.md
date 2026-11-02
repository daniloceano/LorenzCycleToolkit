# Frente E — decisão de arquitetura para aquisição ERA5

**Data das aprovações:** 2026-11-01. **Responsável:** mantenedor do LorenzCycleToolkit. **Estado da decisão:** desenho arquitetural e contrato científico aprovados em conversa; **encerramento operacional não declarado neste registro**. A integração e a retenção exigem comprovação própria. **Base documental e de código inspecionada:** `develop` `835478e101c20d7ed668ac0c692b8bcd626709e8`. Este registro especifica a arquitetura futura; não declara que ela já foi implementada.

## 1. Objetivo e decisão

Preservar a entrada independente por arquivo NetCDF local (`infile`), a CDS API opcional e o adaptador NCAR GDEX experimental, ao mesmo tempo que se distingue a **origem dos dados** do **método de acesso**. Aprova-se a **evolução mínima com contrato explícito do artefato meteorológico e modularização seletiva** (alternativa A+, combinação proporcional de A com requisitos de contrato da B), sem introduzir agora uma arquitetura modular ampla, uma interface universal executável ou um framework novo.

**Separação de conceitos:**

- **Fonte/origem:** NetCDF fornecido localmente, CDS API/Copernicus ou NCAR GDEX. Não tratar `GDEX` como um quarto método concorrente.
- **Método de aquisição na origem NCAR GDEX:** `ncss` (THREDDS NCSS), `subset` (Get a Subset via API assíncrona) ou `kerchunk` (referências e requisições HTTP por faixa). São métodos alternativos de acesso, não produtos meteorológicos cientificamente distintos por definição.
- **Artefato entregue:** arquivo meteorológico NetCDF e, para dados adquiridos programaticamente, evidência de proveniência e integridade correspondente. Um arquivo de entrada ser estruturalmente utilizável não significa equivalência numérica ou física com arquivo de outra fonte.
- **Etapas lógicas:** aquisição → normalização estritamente necessária e declarada → validação de adequação e integridade → preparação do domínio/níveis/tempos → cálculos energéticos. A separação é de responsabilidades; não obriga novos módulos, classes ou conversões nesta frente.

## 2. Interface de linha de comando e compatibilidade

| Operação | Estado atual na base inspecionada | Decisão arquitetural |
|---|---|---|
| Usar `infile` local | Entrada principal, argumento posicional obrigatório | Preservar de modo independente, sem dependência de credenciais ou sidecar obrigatório; manter interface e semântica existentes. |
| Obter via `--cdsapi` | Opção existente, mutuamente exclusiva de `--ncar` | Preservar comportamento existente; continuar opcional. |
| Obter via `--ncar` | Adaptador experimental usando NCSS | Preservar; `--ncar` sem seleção de método significa `ncss` por padrão. |
| Escolher método NCAR na CLI principal | **Não disponível na base atual** | Especificar para implementação futura `--ncar --ncar-method {ncss,subset,kerchunk}`; seletor válido apenas com `--ncar`. Não criar três flags independentes. |
| Selecionar método no cliente independente | `python -m src.utils.ncar_gdex --method {ncss,kerchunk,subset}` | Preservar o caráter experimental, sem declarar os métodos integrados ao comando principal. |

As formas da coluna de decisão são **propostas aprovadas de interface futura**, não comandos operacionais disponíveis para os três métodos na CLI principal. Preservar `--fixed`, `--track` e `--choose`, as configurações atuais e a exclusividade CDS/NCAR; não alterar silenciosamente defaults de resultados, aquisição, precisão ou exclusão de arquivos. Nenhum novo método será considerado suportado oficialmente somente por existir código experimental.

## 3. Contrato científico de entrada — requisitos e limites

A aplicação do contrato deve ser proporcional à origem e à etapa: para arquivo local, verificar adequação ao processamento sem exigir origem conhecida; para aquisições programáticas, documentar também proveniência e integridade. Não declarar que todos os itens abaixo já têm verificadores uniformes implementados.

| Dimensão | Requisito/decisão |
|---|---|
| Variáveis | Temperatura, geopotencial **ou altura geopotencial tratada conforme convenção explícita**, ventos zonal e meridional e velocidade vertical em coordenada de pressão. Identificação dos nomes pelo `namelist`; preservar significado físico e unidades. |
| Unidades | Registrar/validar unidades efetivas e conversões feitas ou presumidas, incluindo pressão hPa/Pa e geopotencial/altura. Unidade ausente não prova unidade correta; não mascarar ambiguidades. |
| Eixos e dimensões | Coordenadas coerentes de tempo, níveis isobáricos, latitude e longitude, com estrutura e ordenação adequadas ao modo solicitado. Não pressupor compatibilidade direta de produto em níveis de modelo. |
| Cobertura | Domínio, tempos e níveis suficientes para cada cálculo, gradientes, bordas e modos fixo/móvel. `850 hPa` é relevante para diagnósticos móveis existentes; níveis completos/comuns dependem do caso. |
| Resolução | **37 níveis e 3 h são referência histórica**, não requisito universal. Diferentes grades, níveis e passos temporais não são automaticamente equivalentes do ponto de vista científico. |
| Precisão | Registrar representação/dtype (`float32`, `float64` etc.) e operações que a alterem. Não promover conversão implícita para igualar fontes; conservar o padrão atual do NCAR (`float32`) e a exceção `float64` do cliente independente até decisão executiva específica. |
| Integridade | Verificar completude e coerência do artefato antes do cálculo conforme o fluxo; verificar hashes e metadados de reaproveitamento quando os mecanismos existirem. Distinguir validade de arquivo, elegibilidade para cálculo, reprodutibilidade e equivalência entre fontes. |
| Proveniência | Para aquisição automatizada: registrar origem, produto/coleção, método, seleção de variáveis/níveis/domínio/tempos, parâmetros relevantes, tamanho/hash da saída quando praticável, versão/configuração, dtype e informações necessárias para auditoria. Para `infile` local: não tornar sidecar ou credenciais pré-requisito; manter metadados observáveis e hash onde já registrados. Não inventar origem externa de arquivo fornecido pelo usuário. |
| Falhas e retomada | Tratar autenticação, disponibilidade, cota, timeout, respostas parciais, duplicação, interrupção e reuso de forma explícita, segura e rastreável. Evitar que saída incompleta seja confundida com entrada válida. |
| Credenciais | Permanecer fora dos arquivos versionados e dos relatórios públicos. Não registrar tokens, cookies ou caminhos privados; evitar requisições remotas implícitas no modo de entrada local. |
| Núcleo científico | Preservar pré-processamento e cálculos do ciclo de Lorenz, equações, convenções e referências. Nenhuma correção ou normalização de dados será introduzida sem caso concreto, critérios e testes autorizados. |

A separação de responsabilidades acima é um **contrato documental**, não autorização para refatorar `src/utils/preprocessing.py`, `src/utils/tools.py`, `src/utils/ncar_gdex.py` ou os frameworks energéticos. Se um adaptador futuro precisar normalizar campos, a transformação deve ser explícita, verificável e não confundir compatibilidade estrutural com equivalência científica.

## 4. Alternativas consideradas

| Alternativa | Benefícios | Custos/riscos | Veredito |
|---|---|---|---|
| **A — evolução mínima** | Preserva a CLI, mantém escopo e risco reduzidos, reaproveita código testado. | Sem contrato explícito, diferenças de validação/proveniência tenderiam a persistir. | **Aceita com contrato e separação de responsabilidades (A+).** |
| **B — interface comum executável de aquisição** | Pode facilitar adaptação de múltiplos provedores e homogeneizar verificações. | Introduz abstrações, dependências e testes para rotas experimentais que ainda não justificam promoção; pode forçar equivalência artificial. | **Adiada**; reconsiderar se houver necessidade comprovada de múltiplas integrações ou duplicação material. |
| **C — modularização ampla de provedores/produtos/armazenamento/execução** | Poderia acomodar armazenamento remoto/Zarr e execução junto aos dados. | Alto custo de manutenção e risco de afetar componentes científicos sem evidência de necessidade. | **Não aprovada** nesta frente; G/H e demais possibilidades permanecem no backlog. |

Não há provedor vencedor previamente designado. O fato de uma rota ter funcionado em uma campanha limitada não autoriza generalizações sobre desempenho, disponibilidade ou equivalência física.

## 5. Estado comprovado e evidência científica

- Na base `835478e`, a CLI principal oferece `--cdsapi` e `--ncar` como opções mutuamente exclusivas; `--ncar` utiliza NCSS. `src/utils/ncar_gdex.py` oferece `--method ncss|kerchunk|subset` apenas no cliente independente; `src/utils/gdex_sources.py` já contém uma abstração interna `_Source` dos métodos NCAR. **Não duplicar por padrão essa abstração na arquitetura principal.**
- A Frente D integrou NCAR NCSS como **experimental** (PR #38) e preservou `infile` e CDS; não demonstrou superioridade sobre CDS nem validação física independente.
- A D-Bench encerrada em 01/11/2026 apresentou comparação primária CDS × NCAR **inconclusiva** (um par válido em dez). Na Etapa B: Get a Subset 9/10, kerchunk 10/10 e NCSS 4/10, sob as limitações de uma campanha; esses números **não** representam comparativo universal de velocidade nem verificação de produção.
- Arquivos NCAR e CDS tiveram diferenças pequenas de temperatura, suficientes para causar desvios energéticos além das tolerâncias pré-registradas em **5 de 9** casos com arquivo CDS. O diagnóstico posterior apoiou a explicação por temperatura nos casos investigados; **não identificou defeito no LECTK** no escopo testado. A origem exata da conversão histórica não foi estabelecida. Não afirmar equivalência irrestrita.
- A Frente C já identificava como **questão a investigar, não defeito confirmado**, o uso de `inputs/namelist_ERA5-cdsapi` em `prepare_data` versus `inputs/namelist` no cálculo. Também merecem atenção futura as diferenças de política de reutilização entre CDS e NCAR. Estas observações não autorizam refatoração corretiva automática.
- Não há exigência de novo baseline NCEP-R2, nem estudo K de sensibilidade vertical, como condição da decisão arquitetural. Preservar os testes e resultados já versionados.

**Fontes permanentes (consultadas na base indicada):** [`front-C-era5-documentary-assessment.md`](front-C-era5-documentary-assessment.md), [`front-D-era5-ncar.md`](front-D-era5-ncar.md), [`front-D-bench-cds-ncar.md`](front-D-bench-cds-ncar.md), [`encerramento D-Bench`](../evidence/records/2026-11-01-dbench-closure.md), [`diagnóstico pós-D-Bench`](../evidence/records/2026-11-01-dbench-postclosure-energy-diagnostic.md), [`documentação dos métodos`](../../source/remote_data_acquisition.rst), [`CLI`](../../../lorenzcycletoolkit.py), [`pré-processamento`](../../../src/utils/preprocessing.py), [`aquisição NCAR`](../../../src/utils/ncar_gdex.py) e [`métodos NCAR`](../../../src/utils/gdex_sources.py).

## 6. Encaminhamento condicional para a Frente F

A **Frente F não está aberta**. Se aprovada em separado, recomenda-se avaliar um **piloto delimitado de Get a Subset** como primeiro candidato, porque operou em nove dos dez cenários da Etapa B, **sem classificá-lo como superior ou já pronto para a CLI**. Kerchunk permanece candidato experimental; NCSS conserva seu papel experimental atual. O piloto deverá decidir, antes de implementar, se é necessária qualquer interface executável comum adicional ou apenas ligação pontual ao contrato de saída.

Critérios a pré-registrar, proporcionais ao escopo do piloto:

1. Compatibilidade da CLI legada: `infile`, `--cdsapi` e `--ncar` isolado; rejeição clara de `--ncar-method` sem `--ncar` se a flag for implementada.
2. Testes offline e controles negativos por método: autenticação/cotas, indisponibilidade, timeout, resultados incompletos, retomada, reutilização, hashes, concorrência e ausência de segredo nos logs.
3. Contrato científico e rastreabilidade: cinco campos e unidades, tempos/níveis/domínio, tratamento explícito de dtype, configurações e proveniência.
4. Casos de execução dos modos pertinentes, com comparação de campos e diagnósticos sob tolerâncias **pré-definidas**, sem pressupor equivalência a CDS ou exigir 37 níveis/3 h universalmente.
5. Custos de manutenção, recursos, disponibilidade e portabilidade com evidência limitada e devidamente contextualizada; nenhuma promessa de ganho de desempenho.

É vedado abrir F automaticamente, alterar o núcleo científico, promover métodos a suportados sem validação, reexecutar a D-Bench, baixar grandes volumes, alterar referências, `main`, tags ou publicações por efeito desta aprovação.

## 7. Estado da Frente E, integração e retenção

**Aprovações substantivas recebidas:** (i) organização `--ncar` como origem e seleção única de método `--ncar-method`, com NCSS padrão para compatibilidade; (ii) contrato científico e limites descritos acima. **Sem implementação de código autorizada ou realizada.**

**Portões de integração:** registrar somente memória/documentação curada por branch nova baseada na `develop`, PR para `develop`, checks GitHub Actions e CircleCI verdes, merge normal, confirmação do SHA final e releitura dos registros. Não alterar o texto protegido de `docs/development_workflow.md`, nem tocar `main`/tags/releases. Uma integração remota não substitui automaticamente sincronização do checkout local.

**Retenção GOV-RET-01 (Marco A): requer auditoria local antes do encerramento.** O material produzido nesta conversa é documental, sem aquisição de ERA5, sem testes executáveis nem downloads; isso **não prova** que não existam arquivos ignorados, ocultos ou externos da própria Frente E no computador do responsável. Antes de encerrar, inventariar esses espaços em modo leitura, classificar cada item, verificar exclusividade e cópias; descarte somente de comprovados dispensáveis dentro da frente com manifesto de caminho, contagem, bytes e hashes, verificação imediata pré-exclusão e relatório durável. Preservar dados D-Bench, `LEC_Results/`, worktrees, branches, stashes e materiais da Frente D. Reabrir revisão das retenções na próxima versão/tag estável.

**Estado de encerramento neste registro:** não declarado; a comprovação da integração documental e da retenção pertence a registro de encerramento separado. Não declarar `ENCERRADA — RETENÇÃO VALIDADA` sem ambos. Após integração e fechamento, fornecer handoff para a conversa-mãe; F somente mediante abertura expressa.
