# Fila de frentes — decisão de 18/10 e atualização de 25/10

A ordem inicial aprovada foi **A → B condicional → C → D → E → F condicional** (D-010). G–J são backlog; K foi acrescentada como hipótese não testada. Fechar uma frente não abre automaticamente a próxima.

| Frente | Objeto | Estado/documento |
|---|---|---|
| A | Auditoria e linha de base científico-computacional | **ENCERRADA 23/10: PASS COM RESSALVAS ACEITAS**; [termo](fronts/front-A-closure.md) |
| B | Lacunas de dados de referência/regressão **do próprio LECTK** | **ENCERRADA 24/10 — triagem concluída; execução de novas referências dispensada por ora**; NCEP-R2 permanece pendente. [Termo](fronts/front-B-closure.md) |
| C | Levantamento de provedores ERA5 e requisitos científicos | **ENCERRADA EM 25/10/2026 — PASS DOCUMENTAL**, aprovação de conteúdo e PR #36 integrados; CI verde e readback no SHA `abaab7add`. [Avaliação homologada e registro de encerramento](fronts/front-C-era5-documentary-assessment.md). Sem benchmarks nem integração de provedores. |
| D | Benchmarks controlados | **Encerrada em 28/10/2026: PASS COM RESSALVAS ACEITAS.** Executada como avaliação de rotas ERA5: adaptador NCAR GDEX **experimental** em `develop` (PR #38, merge `dec13fe`). **Benchmarks controlados e comparação NCAR × CDS não concluídos** (P2-B0 bloqueada); sem superioridade de velocidade nem validação física. Retenção validada: limpeza parcial (837/838 caminhos, ≈ 11 GiB) e retenções temporárias até a próxima versão/tag. [Registro](fronts/front-D-era5-ncar.md) |
| E | Decisão de arquitetura de provedores | **ENCERRADA EM 01/11/2026 — RETENÇÃO VALIDADA.** Decisão arquitetural aprovada e integrada (PR #53, merge `03ed907`): evolução mínima com contrato científico explícito de entrada (A+), `--ncar` como origem, NCSS como padrão e seletor único `--ncar-method` apenas especificado para o futuro; `infile` e CDS preservados. **Sem implementação de código.** [Decisão](fronts/front-E-era5-provider-architecture.md); [encerramento e retenção](evidence/records/2026-11-01-front-E-closure.md) |
| F | Piloto de provedor | Condicional, **não aberto**; abertura só por decisão expressa. Primeiro candidato recomendado, se aprovado: piloto delimitado de Get a Subset ([E §6](fronts/front-E-era5-provider-architecture.md#6-encaminhamento-condicional-para-a-frente-f)) |
| G | Cache e Zarr regional | Backlog |
| H | Execução junto aos dados | Backlog |
| I | Interface visual simples | Backlog, sem tecnologia aprovada |
| J | CI, empacotamento, CLI e docs | Backlog, não manutenção genérica autorizada |
| K | Sensibilidade/convergência versus níveis verticais | Hipótese, **não medida**; backlog separado da elegibilidade documental C |

## Investigação complementar fora da ordem principal

| Investigação | Objeto | Estado/documento |
|---|---|---|
| **D-Bench (complementar, paralela)** | Dez trajetórias reais fornecidas (1979–1980, 438 posições/3 h, modo móvel `--track`) em benchmark CDS × NCAR: aquisição, preparação, LEC, recursos, falhas e latência sob carga | **ENCERRADA EM 01/11/2026 — sem PASS da campanha pré-registrada.** Primária CDS × NCAR INCONCLUSIVA (1 par válido em 10); recuperação NCAR FAIL no diagnóstico; Etapa B (métodos NCAR, sem CDS) executada. Nenhuma superioridade de velocidade NCAR × CDS demonstrada; energética NCAR não equivalente à CDS em 5 de 9 casos. Integração PRs #45, #47–#49. [Encerramento](evidence/records/2026-11-01-dbench-closure.md); [plano e portões](fronts/front-D-bench-cds-ncar.md). **Não reabre D e não altera a fila principal.** Retenção: **ENCERRADA — RETENÇÃO VALIDADA** em 01/11/2026 ([limpeza GOV-RET-01](evidence/records/2026-11-01-dbench-retention-cleanup.md)). |
| **Diagnóstico energético pós-D-Bench (complementar)** | Esclarecimento offline das quatro divergências energéticas remanescentes (19791008, 19800256, 19790166, 19791150) | **CONCLUÍDO, conclusões aceitas pelo responsável em 01/11/2026.** Diferenças de temperatura nos arquivos CDS × NCAR explicam os desvios energéticos, amplificados em estabilidade estática, tendências e resíduos; nenhum defeito do LECTK identificado no escopo testado; correções dispensadas. A origem exata da conversão histórica não foi demonstrada. [Registro pós-encerramento](evidence/records/2026-11-01-dbench-postclosure-energy-diagnostic.md). **Não reabre a D-Bench nem a Frente E.** Retenção: **ENCERRADA — RETENÇÃO VALIDADA** ([limpeza](evidence/records/2026-11-01-dbench-retention-cleanup.md)). |

A proposta de adendo do lote 15 listou **BL-01 a BL-11**, com IDs expressamente **provisórios e não aprovados**. Abrange validação real independente, NCEP-R2, Dz/De, mensagens de domínio inválido, empacotamento, versão, promoções de branch, releases/Pages, documentação, proteção de tags e governança. Não promover esses IDs à fila oficial.

O encerramento A não validou fisicamente o toolkit. ERA5 37 níveis/3 h, NCEP-R2 e comparação física independente permanecem não demonstrados. **B foi encerrada por triagem, não por produção de referência NCEP-R2.** Antes de desenhar benchmarks D, reavaliar se a referência numérica NCEP-R2 será necessária; manter a pendência explícita enquanto não houver protocolo e evidência. **A investigação documental C não exige estudo de sensibilidade a níveis (K) nem abre D automaticamente.**
