# Fila de frentes — replanejamento de 07/11/2026

## Ordem vigente (decisão do responsável, 07/11/2026)

Substitui a fila antiga a partir da Frente F ([D-018](decisions.md); [registro](evidence/records/2026-11-07-replanning-decisions.md)). Fechar uma etapa não abre automaticamente a seguinte.

| Etapa | Objeto | Dependência | Estado |
|---|---|---|---|
| 1 | Registro das decisões de replanejamento (Frente R) | — | Concluída — registro integrado em `develop` (PR #57); encerramento formal da Frente R não documentado |
| 2a | Validação externa: Reg1–Reg3 de Dias Pinto & da Rocha (2011), NCEP/NCAR R1, caixas fixas e fases do artigo | 1 | Em andamento — comparação executada e documentada (D-021, PR #58); FAIL numérico global (53/89, 59,6%); decisão de encerramento pendente |
| 2a.1 | Conferência interna dos termos de fronteira e do fechamento do resíduo (ligada a 2a) | 2a | Concluída em 09/11/2026: implementação conferida (T1–T5 PASS), sem erro de fórmula; ver [registro](evidence/records/2026-11-09-subfrente-2a1-boundary-terms.md) |
| 2b | Benchmark unificado de fontes ERA5 (paralelo a 2a) | 1 | Etapa 0 concluída; pré-registro aprovado; fase 1 em pré-voo — execução no servidor de processamento (emenda a D-017, 09/11/2026); ver [registro](evidence/records/2026-11-09-front-2b-stage0-preregistration.md) |
| 3 | Níveis verticais: medição simples, 10 casos D-Bench, caixa fixa | 2b | Não aberta |
| 4 | Inclusão das fontes escolhidas e documentação de cada uma (absorve F) | 2b | Não aberta |
| 5 | Arrumação final (pendências de D-019) | — | Não aberta |
| 6 | Nova versão oficial única (Marco B de GOV-RET-01) | 1–5 | Exige autorização específica |

G–J seguem como backlog; K continua fora do escopo.

**Backlog (não aberto):**

- Investigar CA e CE 19–25% menores no Reg3 da 2a (não explicáveis pelos termos de fronteira); candidato inicial: σ e SIGMA_FLOOR para CA. Não aberto.

---

# Histórico — fila de 18/10 e atualização de 25/10

*Preservada como escrita; a partir de F, superada pela ordem acima (07/11/2026).*

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
| **ERA5/Zarr e níveis verticais — LZV (complementar, pré-F)** | Um caso real da D-Bench (`19791150`, 27 horários de 6 h): sensibilidade de Ca/Ce/Ck a 32, 29, 17, 14 e 13 níveis efetivos; CDS ARCO Zarr e NCAR NCSS; DestinE bloqueado; correção do alinhamento temporal track × dados | **ENCERRADA EM 06/11/2026 — RETENÇÃO VALIDADA.** Executada em 02–03/11/2026 e integrada em `develop` pelo PR #55 (merge `24ffe5b30dc11ac78c2c0c7c8732733657dff9c5`, com a correção `2101de2` e o ajuste `7e170ac`); GitHub Actions e CircleCI aprovados e releitura concluída. Interpretação aprovada pelo responsável: as diferenças entre conjuntos verticais não alteram substantivamente a interpretação física deste caso (Ca −5 % na média). Sem vantagem geral de velocidade do Zarr; não implementar rota Zarr agora. Correção do alinhamento temporal integrada. [Registro](evidence/records/2026-11-02-lzv-era5-zarr-vertical-levels.md). **Não abre F nem K**; não reabre D, D-Bench ou E. Retenção: **ENCERRADA — RETENÇÃO VALIDADA** em 06/11/2026, após a indisponibilidade do servidor de 03/11 (export redundante descartado; espelho ARCO e NetCDF retidos até o Marco B) ([encerramento e retenção GOV-RET-01](evidence/records/2026-11-06-lzv-closure-retention.md)). |

A proposta de adendo do lote 15 listou **BL-01 a BL-11**, com IDs expressamente **provisórios e não aprovados**. Abrange validação real independente, NCEP-R2, Dz/De, mensagens de domínio inválido, empacotamento, versão, promoções de branch, releases/Pages, documentação, proteção de tags e governança. Não promover esses IDs à fila oficial.

O encerramento A não validou fisicamente o toolkit. ERA5 37 níveis/3 h, NCEP-R2 e comparação física independente permanecem não demonstrados. **B foi encerrada por triagem, não por produção de referência NCEP-R2.** Antes de desenhar benchmarks D, reavaliar se a referência numérica NCEP-R2 será necessária; manter a pendência explícita enquanto não houver protocolo e evidência. **A investigação documental C não exige estudo de sensibilidade a níveis (K) nem abre D automaticamente.**
