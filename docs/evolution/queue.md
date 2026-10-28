# Fila de frentes — decisão de 18/10 e atualização de 25/10

A ordem inicial aprovada foi **A → B condicional → C → D → E → F condicional** (D-010). G–J são backlog; K foi acrescentada como hipótese não testada. Fechar uma frente não abre automaticamente a próxima.

| Frente | Objeto | Estado/documento |
|---|---|---|
| A | Auditoria e linha de base científico-computacional | **ENCERRADA 23/10: PASS COM RESSALVAS ACEITAS**; [termo](fronts/front-A-closure.md) |
| B | Lacunas de dados de referência/regressão **do próprio LECTK** | **ENCERRADA 24/10 — triagem concluída; execução de novas referências dispensada por ora**; NCEP-R2 permanece pendente. [Termo](fronts/front-B-closure.md) |
| C | Levantamento de provedores ERA5 e requisitos científicos | **ENCERRADA EM 25/10/2026 — PASS DOCUMENTAL**, aprovação de conteúdo e PR #36 integrados; CI verde e readback no SHA `abaab7add`. [Avaliação homologada e registro de encerramento](fronts/front-C-era5-documentary-assessment.md). Sem benchmarks nem integração de provedores. |
| D | Benchmarks controlados | **Encerrada em 28/10/2026: PASS COM RESSALVAS ACEITAS.** Executada como avaliação de rotas ERA5: adaptador NCAR GDEX **experimental** em `develop` (PR #38, merge `dec13fe`). **Benchmarks controlados e comparação NCAR × CDS não concluídos** (P2-B0 bloqueada); sem superioridade de velocidade nem validação física. Retenção validada: limpeza parcial (837/838 caminhos, ≈ 11 GiB) e retenções temporárias até a próxima versão/tag. [Registro](fronts/front-D-era5-ncar.md) |
| E | Decisão de arquitetura de provedores | Não aberta; depende de D |
| F | Piloto de provedor | Condicional, não aberto |
| G | Cache e Zarr regional | Backlog |
| H | Execução junto aos dados | Backlog |
| I | Interface visual simples | Backlog, sem tecnologia aprovada |
| J | CI, empacotamento, CLI e docs | Backlog, não manutenção genérica autorizada |
| K | Sensibilidade/convergência versus níveis verticais | Hipótese, **não medida**; backlog separado da elegibilidade documental C |

A proposta de adendo do lote 15 listou **BL-01 a BL-11**, com IDs expressamente **provisórios e não aprovados**. Abrange validação real independente, NCEP-R2, Dz/De, mensagens de domínio inválido, empacotamento, versão, promoções de branch, releases/Pages, documentação, proteção de tags e governança. Não promover esses IDs à fila oficial.

O encerramento A não validou fisicamente o toolkit. ERA5 37 níveis/3 h, NCEP-R2 e comparação física independente permanecem não demonstrados. **B foi encerrada por triagem, não por produção de referência NCEP-R2.** Antes de desenhar benchmarks D, reavaliar se a referência numérica NCEP-R2 será necessária; manter a pendência explícita enquanto não houver protocolo e evidência. **A investigação documental C não exige estudo de sensibilidade a níveis (K) nem abre D automaticamente.**
