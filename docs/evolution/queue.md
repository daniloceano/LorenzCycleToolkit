# Fila de frentes — decisão de 18/10 e atualização de 23/10

A ordem inicial aprovada foi **A → B condicional → C → D → E → F condicional** (D-010). G–J são backlog; K foi acrescentada como hipótese não testada. Fechar A não abre automaticamente B/C.

| Frente | Objeto | Estado/documento |
|---|---|---|
| A | Auditoria e linha de base científico-computacional | **ENCERRADA 23/10: PASS COM RESSALVAS ACEITAS**; [termo](fronts/front-A-closure.md) |
| B | Lacunas de dados de referência/regressão **do próprio LECTK** | **ENCERRADA 24/10 — triagem concluída; execução de novas referências dispensada por ora**; NCEP-R2 permanece pendente. [Termo](fronts/front-B-closure.md) |
| C | Levantamento de provedores ERA5 e requisitos científicos | Aguardar decisão contextual; registro posterior menciona autorização de investigação documental, **não execução por esta migração** |
| D | Benchmarks controlados | Não aberta; exige protocolo e dependências |
| E | Decisão de arquitetura de provedores | Não aberta; depende de D |
| F | Piloto de provedor | Condicional, não aberto |
| G | Cache e Zarr regional | Backlog |
| H | Execução junto aos dados | Backlog |
| I | Interface visual simples | Backlog, sem tecnologia aprovada |
| J | CI, empacotamento, CLI e docs | Backlog, não manutenção genérica autorizada |
| K | Sensibilidade/convergência versus níveis verticais | Hipótese, **não medida**; backlog |

A proposta de adendo do lote 15 listou **BL-01 a BL-11**, com IDs expressamente **provisórios e não aprovados**. Abrange validação real independente, NCEP-R2, Dz/De, mensagens de domínio inválido, empacotamento, versão, promoções de branch, releases/Pages, documentação, proteção de tags e governança. Não promover esses IDs à fila oficial.

O encerramento A não validou fisicamente o toolkit. ERA5 37 níveis/3 h, NCEP-R2 e comparação física independente permanecem não demonstrados. **B foi encerrada por triagem, não por produção de referência NCEP-R2.** Antes de desenhar benchmarks D, reavaliar se a referência numérica NCEP-R2 será necessária; manter a pendência explícita enquanto não houver protocolo e evidência. **C é a próxima frente da fila; sua abertura e seu escopo permanecem sujeitos à decisão específica.**
