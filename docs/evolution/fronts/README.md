# Frentes documentadas

A [Frente A](front-A-audit.md) — auditoria/linha de base do LECTK 2.0.0 — foi **formalmente encerrada em 23/10/2026: PASS COM RESSALVAS ACEITAS**. [Termo](front-A-closure.md).

A [Frente B](front-B-closure.md) foi **encerrada por triagem em 24/10/2026**, sem produção de novas referências; a lacuna NCEP-R2 foi mantida explicitamente como pendência, a reavaliar antes de benchmarks D.

A **Frente C foi encerrada em 25/10/2026 com PASS DOCUMENTAL**, após aprovação do conteúdo, [PR #36](https://github.com/daniloceano/LorenzCycleToolkit/pull/36) integrado em `develop` (`abaab7add04c21c694a4aa90da38c7d1465e5add`), CI verde e releitura pós-merge. [Avaliação, critérios, ressalvas e encerramento](front-C-era5-documentary-assessment.md). Não houve benchmarks nem implementação de provedores; naquela data, a Frente D ainda não havia sido aberta.

A **Frente D** foi executada de 26 a 28/10/2026 como avaliação de rotas alternativas de aquisição ERA5. O adaptador NCAR GDEX foi integrado em `develop` como funcionalidade **experimental** ([PR #38](https://github.com/daniloceano/LorenzCycleToolkit/pull/38), merge `dec13fe6dcacc2e2cd7f5720160ea40d40e75403`). A comparação de tempo com o CDS **não foi realizada**: a P2-B0 foi bloqueada antes de qualquer submissão, um resultado negativo registrado. **Encerrada formalmente em 28/10/2026: PASS COM RESSALVAS ACEITAS**, com retenção validada (limpeza parcial e retenções temporárias até a próxima versão/tag). Benchmarks controlados e comparação NCAR × CDS não concluídos; sem superioridade de velocidade nem validação física demonstradas. A Frente E não foi aberta. [Registro da Frente D](front-D-era5-ncar.md).

A **Frente E** (decisão de arquitetura de provedores) permanece **não aberta** e a **Frente F** (piloto de provedor), **condicional**; G–K são backlog e hipóteses. [Fila](../queue.md).

Após conclusão técnica de cada frente, aplicar [GOV-RET-01 Revisão 2](../retention-policy.md): descartar os resíduos operacionais comprovadamente dispensáveis e justificar somente a retenção mínima necessária. Reavaliar os retidos no fechamento da próxima versão/tag estável, sem exclusão automática. A Frente A receberá auditoria retroativa, sem reabertura da decisão científica.

Frentes futuras devem preservar escopo, SHA de base, decisões e autores, limites de validação, resultados negativos, diff de integração e desfecho formal.

**Comunicação de cada frente:** aplicar [GOV-10](../coordination-protocol.md) em toda resposta da conversa-filha (contexto + checklist curto), e [GOV-09](../coordination-protocol.md) no encerramento (prompt colável, entregue automaticamente para a conversa-mãe).
