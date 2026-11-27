# Visão geral e limites

**Corte: 26/11/2026** (antes, 07/11/2026). O objetivo atual do Projeto Evolução é incluir novas fontes de download de ERA5 no LECTK e documentar cada uma, com validação científica externa, entregando tudo numa única nova versão oficial ao fim (decisões do responsável de 07/11).

O LorenzCycleToolkit (LECTK) realiza diagnósticos energéticos do ciclo de Lorenz. O Projeto Evolução organiza melhorias científicas, computacionais e de obtenção de dados sob decisões explícitas, mantendo o arquivo local `infile` como fluxo principal (D-001) e a CDS API como integração opcional (D-002).

Decisões assentadas: preservar equações e convenções por padrão (D-003); distinguir compatibilidade de ERA5 37 níveis / 3 h de equivalência de resultados (D-004); não adotar nova arquitetura sem diagnóstico/decisão (D-005); interface visual apenas no backlog (D-006). Repositórios de correção de artigo ou reprocessamento são separados do desenvolvimento do toolkit (D-009).

| Data | Marco | Interpretação |
|---|---|---|
| 16/10/2026 | Inventário, propostas e hipóteses | Não eram implementações aprovadas |
| 18/10/2026 | D-010 aprova fila inicial | Aprovação da ordem não inicia as frentes |
| 19/10/2026 | Auditoria original A | Sem PASS global limpo; duas provas parciais |
| 20–23/10/2026 | Ondas e homologações | Unidades/domínio completados e defeito vertical posterior corrigido |
| 23/10/2026 | Encerramento A e D-011 | PASS retrospectivo com ressalvas; fluxo por `develop` |
| 24/10/2026 | Reconciliação documental | 33/33 arquivos indexados conferidos por SHA-256 |
| 25/10–06/11/2026 | Frentes C, D, E, D-Bench e LZV | Ver [fila](queue.md) e [decisões](decisions.md) |
| 07/11/2026 | Replanejamento (corte desta atualização) | Entregável de fontes ERA5 documentadas, versão única ao fim, validação externa Dias Pinto & da Rocha; [registro](evidence/records/2026-11-07-replanning-decisions.md) |
| 08/11/2026 | Etapa 2a encerrada | FAIL pré-registrado (53/89); D-021 |
| 09/11/2026 | Subfrente 2a.1 | Termos de fronteira conferidos, sem erro de fórmula; D-022 |
| 09–19/11/2026 | Etapa 2b (benchmark de fontes) | Concluída em 19/11 (PRs #61–#74, D-023 a D-028); comparação, não eleição de fonte "mais correta" |
| 19/11/2026 | Decisões de rotas e início da etapa 4 | D-029 a D-034; Frente O aberta fora da fila (descontinuada em 26/11, D-046) |
| 23–25/11/2026 | Etapa 4 (4a–4d) | Correções NCAR/CDS, kerchunk padrão, `--wekeo`, `--arco`, `--destine`, documentação, zarr 3; D-035 a D-044 |
| 25/11/2026 | Encerramento da etapa 4 | Retenção validada (PR #90) |
| 26/11/2026 | Atualização de estado (Frente S) | D-045 a D-047; [registro](evidence/records/2026-11-26-frente-s-state-update.md) |
| 27/11/2026 | Etapa 3 (níveis verticais) | Medição simples em 10 casos; 32 → 13 níveis muda os termos em 0,2–5,7 % na mediana; [registro](evidence/records/2026-11-27-stage-3-vertical-levels-results.md). Nota de 27/11/2026: a faixa exclui o termo M (mediana de 14,9 %, classe d); encerrada com retenção validada ([registro](evidence/records/2026-11-27-stage-3-retention.md)) |

**Estado em 27/11/2026:** 1, 2a, 2a.1, 2b, 3 e 4 concluídas; 5 (arrumação) não aberta; 6 (versão única) exige autorização específica.

**Quatro conceitos distintos:** aprovação de uma formulação, reprodução de referência numérica, CI verde e validação física independente. Os três primeiros não provam o quarto. O código posterior ao PyPI 2.0.0 pode continuar registrando o mesmo texto de versão; citar o commit exato nos resultados.
