# Visão geral e limites

**Corte: 07/11/2026.** O objetivo atual do Projeto Evolução é incluir novas fontes de download de ERA5 no LECTK e documentar cada uma, com validação científica externa, entregando tudo numa única nova versão oficial ao fim (decisões do responsável de 07/11).

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

**Quatro conceitos distintos:** aprovação de uma formulação, reprodução de referência numérica, CI verde e validação física independente. Os três primeiros não provam o quarto. O código posterior ao PyPI 2.0.0 pode continuar registrando o mesmo texto de versão; citar o commit exato nos resultados.
