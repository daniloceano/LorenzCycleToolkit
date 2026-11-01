# Decisões e sucessão temporal

| ID | Decisão | Vigência |
|---|---|---|
| D-001 | `infile` local é o fluxo principal | Assentada |
| D-002 | CDS API é opcional | Assentada |
| D-003 | Não alterar equações/convenções/resultados por padrão | Assentada |
| D-004 | ERA5 37 níveis e 3 h é referência histórica; compatibilidade não é equivalência | Assentada |
| D-005 | Não implementar nova arquitetura sem diagnóstico e autorização | Assentada |
| D-006 | GUI permanece no backlog | Assentada |
| D-007 | Coordenação mãe/filha/agente com autorização responsável | Complementada depois |
| D-008 | Não modificar Git durante adoção inicial | **Histórica**, superada por escopos posteriores autorizados |
| D-009 | Repositórios de correção de artigo e rerun fora das frentes do toolkit | Assentada |
| D-010 | Fila A, B condicional, C, D, E, F condicional, backlog | Aprovada 18/10 |
| D-011 | Branch limpa da `develop` → PR → `develop`; promoção a `main` só por PR e decisão separada | Aprovada 23/10; [política](../development_workflow.md) |

**GOV-09 (25/10):** entrega automática, pela filha, de prompt pronto para a mãe ao fechar frente, após verificação de integração, com ressalvas e fila. **GOV-10 (25/10):** checklist visual breve e atualizado em toda resposta de filha, logo após o contexto. [Texto vinculante](coordination-protocol.md).

**Ajustes do modelo mãe–filha–agente (25/10/2026, aprovados):** reservar **todas as mudanças de código e componentes executáveis ao agente**; permitir à filha inspeção, testes existentes e documentação/memória autorizada. Antes de cada prompt técnico, oferecer uma recomendação Codex e uma Claude Code, mais orientação continuar/nova sessão. Respostas da filha com título/subtítulo, resumo humanizado do agente, checklist e evidências/limites; citar o conteúdo das decisões antes de IDs. Padronizar bloco de abertura, Etapa 0 condicional e opções com consequências. [Coordenação](coordination-protocol.md).

**Qualidade científica e de evidências (25/10/2026, aprovada):** critérios de aceitação antes de medir, baselines, controles de falha quando pertinentes, independência proporcional ao risco, números rastreáveis, registro de hipóteses refutadas e erros. [Protocolo técnico](scientific-quality-protocol.md).

**GOV-RET-01 Revisão 2 (24/10):** retenção mínima suficiente: descartar o que perdeu função ao encerrar a frente; exceções temporárias justificadas são reavaliadas no fechamento da próxima versão/tag autorizada. A tag é gatilho de avaliação, não autorização de exclusão. O plano de caminhos específicos continua sujeito à aprovação. [Política vigente](retention-policy.md).

**GOV-05 (21/10):** contextualizar respostas. **GOV-06/07 (21/10):** inspeção e aprovação visual/científica/documental, autonomia não implica silêncio como aprovação. **GOV-08 (22/10):** respostas terminam com `## O que você precisa fazer agora`.

**Emenda de autonomia 23/10:** após aprovação de conteúdo e escopo, o agente pode commit, push, PR e merge em `develop` sujeito a diff, CI e proteções; não libera `main`, tags ou publicação. A regra antiga "agente nunca abre PR" é história, não vigente.

**Exceção pontual expressa de 24/10:** para a migração inicial dos registros `docs/evolution/`, o responsável dispensou inspeção humana individual e autorizou merge direto em `develop` **se verificações passarem**. Continuam obrigatórios saneamento de privacidade, consistência, inspeção automática, CI, escopo exato e pós-merge. A dispensa **não** vale para outros documentos, imagens, código, equações, dados privados, `main` ou exclusões. [Registro](evidence/records/2026-10-24-exception.md).

**Frente B (24/10):** o responsável encerrou a triagem de necessidade sem exigir execução de novas referências. As 49 referências ERA5 versionadas permanecem limitadas aos seus cenários; NCEP-R2 não tem referência numérica aprovada e continua pendente. Reavaliar essa necessidade antes dos benchmarks D; não bloquear a investigação documental da Frente C. [Termo de encerramento B](fronts/front-B-closure.md).

**Frente A:** encerramento formal em 23/10 como **PASS COM RESSALVAS ACEITAS**, referência retrospectiva aos dez critérios §20.B; aceita R1 de `mtime` interno `.git` por autofetch. Preservar diagnóstico original 19/10 e achado novo de Onda 4. [Termo](fronts/front-A-closure.md).

**Cronologia Git:** PR #25 `main` `48fe8e2`; #26–28 `develop`, culminando em `81cc000` (estado 23/10). Textos antigos dizendo “PR #28 aberto” representam uma data anterior ao merge. Sem nova tag ou publicação.

**Frente D — diretrizes aprovadas (27–28/10/2026):** adaptador NCAR GDEX como opção **experimental** `--ncar` mais o comando independente `python -m src.utils.ncar_gdex`; float32 nativo como padrão, float64 opcional só no comando independente; fluxo CDS, NetCDF local (D-001, D-002), núcleo científico e esquema de `lec_metadata.json` preservados; documentação v4 aprovada visualmente; 14 commits preservados, sem squash. Integração: PR #38, merge `dec13fe` em `develop` (28/10), CI verde. **Sem comparação de desempenho com o CDS** (P2-B0 bloqueada), portanto sem recomendação baseada em velocidade. Encerramento formal pendente. [Registro](fronts/front-D-era5-ncar.md).

**Checkout único, artefatos externos e revisão antes do commit (27–28/10/2026):** checkout já existente por padrão, sem worktrees adicionais; artefatos de investigação fora do Git, sob GOV-RET-01; revisão do responsável das alterações não commitadas, com prévia visual quando aplicável, **antes** de `git add`/`git commit`. Subordina a autorização geral da Emenda 23/10; não reescreve histórico. [Regras](coordination-protocol.md).

**Autonomia de ponta a ponta (28/10/2026):** a aprovação de escopo da frente cobre desenvolvimento, testes, staging, commit, push, PR e merge normal em `develop` com checks verdes, além da limpeza de descartáveis comprovados no encerramento, com manifesto verificável e registro durável. Sucede prospectivamente o portão antes de `git add`/`git commit` de 27–28/10, preservado como histórico. `main`, tags, releases, publicações, force-push e destruição fora da frente mantêm autorização específica. [Regras](coordination-protocol.md#autonomia-de-ponta-a-ponta-dentro-do-escopo-aprovado-28102026); [GOV-RET-01 §7](retention-policy.md).

**Encerramento da Frente D (28/10/2026): PASS COM RESSALVAS ACEITAS; ENCERRADA — RETENÇÃO VALIDADA.** Adaptador NCAR GDEX experimental integrado e testado no escopo descrito. Ressalvas: benchmarks controlados não concluídos; comparação NCAR × CDS não realizada; nenhuma superioridade de velocidade nem validação física demonstradas; limitações de precisão e cobertura mantidas. Limpeza parcial (G6 residual recriável) e retenções temporárias até a próxima versão/tag. A Frente E não foi aberta. [Registro](fronts/front-D-era5-ncar.md).

**D-Bench — dispensa do cancelamento seletivo CDS (29/10/2026, decisão do responsável):** credencial CDS designada pelo responsável, mantida fora do repositório, declarada por ele válida e exclusiva da D-Bench (não verificada pelo agente); cancelamento seletivo declarado inviável. O portão "nenhuma submissão CDS sem cancelamento seletivo" é substituído por mitigações obrigatórias: um pedido CDS por vez, ledger privado antes de cada submissão, limites de tempo/tentativas no cliente (mudança de código ainda a autorizar e implementar), parada após erro de cota/autenticação/política/timeout e aceitação explícita do risco de job órfão consumir cota. **Nenhuma submissão CDS autorizada por este registro**; os demais portões continuam. [Registro](evidence/records/2026-10-29-dbench-cds-cancellation-waiver.md).
