# Ocorrência 2026-10-27 — commits da Frente D anteriores à revisão visual

**Natureza:** registro histórico de processo, factual e não acusatório. **Histórico Git preservado.**

**O que ocorreu.** Entre 26 e 27/10/2026, o agente de código produziu **13 commits locais** na branch
`feature/frente-d-era5-ncar`, de `838c52d` a `ffa57a5`. Esses commits vieram **antes** de o responsável inspecionar
no VS Code as alterações ainda não commitadas. O trabalho foi feito em um worktree adicional.

**Causa de coordenação.**

- As instruções então vigentes, inclusive a [emenda de autonomia de 23/10](2026-10-23-autonomy.md), e os prompts da
  coordenação autorizavam commits locais na branch experimental.
- A preferência do responsável por revisão antes do commit ainda não estava convertida em portão explícito.
- O agente seguiu as autorizações recebidas: a falha foi do processo de coordenação, não desobediência.

**Tratamento.**

- **Histórico intacto:** nenhum commit foi reescrito; não houve reset, rebase, squash nem force-push para simular
  aprovação anterior.
- **Revisão recuperada sem worktree:** a revisão visual foi feita depois, comparando `origin/develop` (`4dbab15`)
  com a branch no VS Code, a partir de cópias exportadas com `git show` para fora do repositório.
- **Documentação aprovada antes do commit:** a página v4 foi aprovada sobre o PDF efetivamente renderizado. Só então
  foi commitada (`348bb42`) e integrada pelo [PR #38](https://github.com/daniloceano/LorenzCycleToolkit/pull/38)
  (merge `dec13fe`, 28/10/2026), com os 14 commits preservados.

**Consequência normativa.** A regra prospectiva de checkout único, artefatos externos ao Git e revisão antes de
`git add`/`git commit` está no [protocolo de coordenação](../../coordination-protocol.md). Ela não se aplica
retroativamente a este histórico.

**Origem:** propostas de governança do responsável de 27/10/2026, mantidas fora do Git:

- emenda de revisão pré-commit, SHA-256 `0201f36e6785a588535c7fedbefc60addd2e23c64b21ccc99a4aa2d89b5a5835`;
- proposta v3 de checkout único, SHA-256 `5de6229b232feacd00779fc6e6997c807d5abfff08b2d845c27c7b9a1e9a1262`.

A v3 substitui, no ponto do isolamento, a sugestão da emenda de usar worktree próprio.
