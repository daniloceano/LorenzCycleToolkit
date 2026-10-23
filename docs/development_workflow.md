# Development workflow

<!-- EDITORIAL PROPOSAL (2026-10-23): this introduction is not part of the approved policy text and awaits the project lead's review. -->

This page records the branch, integration and release policy adopted by the maintainers of LorenzCycleToolkit on 2026-10-23. In
short:

- new work starts on a new branch created from the latest `develop`;
- pull requests target `develop`;
- `main` receives deliberate promotions from `develop`, each by a pull request of its own;
- tags are created on `main`, and publishing to PyPI, TestPyPI or GitHub Pages is a separate, manual step.

The contributing guide (`docs/source/contributing.rst`) summarizes what this means for contributors.

The policy follows in its approved wording, in Portuguese, without changes. The section "Adoção imediata — PR #26" describes the
transition of 2026-10-23 and its SHAs are a record of that day.

<!-- BEGIN APPROVED POLICY TEXT: verbatim copy of LECTK_politica_branches_2026-10-23.md (sha256 e9f4e6bb9f961d6a4a0ba1ec37dd6b98728837cfea04ac2e0af44ba6823ac0ab). Changes inside these markers need a new approval. -->

# LECTK — Política de branches, integração e releases

**Data da decisão:** 2026-10-23  
**Estado:** princípio de workflow aprovado pelo responsável nesta conversa; texto preparado para incorporação à documentação permanente, não sincronizado automaticamente com o GitHub ou os arquivos do Project.

## Decisão operacional

O fluxo obrigatório para novas frentes de desenvolvimento é:

`feature/<frente>` → Pull Request → `develop` → Pull Request de promoção → `main` → tag anotada → publicação manual autorizada.

1. **`main`** é o estado integrado estável e a origem das tags de release; não é o destino rotineiro de PRs de feature.
2. **`develop`** integra trabalhos concluídos e homologados que ainda não foram promovidos à `main`. Criá-la inicialmente no mesmo SHA de `origin/main`. Ela pode ficar à frente de `main`, mas não deve ficar atrás: quando `main` avançar por uma correção urgente ou uma promoção, incorporar essas mudanças a `develop` de modo verificável.
3. **Nova frente** começa sempre em branch de trabalho nova e limpa, criada a partir da versão mais recente de `origin/develop`, em worktree isolada quando houver modificações pendentes em outra árvore. Antes de trabalhar, registrar SHA-base, estado limpo, ausência de divergência relevante e expectativa de testes. Não começar a partir de uma branch experimental ou antiga.
4. **Antes do PR para `develop`**, verificar que a branch incorpora as atualizações de `develop`; diferenças de CI entre head e base não podem ser ignoradas. O PR deve conter somente a frente autorizada, com testes e evidências exigidos.
5. **CI** deve testar tanto o checkout da branch quanto a integração com sua base (GitHub Actions e CircleCI conforme aplicável); resultados verdes, sem skips críticos e sem mascarar falhas, são condição para merge. Verificar diferenças de plataforma, instalação efetiva do pacote, dados auxiliares verificados e publicação desabilitada por padrão.
6. **Governança**: aprovação de desenho científico, apresentações, figuras e documentação antes dos commits pertinentes (GOV-06/07). Commits, push, PR e merge requerem autorização no escopo aplicável; a automação de publicação não é autorizada por merge.
7. **Promoção `develop` → `main`** ocorre em lote deliberado, por PR próprio, quando o conjunto incluído em `develop` estiver pronto e homologado. Registrar escopo, versão, testes, regressões e impacto científico. Uma vez integrada, criar tag anotada no commit exato de `main` **somente após autorização explícita da versão e da tag**. Publicar pacote e documentação somente por passos manuais explicitamente autorizados.
8. **Correções urgentes na `main`**: fazer branch limpa da `main`, PR urgente para `main` com aprovação própria; em seguida sincronizar a correção com `develop` antes da próxima feature. Nunca redefinir ou descartar trabalho acumulado em `develop`.
9. **Preservação**: sem force-push, reset destrutivo, limpeza de branches históricas, exclusão de worktrees, releases ou mudança de versões implícita. Preferir merges que preservem histórico; documentar exceções de transição.

## Adoção imediata — PR #26

**Constatação feita no GitHub em 2026-10-23:** `origin/main` aponta para `48fe8e2e75d0c1527214cedda145924ab13ae9ea`; `develop` ainda não existe no remoto; PR #26 continua aberto contra `main`, head `ecc79ed4e63998f287dc2748719f30a13817337c`. Esses dados exigem reconfirmação imediatamente antes de qualquer escrita remota.

1. Suspender/substituir a orientação anterior de merge do PR #26 diretamente em `main`.
2. Conferir `origin/main`, worktrees e a ausência de uma branch `develop` divergente em outros remotos/locais. Criar `develop` **no mesmo commit atual da `main`**, sem alterar a `main` nem o worktree que contém DEC-09/DEC-11. Publicar essa branch somente após preflight e com autorização aplicável.
3. **Exceção única de migração:** como a branch da DEC-09/DEC-11 já foi iniciada antes desta política, preservá-la e fazer merge normal de `origin/develop` nela, sem rebase destrutivo nem force-push. Não refazer os commits aprovados. Conferir que seu diff relativo a `develop` continua contendo apenas DEC-09/DEC-11 e que o comportamento numérico permanece igual.
4. Alterar a base do PR #26 de `main` para `develop`, registrar no texto do PR a mudança do fluxo e atualizar a branch. Não abrir PR duplicado desnecessariamente.
5. Exigir CI verde tanto no push da branch quanto no PR contra `develop`, além da verificação de segurança de publicação. Depois de confirmar os checks e a integridade, o merge do PR #26 deve ser **somente em `develop`**, respeitando a autorização e proteções existentes. `main` permanece em `48fe8e2` até uma promoção futura deliberada.
6. Registrar evidência: SHA inicial de `develop`, SHA da atualização do head, diff final do PR, checks completos, SHA do merge em `develop`, `main` intacta e ausência de publicação.
7. Preparar atualização dos arquivos permanentes de governança e orientação aos colaboradores para inspeção do responsável antes de seus commits. Não declarar que o texto deste arquivo foi automaticamente salvo no Project ou no repositório.

## Interpretação de “sincronizado com main”

No momento inicial, `develop == main`. Após integrar features, `develop` estará normalmente **à frente de `main`**. A invariante correta é que `develop` contenha o commit mais recente de `main` em seu histórico (`git merge-base --is-ancestor origin/main origin/develop`) e mantenha os trabalhos ainda não promovidos. Exigir igualdade permanente eliminaria o propósito de `develop`.

## Limites

Esta política não determina um número de versão, uma tag, uma data de release ou a publicação do pacote. Não altera metodologia científica, testes, dados de referência ou os resultados da DEC-09/11. As falhas preexistentes de validação de domínios geográficos permanecem em frente separada.

<!-- END APPROVED POLICY TEXT -->
