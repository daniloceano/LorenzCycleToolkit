# Governança da memória e execução

O coordenador organiza decisões, a conversa-filha investiga e o agente executa tarefas aprovadas. O responsável autoriza mudanças, define escopo e aceita ressalvas. O workflow vinculante do repositório é [docs/development_workflow.md](../development_workflow.md); não reproduzir seu texto protegido.

**Regras gerais:** nova branch de `develop` atualizado; diff limitado; visualização e aprovação de mudanças documentais/figuras/formulações conforme GOV-06/07; para formulação, comparação científica explícita; testes de push/PR e CI de Actions/CircleCI verdes; preservar o histórico; sem bypass, force-push, publicação, promoção à `main` ou exclusão destrutiva sem decisões específicas.

**Emenda 23/10:** uma aprovação expressa de modificações e escopo autoriza atos posteriores de commit, push, PR e merge para a branch pertinente **apenas com portões satisfeitos**. Veja [registro da emenda](evidence/records/2026-10-23-autonomy.md).

**Exceção 24/10:** a pedido expresso do responsável, somente os documentos gerados por agentes para a migração inicial `docs/evolution/` podem ser integrados em `develop` **sem inspeção humana individual**, desde que verificados automaticamente e com CI verde. Não é revogação geral GOV-06/07. Veja [registro](evidence/records/2026-10-24-exception.md).

**Repositório público:** não publicar segredos, paths pessoais, datasets não redistribuíveis, logs/arquivos científicos volumosos, PDFs protegidos ou instruções internas integrais do Project. Publicar sínteses e índices saneados. Os **33 hashes** comprovam igualdade dos documentos selecionados, não backup integral. Manter as pastas originais `LECTK_Dev`, `LECTK_Auditoria`, `LECTK_Retomada` até inventário completo, backups independentes, teste de restauração e autorização específica para exclusão.

**GOV-RET-01 (24/10/2026):** a [auditoria de retenção](retention-policy.md) é o último portão operacional de cada frente. Inventário inicialmente somente leitura; qualquer exclusão exige manifesto específico aprovado. Para a Frente A, aplicação retroativa sem reabrir o encerramento científico.

Após merge, conferir SHA final de `develop`, reler o material versionado e registrar a situação real; sem readback não declarar sucesso.
