# GOV-RET-01 — Retenção mínima e descarte no encerramento de frentes

**Revisão 2 — aprovada pelo responsável em 2026-10-24.** Substitui prospectivamente a redação inicial da GOV-RET-01, preservando-a no histórico Git. Aplica-se ao LorenzCycleToolkit | Evolução e complementa a [política de desenvolvimento](../development_workflow.md), sem modificá-la.

## 1. Princípio: retenção mínima, temporária e vinculada à finalidade

**O descarte é a regra para os artefatos operacionais de frentes encerradas.** O objetivo de encerrar uma frente é preservar **conhecimento, código e evidências suficientes**, não manter indefinidamente seus diretórios de trabalho.

A memória versionada em `docs/evolution/` preserva decisões, achados e limitações; o Git preserva o código e os arquivos efetivamente versionados. **Nenhum deles, sozinho, preserva os bytes de dados, scripts ou resultados científicos exclusivos que não estejam no repositório.** Portanto, exceções à limpeza precisam de razão concreta, dependência identificada e um marco de reavaliação.

**Não é exigido backup de arquivos comprovadamente descartáveis.** Backup e teste de restauração são necessários para aquilo cuja retenção seja justificada e cuja única cópia seria perdida.

## 2. Dois marcos de limpeza

### Marco A — encerramento de cada frente

Após homologação técnica, integração autorizada, verificação de resultados e consolidação da memória documental:

1. Auditar os artefatos locais relevantes, inicialmente em modo **somente leitura**.
2. Identificar o que já está de fato preservado no Git ou em acervo verificável e o que permanece exclusivo.
3. **Preparar para descarte** caches, temporários, builds, ambientes comprovadamente recriáveis, execuções intermediárias desnecessárias e duplicatas verificadas.
4. **Reter excepcional e temporariamente** apenas os itens que ainda sejam necessários à reprodução justificável, validação, integração ou continuidade. Registrar **por que**, **até quando/qual evento** e **onde** estão.
5. Executar apenas as exclusões especificamente aprovadas pelo responsável em manifesto de caminhos.

### Marco B — fechamento de nova versão estável / tag

Depois que uma nova versão tiver sido **deliberadamente homologada, promovida conforme a governança e marcada por tag autorizada**, reavaliar os materiais **retidos de frentes anteriores**.

A criação da versão/tag é um **gatilho de revisão para descarte**, **não** prova automática de que arquivos externos ao Git foram preservados e **não** concede autorização automática de exclusão. Eliminar o que deixou de ser necessário, mediante manifesto aprovado. Permanências além desse marco exigem justificativa renovada e próximo marco de revisão.

A auditoria de Marco B não autoriza, não cria e não publica versão, tag, `main`, PyPI ou Pages: esses atos mantêm aprovações próprias.

## 3. Classes e critério de decisão

| Classe | Regra | Necessidade de backup |
|---|---|---|
| **DESCARTAR — CANDIDATO** | Operacional, recriável, redundante ou sem finalidade remanescente **demonstrada**; registrar caminho e justificativa | **Não**, quando a dispensabilidade estiver comprovada |
| **RETER TEMPORARIAMENTE** | Exclusivo e ainda necessário; explicitar dependência, responsável, destino e marco de reavaliação (preferencialmente próxima versão/tag) | **Sim**, quando a perda da cópia impediria a finalidade alegada |
| **ARQUIVAR EXCEPCIONALMENTE** | Relevância histórica/científica não suprida por memória/Git; retenção mais longa com justificativa | **Sim**, com integridade e restauração verificadas |
| **INVESTIGAR** | Evidência insuficiente para avaliar exclusividade, reconstrução ou dependência | **Não apagar** até resolver; priorizar investigação |

**Ônus:** justificar cada retenção. **Segurança:** incerteza sobre a função de um arquivo nunca constitui demonstração de dispensabilidade.

## 4. Procedimento verificável

1. **Inventariar em modo somente leitura:** totais, tamanhos, arquivos ocultos e ignorados pelo Git, origens, hashes seletivos, symlinks/hardlinks, dependências entre frentes, entradas e outputs científicos, artefatos de testes negativos, versões locais não publicadas e histórico Git exclusivo. Não tratar uma pasta inteira como descartável sem seu inventário.
2. **Classificar e justificar:** primeiro separar materiais ainda necessários; depois identificar descartáveis. Para retenções, registrar data/marco de revisão; para descarte, explicar a prova de recriabilidade ou redundância. Conferir se os objetos supostamente presentes no Git estão **efetivamente rastreados no commit/na tag relevante**, e não apenas mencionados na memória.
3. **Proteger apenas o necessário:** para os artefatos retidos, confirmar cópia independente e integridade/restauração quando cabível; **não copiar lixo operacional por precaução genérica**.
4. **Apresentar manifesto de exclusão:** caminhos/IDs exatos, quantidade de bytes potenciais, dependências, avaliação de risco, justificativa e evidência de descarte seguro. Marcar itens não decididos como INVESTIGAR, fora da lista autorizável.
5. **Executar somente após aprovação expressa do manifesto pelo responsável**, sem autorizações arquivo a arquivo, mas **sem ampliar o escopo**. Não usar exclusões destrutivas em Git, não seguir links para fora da raiz autorizada; worktrees, refs e stashes exigem procedimento e autorização próprios.
6. **Conferir e registrar:** listar o que realmente foi excluído, medir espaço efetivamente recuperado, validar que itens retidos permanecem acessíveis e registrar pendências. Na próxima versão/tag repetir a revisão dos remanescentes.

## 5. Proteções que permanecem

- Achados científicos já documentados **não implicam** que todos os dados subjacentes precisem ser mantidos; avaliar se seus bytes são necessários ou se parâmetros, códigos e resultados versionados bastam para a finalidade aceita.
- Também **não implicam** que todos os dados externos sejam dispensáveis: entradas únicas, resultados negativos exclusivos, scripts não versionados, refs/stashes/commits só locais e dados de artigos não podem desaparecer sem análise de dependência.
- Arquivo grande, antigo, gerado por agente ou aparentemente duplicado **não** equivale a arquivo comprovadamente recriável.
- GitHub público não é destino apropriado para material restrito ou segredo; uma tag não contém arquivos externos que nunca foram commitados.
- Esta política **não é** autorização de exclusão de Project Sources, diretórios pessoais, dados, worktrees ou histórias Git. Cada limpeza depende de autorização específica para seus caminhos.

## 6. Estados operacionais e aplicação retroativa

- **ENCERRADA — RETENÇÃO VALIDADA:** conhecimento/produto integrados, descartáveis aprovados removidos (ou ausência de descarte justificada), e retenções excepcionais documentadas com finalidade e marco. Registrar o espaço liberado.
- **CONCLUÍDA TECNICAMENTE — RETENÇÃO PENDENTE:** resta manifesto não aprovado, investigação essencial, preservação excepcional sem destino ou verificação incompleta. A parte científica pode estar encerrada sem que o portão operacional esteja finalizado.

A Frente A **permanece formalmente encerrada em 2026-10-23, PASS COM RESSALVAS ACEITAS**; não reabrir seu julgamento científico. Aplicar este procedimento **retroativamente** a `LECTK_Retomada` (notebook) e `LECTK_Dev`/`LECTK_Auditoria` (Mac). O inventário do notebook relatou aproximadamente 45,92 GiB de candidatos; esse volume é **potencial**, ainda não é exclusão autorizada nem espaço confirmado como recuperável. As faixas envolvendo worktrees, refs e experimentos exclusivos exigem triagem específica. A falta de acesso ao Mac não impede a análise das faixas comprovadamente independentes do notebook.

**Em qualquer marco, aprovação desta regra ≠ aprovação de exclusão específica.**

## 7. Emenda de 28/10/2026 — descarte sob a autorização geral da frente

Decisão do responsável que altera prospectivamente o passo 5 do Marco A e o passo 5 do §4, sem mudar as classes do §3.

- **Descartáveis comprovados** (caches, builds, ambientes recriáveis, temporários de testes, duplicatas verificadas
  byte a byte) podem ser eliminados no encerramento da frente **sob a autorização geral da frente**, sem aprovação
  caminho a caminho, desde que haja:
  - manifesto verificável, com caminhos exatos, contagens, bytes e hash de conteúdo, fixado por sha256;
  - conferência imediatamente antes de cada exclusão, sem seguir links e sem ampliar o escopo;
  - registro durável da execução, com espaço estimado e recuperado.
- **Não entram nesse descarte automático:** RETER TEMPORARIAMENTE, ARQUIVAR, INVESTIGAR, materiais exclusivos,
  resultados científicos ainda necessários, worktrees, branches, stashes e qualquer item fora da frente.
- Falha de validação ou de execução interrompe a limpeza e é registrada; não se contorna a falha nem se apaga à mão.

O Marco B (versão/tag) segue o mesmo regime para retidos que tenham se tornado descartáveis comprovados.
