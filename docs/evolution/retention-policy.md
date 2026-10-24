# GOV-RET-01 — Retenção, arquivamento e limpeza no encerramento de frentes

**Estado:** APROVADA pelo responsável em 2026-10-24. **Âmbito:** LorenzCycleToolkit | Evolução. Complementa a [política de desenvolvimento](../development_workflow.md), sem alterá-la.

## 1. Regra central

Toda frente terá uma **auditoria de retenção e limpeza como último portão de encerramento operacional**, após conclusão técnica, homologação, integração autorizada e documentação verificada. Isso não é autorização automática de apagar arquivos.

## 2. Classificação obrigatória

| Classe | Tratamento |
|---|---|
| **PRESERVAR** | Fontes originais, evidências científicas indispensáveis, dados/entradas de referência, resultados de validação, decisões e arquivos exclusivos; guardar integridade e proveniência |
| **ARQUIVAR** | Documentos/resultados úteis a auditoria e reprodução, mas dispensáveis no armazenamento ativo; manter cópia recuperável |
| **EXCLUIR** | Caches, ambientes recriáveis, builds, temporários e duplicatas cuja redundância esteja **comprovada**; executar somente após aprovação de manifesto específico |
| **INVESTIGAR** | Itens de origem, reconstrução ou relevância incertas; **não excluir** até esclarecer |

Nunca classificar um diretório inteiro como descartável sem conhecer seu conteúdo, inclusive ocultos e arquivos ignorados pelo Git.

## 3. Procedimento

1. **Inventário somente leitura:** relacionar caminhos exatos, tamanhos, origem, SHA quando pertinente, dependências, arquivos não versionados/ocultos, links, duplicatas por conteúdo, dados de entrada, resultados e itens necessários a trabalhos ainda ativos. Não modificar originais.
2. **Preservação verificável:** confirmar armazenamento independente para **PRESERVAR** e **ARQUIVAR**, integridade e capacidade de **restaurar e ler**. Hash ou resumo no Git não substitui os bytes não versionados.
3. **Manifesto de limpeza:** documentar caminho exato, classe, justificativa, risco, tamanho recuperável, destino de backup e prova de recuperabilidade. Separar espaço potencial de espaço efetivamente liberado.
4. **Execução autorizada:** **após aprovação expressa do responsável sobre o manifesto de exclusões**, agentes podem remover somente os caminhos autorizados sem solicitar aprovação individual por arquivo. Proibido ampliar escopo, apagar diretórios-pai genericamente, seguir links para fora do escopo, destruir históricos Git, ou apagar arquivos duvidosos. Diante de divergência, **interromper** a exclusão afetada.
5. **Conferência final:** verificar integridade/acesso dos preservados, registrar cada exclusão, espaço recuperado e pendências, e produzir relatório com data, responsável e evidências.

## 4. Proteções

- **Não eliminar** datasets meteorológicos únicos, entradas de análises, referências numéricas congeladas, scripts irreconstruíveis, experimentos negativos, configurações essenciais e arquivos ligados a artigos sem evidência de dispensabilidade.
- Idade, volume, origem por agente e aparência de duplicata **não bastam** para exclusão.
- Incluir na auditoria links simbólicos, hardlinks, `.gitignore`, arquivos ocultos, árvores não versionadas e dependências entre frentes.
- Não publicar dados privados no GitHub público como suposto backup. Cópias no mesmo volume e meros índices de hash não são backups independentes; testar restauração representativa.
- Nenhuma autorização implícita para apagar Project Sources, `.git`, refs, worktrees, diretórios científicos ou dados dos computadores.

## 5. Estados de encerramento

- **ENCERRADA — RETENÇÃO VALIDADA:** destinos definidos e verificáveis; limpeza autorizada realizada/documentada, **ou** ausência justificada de exclusões possíveis.
- **CONCLUÍDA TECNICAMENTE — RETENÇÃO PENDENTE:** há material a investigar, backup a criar ou plano de exclusão pendente de autorização. A frente não deve ser declarada com retenção validada.
- Esta auditoria não substitui homologação, testes, CI ou validação científica independente.

## 6. Aplicação retroativa

A **Frente A permanece formalmente ENCERRADA em 2026-10-23, PASS COM RESSALVAS ACEITAS**. Sua decisão científica não será reaberta. Aplicar auditoria de retenção **retroativa e complementar** aos diretórios históricos `LECTK_Dev`, `LECTK_Auditoria` (Mac) e `LECTK_Retomada` (notebook); a primeira etapa será **somente leitura**, sem exclusões.

**Aprovar esta política não é aprovar qualquer exclusão específica.** O plano de remoção de arquivos exige decisão separada.
