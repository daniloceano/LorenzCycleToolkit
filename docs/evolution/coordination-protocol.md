# GOV-09 e GOV-10 — Handoff à conversa-mãe e checklist da conversa-filha

**Decisão do responsável: 25/10/2026. Estado: APROVADAS.** Complementam GOV-05 (contexto no início), GOV-08 (ação final), o modelo mãe/filha e a [governança](governance.md). Não alteram competências científicas, autorização de execução, critérios de merge ou a [política Git](../development_workflow.md).

## GOV-09 — Prompt automático de encerramento para a conversa-mãe

**Toda conversa-filha, assim que confirmar o encerramento de sua frente, entrega na mesma resposta, sem solicitação adicional do responsável, um texto autocontido e pronto para copiar/colar na conversa-mãe.** A entrega é parte obrigatória do fechamento, além do registro versionado e da consolidação da memória.

O texto deve ser um **prompt de transferência**, não apenas "frente concluída". Deve conter concisamente:
1. Identidade da frente, escopo, data, estado e tipo de encerramento (PASS, PASS com ressalvas, dispensa de execução por triagem etc.).
2. O que foi feito e verificado; principais evidências, documentos/PRs, branch e SHA final confirmado. Distinguir fatos verificados de afirmações do agente.
3. Achados, limitações, resultados negativos e pendências aceitas; indicar explicitamente o que **não** foi demonstrado.
4. Próxima frente conforme a fila aprovada, dependências e eventual decisão ainda necessária.
5. Instrução à conversa-mãe para **reler o registro integrado**, conferir a fila e atualizar o estado, **sem iniciar automaticamente outra frente ou modificar a ordem**.

**Momento e ressalvas:** só rotular o handoff como "fechamento confirmado" depois de satisfeitos os portões da frente, incluindo confirmação do merge/readback quando aplicável. Se decisão de fechamento já existe mas a integração documental está pendente, entregar um **handoff provisório** claramente identificado e depois entregar o **definitivo automaticamente** quando os portões passarem. Não criar ciclos de confirmação desnecessários.

**Limite de capacidade:** uma conversa-filha não envia mensagens diretamente para outra conversa do ChatGPT. "Automaticamente" significa **gerar e exibir o prompt na resposta de encerramento**, sem precisar que o responsável peça. O responsável o cola na mãe, salvo se existir futuramente uma integração explícita e efetivamente disponível.

O formato do prompt deve ser de fácil cópia, curto o bastante para transferência (preferência por 5–10 linhas substantivas, mais links/SHAs essenciais) e **sem exigir que a mãe leia relatórios privados completos** para entender a decisão. A mãe mantém autoridade para verificar a documentação e atualizar a fila.

## GOV-10 — Checklist visual obrigatório em toda resposta da filha

**Toda resposta substantiva de conversa-filha deve começar com contexto autossuficiente breve** (GOV-05) e incluir, logo após o contexto/veredito, um **checklist visual simples da frente em curso**, atualizado a cada resposta.

Regras de apresentação:
- Preferir **3 a 6 linhas curtas**, uma por etapa principal; não reproduzir subtarefas técnicas, logs ou evidências extensas.
- Usar caixas Markdown `- [x]` (concluída e verificada), `- [ ]` (pendente) e, quando necessário, qualificador legível **em andamento** ou **bloqueado**; nunca marcar concluído algo apenas planejado, aprovado mas não executado, ou CI ainda pendente.
- As etapas dependem da frente: exemplo **diagnóstico → decisão do desenho → execução/testes → homologação/PR → retenção/fechamento → handoff à mãe**. Em frentes documentais, substituir pelas etapas realmente aplicáveis.
- Manter **o mesmo conjunto de marcos principais** ao longo da frente, atualizando estado para permitir orientação rápida. Acrescentar/remover um marco somente se o escopo mudar comprovadamente.
- Distinguir **fechamento científico/técnico**, **integração documental**, **retenção operacional** e **transferência à mãe** quando forem estados diferentes.
- Em respostas após fechamento, checklist pode mostrar todos os marcos concluídos e apenas o handoff pendente, se for o caso. **Não gerar um plano extenso.**

Ordem de resposta da filha: **título da frente → contexto em 1–3 frases → veredito, quando houver → checklist curto → resposta à questão atual → `## O que você precisa fazer agora`**. A regra **não** exige checklist na conversa-mãe, embora a mãe possa usá-lo quando for útil.

## Exemplo mínimo de resposta de filha

**Contexto:** A triagem de referências da Frente B terminou. ERA5 tem 49 referências versionadas, e NCEP-R2 permanece sem referência numérica aprovada.

**Checklist — Frente B**
- [x] Inspecionar referências e testes atuais
- [x] Decidir tratamento de NCEP-R2 (pendência documentada)
- [x] Homologar encerramento e integrar registro em `develop`
- [ ] Encaminhar prompt final à conversa-mãe

**Conclusão:** B encerrada por triagem; não declarar regressão numérica NCEP-R2 validada.

Ao entregar o prompt colável de GOV-09 na mesma resposta, marcar o último item como concluído (a transferência efetiva para a mãe depende de o responsável colar o prompt).
