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

## Ampliação aprovada em 25/10/2026 — formato e divisão de trabalho

### Limite de edição de código pela conversa-filha

A filha **pode** consultar GitHub, inspecionar código, executar testes e verificações seguras sem alterar fontes, escrever e integrar documentação/memória autorizada e coordenar o agente. A filha **nunca cria, edita ou remove código-fonte, testes-fonte, scripts, notebooks executáveis, dependências, configurações com efeito executável ou CI**, nem por automações ou APIs. Tais mudanças são **exclusivas do agente de código** após autorização de escopo e aprovação de desenho quando aplicável. A filha pode produzir especificações e prompts, nunca aplicar a alteração executável por conta própria.

### Recomendação dual de modelo e de sessão, em cada prompt técnico

Imediatamente antes do prompt para o agente, exibir **somente**: `Codex — [modelo e esforço]`, `Claude Code — [modelo e esforço]`, `Sessão — continuar/nova`. Recomendar as melhores opções atuais e efetivamente disponíveis/informadas, sem justificativa; não inventar acesso, nomes ou esforço. Decidir continuar/nova por qualidade e tamanho do contexto, continuidade da frente, risco de confusão e economia de tokens. Quando indicar sessão nova, incluir no prompt um resumo mínimo de contexto com fonte e SHA. Não exigir que o usuário use as duas ferramentas.

### Forma humanizada, com títulos e evidências

Em **toda resposta substantiva da filha**, usar a ordem: **título da frente → subtítulo da etapa → contexto breve → resumo em linguagem humana do passo anterior do agente (o que fez, descobriu e por que importa; ou indicar que não houve agente) → veredito se aplicável → checklist visual de 3–6 marcos → explicação e escolhas → evidências e riscos/limites → `## O que você precisa fazer agora`**. Presumir que o responsável não leu logs/relatórios. Códigos internos de decisão, ação ou evidência servem à rastreabilidade, mas **sempre nomear primeiro o assunto em português**; não escrever somente “DEC-01” ou similar. Se não houver risco novo, indicar sucintamente, sem seção longa.

### Bloco de abertura e decisões

A conversa-mãe fornece à filha um **bloco de abertura autocontido** com problema, objetivo, fatos já verificados e respectivas fontes/SHAs, incertezas, limites, pré-requisitos, portões de aceitação definidos antes de testes e a exclusão explícita do que não pertence à frente. Não reabrir diagnóstico já comprovado sem motivo. Se a causa estiver incerta, autorizar **Etapa 0 de leitura/medição sem implementação**; a hipótese pode ser descartada legitimamente. Mudanças de ordem da fila e abertura de outra frente exigem decisão do responsável. Em escolhas substantivas, a filha apresenta alternativas, consequências e recomendação; não substituir autorização por ambiguidade.

As regras de portões científicos, baselines, independência de validação, números e erros estão no [protocolo de rigor científico e verificação computacional](scientific-quality-protocol.md) e são aplicadas conforme o risco, sem reescrever a política de branches ou autorizar mudanças não aprovadas.

## Checkout único, artefatos externos e revisão antes do commit (27–28/10/2026)

> **Sucedida em parte em 28/10/2026.** O portão de autorização antes de `git add`/`git commit` (passo 4 do fluxo
> abaixo), a frase obrigatória e a exigência de manifesto aprovado caminho a caminho para descartáveis foram
> substituídos prospectivamente pela [autonomia de ponta a ponta dentro do escopo aprovado](#autonomia-de-ponta-a-ponta-dentro-do-escopo-aprovado-28102026).
> Continuam valendo o checkout único, a separação entre Git e área externa e a triagem dos artefatos. O texto abaixo
> fica preservado como registro histórico.

**Origem:** preferência expressa do responsável em 27/10/2026 (proposta v3 de checkout único e emenda de revisão
pré-commit). Redação submetida à revisão do responsável antes do commit, em 28/10/2026. Vale para todas as frentes
futuras e para qualquer trabalho ainda não commitado, seja qual for o agente de código. Não altera retroativamente o
histórico (ver o [registro da Frente D](evidence/records/2026-10-27-front-D-precommit.md)).

### Checkout único por padrão

O agente trabalha numa branch do checkout já usado pelo responsável e aberto no VS Code. Não cria worktree, clone ou
checkout adicional por conveniência.

- **Antes de mudar de branch:** verificar o estado Git e possíveis conflitos.
- **Se a branch estiver ocupada noutro worktree, ou houver alterações pendentes incompatíveis:** parar e pedir
  decisão, sem `stash`, `reset`, `clean`, cópia ou troca automática de branch.
- **Exceções:** worktree, commit temporário ou outro atalho só com autorização expressa e anterior, que nomeie a
  necessidade técnica, o escopo e o plano de retorno. Isso é compatível com a
  [política de desenvolvimento](../development_workflow.md), que admite worktree isolado quando houver modificações
  pendentes noutra árvore.

### Dois espaços com finalidades distintas

1. **Repositório Git, na branch de trabalho:** só o que se destina ao produto ou à memória:
   - código-fonte e testes versionáveis;
   - configuração e CI, quando autorizadas;
   - exemplos e dados auxiliares mínimos e redistribuíveis;
   - documentação do usuário aprovada (`docs/source/`);
   - sínteses permanentes curadas (`docs/evolution/`).
2. **Área de trabalho externa ao Git:** tudo o que é de investigação, preservando proveniência e hashes:
   - downloads e respostas brutas, NetCDF intermediários;
   - experimentos e harnesses descartáveis, logs extensos;
   - ambientes, builds e temporários;
   - pré-registros e selos, relatórios completos de investigação, patches de apoio;
   - capturas, prévias PDF/HTML e figuras antes × depois, inventários.

   Estar fora do Git **não** significa reter indefinidamente: aplica-se a [GOV-RET-01](retention-policy.md).

**Triagem antes de integrar.** Cada item gerado recebe uma classe explícita:

- **produto versionável:** vai ao Git após revisão visual, aprovação de conteúdo e autorização de commit;
- **memória permanente curada:** decisões, escopo, critérios, resumo de evidências inclusive negativas, ressalvas,
  SHA e proveniência; vai a `docs/evolution/` sem logs integrais, dados brutos, caminhos pessoais ou segredos;
- **evidência temporária externa:** fica fora do Git até o portão de retenção;
- **descartável recriável:** identificado em inventário, sem exclusão por iniciativa própria.

Um hash ou um caminho na máquina não é, sozinho, backup.

### Fluxo com portão antes de `git add` e `git commit`

1. **Preflight:** checkout único, branch, base e estado Git; delimitar tarefas e destinos das saídas.
2. **Trabalho:** investigar e testar. As evidências vão para a área externa. As mudanças propostas ao produto ficam
   na branch **sem staging e sem commit**, e os arquivos novos ficam visíveis como *untracked*: não escondê-los com
   `.gitignore`.
3. **Apresentação:** mostrar ao responsável:
   - `git status --short`, o resumo por arquivo e o diff não commitado no VS Code;
   - os testes, com fontes e hashes;
   - para documentação, interfaces e figuras, a renderização visual e, quando houver base, o **antes × depois**,
     com a prévia mantida fora do repositório.

   A aprovação técnica do código não substitui a aprovação visual/documental.
4. **Portão:** **aguardar autorização expressa antes de `git add` ou `git commit`.** Aprovar escopo ou execução
   autônoma não autoriza commit. Push, PR, merge, `main`, tag, release, publicação e exclusões têm portões próprios.
5. **Integração:** incluir no repositório **somente os arquivos aprovados** (produto e memória curada), nunca a
   pasta inteira do experimento.
6. **Retenção:** reavaliar a área externa no encerramento da frente e após versão/tag estável, com manifesto
   aprovado antes de exclusões concretas.

**Ausência do responsável.** O agente pode experimentar, corrigir e testar, preservando checkpoints não destrutivos
fora do Git. Ao chegar ao ponto de commit sem aprovação, não faz commit: deixa as mudanças intactas e visíveis.

**Frase obrigatória nos prompts técnicos ao agente:**

> Trabalhe na branch do checkout já existente, sem criar worktree ou clone adicional. Não faça commits nem stage
> antes de eu examinar no VS Code as alterações não commitadas e autorizar o commit explicitamente. Guarde
> downloads, resultados de experimento, logs completos e prévias visuais fora do repositório, com proveniência, e
> proponha ao Git apenas o produto versionável e a memória curada. Se o checkout estiver bloqueado ou em conflito,
> pare e peça instrução.

Mostrar prévia visual de documentação e figuras também é dever explícito do agente.

## Autonomia de ponta a ponta dentro do escopo aprovado (28/10/2026)

**Origem:** decisão expressa do responsável em 28/10/2026. Sucede prospectivamente o portão antes de
`git add`/`git commit` de 27–28/10, sem reescrever o histórico nem a ocorrência registrada da Frente D.

1. **Autorização inicial da frente.** Delimita objetivos, mudanças permitidas, critérios científicos pré-declarados
   e categorias de artefatos (produto versionável, memória curada, evidência externa, descartável recriável).
2. **Execução dentro do escopo.** O agente pode desenvolver, testar, revisar, fazer staging, commit e push, abrir PR
   e executar **merge normal em `develop`**, desde que os checks obrigatórios (GitHub Actions e CircleCI) estejam
   verdes e as proteções da branch satisfeitas, sem bypass.
3. **Retenção no encerramento.** A retenção é auditada automaticamente. Caches, builds, temporários e duplicatas
   **comprovadamente** descartáveis podem ser eliminados como parte do encerramento, sob a autorização geral da
   frente, com manifesto verificável (caminhos exatos e hashes) e registro durável da execução
   ([GOV-RET-01, §7](retention-policy.md)).
4. **O que nunca é apagado automaticamente:** materiais exclusivos, resultados científicos ainda necessários e itens
   de função incerta (classe INVESTIGAR).
5. **Sem portões repetidos.** Não se exige supervisão, inspeção visual individual nem aprovações repetidas para
   passos cobertos pela autorização inicial. Relatórios, diffs e prévias continuam sendo produzidos e ficam
   inspecionáveis, sem se tornarem bloqueios.
6. **Quando parar e consultar:** alteração substancial de escopo, risco científico novo, conflito relevante, falha
   de segurança ou possível perda de evidência não descartável.
7. **Continuam exigindo autorização específica:** `main`, tags, releases, publicações (PyPI, TestPyPI, Pages),
   force-push e alterações destrutivas fora da frente. Worktrees, branches e stashes de outras frentes também.

O responsável pode revisar o trabalho a qualquer momento; a revisão posterior não invalida a autonomia. A regra de
releitura após o merge (SHA final de `develop` e conteúdo versionado) permanece.
