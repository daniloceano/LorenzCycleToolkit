# Frente E — encerramento: integração da decisão arquitetural e retenção GOV-RET-01 (01/11/2026)

**Estado operacional: ENCERRADA — RETENÇÃO VALIDADA**, com a decisão arquitetural aprovada pelo responsável em
01/11/2026 integrada em `develop` e a auditoria de retenção (Marco A) concluída. A confirmação do merge deste
próprio registro e a sua releitura não podem constar dele: ficam no relatório de encerramento e no handoff à
conversa-mãe. **Nenhum código foi implementado nem alterado.** A Frente F continua condicional e **não aberta**.

## Escopo autorizado

O responsável aprovou em 01/11/2026 o desenho e o contrato descritos na
[decisão da Frente E](../../fronts/front-E-era5-provider-architecture.md) e autorizou, para esta frente: integração
documental por branch, PR e merge normal em `develop` com os portões satisfeitos; inventário de retenção; e descarte
apenas de descartáveis comprovados da própria frente ([GOV-RET-01 §7](../../retention-policy.md)). Fora do escopo:
código, testes, scripts, dependências, CI, dados, flag nova (`--ncar-method` fica só especificada), Frente F,
`main`, tags e releases.

## Integração documental

- **Documento aprovado:** recebido com SHA-256
  `d67f987722f1b06f4d1db123eb9c82f01acc456d27d660f68fd1018edff61133` (conferido antes de editar). Integrado como
  [`fronts/front-E-era5-provider-architecture.md`](../../fronts/front-E-era5-provider-architecture.md) **byte a
  byte**, sem ajuste editorial: seus links relativos já apontavam para esse destino e resolvem para arquivos
  rastreados. O estado da frente fica nos índices e neste registro, para não alterar o texto aprovado.
- **Índices:** [`queue.md`](../../queue.md), [`fronts/README.md`](../../fronts/README.md) e
  [`decisions.md`](../../decisions.md) passaram a registrar a frente aberta, com decisão aprovada, em
  integração/encerramento. Passagens datadas de C, D e D-Bench ficaram inalteradas.
- **PR [#53](https://github.com/daniloceano/LorenzCycleToolkit/pull/53):** base `develop` em `835478e`, head
  `7ef84af`, merge normal `03ed907dd18b1bcc3fca67c78c2a110be449cbdb` (pais `835478e` e `7ef84af`). Checks
  obrigatórios verdes: `build` do GitHub Actions nos eventos push e pull_request (573 testes aprovados, nenhum
  pulado, em cada um) e `ci/circleci: build_test`.
- **Releitura em `develop`:** o merge altera só os quatro arquivos Markdown acima; o documento integrado tem o mesmo
  SHA-256 do aprovado, também no conteúdo servido pelo GitHub; 47 de 47 links relativos e a âncora do §6 resolvem;
  `main` permanece em `48fe8e2` e é ancestral de `develop`. O checkout local foi sincronizado por fast-forward,
  sem alterações pendentes. O GitHub Actions do push de `03ed907` em `develop` também passou (573 aprovados); o
  CircleCI desse commit ainda estava em execução quando este registro foi escrito.

## Retenção GOV-RET-01 — Marco A

**Raízes inventariadas, em modo somente leitura e delimitadas:** o checkout (arquivos não rastreados e ignorados,
com tamanho e data, antes de qualquer escrita e no fim); a área temporária do agente nesta sessão; os dois arquivos
de entrada fornecidos pelo responsável; a área externa criada para a Frente E; o armazenamento interno da ferramenta
do agente; e uma busca rasa por nome nos temporários do sistema. Diretórios pessoais não foram varridos.

| Item | Classe | Resultado |
|---|---|---|
| Decisão e índices integrados (PR #53) e este registro | Produto versionado e memória permanente | no Git |
| Checkout | — | **sem resíduos da frente**: nenhum arquivo não rastreado; os 17 515 arquivos ignorados pré-existentes são idênticos à linha de base (caminho, tamanho e data) |
| Área temporária do agente: 26 arquivos, 1 link simbólico, 5 diretórios | DESCARTAR — **executado** | 21 duplicatas byte a byte das evidências retidas; 4 prévias HTML reproduzidas bit a bit a partir de `7ef84af`; 2 itens triviais do teste de controle |
| Inventário, manifesto, executor, registro da execução, pós-verificação, linhas de base, verificador de links, controles do executor, capturas de CI e releitura | RETER TEMPORARIAMENTE, fora do Git | auditoria desta limpeza; revisão na próxima versão/tag estável (Marco B) |
| Entradas do responsável (decisão aprovada e instrução da tarefa) | Fora do descarte | preservadas intactas; a decisão também está no Git, byte a byte |
| Armazenamento interno da ferramenta do agente | Fora do escopo | não tocado |
| Branches da frente (local e remota) | Refs Git | retidas; removê-las não fazia parte da autorização |
| INVESTIGAR | — | nenhum item |

**Execução:** manifesto fixado em SHA-256 `2c7bf22fc082e5c171a034d9347b8fb2b8d1283dec5624206c1abecab5764643`
(inventário `0f5c8825f51fccf1dd30f582e1b93ff7aef90831c4d8efaed5eaa6d53f55e2ad`); executor
`4d95d3798ab8ec4c35934ba836d3961647a130c2e9eb881473ade8f7d9ace173`, que confere imediatamente antes de cada
exclusão a raiz autorizada, a ausência de links no caminho, o tipo de arquivo, a contagem de hardlinks, o tamanho e o
SHA-256. Resultado: **26 arquivos (3 316 288 bytes, ≈ 3,2 MiB), 1 link e 5 diretórios removidos, 0 falhas**; o
espaço livre do volume subiu cerca de 3,2 MiB. Registro durável fora do Git: execução
`889c2ec0414207a9cd88a134b35bf539b19a69a6c10a547082e6a545102b3899`, resumo
`73b4c0328ddaee135b0946ac2b5258b6143d4db758e09055b49813b8104b57d1`, inventário final privado
`4b1ec8c5256c24d2b016d607bfe44fbd7bda2115a7a2d8be04ce36b4bf4e6f7d`.

**Pós-verificação:** 21 de 21 cópias retidas íntegras; os 32 itens do manifesto ausentes; entradas do responsável
com os SHA-256 originais; checkout sem resíduos; worktrees e stash de outras frentes inalterados.

**Controle do executor e erro registrado:** antes do uso real, o executor bloqueou um hash errado, um caminho por
link simbólico, um manifesto com hash divergente e um link com alvo divergente, e removeu corretamente os itens
válidos, sem tocar o alvo do link. O controle revelou um defeito do próprio executor: execuções no mesmo segundo
sobrescreviam o registro umas das outras. Foi corrigido antes da execução real e os controles foram repetidos;
nenhum material preservado foi afetado.

**Não tocados:** D-Bench, `LEC_Results/`, materiais e worktrees da Frente D, outros worktrees, stash, dados
científicos, `main`, tags e releases. Nenhum download, aquisição ou execução remota.

## Ressalvas

- A decisão é de **desenho**: `--ncar-method`, a verificação uniforme do contrato e qualquer piloto dependem de
  autorização própria; os métodos NCAR continuam experimentais.
- O documento integrado preserva a redação aprovada, inclusive o §7, que remete o encerramento a este registro.
- A revisão das retenções deve ser reaberta na próxima versão/tag estável; isso não autoriza criar tag.
