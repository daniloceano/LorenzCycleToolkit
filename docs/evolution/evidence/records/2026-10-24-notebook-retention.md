# Auditoria de retenção retroativa — notebook (24/10/2026)

**Status:** limpeza F1A e limpeza consolidada **HOMOLOGADAS DOCUMENTALMENTE** pela coordenação em 24/10/2026, com base em pacotes de evidências inspecionados. **A Frente A mantém seu encerramento científico de 23/10, PASS COM RESSALVAS ACEITAS;** esta é uma conclusão operacional parcial da auditoria de retenção, não uma nova validação científica. Política aplicável: [GOV-RET-01 Revisão 2](../../retention-policy.md).

## Escopo e resultados

| Etapa | Material | Resultado relatado e auditado |
|---|---|---|
| F1A | 78 grupos de temporários pytest; excluído o grupo de traceback exclusivo NB-1528 | **126.593 arquivos removidos**, sem falhas; **26,30 GiB** adicionais de espaço livre no volume, conforme medição local |
| Limpeza consolidada | Categorias F1B–F1D, F2, F3 e EXP; exatamente 756 grupos autorizados | **158.788 itens removidos**: 139.580 arquivos, 1.387 symlinks e 17.821 diretórios; **25,46 GiB** adicionais de espaço livre, conforme medição local |
| Total das duas etapas | Resíduos operacionais de diretório histórico da Frente A | **285.381 itens removidos**; **aproximadamente 51,76 GiB** somados de espaço livre; ocupação do diretório histórico reduzida de **56,14 GiB** no inventário inicial para **4,42 GiB** após a limpeza consolidada |

**Documento de evidência consolidada recebido e conferido:** arquivo ZIP local, SHA-256 `9bf6ff5782b427e12eb7e76e17f0f5b8a466ef34873f7cf8e8e2a16db33dc8f0`. Manifesto aprovado `MANIFESTO_CONSOLIDADO_v1.csv`, SHA-256 `37397b95d697a74e36ee958bed6eda72a52837b4dcec2c08779255c8259cb8c5`. O ZIP e os logs detalhados **não foram publicados** no repositório público.

## Auditoria documental independente do ZIP consolidado

- Arquivo ZIP legível e íntegro; **23/23 checksums de conteúdo** constantes no `SHA256SUMS.txt` interno correspondem exatamente aos respectivos arquivos.
- Manifesto recuperado com hash idêntico ao aprovado; **756 grupos** pertencem somente às sete categorias autorizadas.
- Log de execução: **158.788 linhas de exclusão bem-sucedida**, sem outra ação; exatamente os mesmos caminhos da simulação.
- Inventário antes/depois: **158.788 caminhos removidos**, sem inclusões, e conjunto removido exatamente igual ao log; nenhuma remoção fora das categorias selecionadas; tipo e inode de cada entrada removida conferidos contra o inventário anterior e a lista original do manifesto.
- Arquivos preservados com mudanças de metadados: **12.922 hardlinks** compartilhados com o cache removido passaram de `nlink=2` para `nlink=1`; o conjunto de bytes alocados correspondente é **352.874.496**, exatamente a diferença entre a soma por arquivo removido e a estimativa de bytes recuperáveis. Os demais atributos de conteúdo e identificação desses arquivos não divergiram no inventário.
- As medições locais registradas mostram aumento de **25,46 GiB** em `df` e redução de **25,44 GiB** em `du` na execução consolidada; essas medidas não são equivalentes e refletem diferenças de contabilização e metadados.
- As verificações locais pós-execução relatam como intactos: **4.492 evidências F1A**, **NB-1528**, **26 entradas canônicas dos experimentos**, ambiente Python-base, wheels, referências e worktrees Git.

**Limite da homologação:** a inspeção independente validou **os bytes do pacote fornecido, consistência interna do manifesto, logs e inventários**; o disco físico, os repositórios locais e as medições `df`/`du` foram observados exclusivamente pelo agente local e constam dos registros anexados ao ZIP privado. Não declarar verificação física independente. A afirmação sobre disponibilidade pública das dependências foi verificada pelo agente local, não reexecutada nesta homologação.

## Pendências protegidas — não reabrir como falha desta limpeza

1. **F4 (≈0,58 GiB):** dez worktrees e um clone auxiliar não estavam autorizados; continuam preservados, inclusive operações Git potencialmente destrutivas.
2. **Resíduo F1A:** symlinks e diretórios vazios, sem ganho relevante de espaço; não autorizado para exclusão.
3. **Retenção temporária R-01 a R-09 e C-01 (ensaios com commits Git não publicados):** reavaliar no **Marco B**, no fechamento de nova versão/tag; nenhum descarte automático.
4. **Backup independente dos itens efetivamente retidos:** ainda não configurado; avaliar conforme finalidade, sem criar cópias de arquivos descartáveis.
5. **Outra máquina:** auditoria dos seus diretórios de trabalho permanece separada e não foi executada por este registro.

**Decisão:** homologar a limpeza executada **sem ampliar autorização para qualquer outro caminho ou categoria**. O portão operacional de retenção do notebook ainda mantém pendências documentadas; não declarar GOV-RET-01 integralmente encerrada para todos os equipamentos.
