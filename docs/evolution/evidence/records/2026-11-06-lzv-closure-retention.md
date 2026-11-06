# LZV — encerramento: integração, auditoria do servidor e retenção GOV-RET-01 final (06/11/2026)

**Estado operacional: ENCERRADA — RETENÇÃO VALIDADA.**
- **Integração:** a investigação complementar ERA5/Zarr e níveis verticais ([registro de 02/11](2026-11-02-lzv-era5-zarr-vertical-levels.md)) está integrada em `develop`.
- **Retenção final:** executada e verificada em 06/11/2026, depois da indisponibilidade do servidor de processamento em 03/11.
- **Merge e releitura deste registro:** não podem constar dele. Ficam no relatório de encerramento e no handoff à conversa-mãe.
- **Nada executável mudou:** nenhum código, teste, script, dependência ou CI foi alterado.
- **Sem abertura de frentes:** **não abre a Frente F nem a hipótese K**, não reabre D, D-Bench ou E, e `main`, tags e releases estão inalterados.

## Escopo autorizado

O responsável autorizou em 06/11/2026:
- preservar fora do servidor os registros necessários à auditoria e à reprodução;
- descartar **exclusivamente** o export redundante da correção, conforme [GOV-RET-01 §7](../../retention-policy.md);
- sincronizar com `develop`, só por fast-forward, o checkout do servidor;
- registrar o encerramento por branch, PR e merge normal em `develop`.

Ficaram fora do escopo:
- aquisições ERA5 e medições de desempenho;
- a Frente F e outras investigações;
- os dados retidos da D-Bench e `LEC_Results/`;
- `main`, tags e releases.

## Integração (verificada)

- **Merge:** o PR [#55](https://github.com/daniloceano/LorenzCycleToolkit/pull/55) foi integrado em `develop` em 03/11/2026 por merge normal `24ffe5b30dc11ac78c2c0c7c8732733657dff9c5` (pais `b4a5f77` e `7e170ac`). Contém a correção `2101de2` e o ajuste de revisão `7e170ac`, e a árvore do merge é idêntica à de `7e170ac`.
- **CI:** os checks obrigatórios de `develop` no merge aprovaram: GitHub Actions `build` e CircleCI `build_test`. `main` permanece em `48fe8e2`, ancestral de `develop`.
- **Releitura de 03/11:** o conteúdo integrado é igual ao aprovado. A suíte local deu 595 aprovados e 1 falha preexistente (`test_packaging`, metadados do pacote instalado no ambiente local). O Sphinx em modo estrito passou.
- **Servidor, 06/11:** depois da sincronização descrita abaixo, os 23 testes de `tests/test_track_alignment.py` passaram no ambiente científico do grupo.

## Indisponibilidade e auditoria do servidor

- **03/11/2026:** a validação final da retenção ficou **BLOQUEADA** porque o servidor não respondia a partir do host de acesso. Os itens retidos não foram tocados.
- **06/11/2026, auditoria somente leitura:**
  - o servidor tinha voltado depois de uma reinicialização, sem processos antigos;
  - as evidências da LZV estavam íntegras;
  - o export da correção foi classificado como **DESCARTÁVEL COMPROVADO**;
  - o checkout do servidor ainda estava em `835478e`, sem a correção.

## Cópia independente fora do servidor

**Onde e quando:** antes de qualquer exclusão, na área externa ao Git do notebook do responsável.

**Manifesto:** o manifesto próprio tem 675 entradas SHA-256. Todas foram conferidas contra a referência do momento da cópia: a lista tirada no servidor imediatamente antes da transferência, o inventário v2 selado, o manifesto da D-Bench ou a lista do download.

| Conjunto | Arquivos | Tamanho | Verificação |
|---|---:|---:|---|
| Memória permanente da LZV: relatório, protocolo selado, adendo v2, inventários, manifestos v3/v4, logs (inclusive `deviations.log`), ledger HTTP, resultados (inclusive negativos), 15 execuções LECTK, ferramentas, figuras, receitas de ambiente e registros de retenção | 618 | 19,1 MB | 618/618 |
| NetCDF ARCO A1 e A3 e NetCDF NCAR, retidos até o Marco B | 3 | 1,18 GB | 3/3; cópia retomável, renomeada só depois de o hash conferir |
| Os dez tracks originais da D-Bench, copiados sem alterar a origem | 10 | 15,6 kB | 10/10 contra o `tracks_sha256` do manifesto da D-Bench |
| As dez referências kerchunk da GDEX usadas na Etapa B (T, Z, W, U, V; décadas de 1970 e 1980) | 10 | 170 MB | 10/10 contra o servidor e contra a lista do download de 30/10 |
| Receita do ambiente científico e recursos da máquina | 8 | 18 kB | capturados nesta data |

**Teste de restauração:**
- **Resultado:** as cópias abrem normalmente. Com as coordenadas alinhadas, ARCO A1 = A3 = NCAR bit a bit nos cinco campos (13 níveis × 27 horários), como registrado em 02/11.
- **Erro do instrumento de teste, corrigido:** uma primeira comparação, feita sem alinhar as coordenadas, indicou diferenças espúrias. As latitudes estavam em ordem inversa e as longitudes em convenções diferentes (0–360 × −180–180).

**Continua só no servidor:**
- **O espelho ARCO A1** (1788 objetos, 2,6 GiB), retido até o Marco B. É recuperável byte a byte da fonte pública, e o SHA-256 de cada objeto está no inventário v2, que tem cópia independente.
- **Os dados retidos da D-Bench**, que seguem o [registro de retenção](2026-11-01-dbench-retention-cleanup.md) próprio. A cópia não os incluiu, exceto os tracks e as referências kerchunk.

## Retenção GOV-RET-01 — rodada 4

- **Escopo exclusivo:** o export da correção (`code/lectk_fix_20261102/`), única retenção condicionada ao commit.
- **Prova, reverificada imediatamente antes:**
  - 410/410 arquivos com SHA-256 idêntico ao blob Git do mesmo caminho em commits publicados e alcançáveis a partir de `develop` `24ffe5b`: 395 em `24ffe5b`, 13 em `b4a5f77` e 2 em `2101de2`;
  - a lista própria do export confere 238/238;
  - 50 diretórios, sem links, hardlinks ou entradas especiais;
  - nada modificado depois de 02/11.
- **Manifesto e executor:**
  - manifesto fixado em SHA-256 `d3dd4b120db7e9723e482905076d2cf5b696adece6609636c26e6526f277fb6b`;
  - executado pelo executor da rodada 3, sem alteração (SHA-256 `ed680f3defe671a973be51486addb362344e2892ba1e6d4b27839697d8f8ace1`);
  - o executor confere tipo, tamanho e SHA-256 antes de cada exclusão, não sai da área e para na primeira falha.
- **Resultado:**
  - 410 arquivos excluídos, **0 falhas**, 49 370 918 bytes de conteúdo;
  - 50 diretórios vazios removidos;
  - **50 540 544 bytes recuperados em disco** (≈ 48,2 MiB);
  - no sistema de arquivos, a variação foi de −50 466 816 bytes, porque os novos registros da rodada ocupam 73 728 bytes;
  - a área foi de 3 965 112 320 para 3 914 645 504 bytes em disco.
- **Registro durável na própria área:** `report/RETENCAO_RODADA4_2026-11-06.md` (SHA-256 `347cd5bf1ec14b6b9804ca4f3d403a33d5f7f19685fc379130ae4d9e8533e0fb`), com log por arquivo e verificação posterior.
- **Verificação posterior:**
  - os selos do pré-registro, do adendo, do inventário v2 e dos manifestos v3 e v4 conferem;
  - o índice de entregáveis confere 85/86 (ver o `deviations.log` abaixo);
  - o inventário v2 não tem ausência sem explicação: os 1795 itens retidos estão íntegros, e os 1837 ausentes estão todos listados no manifesto v3;
  - nenhum item dos manifestos v3 e v4 permanece.
- **Desvio registrado:** a primeira execução do script de verificação posterior parou antes de gravar qualquer coisa, porque um selo usa caminho relativo. Foi reexecutada depois de ajuste no próprio script, sem alterar a área.

**Histórico do `deviations.log`:**
- **Selo:** o índice de entregáveis selou o arquivo com cinco linhas, e as cinco primeiras linhas atuais reproduzem o SHA-256 selado (`665573c6…`).
- **Sexta linha:** foi acrescentada em 02/11/2026, às 18:02:45 UTC, durante a rodada 3. Ela registra o `chmod u+w` nos diretórios do export somente leitura de `835478e`, necessário à exclusão; os arquivos não foram alterados.
- **Situação atual:** é por isso que o arquivo diverge do índice. A ocorrência foi preservada e nada foi acrescentado nesta rodada; os desvios da rodada 4 estão no registro dela.

| Retido até o Marco B | Cópia independente | Finalidade |
|---|---|---|
| Espelho ARCO A1 | Não. É recuperável da fonte pública e tem hash por objeto | Bytes brutos da aquisição a frio |
| NetCDF ARCO A1 e A3 | Sim, verificada | Entradas do LECTK e prova da igualdade ARCO = NCAR |
| NetCDF NCAR e proveniência | Sim, verificada | Entrada da verificação da correção |
| Track de 3 h, track de 6 h, horários e namelist | Sim, verificada | Entradas pequenas dos experimentos |

**Memória permanente da área:** relatório, protocolo, adendo, ledger, logs, resultados, execuções, ferramentas, figuras, manifestos e registros de retenção ficam no servidor, com cópia independente verificada. Reproduzir C0–C3 continua exigindo o arquivo CDS de referência retido pela D-Bench.

## Checkout do servidor

**Preflight:**
- `develop` estava em `835478e`, sem commits exclusivos, stash, hooks ou alterações rastreadas;
- havia um arquivo não rastreado e dois ignorados, de outras atividades, cujos SHA-256 foram registrados;
- a worktree de outra investigação estava em `d38cda7`;
- o avanço até `24ffe5b` não adicionava caminhos em `inputs/`, `LEC_Results/` ou caches, nem mudava dependências.

**Execução:** `git fetch` e merge **somente fast-forward** (`835478e` → `24ffe5b`, 28 arquivos).

**Depois:**
- os 416 arquivos rastreados conferem com os blobs de `HEAD`;
- os arquivos locais mantêm os mesmos SHA-256 e datas;
- a worktree da outra investigação e os arquivos ignorados ficaram inalterados;
- nenhuma dependência foi instalada;
- não houve `reset`, `clean`, `stash` nem force-push.

## Portões do encerramento

| Portão | Situação |
|---|---|
| Científico | Interpretação aprovada pelo responsável em 02/11 (§5 do [registro](2026-11-02-lzv-era5-zarr-vertical-levels.md)) |
| Integração | PR #55 em `develop` (`24ffe5b`), CI aprovada, releitura concluída |
| Retenção operacional | Rodada 4 executada e verificada; retenções documentadas com finalidade, cópia e marco; nenhum item em INVESTIGAR |
| Documental | Este registro, a [fila](../../queue.md) e o registro de 02/11 atualizados, preservando a redação original |

## Ressalvas

- **O espelho ARCO A1 continua só no servidor.** Perdê-lo exigiria nova aquisição da fonte pública.
- **Marco B:** a revisão das retenções deve ser reaberta na próxima versão/tag estável. Isso não autoriza criar tag.
- **Ressalvas científicas do §6 do registro de 02/11, inalteradas:**
  - vantagem geral do Zarr não demonstrada;
  - GRIB não verificado;
  - um único caso examinado;
  - DestinE bloqueado;
  - nenhuma validação física independente.
- **D-Bench:** a cópia externa desta investigação não cobre os dados retidos da D-Bench.
