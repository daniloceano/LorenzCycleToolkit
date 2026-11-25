# Etapa 4 — retenção dos dados na Swell (GOV-RET-01, D-028) e D-044 (25/11/2026)

**Estado operacional: ENCERRADA — RETENÇÃO VALIDADA.**

**Natureza:** só documentação no repositório; nenhum código, teste, dependência ou CI mudou. Na Swell, as exclusões
foram feitas exclusivamente pelo manifesto aprovado. `main`, tags e releases ficam inalterados.

**Base:** `develop` em `b4c68bd` (PR #89); `main` em `48fe8e2`. Branch única do encerramento (D-042):
`docs/stage-4-closure-retention-2026-11-25`.

## Resumo em linguagem simples

As campanhas da 2b e os testes da etapa 4 tinham deixado **80,9 GB** em dez pastas na Swell. Foram **apagados
56,84 GB (52,93 GiB)**: downloads brutos já cobertos pelos arquivos de entrada guardados, cópias idênticas dos mesmos
arquivos, peças de pedidos que falharam, cópias de código idênticas ao Git, ambientes recriáveis e temporários. Ficam
**24,04 GB**, que são as entradas das 78 comparações de regressão (e os resultados de referência), os manifestos, selos
e relatórios da 2b, o código exclusivo (harness e scripts de comparação) e poucos extras. Esses retidos continuam
**só na Swell**: o responsável **dispensou a cópia independente** (§4). Nada foi baixado.

## 1. Escopo e autorização

- **Política:** [GOV-RET-01](../../retention-policy.md), incluindo a emenda do §7; D-028 (dados guardados até o fim da
  etapa 4), D-034 (Swell com consulta prévia), D-042 (branch única) e D-044 (§9).
- **Parte 1:** inventário só de leitura (listagem completa, sha256 de todos os arquivos de 1 MB ou mais, links,
  hardlinks, manifestos da 2b, blobs de cada cópia de código contra o Git). Ver §2.
- **Parte 2:** manifesto apresentado e **aprovado expressamente pelo responsável em 25/11/2026** (todo o manifesto,
  inclusive D-VALORES e a retenção de R-DIAG). Ver §3.
- **Parte 3:** cópia independente dispensada por escolha do responsável (§4); execução em §5.

## 2. Inventário (somente leitura)

Pastas em `/p1-swell/danilocs/`: `lectk_bench2b_20261109` (fase 1), `lectk_bench2b_round2_20261113` (rodada
complementar), `lectk_bench2b_phase2_20261119` (fase 2), `lectk_stage4a_diag_20261122`, `lectk_stage4a_kerchunk_20261123`,
`lectk_stage4a1_default_cap_20261123`, `lectk_stage4b1_wekeo_20261123`, `lectk_stage4b2_arco_20261124`,
`lectk_stage4b3_destine_20261124` e `lectk_stage4d_zarr3_20261125`: **22.954 arquivos e links** (26.439 entradas com
diretórios), 80,9 GB, sem hardlinks.

Achados que sustentam as classes:

- **Dependência entre pastas:** as 68 entradas da fase 2 são links (`lectk_bench2b_phase2_20261119/inputs/*`) para
  arquivos `common_*` da fase 1 e da rodada complementar. Das 68 execuções, os 57 destinos distintos têm só **30
  conteúdos distintos**, porque rotas diferentes dão o mesmo arquivo (ARCO = kerchunk = Get a Subset; CDS = WEkEO).
  Os 68 sha256 de entrada conferem com os selados nos manifestos.
- **O que as 78 comparações leem** (`regress_4d.py` da 4d, sobre `regress_f923d80`): 10 arquivos de entrada da etapa 4
  (caso 19790205), mais `inputs/`, o CSV de referência e a proveniência de cada pasta de referência da etapa 4; e, da
  fase 2, `runs/*/*/run_manifest.json`, `inputs/` e as 57 entradas. Os prefixos de sha256 dos 10 da etapa 4 conferem com
  `regression.tsv` da 4d.
- **Manifestos da 2b:** 223 manifestos de pedido; 159 arquivos `common_*` e 159 de brutos listados, todos com sha256 e
  bytes iguais aos do disco. Brutos de pedidos com sucesso: 111, todos com `common_*` de conteúdo retido. Peças T/U/V/W/Z e
  `JOB.*`: 651, todas de pedidos que falharam e sem `common_*`. ARCO e DestinE não têm bruto (só o `common_*`).
- **Código:** os `code/lectk_develop` da fase 1, rodada 2 e fase 2 têm árvore idêntica, arquivo a arquivo, à dos commits
  `cac34e96`, `20aa986b` e `068c3a99`; `code_f923d80` ao `f923d80`; `4b3/code` ao `cc9729f`. Os tar da fase 1 e da
  rodada 2 têm o mesmo sha256 que `git archive` desses commits. O clone aninhado `4b3/code/.git` tem um único ref local
  (`feat/stage-4b3-destine` em `cc9729f`, que está em `origin/develop`); as tags v1.0.0 e v1.0.7 também estão no
  GitHub; sem stash.
- **Fora do GitHub:** 19 arquivos de cópias de trabalho (rascunhos e versões locais, entre eles `arco_zarr.py` e
  `test_arco_zarr.py` da versão 1 da 4b·2 e `inputs/namelist`, `track`, `track-Melissa` modificados) e o `dev.tgz`.
- **Referências kerchunk:** as da 2b e as da 4a têm sha256 diferentes (a 4a corrigiu o endereço); as da rodada 2 são
  idênticas às da fase 1.
- **Valores da 4d:** por xarray (todas as variáveis e coordenadas), `real_arco` = `4b2/run` e `real_destine_v3` =
  `4b3/run` em valores; só `history` (e `source` e a nota de proveniência no v3) diferem. `4b2/run2` = `4b2/run` em valores.
- **Ambientes:** os cinco `venv` reaproveitam o Python da env `lorenz` (`--system-site-packages`). Pacotes extras:
  aiohappyeyeballs 2.7.1, aiohttp 3.14.5, aiosignal 1.4.0, asciitree 1.0.0, deprecated 3.0.1, fasteners 0.20,
  frozenlist 1.8.0, hda 2.40 (fase 1 e 4b·1), multidict 7.2.1, numcodecs 0.15.1, propcache 0.5.4, wrapt 2.5.2,
  yarl 1.25.1, zarr 2.18.4.

## 3. Classes e manifesto aprovado

Manifesto fixado por sha256 (`MANIFESTO_EXCLUSAO_arquivos.tsv`: caminho, tipo, bytes, sha256, regra e evidência de cada
entrada; `MANIFESTO_EXCLUSAO_diretorios.tsv`: diretórios que ficam vazios):

- arquivos: `246e4cf9150671a5d37c450136d6e2cbc89a54e6b6eb3ce8f5e7bedadc080d05` (15.776 entradas);
- diretórios: `2252064e44f2827616520d2826c00d8a92fc55d8dd0703f52c95da35b4c29ef3` (1.899 entradas).

| Classe | Regra | Entradas | Bytes | Prova |
|---|---|---:|---:|---|
| DESCARTAR | D-DUP-COMUM | 103 | 22.268.330.677 | sha256 igual ao de arquivo retido (gêmeo citado) |
| DESCARTAR | D-BRUTO | 112 | 19.970.745.066 | bruto de pedido com sucesso; manifesto liga ao `common_*` de conteúdo retido; mais o bruto da 4a (gerou o common retido `567ccdce`) |
| DESCARTAR | D-PARCIAL | 651 | 8.275.398.318 | peças de pedidos que falharam, sem `common_*` |
| DESCARTAR | D-VALORES | 2 | 3.099.216.654 | valores idênticos (xarray) aos retidos `4b2/run` e `4b3/run`; sha256 difere |
| DESCARTAR | D-DUP-4D | 2 | 1.380.845.902 | sha256 igual ao de `4b3/cmp/destine_new_19790205.nc` |
| DESCARTAR | D-CODIGO | 8.812 | 1.123.763.162 | cada arquivo é blob do GitHub, ou lixo, bytecode ou artefato de build; tars e `.tgz` |
| DESCARTAR | D-REFS-DUP | 20 | 339.728.138 | sha256 igual ao de cópia retida |
| DESCARTAR | D-ENV | 6.036 | 212.435.475 | `venv` recriável (pacotes na §2) |
| DESCARTAR | D-CODIGO-GIT | 28 | 107.773.534 | clone aninhado inteiramente no GitHub |
| DESCARTAR | D-TEMP | 2 | 59.601.534 | cópias perturbadas do autoteste C3, recriáveis |
| DESCARTAR | D-LINK | 8 | 0 | links que ficariam pendentes, apagados sem seguir |
| RETER | R-ENTRADA | 67 | 23.606.544.931 | entradas das 78 comparações |
| RETER | R-REFS | 15 | 254.806.663 | referências kerchunk |
| RETER | R-DIAG | 2 | 105.651.172 | captura bruta do diagnóstico da 4a |
| RETER | R-EXTRATO | 20 | 24.983.847 | arquivos fora do GitHub, no lugar |
| RETER | R-EVIDENCIA | 7.074 | 50.001.600 | manifestos, selos, relatórios, logs, scripts, CSVs de referência |
| INVESTIGAR | — | 0 | 0 | nenhum arquivo de 1 MB ou mais ficou sem classe |

**Fora do escopo e não tocados:** `/p1-swell/danilocs/LorenzCycleToolkit` (clone em `develop` `c0e3942`, mais um arquivo
não rastreado `inputs/namelist_WRF_bigua`; **INVESTIGAR** em outra ocasião) e `lectk_zarr_vertical_20261102` (LZV, outra
frente já encerrada, retido até o Marco B).

**Ambientes conda:** `lorenz-4d` (744.745.439 B) é o ambiente do projeto (D-043) e fica. `lorenz` (824.371.471 B) **não
foi apagado** (proposta do responsável; pode servir fora do projeto); reavaliar no Marco B.

## 4. Cópia independente: dispensada por escolha do responsável

- **Proposta:** camadas 0 e 1 (evidência, código e scripts exclusivos; os arquivos da DestinE v2, 4,52 GiB), em pasta local
  fora do Git; camada 2 (demais entradas, 11,45 GiB) dependia de decisão, pelo espaço no Mac Mini.
- **Decisão do responsável (25/11/2026):** pular a cópia. Motivos dados: é muito mais rápido baixar de novo na Swell que
  transferir, e ele não quer muitos dados no computador. **A escolha é dele**; ficou registrada aqui.
- **O que aconteceu:** a camada 0 (7.053 arquivos, 435 MB) foi copiada e conferida; a camada 1 foi iniciada e
  interrompida em cerca de 2,1 GB de 4,85 GB. **A cópia local foi apagada por inteiro** a pedido do responsável (a
  primeira tentativa de `rm` falhou em arquivos selados somente leitura; liberei a escrita só dentro dessa pasta e apaguei).
  Nada dessa cópia resta no Mac Mini.
- **Consequência:** os retidos têm **uma só cópia, na Swell**. Em particular, os 12 conteúdos da DestinE v2 (4,52 GiB: 10
  `common_destine_*` da fase 2, `4b3/run/case_19790205.nc` e `4b3/cmp/destine_new_19790205.nc`) só se recuperariam de
  novo se o store v2 continuar sendo servido (D-043 o mantém como alternativa enquanto for servido) e, mesmo então, a
  igualdade de bytes não está garantida. O responsável decidiu assumir isso.

## 5. Execução

- **Executor:** conferência por entrada imediatamente antes de apagar (tipo, tamanho, `nlink` = 1, sha256 igual ao do
  manifesto, gêmeo retido presente, nenhum componente do caminho é link); não segue links; usa só `unlink` e `rmdir`; para na
  primeira falha. Antes de apagar, confere o sha256 do manifesto, a existência de todas as entradas e o sha256 dos 82
  retidos de dados e referências.
- **Tentativa 1 (22:29 UTC):** falhou na primeira entrada com `PermissionError`; só o tar
  `code/lectk_cac34e96…tar` (50.278.400 B, sha256 conferido) foi apagado. Causa: as árvores `code/lectk_develop` da fase 1
  e da rodada 2 estavam somente leitura (diretórios `dr-x------`). A execução parou, como a regra exige.
- **Autorização do responsável (25/11/2026):** liberar `u+w` apenas nos diretórios sem escrita que estão na lista de
  diretórios do manifesto, e retomar. O executor foi alterado só para isso e para aceitar o tar já apagado.
- **Tentativa 2 (22:33 UTC, 1 min):** `chmod u+w` em **104 diretórios do manifesto** (todos dentro dessas duas árvores);
  **15.775 arquivos e links apagados (56.787.560.060 B) e 1.899 diretórios removidos; 0 falhas.** Com o tar da tentativa 1:
  15.776 entradas, 56.837.838.460 B.
- **Registro durável:** `/p1-swell/danilocs/lectk_stage4d_zarr3_20261125/retention_stage4_20261125/` (7,0 MB;
  `SHA256SUMS` conferido): manifestos, o documento de aprovação, os dois logs por arquivo, a conferência final e os
  três scripts. sha256: `exec_log.txt` `60e4ea3e…7c42`, `exec_log_tentativa1.txt` `a97acf44…0490`, `verify_after.txt`
  `2a5fc0a9…4cec6`, `exec_manifest.py` `dae45898…3c36`, `exec_manifest2.py` `122b623f…2c6`, `verify_after.py` `f40914f0…6666`.
  Esta pasta foi criada depois da execução e não estava no manifesto; é só registro.

## 6. Conferência posterior (somente leitura)

- **Apagado:** 0 dos 15.776 arquivos/links e dos 1.899 diretórios do manifesto continuam presentes.
- **Escopo:** a pasta tem exatamente as 8.764 entradas esperadas (inventário inicial menos o manifesto): **0 a mais, 0 a
  menos**.
- **Retidos íntegros:** 7.120 arquivos retidos (67 entradas, 15 referências e demais de evidência) com o **mesmo sha256**
  de antes; 58 links retidos, **nenhum quebrado**.
- **Espaço recuperado** (`du -sb`, tamanho aparente; o sistema de arquivos é compartilhado e o `df` mudaria por outros
  usuários):

| Pasta | Antes (B) | Depois (B) | Recuperado (B) |
|---|---:|---:|---:|
| `lectk_bench2b_20261109` | 33.164.367.524 | 8.066.218.553 | 25.098.148.971 |
| `lectk_bench2b_round2_20261113` | 29.804.646.344 | 4.825.795.756 | 24.978.850.588 |
| `lectk_bench2b_phase2_20261119` | 175.322.107 | 14.864.322 | 160.457.785 |
| `lectk_stage4a_diag_20261122` | 105.671.662 | 105.671.662 | 0 |
| `lectk_stage4a_kerchunk_20261123` | 1.729.278.229 | 711.312.086 | 1.017.966.143 |
| `lectk_stage4a1_default_cap_20261123` | 1.705.346.010 | 1.654.718.930 | 50.627.080 |
| `lectk_stage4b1_wekeo_20261123` | 2.941.351.002 | 2.879.389.824 | 61.961.178 |
| `lectk_stage4b2_arco_20261124` | 3.396.461.622 | 3.209.527.768 | 186.933.854 |
| `lectk_stage4b3_destine_20261124` | 2.966.222.301 | 2.531.400.480 | 434.821.821 |
| `lectk_stage4d_zarr3_20261125` | 4.891.167.137 | 43.095.266 | 4.848.071.871 |
| **Total** | **80.879.833.938** | **24.041.994.647** | **56.837.839.291 (52,93 GiB)** |

(O `4d` depois da conferência ganhou 7,0 MB de registro, §5, fora desta medida.)

## 7. Retenções (todas até o Marco B, D-013, salvo nota)

| O que | Onde (`/p1-swell/danilocs/`) | Finalidade | Reavaliação |
|---|---|---|---|
| 57 destinos de entrada da fase 2 (30 conteúdos) | `lectk_bench2b_20261109/runs/campaign/…`, `lectk_bench2b_round2_20261113/runs/campaign/…` (links em `lectk_bench2b_phase2_20261119/inputs/`) | base de regressão do ambiente (78 comparações) | Marco B |
| 5 entradas principais e 5 de grade comum da etapa 4 (caso 19790205) | `lectk_stage4{a1,b1,b2,b3}_*/run{,2}/case_19790205.nc`, `lectk_stage4a_kerchunk_*/out/common_kerchunk_*`, `lectk_stage4b1_*/lec34/{wekeo_new,cds_2b}_*`, `lectk_stage4b2_*/cmp/arco_new_*`, `lectk_stage4b3_*/cmp/destine_new_*` | idem | Marco B |
| Pastas de referência da etapa 4 (`inputs/`, CSV de referência, proveniência) e `runs/` da fase 2 (CSVs de referência) | nas mesmas pastas | idem | Marco B |
| Manifestos, selos, tabelas por pedido, relatórios, logs, `state/` da 2b | `lectk_bench2b_*/{manifests,state,report,logs,comparisons}` | rastreabilidade dos números da documentação | Marco B |
| Harness da 2b e scripts de comparação (exclusivos, fora do Git) | `lectk_bench2b_*/harness/`, scripts de `lectk_stage4*_*/` (inclui `regress_4d.py`; `run_stage4a_kerchunk.py` importa o harness da rodada 2) | rodar de novo a regressão e as comparações | Marco B |
| Referências kerchunk | `lectk_bench2b_20261109/ref/kerchunk_refs`, `lectk_stage4a1_default_cap_20261123/run/inputs/kerchunk_refs` | a regressão da 4a1 as lê; os bytes da 2b e da 4a diferem | Marco B |
| Captura bruta do diagnóstico da 4a (R-DIAG, 105 MB) | `lectk_stage4a_diag_20261122/data/` | o pedido 989514 foi apagado no NCAR; não se recria | Marco B (decidir descartar ou não) |
| 19 arquivos fora do GitHub e `dev.tgz` | em `lectk_stage4b2_*/code{,_v1}`, `lectk_stage4b3_*/work{,_z}`, `lectk_stage4d_*/dev{,2}`, `lectk_stage4a1_*/code` | rascunhos e versões locais não publicadas | Marco B |
| Ambiente `lorenz-4d` | `/home/danilocs/.conda/envs/lorenz-4d` | ambiente do projeto (D-043) | Marco B |
| Ambiente `lorenz` | `/home/danilocs/.conda/envs/lorenz` | pode servir fora do projeto | Marco B |

**Recriar código das campanhas:** `git archive <commit>` com os commits da §2.

## 8. Pendências

- **INVESTIGAR:** o clone `/p1-swell/danilocs/LorenzCycleToolkit` (fora da lista).
- **Backlog de código** (registrado em [queue.md](../../queue.md)): trocar `drop` por `drop_vars` em
  `src/utils/preprocessing.py` (linhas 428 e 430). É o aviso de depreciação do xarray; não muda resultados.
- **Resumo de transferência:** fora do repositório (D-044).

## 9. Fechamento (D-044)

Aprovado pelo responsável, PR #90.
