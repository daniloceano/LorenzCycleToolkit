# D-Bench e diagnóstico pós-encerramento — limpeza GOV-RET-01 executada (01/11/2026)

**Estado operacional:** D-Bench e diagnóstico energético pós-encerramento **ENCERRADOS — RETENÇÃO VALIDADA**.
Frente D original: **manutenção de retenção pendente no ambiente original** (inacessível a partir do servidor; ver
abaixo). Este registro é posterior aos de [encerramento](2026-11-01-dbench-closure.md) e do
[diagnóstico](2026-11-01-dbench-postclosure-energy-diagnostic.md), que permanecem inalterados. Nenhum resultado
científico foi modificado; a D-Bench não foi reaberta e a Frente E não foi aberta.

## Base e procedimento

- Auditoria somente leitura com SHA-256 de todos os arquivos legíveis (5 365 nos conjuntos D-Bench, 21 733
  temporários de pytest) e manifesto de exclusão fixado (sha256
  `379896b554a96d87b1d8f0180409e0563f0fa77d51f7840b8718175994bacb37`). Autorização do responsável para S1–S3 e
  para C1–C6 condicionados à demonstração de dispensabilidade.
- Executor com conferência imediatamente antes de cada exclusão: raiz autorizada, sem seguir links, arquivo
  regular, sem hardlink, tamanho e SHA-256. Registro durável por arquivo, fora do Git.
- Nenhum download, aquisição CDS/NCAR ou benchmark. As únicas computações foram reconstruções determinísticas de
  híbridos para provar recriabilidade, feitas fora das áreas preservadas.

## Execução

| Grupo | Conteúdo | Decisão | Arquivos | Bytes | Prova |
|---|---|---|---:|---:|---|
| S1 | partes Get a Subset duplicadas (diagnóstico de 30/10) | **excluído** | 15 | 682 598 311 | idênticas byte a byte à cópia mantida na Etapa A, reconferida antes |
| S2 | cache Python | **excluído** | 1 | 10 085 | recriável |
| S3 | temporários de pytest das suítes D-Bench | **excluído** | 21 733 | 4 454 684 081 | recriáveis; suítes integradas nos PRs #47–#49 com CI verde |
| C1 | 22 híbridos do diagnóstico pós-encerramento | **excluído** | 22 | 22 963 902 178 | 3/3 reconstruções (base kerchunk, base CDS, troca t+z) com SHA-256 idêntico; entradas retidas; saídas LEC e logs fixados em manifesto complementar; 22 links registrados e removidos |
| C2 | partes brutas Get a Subset da Etapa B | **retido** | 305 | 33 952 311 748 | proveniência exclusiva; pedidos GDEX purgados; o NetCDF montado deriva delas, mas não as substitui |
| C3 | híbridos t/z de 19790644 (Etapa B) | **excluído** | 2 | 3 186 924 805 | 2/2 reconstruções com SHA-256 idêntico e estatísticas idênticas; entradas retidas; 2 links registrados |
| C4 | partes brutas da Etapa A | **retido** | 65 | 7 399 330 289 | mesmo critério de C2 |
| C5 | diagnóstico kerchunk de 30/10 (2 NetCDF, 380 partes) | **excluído** | 382 | 8 545 088 482 | campos bit a bit iguais aos da Etapa B kerchunk (tempo igual, só a unidade de codificação difere); 380/380 partes iguais às fatias do NetCDF; ledgers por chunk, `refs/` e README mantidos |
| C6 | piloto NCAR da tentativa a | **excluído** (só o NetCDF) | 1 | 235 588 133 | campos bit a bit iguais ao piloto da tentativa b; registros e ledger mantidos |
| C6 (erro do manifesto) | 4 amostras do código congelado | **excluído do descarte** | 4 | 35 359 062 | pertencem ao código congelado verificado; classificação incorreta corrigida antes da execução |
| Scratchpads de sessões anteriores | patches, árvores, builds, sondas | **excluído** | 9 184 | 1 442 370 879 | 8 002 cópias de blobs Git ou de arquivos preservados, 604 temporários; os 578 exclusivos foram arquivados em pacote verificado (578/578) antes da exclusão |

**Total:** 31 340 arquivos, 41 511 166 954 bytes (38,66 GiB), **0 falhas**. O espaço medido no volume de dados e no
disco de sistema caiu cerca de 40,3 GB no total; a medida inclui atividade concorrente de outros usuários.

| Conjunto | Inicial | Excluído | Remanescente |
|---|---:|---:|---:|
| Área privada D-Bench | 90,26 GiB | 11,78 GiB | 78,48 GiB |
| Diagnóstico pós-encerramento | 21,42 GiB | 21,39 GiB | 0,04 GiB |
| Temporários de pytest | 4,15 GiB | 4,15 GiB | 0 |
| Scratchpads anteriores | 1,34 GiB | 1,34 GiB | 0 (exclusivos arquivados, 63 MiB) |

## Retenções e cópia independente

Tudo o que foi retido tem **cópia independente num dispositivo físico distinto do mesmo servidor** (não é cópia
externa). São 4 729 arquivos, 78,5 GiB, com SHA-256 conferido arquivo a arquivo, sem divergências, e leitura de
restauração testada num NetCDF CDS. Credenciais não foram copiadas.

| Retido | Tamanho | Classe | Motivo |
|---|---:|---|---|
| 9 NetCDF CDS de referência | 9,46 GiB | ARQUIVAR | base de todas as equivalências; só se refazem com nova aquisição |
| NetCDF NCAR montados (Etapa B: Get a Subset, kerchunk, NCSS; Etapa A ×3) | 27,9 GiB | RETER | sustentam a equivalência de métodos, a tabela LEC 4/9 × 5/9 e as entradas dos híbridos |
| Partes brutas Get a Subset (C2, C4) | 38,5 GiB | RETER | proveniência de aquisição não reproduzível (pedidos purgados) |
| Piloto NCAR (tentativa b), partes do recovery1, diagnóstico de tamanho NCSS | 1,6 GiB | RETER | resultados negativos |
| Protocolos, selos, ledgers, registros, auditorias, CSV LEC, referências kerchunk, código congelado com diff local | ≈0,6 GiB | ARQUIVAR | história das decisões e código efetivamente executado |
| Diagnóstico pós-encerramento: scripts, resultados, logs, saídas LEC, manifestos | 0,04 GiB | ARQUIVAR | sustenta a conclusão aprovada; manifesto original 112/112 e complementar de 1 343 entradas |

**Próximo marco de revisão:** próxima versão/tag estável (Marco B). Candidatos naturais nessa revisão: NetCDF NCAR
redundantes entre métodos e partes brutas, se a montagem for então considerada auditada o bastante. Nenhuma tag foi
criada para disparar o marco.

## Fora do descarte

- **Frente D original:** os grupos G2, G5, G6, P2-A, P3-L e P3-L2 não estão no servidor da D-Bench, e não há acesso
  configurado a partir dele ao ambiente original da Frente D (nenhum host remoto configurado). As retenções de 28/10
  continuam como registradas na [Frente D](../../fronts/front-D-era5-ncar.md). G2-idênticas e o resto de G6 seguem
  como candidatos à próxima manutenção, a executar naquele ambiente. A limpeza da Frente D original **não** é
  declarada concluída.
- **Resultados de outra atividade no checkout:** `LEC_Results/` (anterior à D-Bench; diverge da climatologia
  corrigida) e `inputs/namelist_WRF_bigua` ficaram intocados. A propriedade de `LEC_Results/` está por definir.
- **Dois caches `.pyc` ilegíveis** na ferramenta da Etapa B: não copiados; permissões não alteradas.
- **Branches locais:** cinco branches D-Bench, todas ancestrais de `origin/develop` sem commits exclusivos, foram
  removidas do checkout. Branches remotas, o worktree de outra atividade e stashes (nenhum) não foram tocados.
  O clone congelado da campanha foi preservado.
