# D-Bench — registro de encerramento técnico (01/11/2026)

**Estado:** ENCERRADA EM 01/11/2026 — **sem PASS da campanha pré-registrada**. Campanha primária CDS × NCAR **INCONCLUSIVA** (1 par válido em 10); fase de recuperação NCAR **FAIL no diagnóstico**; Etapa B (comparação entre métodos NCAR, sem CDS) executada e documentada. **Nenhuma superioridade de velocidade NCAR × CDS demonstrada.** Complementa a [Frente D](../../fronts/front-D-era5-ncar.md), que continua encerrada; a Frente E não foi aberta.

Este registro resume resultados já selados. Os registros de 29/10 ([pré-registro](2026-10-29-dbench-preregistration.md), [dispensa do cancelamento seletivo](2026-10-29-dbench-cds-cancellation-waiver.md)) e de 30/10 ([adendo de recuperação](2026-10-30-dbench-ncar-recovery-addendum.md), [revisão 1](2026-10-30-dbench-ncar-recovery-addendum-rev1.md)) permanecem como foram escritos e refletem as decisões vigentes nas respectivas datas.

## Integração

- PRs #47 (métodos de aquisição, merge `1c8f475`), #48 (runner e relatório da Etapa B, merge `1073cac`) e #49 (documentação e registros selados, merge `642a21d`) em `develop`.
- PR #45 (dispensa de 29/10) integrado depois, sem alteração do conteúdo histórico: merge `013a6fd`.
- `main`, tags, releases e publicações não foram alterados.

## Resultados

### Campanha primária `campaign_20261030b` (30/10/2026, r1)

- **CDS:** 9 aquisições frias com sucesso (pedidos sequenciais, 24–283 min de relógio externo por caso). O caso 19800318 não tem execução CDS registrada na r1.
- **NCAR (NCSS, `--ncar`):** piloto 19800296 com sucesso; 19790644 e 19790166 falharam (timeouts de leitura de 300 s após 4 tentativas, HTTP 503/504 do front end do GDEX); a regra de duas falhas consecutivas parou o NCAR nos demais casos.
- **Resultado:** 1 par válido em 10, abaixo do mínimo de análise (n ≥ 6). Pela regra pré-registrada: **INCONCLUSIVO**, sem comparação de velocidade. A r2 não foi executada (o controlador exigia dez pares equivalentes).

### Fase suplementar `recovery1` (30/10/2026)

- Diagnóstico 19790644 (aquisição NCAR fria, iniciado às 15:42 UTC conforme a revisão 1): **failed** após 1 922 s. Pela regra do adendo, a fase parou sem outras aquisições. Nenhum par suplementar foi formado.

### Etapa B — métodos NCAR (31/10/2026)

Três métodos de acesso ao GDEX `d633000`, uma execução fria por caso e método (30 planejadas). Detalhes, tabela pública e figuras em `docs/source/remote_data_acquisition.rst`.

- **Sucessos frios:** Get a Subset 9/10 (1 não submetido por deduplicação), kerchunk 10/10, NCSS 4/10 (2 falhas e 4 não executados após a regra de parada).
- **Pareado nos mesmos casos:** Get a Subset / kerchunk, razão mediana 0,42 (9 de 9 mais rápidos); Get a Subset / NCSS, 0,33 (3 de 3); kerchunk / NCSS, 0,74 (3 de 4).
- **Equivalência com o arquivo CDS da campanha primária:** campos dentro das tolerâncias em todos os casos com referência; energética LEC (`-t -r`, tolerância por coluna 10⁻⁴ × escala) **equivalente em 4 de 9 casos e não equivalente em 5** (excesso máximo de 1,24 a 6,66 × a tolerância, concentrado em termos de tendência e residuais: ∂Az/∂t, RGz, BAz e outros).
- **Diagnóstico 19790644 (não é resultado da campanha):** a substituição da temperatura GDEX pela CDS torna a energética equivalente; a do geopotencial não. As diferenças de temperatura (≤ 2,4 × 10⁻⁴ K, 1000–800 hPa) são amplificadas em Az e em sua tendência.

## Resultados negativos

1. A hipótese de vantagem prática NCAR/CDS (mediana ≤ 0,80 e IC 95 % < 1) **não pôde ser testada**: não há amostra pareada válida.
2. O acesso NCSS pelo `--ncar` falhou em 4 de 7 aquisições frias de domínios de 1,9 a 3,9 GiB em dois dias distintos.
3. Arquivos NCAR e CDS **não produzem energética equivalente** dentro da tolerância de coluna pré-registrada em 5 de 9 casos, embora os campos estejam dentro das tolerâncias de campo.
4. A hipótese de concorrência (k = 2) não foi avaliada.

## Ressalvas científicas

- Uma execução por caso e método, em um dia e de um local; sem repetições em outras janelas. Contagens de sucesso não estimam taxas gerais de falha.
- População de dez trajetórias de 1979–1980 no Atlântico Sul; não representa a disponibilidade ERA5 contemporânea nem outros domínios.
- A Etapa B não incluiu o CDS; os tempos CDS da campanha primária são de outro dia e outra fila e não foram comparados.
- A atribuição do mecanismo de temperatura aos outros quatro casos não equivalentes é **hipótese não testada**. A magnitude das diferenças fica muito abaixo dos termos de energia discutidos em estudos de ciclones, mas não se deve esperar concordância NCAR × CDS a 10⁻⁴ da escala de cada coluna.
- Os clientes kerchunk e Get a Subset da Etapa B estão fora do toolkit; a integração no toolkit só foi testada offline e seus tempos não foram medidos.
- Tolerâncias, população e critérios não foram alterados após a observação dos resultados.

## Retenções temporárias e pendências

- **Área privada `dbench_private/` (fora do Git):** protocolos selados, ledgers, registros, auditoria, NetCDF, logs e chaves. Retida integralmente até nova decisão GOV-RET-01; contém detalhes de conta e não deve ser publicada.
- **Checkout:** `LEC_Results/` retido integralmente até se comprovar cópia íntegra ou dispensabilidade de cada resultado; `inputs/namelist`, `inputs/track` e `inputs/namelist_WRF_bigua` (não versionado, fora do escopo da D-Bench) preservados sem alteração.
- **Descartes desta etapa:** apenas 17 arquivos não rastreados e 11 modificações do branch local, verificados idênticos byte a byte a `origin/develop` antes da remoção; branch local experimental sem commits exclusivos removido. Nenhum artefato científico apagado.
- **Pendências sem frente aberta:** teste do mecanismo de temperatura nos outros quatro casos; repetições em outras janelas; decisão sobre a integração dos métodos kerchunk e Get a Subset (experimentais, não lançados); referência NCEP-R2 (herdada da Frente B). Nenhuma nova aquisição CDS ou NCAR autorizada.
