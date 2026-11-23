# Etapa 4b·x — namelist do LEC com `--cdsapi`/`--wekeo` e conferência dos 34 × 35 horários (23/11/2026)

**Natureza:** registro de implementação (Tarefa A, correção de código) e de investigação (Tarefa B, só leitura).
Nenhum download real. Nenhum valor de saída do LEC muda. `main`, tags e releases continuam inalterados.

**Base:** `develop` em `b4d5d17` (PR #79, conferido antes da branch); `main` em `48fe8e2`.

## Resumo em linguagem simples

1. **Namelist.** Com `--cdsapi` ou `--wekeo`, o programa baixava e preparava os dados com o namelist ERA5 próprio da
   rota (`inputs/namelist_ERA5-cdsapi`), mas, na hora de calcular o LEC, lia de novo o `inputs/namelist`. Se os dois
   fossem diferentes, o cálculo falhava depois do download. Agora o namelist é escolhido uma vez e o mesmo arquivo
   serve para preparar os dados e para o LEC. Com arquivo local e com `--ncar` nada mudou.
2. **34 × 35 horários.** O comando principal pede 34 horários e o pedido selado da 2b tinha 35. O 35º (27/03 00 UTC)
   está **depois** da última posição do track (26/03 21 UTC). O comando principal não perde nenhuma posição que o
   track tem. Nada foi mudado.

## 1. Tarefa A — namelist do LEC

**Defeito (registro da 4b·1, §7).** `prepare_data` trocava para `inputs/namelist_ERA5-cdsapi` com `--cdsapi` ou
`--wekeo`; `run_lec_analysis` relia sempre `inputs/namelist`. Com os dois diferentes, o LEC falhava depois do download
(`KeyError: 'time'`, ou `KeyError` do nome da variável, conforme o que difere).

**`--ncar`:** não tem o problema. `prepare_data` e `get_ncar_data` recebem `inputs/namelist` e o LEC o relê: é o mesmo
arquivo. Mantido como estava.

**Correção (mínima).**

- `preprocessing.namelist_for_run(args)`: devolve `inputs/namelist_ERA5-cdsapi` com `--cdsapi` ou `--wekeo`; nos
  demais casos, `inputs/namelist` (ou seu `.default`, como antes). `prepare_data` usa a mesma constante.
- `main` escolhe o namelist uma vez, antes de qualquer pedido, e o passa a `prepare_data` e a `run_lec_analysis` (novo
  parâmetro opcional `namelist`; sem ele, `inputs/namelist`, como antes).
- Namelist ausente: `prepare_data` valida o arquivo antes de `get_data`, portanto a mensagem já aparecia antes de
  qualquer pedido; isso foi mantido e agora está protegido por teste.
- Documentação: só `docs/source/configuration.rst` (a seção diz que `inputs/namelist` não é lido nesse modo e inclui
  `--wekeo`). A reescrita completa é a 4c.

**Arquivos:** `lorenzcycletoolkit.py`, `src/utils/preprocessing.py`, `tests/test_run_namelist.py` (novo),
`tests/test_track_alignment.py` (um teste novo, Tarefa B), `docs/source/configuration.rst`, `CHANGELOG.md`.

### 1.1 Testes (Mac Mini, ambiente `lorenz`, offline)

`tests/test_run_namelist.py`: o arquivo de entrada já existe (as rotas o reaproveitam; nada é baixado e qualquer
chamada aos clientes falha o teste), e `inputs/namelist` é a variante ERA5 com outros nomes (`time`, `level`, `T`…).

| Teste | Antes (`b4d5d17`) | Depois |
|---|---|---|
| `--cdsapi` com `inputs/namelist` diferente: o LEC roda | falha (`KeyError`) | passa |
| `--wekeo` com `inputs/namelist` diferente: o LEC roda | falha (`KeyError`) | passa |
| `namelist_for_run` por rota (`--cdsapi`, `--wekeo`, `--ncar`, local) | falha (função inexistente) | passa |
| Preset ausente com `--cdsapi` / `--wekeo`: erro antes do download, sem arquivo criado | passa | passa |
| Arquivo local (`infile`): continua usando `inputs/namelist` | passa | passa |

Os testes de namelist ausente e de arquivo local já passavam antes: protegem o comportamento que a correção não pode
mudar. Suíte completa: 1 falha (`test_installed_metadata_matches_pyproject`, instalação local desatualizada, igual na
`develop` limpa) e 16 erros (fixture do Natural Earth, ficam para o CI), os mesmos de antes; o restante passa, com os
testes novos. Nenhum teste existente mudou de resultado, e o LEC não foi alterado.

## 2. Tarefa B — 34 × 35 horários (só leitura)

**Autorização:** o responsável autorizou, em 23/11/2026, a leitura na Swell (D-034) do track e do código do harness.
Foram lidos: o track `track_19790205.txt` (copiado para fora do repositório; sha256 conferido), e `bench2b.py` e
`boxes/boxes.json` (selado) da 2b. Nada foi alterado, copiado além do track, nem baixado.

**O track (sha256 `9da9d809…`, o selado):** 68 posições a cada 3 h, de 18/03/1979 12 UTC a 26/03/1979 21 UTC.

**Comando principal** (`daily_hours` → `planned_request_times`, `--time-resolution 6`): a cada dia do período, horários
múltiplos de 6 h; o primeiro dia começa na hora do primeiro ponto arredondada para baixo; o último dia termina no
último horário de 6 h que não passa do último ponto. Resultado: **34 horários, de 18/03 12 a 26/03 18 UTC**.

**Harness da 2b:** regra gravada em `boxes.json` (`rule`): "todos os 00/06/12/18 UTC do primeiro horário do track
arredondado para baixo ao último arredondado para cima". O último ponto, 26/03 21, sobe para 27/03 00. Resultado:
**35 horários, de 18/03 12 a 27/03 00 UTC**.

| | Comando principal | Harness da 2b |
|---|---|---|
| Fim do pedido | último horário de 6 h ≤ último ponto (arredonda para baixo) | último ponto arredondado para cima |
| Horários | 34 | 35 |
| Só neste | nenhum | 27/03 00 UTC |
| Posições do track em horários de 6 h cobertas | 34 de 34 | 34 de 34 |

**Evidência (offline, sobre o track real):** os dois conjuntos diferem só em 27/03 00; esse horário não é horário do
track (nenhuma posição do track cai nele). As 34 posições do track em horários de 6 h estão nos dois. No modo
`--track`, com os dados de 35 horários, o alinhamento usa as mesmas 34 posições (`data.sel` nos horários do track), de
modo que o 35º horário não entra no cálculo. Um teste offline (só horários, sem coordenadas) guarda isso.

**Conclusão: 1. O comando principal está certo; o 35º horário do harness estava fora do track. Nada foi mudado.**

**Nuance (informativa, sem defeito):**

- Em `-f` (caixa fixa, usada nas medições da 2b e da 4a.1/4b·1) as linhas do CSV são os horários dos dados, e o
  comando principal termina em 26/03 18, três horas antes do último ponto (21 UTC, que não está na grade de 6 h). O
  harness fecha o período por cima, com um horário depois do track. É uma escolha de arredondamento (para baixo ou
  para cima), não perda de dado do track. Isso explica por que o CSV da fase 2 tinha 35 linhas e o do comando
  principal tem 34.
- Com `--time-resolution 3`, o horário 21 UTC entra e os dois critérios coincidem.

## 3. Pendências para a ARCO (4b·2)

- A rota ARCO deve usar a mesma regra de horários do comando principal (`planned_request_times`), para que `--wekeo`,
  `--cdsapi`, `--ncar` e ARCO peçam os mesmos horários.
- Nenhum namelist próprio: se a ARCO usar nomes diferentes dos do `namelist_ERA5-cdsapi`, ela deve passar pelo mesmo
  mecanismo (`namelist_for_run`); caso contrário a falha corrigida aqui voltaria.
- Em testes comparativos com a 2b, lembrar que a caixa (15° contra 5°), os níveis (37 contra 13) e os horários (34
  contra 35) do comando principal são outros.

## Fora deste registro

ARCO, DestinE, robustez do kerchunk, reescrita da documentação (4c). Aprovação, merge, CI e releitura entram por nota
posterior.
