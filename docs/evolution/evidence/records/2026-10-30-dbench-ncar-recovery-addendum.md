# D-Bench — adendo ao pré-registro: fase suplementar de recuperação NCAR (30/10/2026)

**Estado:** escrito e selado **antes de qualquer aquisição da fase de recuperação**, com autorização expressa do responsável (30/10/2026). O selo é o SHA-256 deste arquivo gravado em `ROOT_REC/manifest/addendum.json` antes do primeiro pedido. Complementa o [pré-registro de 29/10/2026](2026-10-29-dbench-preregistration.md) e **não altera** nenhum critério dele. Pela regra do pré-registro, um desvio necessário é registrado como desvio, com motivo, e os resultados afetados são reportados separadamente. Este adendo é esse registro.

## 1. Falhas NCAR originais e regra de parada

Campanha primária: `campaign_20261030b`, r1, código congelado em `dc4ed748`.

| Ordem | Caso | Janela UTC | Resultado | Causa observada |
|---|---|---|---|---|
| 0 | 19800296 (piloto) | 02:03–02:19 | `success` | 15 respostas, 532 MB |
| 1 | 19790644 | 02:19–03:00 | `failed` | `NcarTransferError`: `T_19790724_1` com 4 × `TimeoutError` de leitura (300 s), antes houve HTTP 504/503 do nginx do GDEX em `W_19790731_0`; 11/80 respostas verificadas |
| 2 | 19790166 | 05:22–06:09 | `failed` | `NcarTransferError`: `T_19790304_0` com 4 × `TimeoutError` de leitura (300 s); 10/90 respostas verificadas |

- **Origem das respostas:** o servidor THREDDS NCSS `tds.gdex.ucar.edu` (dataset `d633000`), recortando server-side arquivos diários completos para domínios grandes (177 × 501 e 177 × 525 pontos, cerca de 3,6–3,9 GiB por caso). Uso de CPU desprezível e respostas recebidas íntegras: nada indica falha do cliente. A causa provável é lentidão ou carga do serviço, mas **isso não foi demonstrado**.
- **Regra de parada:** as duas falhas consecutivas acionaram a regra `NCAR para após duas falhas consecutivas` (§6) às 06:09:58 UTC. A partir daí o controlador pulou o NCAR dos demais casos e seguiu só com o CDS, como pré-registrado.
- **Consequência:** sem esta fase, a r1 termina com no máximo **1 par válido em 10**, insuficiente para a análise de §8 (n ≥ 6).

## 2. Natureza da fase

- **Fase suplementar `recovery1`.** Não substitui a r1 e não é a r2.
- **A r1 primária fica intacta:** registros, `campaign.jsonl`, `compare/`, diretórios `runs/`, `.parts` e ledger CDS de `campaign_20261030b` não são alterados, movidos nem reprocessados. As falhas originais continuam contando como falhas da r1.
- **Sem novas aquisições CDS:** cada par de recuperação usa a aquisição CDS fria primária da r1, sem modificação.
- **Código:** o mesmo snapshot congelado (`campaign_20261030/code`, commit `dc4ed748`). A função `run()` do runner é chamada sem alteração. Só um controlador externo novo (`launcher/run_recovery.py`) define a raiz, a ordem e a parada desta fase.

## 3. Casos, diretórios e aquisição fria

- **Casos elegíveis:** todo caso da r1 com CDS `success` e sem NCAR `success`. São as duas falhas (19790644, 19790166) e os casos cujo NCAR foi pulado pela regra de parada, contanto que o CDS deles tenha terminado em `success` ao fim da campanha primária. Um caso com CDS sem sucesso não entra, porque não haveria par.
- **Raiz nova:** `ROOT_REC = dbench_private/campaign_20261030b_ncar_recovery1`.
  - Aquisição: `ROOT_REC/runs/<caso>/ncar/r1/`, criado de forma exclusiva pelo `run()` congelado e recusado se já existir.
  - Registro do runner: `ROOT_REC/records/<caso>_ncar_r1.json`.
  - Registro suplementar: `ROOT_REC/records/<caso>_ncar_recovery1.json`, com fase, tentativa, janela, hora em Boulder e referência ao registro original falho com SHA-256.
- **Aquisição fria:** diretório novo e vazio. **Nenhum `.parts` da r1 é copiado, ligado ou lido.** O contador `reused_responses` da proveniência deve ser 0, e uma recuperação com valor diferente é rotulada como não fria e excluída das medidas de tempo.
- **Ordem:**
  1. **Diagnóstico:** um caso só, **19790644**, a primeira falha e o menor dos dois domínios que falharam.
  2. **Continuação:** só se o diagnóstico terminar em `success`. Depois os demais casos elegíveis, em sequência, na ordem pré-aleatorizada da r1 (`dbench.py plan`).
- **Janela operacional:** o diagnóstico começa num período de menor demanda presumida, sábado 31/10/2026 às 07:00 UTC (01:00 em Boulder, onde fica o GDEX). O horário é registrado como **condição operacional**, não como garantia de sucesso. A continuação roda em sequência sem escolha de janela, com o horário registrado em cada caso.
- **Concorrência:** a fase só começa depois do fim do controlador primário, verificado por `controller_end.json` e pela ausência de processos. Um caso por vez.

## 4. Configuração científica e limites mantidos

Valem sem alteração:

- namelist (SHA-256 `bc24c204…9d0c7`), campos, níveis, grade, passo de 3 h, regra de área e `-t -r`;
- tolerâncias de §5 e regra das colunas do CSV;
- limites NCAR nativos: 8 GiB por execução, 4 tentativas, timeout de 300 s;
- prazo de 12 h por caso e parada externa com margem de 30 min;
- regra de espaço livre e definição de `success`, `censored` e `failed`.

Timeouts, tentativas e prazos **não são aumentados**, nem automaticamente nem manualmente.

## 5. Teto de 96 h

O teto é o da campanha primária, contado de `campaign_20261030b/records/campaign_start.json`, **sem reinício**. Um caso de recuperação só começa se o tempo decorrido mais 12 h 30 min ainda couber nas 96 h. Caso contrário a fase para e os pares faltantes ficam registrados como não concluídos.

## 6. Paradas da fase

- **Diagnóstico:** se falhar ou for censurado, por qualquer motivo, a fase para. Nenhuma outra aquisição de grande volume é feita, e o resultado é apresentado ao responsável.
- **Continuação:** para após **duas falhas NCAR consecutivas**, a mesma regra de §6.
- **Outras paradas:** teto de 96 h, espaço insuficiente ou interrupção manual.
- **Sem repetição:** um caso de recuperação falho não é repetido nesta fase.

## 7. Comparação e análise

- **Pares de recuperação:** CDS r1 primário × NCAR `recovery1`. São avaliados pelo `compare()` congelado, com as mesmas tolerâncias, numa visão própria `ROOT_REC/pairs/` (ligações simbólicas somente leitura para os dois diretórios). O resultado vai para `ROOT_REC/pairs/compare/` e não toca `campaign_20261030b/compare/`.
- **Separação na análise final:**
  1. **Primária (r1):** o `summary` pré-registrado sobre `campaign_20261030b`, sem alterações. Hoje dá 1 par válido e "nenhuma comparação de velocidade".
  2. **Suplementar:** pares formados com aquisição NCAR de recuperação, sempre rotulados assim. Equivalência por caso e razões de tempo de aquisição reportadas **à parte**.
  3. **Conjunta:** razão mediana sobre pares primários e de recuperação juntos, só como análise **de sensibilidade**, explicitamente rotulada, nunca como resultado primário.
- **Limites de comparabilidade:** as medições NCAR de recuperação foram feitas em outro dia, em outra janela e após uma falha conhecida do serviço. O par mistura condições de serviço de momentos diferentes. A razão NCAR/CDS desses pares descreve essas condições e não equivale a uma medição pareada simultânea.
- **Falhas da recuperação:** reportadas como tais, sem imputação.
