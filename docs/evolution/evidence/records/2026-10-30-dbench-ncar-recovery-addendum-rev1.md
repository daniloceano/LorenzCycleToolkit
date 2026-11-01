# D-Bench — revisão 1 do adendo de recuperação NCAR: horário de início (30/10/2026)

**Estado:** registrada **antes de qualquer aquisição da fase `recovery1`**, por decisão expressa do responsável (30/10/2026, cerca de 11:10 UTC). Revisa um único ponto do [adendo de recuperação](2026-10-30-dbench-ncar-recovery-addendum.md), que **permanece inalterado** (SHA-256 `abcfebee073471bd06bab5eaabfb53a02bf0b974005788a035f3921f5cad97a8`, selado às 10:53:32 UTC). O selo desta revisão é o SHA-256 dela em `ROOT_REC/manifest/revision1.json`.

## 1. Alteração

- **Antes (adendo §3, "Janela operacional"):** o diagnóstico 19790644 começaria no sábado, 31/10/2026, às 07:00 UTC, num período de menor demanda presumida.
- **Agora:** o diagnóstico 19790644 começa **assim que a campanha primária `campaign_20261030b` terminar**. O fim é verificado como antes: existe `controller_end.json` e nenhum processo da campanha primária está rodando.
- **Motivo:** não há evidência de que a janela noturna nos EUA seja melhor para o GDEX. A espera era só uma hipótese e atrasava a fase em cerca de 18 h. As duas falhas originais aconteceram entre 02:19 e 06:09 UTC (20:19–00:09 em Boulder), já no início da noite local.
- **Registro:** o horário de início de cada caso continua registrado em UTC e na hora de Boulder, como condição operacional.

## 2. O que não muda

Todo o resto do adendo vale sem alteração:

- aquisição fria sem nenhum `.parts` da r1, em diretórios novos e separados;
- diagnóstico de um único caso, 19790644; qualquer resultado diferente de `success` para a fase;
- continuação sequencial pelos demais casos elegíveis, na ordem da r1, só após o sucesso do diagnóstico, com parada após duas falhas NCAR consecutivas;
- prazo de 12 h por caso com margem de parada de 30 min, limites NCAR nativos e tolerâncias;
- teto de 96 h contado do início original da campanha primária, sem reinício;
- código congelado `dc4ed748` sem alteração e análise primária, suplementar e de sensibilidade separadas.

## 3. Procedimento operacional da troca

- **Primeiro lançamento** (`--start-at 2026-10-31T07:00:00Z`): o controlador, que só estava esperando o fim da primária, sem nenhuma aquisição e sem diretório `runs/`, é parado por sinal.
- **Preservação:** os arquivos de estado e os logs desse lançamento vão para `ROOT_REC/logs/launch1_superseded/`, sem edição.
- **Novo lançamento:** o **mesmo** controlador, com o mesmo SHA-256 selado no adendo, é iniciado com `--start-at` igual ao horário de registro desta revisão, ou seja, sem espera adicional após o fim da primária.
