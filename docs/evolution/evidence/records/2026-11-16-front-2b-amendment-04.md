# Frente 2b — emenda 4 ao pré-registro: passada extra só do CDS (16/11/2026)

**Natureza:** registro documental. Torna permanente a emenda 4 ao
[pré-registro v1](2026-11-09-front-2b-stage0-preregistration.md), que acrescenta uma passada à rodada complementar
da [emenda 2](2026-11-13-front-2b-amendment-02.md), já ajustada pela [emenda 3](2026-11-13-front-2b-amendment-03.md).
**Nada executável muda no repositório:** nenhum código, teste, script, dependência ou CI. O harness, os manifestos e
os dados da rodada ficam fora do Git, no diretório de trabalho do servidor de processamento. `main`, tags e releases
continuam inalterados.

**Base:** `develop` em `b13f422` (emenda 3).

O pré-registro v1 e as emendas 1 a 3 não são reescritos: o que segue é uma emenda datada, aprovada durante a rodada
complementar e **antes** de qualquer pedido desta passada.

## Emenda 4

Aprovada pelo responsável em 16/11/2026, durante a rodada complementar e antes de qualquer pedido desta passada.

**Motivo:** a rota CDS não obteve medições válidas na repetição 1 (um download de 0 bytes travou a proteção do CDS
até 16/11) e é a última da repetição 2; a repetição 3 deve ser cortada pelo limite de 5 dias. Sem esta passada, a
rota usada hoje pelo LECTK ficaria quase sem estimativa de tempo.

- **Passada:** depois do fim da rodada (limite de 5 dias ou conclusão), uma repetição adicional **somente** do CDS:
  - os 10 casos, mesmas caixas seladas e pedido comum;
  - mesma chave (nº 9) via `CDSAPI_RC`;
  - ordem sorteada com nova semente registrada;
  - um pedido por vez.
- **Início e limites:** começa pelo menos 12 h depois do último pedido CDS da rodada; duração máxima de 24 h; conta
  no teto de 300 GB da rodada.
- **Proteção do CDS:** se travar por download vazio, o agente pode liberá-la **uma** vez, depois de conferir,
  somente leitura, que a conta não tem pedidos pendentes, com registro selado. Se travar de novo, a passada para.
- **Demais regras:** classes de falha e demais regras como na emenda 2.
- **Resultados:** aparecem como "passada extra" e com a data, e não substituem as repetições da rodada.

## Fora deste registro

A execução da passada, seus resultados, os resultados da rodada e a fase 2 ficam para registros posteriores. Merge,
CI e releitura deste documento não podem constar dele.
