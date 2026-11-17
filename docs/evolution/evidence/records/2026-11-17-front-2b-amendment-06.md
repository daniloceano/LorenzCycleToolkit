# Frente 2b — emenda 6 ao pré-registro: cancelamento da passada extra do CDS (17/11/2026)

**Natureza:** registro documental. Torna permanente a emenda 6 ao
[pré-registro v1](2026-11-09-front-2b-stage0-preregistration.md), que cancela a passada extra da
[emenda 4](2026-11-16-front-2b-amendment-04.md). **Nada executável muda no repositório:** nenhum código, teste,
script, dependência ou CI. O harness, os manifestos e os dados da rodada ficam fora do Git, no diretório de trabalho
do servidor de processamento. `main`, tags e releases continuam inalterados.

**Base:** `develop` em `1b0ab74` (emenda 5).

O pré-registro v1 e as emendas 1 a 5 não são reescritos: o que segue é uma emenda datada, aprovada durante a rodada
complementar e **antes** de qualquer pedido da passada extra.

## Emenda 6

**Decisão do responsável** em 17/11/2026, durante a rodada complementar (repetição 2 em curso).

**Motivo:** a passada extra existia porque o CDS ficaria quase sem medições: a repetição 1 foi bloqueada pela
proteção do CDS e a repetição 3 seria cortada pelo limite de 5 dias. Com a [emenda 5](2026-11-17-front-2b-amendment-05.md),
o CDS da repetição 3 roda durante as esperas do Get a Subset e deve completar os 10 casos; a repetição 2 do CDS já
tem medições válidas. A passada extra acrescentaria de 12 a 36 h depois do fim da rodada sem ganho proporcional.

- **O que muda:** a passada extra só do CDS **não será executada**. Nenhum pedido dela foi feito.
- **Fechamento da rodada:** o fechamento previsto na emenda 4 (remover as linhas da rodada do agendador, selar
  manifestos e a tabela por pedido, gerar o resumo) continua o mesmo, mas passa a ser feito pelo agente depois do
  fim da rodada e do executor da emenda 5, em vez de disparado pelo fim da passada. O registro da passada indicará
  "cancelada pelo responsável (emenda 6)".
- **Ressalva:** se a repetição 3 do CDS não completar (falhas ou limite de 5 dias), o CDS ficará com menos medições
  do que a passada extra daria. Isso aparecerá no relatório; reabrir uma passada exigiria nova decisão do
  responsável.

## Fora deste registro

A execução, os resultados da rodada e a fase 2 ficam para registros posteriores. Merge, CI e releitura deste
documento não podem constar dele.
