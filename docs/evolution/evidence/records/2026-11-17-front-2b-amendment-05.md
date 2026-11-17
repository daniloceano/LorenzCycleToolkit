# Frente 2b — emenda 5 ao pré-registro: CDS da repetição 3 durante a espera do Get a Subset (17/11/2026)

**Natureza:** registro documental. Torna permanente a emenda 5 ao
[pré-registro v1](2026-11-09-front-2b-stage0-preregistration.md), que altera o agendamento da rodada complementar
da [emenda 2](2026-11-13-front-2b-amendment-02.md), já ajustada pelas emendas [3](2026-11-13-front-2b-amendment-03.md)
e [4](2026-11-16-front-2b-amendment-04.md). **Nada executável muda no repositório:** nenhum código, teste, script,
dependência ou CI. O harness, os manifestos e os dados da rodada ficam fora do Git, no diretório de trabalho do
servidor de processamento. `main`, tags e releases continuam inalterados.

**Base:** `develop` em `5145f7b` (emenda 4).

O pré-registro v1 e as emendas 1 a 4 não são reescritos: o que segue é uma emenda datada, aprovada durante a rodada
complementar e **antes** de qualquer pedido da repetição 3.

## Emenda 5

**Decisão do responsável** em 17/11/2026, durante a rodada complementar (repetição 2 em curso), para reduzir a espera.

**Motivo:** na repetição 3 a ordem sorteada é kerchunk → Get a Subset → NCSS → CDS. O Get a Subset passa a maior
parte do tempo esperando vaga (emenda 3; cerca de 23 h para os 10 casos na repetição 2), e o CDS, último da fila,
seria cortado pelo limite de 5 dias depois de 1 ou 2 casos. O kerchunk não precisa de mudança: a repetição 2 dele
já terminou e ele é o primeiro da repetição 3.

- **O que muda:** na repetição 3, a rota CDS sai da fila do controlador e roda num executor próprio, ao mesmo tempo
  que as demais rotas, em especial durante as esperas de vaga do Get a Subset.
- **O que não muda:**
  - **um download por vez no servidor** (trava global da rodada): os pedidos das rotas não se sobrepõem, e os
    tempos medidos continuam comparáveis;
  - mesmos 10 casos, mesma ordem sorteada pela semente da rodada, mesma chave (nº 9) via `CDSAPI_RC`;
  - mesmo código de pedido, verificação, classes de falha e contagem de tráfego da rodada;
  - mesmas regras de parada: 2 falhas da fonte ou 2 falhas do cliente seguidas param a rota na repetição; arquivo
    STOP; limite de 5 dias da rodada; teto de 300 GB; disco livre mínimo;
  - intervalo de 12 h: o executor só começa 12 h depois do fim da repetição 2 (mesmo início da repetição 3 do
    controlador);
  - as repetições 1 e 2 e as demais rotas da repetição 3 seguem como estavam.
- **Ressalvas:**
  - quando o Get a Subset ganhar vaga enquanto um pedido CDS está em curso, ele espera o fim desse pedido
    (até cerca de 2 h), o que pode atrasar o Get a Subset e o NCSS da repetição 3;
  - a ordem real em que as rotas usaram o servidor na repetição 3 será reconstruída pelos horários dos manifestos e
    aparecerá no relatório.
- **Emenda 4:** sem mudança neste registro. A passada extra continua prevista para depois do fim da rodada; manter
  ou não essa passada é decisão separada do responsável.

## Fora deste registro

A execução, os resultados da rodada e da passada extra e a fase 2 ficam para registros posteriores. Merge, CI e
releitura deste documento não podem constar dele.
