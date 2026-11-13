# Frente 2b — emenda 2 ao pré-registro: rodada complementar (13/11/2026)

**Natureza:** registro documental. Torna permanente a emenda 2 ao
[pré-registro v1](2026-11-09-front-2b-stage0-preregistration.md), que executa o passo 3 de
[D-026](../../decisions.md) ([notas de coordenação](2026-11-13-coordination-notes-2a-2b.md)). **Nada executável
muda no repositório:** nenhum código, teste, script, dependência ou CI. O harness, os manifestos e os dados da
rodada ficam fora do Git, no diretório de trabalho do servidor de processamento. `main`, tags e releases
continuam inalterados.

**Base:** `develop` em `20aa986` (PR #65, correção dos clientes NCAR); [pré-voo e emenda 1](2026-11-10-front-2b-preflight-amendment.md);
[resultados da fase 1](2026-11-13-front-2b-phase1-results.md).

O pré-registro v1 e a emenda 1 não são reescritos: o que segue é uma emenda datada, aprovada **antes** de
qualquer download desta rodada.

## Emenda 2

Aprovada pelo responsável em 13/11/2026, antes de qualquer download desta rodada. Executa o passo 3 de D-026.

- **Rotas:**
  - CDS API, com outra chave da conta do responsável, livre de outros trabalhos;
  - NCAR Get a Subset e NCAR kerchunk, com os clientes corrigidos no PR #65;
  - NCAR NCSS: nova medição; falhou na fase 1 por instabilidade do servidor NCAR.

  ARCO, WEkEO e DestinE não rodam (completas na fase 1).
- **Pedido:**
  - as mesmas caixas seladas da fase 1 (`boxes.json`, sha256 `791a06e5…be89c`);
  - 00/06/12/18 UTC; 13 níveis, de 1000 a 50 hPa; 0,25°; T, z, u, v, ω;
  - 3 repetições separadas por ≥ 12 h; um pedido por vez;
  - ordem sorteada com nova semente registrada.
- **Teste mínimo de Greenwich primeiro:** 1 pedido de poucos MB no Get a Subset (purgado depois) e 1 no
  kerchunk, numa caixa que cruza 0°. Se falhar, a rota sai da rodada, com o motivo registrado; as demais
  seguem.
- **Classes de falha, fixadas antes:**
  - falha da fonte (resposta ou erro do serviço);
  - falha de rede;
  - interrupção externa;
  - **falha do cliente:** erro produzido pela lógica do código LECTK ou do harness, sem erro do serviço.

  Falha do cliente não conta para a regra das 2 falhas da fonte, mas 2 seguidas param a rota naquela
  repetição.
- **Limites:** tráfego ≤ 300 GB; duração ≤ 5 dias; disco livre ≥ 1 TB. As demais regras de parada, o modo
  protegido CDS, a retomada automática e a marcação de repetições afetadas seguem o pré-registro v1.
- **Código:** export de `develop` em `20aa986`; o adaptador de barra final da fase 1 é retirado.
- **Ressalva de época:** a rodada é medida em data diferente da fase 1; os números descrevem cada rota e não
  comparam rotas como se fossem simultâneas.
- **Uso de várias chaves CDS:** chaves sob controle do responsável, que confirmou o uso dentro das regras das
  contas.
- **Resultado esperado:** tabela única com tempo típico (mediana) e faixa por caso para as 7 rotas (fase 1 +
  rodada complementar), com a data de cada medição.

## Fora deste registro

A execução da rodada, seus resultados e a fase 2 ficam para registros posteriores. Merge, CI e releitura deste
documento não podem constar dele.
