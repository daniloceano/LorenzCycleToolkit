# Frente 2b — fase 1: resultados da campanha e conclusões homologadas (13/11/2026)

**Natureza:** registro documental. Torna permanentes o resultado da campanha da fase 1 da Frente 2b (10–12/11/2026)
e as conclusões por rota que o responsável homologou em 13/11/2026. **Nada executável muda no repositório:** nenhum
código, teste, script, dependência ou CI. O harness, os manifestos, a tabela por pedido e os dados ficam fora do
Git, em `<diretório de trabalho no servidor>`. `main`, tags e releases continuam inalterados.

**Base:**
- `develop` em `8911708327e63f624b4d73c280d595a2528a45c3` (PR #62);
- [pré-registro v1](2026-11-09-front-2b-stage0-preregistration.md);
- [pré-voo e emenda 1](2026-11-10-front-2b-preflight-amendment.md).

**Fora do Git:**
- relatório da fase 1, sha256 `69c64bc4d225bb58c24a786c780accf38ea517e323e392e1d709caa55e21cada`;
- tabela com os 134 pedidos, sha256 `925ac18c3fb4b2a0126e55d1e0eb2e421b195ef2587ae953e9505af7f6913447`.

## 1. Execução

**Identidade da campanha:**
- **Harness:** v5, sha256 `cc1538a0…`. Os 134 pedidos rodaram nessa versão.
- **Código do LECTK:** export de `develop` em `cac34e9`.
- **Caixas:** `791a06e5…`.
- **Semente:** 20261109.

| Repetição | Início → fim (UTC) | Duração |
|---|---|---|
| 1 | 10/11 10:30 → 19:34 | 9,1 h |
| 2 | 11/11 07:34 → 14:11 | 6,6 h |
| 3 | 12/11 02:11 → 06:41 | 4,5 h |

**Condições de execução:**
- **Intervalo entre repetições:** 12 h exatas.
- **Paradas globais:** nenhuma.
- **Interrupções:** nenhuma interrupção externa e nenhuma repetição afetada.
- **Fim da retomada automática:** o vigia agendado foi removido depois do fim da campanha.

**Tráfego:**
- **Sem a DestinE:** 200,9 GB (teto da emenda 1: 420 GB).
- **DestinE:** 133,4 GB (teto próprio: 500 GB), em 22.148 requisições HTTP, abaixo de 150.000 no mês. A campanha
  não cruzou a virada do mês.

## 2. Resultado pré-registrado por rota

As células mostram sucessos / pedidos feitos em cada repetição. As classes são as pré-registradas.

| Rota | Rep. 1 | Rep. 2 | Rep. 3 | Total | Falhas |
|---|---|---|---|---|---|
| CDS ARCO Zarr | 10/10 | 10/10 | 10/10 | 30/30 | — |
| WEkEO HDA | 10/10 | 10/10 | 10/10 | 30/30 | — |
| DestinE EDH | 7/10 | 10/10 | 10/10 | 27/30 | 2 de rede, 1 da fonte |
| CDS API | 4/6 | 0/2 | 0/2 | 4/10 | 6 da fonte |
| NCAR kerchunk | 0/2 | 3/8 | 0/2 | 3/12 | 9 da fonte |
| NCAR Get a Subset | 0/2 | 1/4 | 1/4 | 2/10 | 8 da fonte |
| NCAR NCSS | 0/3 | 0/4 | 0/5 | 0/12 | 6 da fonte, 6 de rede |

CDS, kerchunk, Get a Subset e NCSS pararam nas 3 repetições pela regra das 2 falhas consecutivas da fonte. Como
nenhuma repetição teve parada global, as 7 rotas são documentáveis pela regra de validade (§2.10 do pré-registro).

## 3. Métricas das rotas com sucesso

Excesso = bytes recebidos ÷ bytes úteis (grade comum em float32).

| Rota | Tempo total da rota (rep. 1 / 2 / 3) | Recebido por repetição | Excesso |
|---|---|---|---|
| ARCO | 76 / 87 / 55 min | 23,4 GB | 10,1× |
| WEkEO | 48 / 13 / 8 min | 1,28 GB | 0,55× (arquivo comprimido) |
| DestinE | 110 / 73 / 78 min | 40,9 / 46,3 / 46,3 GB | 18–20× |
| CDS (4 pedidos) | 4,9 min | 0,36 GB | 0,58×; fila ≈ 0,3 s |

**Degradação:** calculei o τ de Kendall entre a posição do pedido na sequência e os segundos por GB útil, para
ARCO, WEkEO e DestinE em cada repetição. Não houve tendência significativa: todos com p ≥ 0,11, e os sinais variam
entre repetições.

## 4. Causas das falhas (análise posterior, feita depois de ver os resultados)

Esta seção explica as falhas e não substitui as classes pré-registradas da §2.

- **Kerchunk e Get a Subset, caixas que cruzam Greenwich:** todos os casos que cruzam falharam e todos os outros
  passaram (verificado). 7 dos 10 casos cruzam Greenwich.

  | Rota | Cruza Greenwich | Não cruza |
  |---|---|---|
  | kerchunk | 9 de 9 falharam | 3 de 3 passaram |
  | Get a Subset | 8 de 8 falharam | 2 de 2 passaram |

  - **Kerchunk (mecanismo verificado pela contagem de requisições):** o cliente divide o pedido em dois segmentos de
    longitude, e cada um baixa o mesmo bloco global. O limite de bytes do harness, de 1,5× a estimativa, que conta
    cada bloco uma vez, se esgota.
  - **Get a Subset:** um arquivo excede o limite de tamanho por arquivo do cliente. O mecanismo exato é inferido.
- **NCSS:** 4 HTTP 503 e 6 timeouts do servidor THREDDS. Em 2 pedidos, as retentativas por timeout esgotaram o
  limite de bytes.
- **DestinE:**
  - 2 respostas truncadas na repetição 1;
  - 1 caso com valores não finitos em z e w na repetição 1, que passou nas repetições 2 e 3 (passageiro).
- **CDS:**
  - 1 rejeição por fila cheia na conta, ocupada por outro trabalho do responsável;
  - depois, 5 bloqueios da guarda, mantida bloqueada por decisão do responsável de 11/11. Nenhum pedido foi
    enviado e nenhum byte foi recebido.

## 5. Conclusões por rota (homologadas pelo responsável em 13/11/2026)

| Rota | Conclusão da fase 1 | Situação |
|---|---|---|
| CDS ARCO Zarr | Confiável: 30/30, ≈ 1–1,5 h por repetição, excesso de 10× | documentada |
| WEkEO HDA | Confiável: 30/30, a mais leve em bytes; tempos das rep. 2 e 3 talvez sob cache | documentada, com ressalva de cache |
| DestinE EDH | Confiável a partir da rep. 2: 27/30, excesso de ≈ 20×; 1 defeito passageiro de dado | documentada |
| CDS API | Rápida quando a fila da conta está livre: 4/4 efetivos | **parcial**: refazer com outra chave |
| NCAR kerchunk | Não funciona com o cliente atual para caixas que cruzam Greenwich; nas outras, excesso de 69× | resultado negativo, com causa no cliente e no harness |
| NCAR Get a Subset | Não funciona com o cliente atual para caixas que cruzam Greenwich | resultado negativo, com causa no cliente |
| NCAR NCSS | Falhou em tudo por instabilidade do servidor (503 e timeouts) | resultado negativo, com causa na fonte |

Nenhuma conclusão elege "a melhor fonte" nem compara produtos ([D-015](../../decisions.md)).

## 6. Ressalvas

1. O limite de 1,5× sobre a estimativa de bytes é uma escolha do harness, não do pré-registro. Ele transformou a
   duplicação do kerchunk em falha da fonte.
2. A classe pré-registrada "falha da fonte" inclui limites do harness e do cliente. A §4 separa as causas, mas como
   análise posterior.
3. Os pedidos foram idênticos nas 3 repetições, então a WEkEO e a CDS podem ter servido resultados guardados em
   cache.
4. Só houve um ponto de medida (um servidor, uma rede), e as repetições ocorreram em horários diferentes.

## 7. Defeitos do LECTK registrados (não corrigidos nesta frente)

1. **Get a Subset (`src/utils/gdex_sources.py`), endereços:** a API do GDEX responde HTTP 301 para endereços com
   barra final, e o cliente recusa redirecionamentos (registrado no [pré-voo](2026-11-10-front-2b-preflight-amendment.md)).
2. **Get a Subset, Greenwich:** o limite de tamanho por arquivo falha para caixas que cruzam Greenwich.
3. **Kerchunk (`src/utils/gdex_sources.py`), Greenwich:** com dois segmentos de longitude, o mesmo bloco global é
   baixado duas vezes, e a estimativa de bytes do plano não conta essa duplicação.

## 8. Pendências

- **Rodada da CDS com outra chave:** antes dela, fixar o pedido (as mesmas caixas seladas), as regras de parada e
  como os resultados se juntam a esta campanha.
- **Correção dos defeitos da §7:** em frente própria; depois dela, possível rodada complementar de kerchunk e
  Get a Subset.
- **Retenção GOV-RET-01 dos dados da campanha:** ≈ 335 GB baixados, no servidor.
- **Fase 2 (cálculos LEC e comparação de produtos):** o pré-requisito (2a.1) está concluído desde 09/11; a
  abertura depende de decisão do responsável.

## Fora deste registro

Merge, CI e releitura deste documento não podem constar dele.
