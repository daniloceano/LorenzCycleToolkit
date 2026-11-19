# Frente 2b — resultados da rodada complementar (emendas 2 a 6) (19/11/2026)

**Natureza:** registro documental. Torna permanentes os resultados da rodada complementar da
[emenda 2](2026-11-13-front-2b-amendment-02.md), com as emendas [3](2026-11-13-front-2b-amendment-03.md),
[4](2026-11-16-front-2b-amendment-04.md), [5](2026-11-17-front-2b-amendment-05.md) e
[6](2026-11-17-front-2b-amendment-06.md). O relatório da rodada foi **aprovado pelo responsável em 19/11/2026**.
**Nada executável muda no repositório:** nenhum código, teste, script, dependência ou CI. O harness, os manifestos,
os dados e o relatório original ficam fora do Git, no servidor de processamento. `main`, tags e releases continuam
inalterados.

**Base:** `develop` em `ee53db6` (emenda 6).

**Estado:** rodada CONCLUÍDA em 18/11/2026 21:14 UTC, antes do limite de 5 dias; fechamento e selagem em 19/11.
Fase 2 não aberta.

## 1. Execução

- **Rotas da rodada:** CDS API, NCAR Get a Subset, NCAR kerchunk e NCAR NCSS. Foram 10 casos selados, 3 repetições,
  chave CDS nº 9 e um download por vez no servidor.
- **Período:** de 14/11 02:30 a 18/11 21:14 UTC.

| Repetição | Início → fim (UTC) |
|---|---|
| 1 | 14/11 02:30 → 15/11 04:33 |
| 2 | 15/11 16:33 → 17/11 20:05 |
| 3 | 18/11 08:05 → 18/11 21:14 |

- **Mudanças de agendamento aprovadas durante a rodada:**
  - **emenda 3:** o Get a Subset espera vaga na lista do NCAR;
  - **emenda 5:** o CDS da repetição 3 rodou num executor próprio, das 08:05 às 11:08 do dia 18, durante as esperas
    do Get a Subset;
  - **emenda 6:** a passada extra só do CDS (emenda 4) foi cancelada sem nenhum pedido.
- **Tráfego total:** 236,7 GB, dentro do teto de 300 GB.
- **Selos (fora do Git, no servidor de processamento):** selo final f0b5bfc2…; relatório d1881585…. O selo final cobre a tabela por pedido (91 linhas), o selo dos
  manifestos, o resumo, o estado da campanha e as cópias do agendador.

## 2. Tabela das 7 rotas (tempo por caso bem-sucedido)

As rotas ARCO, WEkEO e DestinE foram medidas só na fase 1; as outras quatro, na rodada 2. O tempo é o total do
pedido, incluindo fila e verificação. No Get a Subset, o tempo **não** inclui a espera de vaga.

| Rota | Medição | Sucessos | Mediana por caso | Faixa | Recebido por caso (mediana) |
|---|---|---|---|---|---|
| CDS ARCO Zarr | fase 1, 10–12/11 | 30/30 | 6,8 min | 2,6–20 min | 2,5 GB |
| WEkEO HDA | fase 1, 10–12/11 | 30/30 | 1,2 min | 0,4–7,3 min | 0,11 GB |
| DestinE EDH | fase 1, 10–12/11 | 27/30 | 7,8 min | 2,6–19 min | 3,1 GB |
| NCAR kerchunk | rodada 2, 14–18/11 | **30/30** | 16 min | 6,9–33 min | 7,5 GB |
| CDS API, rep. 2 | rodada 2, 16–17/11 | 10/10 | **2,7 h** | 6,5 min–4,4 h | 0,09 GB |
| CDS API, rep. 3 | rodada 2, 18/11 | 10/10 | **30 s** | 15–73 s | 0,09 GB |
| NCAR Get a Subset | rodada 2, 14–18/11 | 13/20 rodados | 6,6 min | 2,6–12,6 min | 1,1 GB |
| NCAR NCSS | rodada 2, 14–18/11 | **0/16** | — | — | — |

Na fase 1, o CDS teve 4/10 sucessos, com mediana de 46 s, o kerchunk 3/12 e o Get a Subset 2/10.

## 3. Resultados e falhas por rota

| Rota | Rep. 1 | Rep. 2 | Rep. 3 | Falhas por classe |
|---|---|---|---|---|
| kerchunk | 10/10 | 10/10 | 10/10 | nenhuma |
| CDS | 0/3 (rota parada) | 10/10 | 10/10 | 1 da fonte (download de 0 B); 2 do cliente (proteção do CDS travada) |
| Get a Subset | 6/9 (1 não rodado) | 6/8 (2 não rodados) | 1/3 (7 não rodados) | 7 da fonte |
| NCSS | 0/10 | 0/3 (rota parada) | 0/3 (rota parada) | 6 da fonte (HTTP 504, resposta incompleta); 10 + 2 de rede (tempo esgotado, DNS) |

- **kerchunk:** nas três repetições, todos os 10 casos passaram, inclusive os que cruzam Greenwich, com tempos
  estáveis entre as repetições. A correção do cliente feita na fase 1 (PR #65) resolveu o problema.
- **CDS:**
  - o tempo é quase todo fila;
  - na repetição 2 a fila ficou entre 1,5 e 4,3 h por caso;
  - na repetição 3, os mesmos pedidos ficaram prontos quase sem fila (cerca de 0,3 s). A explicação mais provável é
    que o CDS reaproveitou resultados guardados dos pedidos idênticos de 1 a 2 dias antes. **[NÃO CONFIRMADO]**
    Portanto a repetição 3 mede o melhor caso, com resultado em cache, e não o tempo típico de um pedido novo. Para
    estimar um pedido novo, vale a repetição 2.
- **Get a Subset:**
  - **onde falha:** as 7 falhas da fonte se concentram em 3 casos. Os casos 19790205 e 19790585 falharam nas três
    repetições com o erro "arquivo maior que o limite". Esse limite por arquivo é calculado pelo próprio cliente
    LECTK; pelo mapeamento selado a falha é classificada como da fonte, mas a causa é a estimativa do cliente. O caso
    19800318 falhou na repetição 1 porque a submissão ficou sem resposta;
  - **por que a rota parou:** essas falhas, duas seguidas, pararam a rota nas três repetições, e por isso 10 casos
    não rodaram;
  - **esperas de vaga:** houve 20 esperas, somando 46 h; 7 delas passaram de 1 h, e as maiores ficaram entre 6 e
    8 h, esperando a limpeza da lista no NCAR.
- **NCSS:** não houve nenhum sucesso. As falhas foram tempo esgotado na repetição 1 e erros do servidor do NCAR nas
  repetições 2 e 3 (HTTP 504 e resposta incompleta).

## 4. Intervenções durante a rodada (todas autorizadas e seladas)

1. **13/11:** purga de 6 pedidos do Get a Subset deixados pela fase 1.
2. **16/11:**
   - liberação da proteção do CDS, travada desde 14/11 por um download de 0 B, depois de conferir que a conta não
     tinha pedidos pendentes;
   - purga de 5 pedidos de casos com falha da repetição 1, porque o harness só purgava casos bem-sucedidos.
3. **17/11:**
   - **emenda 5:** executor próprio para o CDS da repetição 3;
   - **emenda 6:** cancelamento da passada extra. Na ausência da passada, o fechamento foi feito à mão pelo comando
     selado de fechamento da emenda 4, com o estado da passada marcado como "cancelada pelo responsável (emenda 6)".

## 5. Ressalvas

- **Épocas diferentes:** as rotas foram medidas em dias diferentes (fase 1 em 10–12/11, rodada 2 em 14–18/11), e
  a carga dos servidores muda; a comparação entre as 7 rotas é indicativa.
- **CDS:**
  - a repetição 3 provavelmente foi atendida por cache **[NÃO CONFIRMADO]**;
  - a repetição 1 não tem medição válida;
  - uma única falha de download trava a proteção do CDS e, com ela, a rota inteira.
- **Get a Subset:**
  - **limites de pedidos:** com o limite de 7 pedidos listados e a limpeza do NCAR em lote, o tempo real de uso é
    dominado pela espera. Somando as esperas, a rota levou dias por repetição, contra minutos por pedido;
  - **estimativa do cliente:** as falhas "arquivo maior que o limite" vêm da estimativa do cliente LECTK.
- **NCSS:** a rota mostrou-se inviável nas duas campanhas.
- **Teto de tráfego:** a estimativa inicial (298 GB) estava perto do teto de 300 GB; o uso real ficou em 236,7 GB
  porque o Get a Subset e o NCSS pararam cedo.
- **Desvios de procedimento:**
  - a chave 8 foi pulada (conta sem a licença exigida pelo ERA5);
  - a preparação da rodada começou antes do merge do PR #66;
  - houve as duas intervenções de 16/11.
- **Uso do servidor na repetição 3:** com a emenda 5, as rotas da repetição 3 se alternaram no servidor, sempre com
  um download por vez. A ordem real pode ser reconstruída pelos horários dos manifestos.

## 6. Achados para decisão futura (sem ação agora)

- **Cliente Get a Subset:**
  - rever a estimativa do limite por arquivo;
  - não contar "Set for Purge" no limite de 7.
- **Harness:** purgar também os pedidos de casos com falha.
- **Proteção do CDS:** rever o bloqueio total da rota depois de um único download vazio.
- **Rota CDS:** decidir se a fila típica de cerca de 2,7 h por pedido novo é aceitável para o LECTK.

## Fora deste registro

Ficam de fora a fase 2 e a retenção ou limpeza dos dados (cerca de 237 GB no diretório da rodada). Nada foi apagado. Merge, CI e releitura deste documento não podem constar dele.
