# Frente 2b — emenda 3 ao pré-registro: espera de vaga no Get a Subset (13/11/2026)

**Natureza:** registro documental. Torna permanente a emenda 3 ao
[pré-registro v1](2026-11-09-front-2b-stage0-preregistration.md), que ajusta a
[emenda 2](2026-11-13-front-2b-amendment-02.md) (rodada complementar). **Nada executável muda no repositório:**
nenhum código, teste, script, dependência ou CI. O harness, os manifestos e os dados da rodada ficam fora do Git,
no diretório de trabalho do servidor de processamento. `main`, tags e releases continuam inalterados.

**Base:** `develop` em `be1675e` (emenda 2); clientes NCAR corrigidos no PR #65 (`20aa986`).

O pré-registro v1 e as emendas 1 e 2 não são reescritos: o que segue é uma emenda datada, aprovada **antes** de
qualquer download do Get a Subset nesta rodada e antes do início da campanha.

## Motivo

- A conta NCAR aceita no máximo 8 pedidos Get a Subset abertos. O cliente LECTK corrigido no PR #65 não submete
  quando os pedidos listados mais os novos passariam de 7, e conta também os pedidos em estado "Set for Purge".
- Um pedido purgado continua listado até o NCAR processar a purga, de forma assíncrona: em 13/11, seis pedidos
  purgados continuavam listados mais de 2 h depois; na fase 1, pedidos purgados saíram da lista em menos de ~39 h.
- 7 dos 10 casos selados cruzam 0° e precisam de 2 pedidos cada: são 17 pedidos por repetição. Sem espera, cada
  recusa contaria como falha da fonte, e a rota pararia depois de cerca de 3 casos em cada repetição.
- Os seis pedidos que a fase 1 deixou na conta foram purgados em 13/11, com autorização do responsável; a purga
  está registrada no diretório de trabalho da rodada.

## Emenda 3

Aprovada pelo responsável em 13/11/2026 (opção "esperar vaga"; limite de duração mantido em 5 dias).

- **Espera de vaga:** antes de cada caso do Get a Subset, o harness consulta a lista de pedidos da conta. Se os
  listados mais os pedidos do caso (1, ou 2 quando a caixa cruza 0°) passarem de 7, espera e consulta de novo a
  cada 10 min.
- **Registro:** o tempo de espera fica registrado à parte, por caso. Não entra no tempo medido (fila, download e
  total) e não é falha.
- **Limite da espera:** 36 h por caso. Esgotado o limite, o caso fica **incompleto**: não conta para nenhuma regra
  de falha, não é refeito na retomada, e a rota segue para o caso seguinte.
- **Ordem:** durante a espera, a campanha fica parada; a ordem sorteada e a regra de um pedido por vez não mudam.
  O arquivo de parada e o limite de duração continuam valendo durante a espera.
- **Duração:** o limite de 5 dias da emenda 2 é mantido. Se esgotar, a campanha para como no pré-registro, e o
  que não foi medido fica registrado como não realizado.
- **Teste de Greenwich do Get a Subset:** também espera vaga, com o mesmo limite; a campanha só começa depois
  dele. Se falhar, a rota sai da rodada, com o motivo registrado, e as demais seguem (como na emenda 2).
- **Consultas à lista:** não baixam dados; o tráfego é desprezível e não entra no teto de 300 GB.
- **Demais regras:** inalteradas (emenda 2 e pré-registro v1).

**Ressalva:** a medição do Get a Subset passa a depender do ritmo de limpeza da lista no NCAR; o número de casos
medidos por repetição pode ficar abaixo de 10, e a rodada pode terminar pelo limite de 5 dias.

## Fora deste registro

A execução da rodada, seus resultados e a fase 2 ficam para registros posteriores. Merge, CI e releitura deste
documento não podem constar dele.
