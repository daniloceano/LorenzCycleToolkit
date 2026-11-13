# Notas de coordenação — encerramento da 2a e decisões da 2b (10–13/11/2026)

Decisões do responsável transmitidas pela coordenação em 10/11/2026 e tomadas em 13/11/2026. Este registro não altera dados nem medições.

## 1. Frente 2a — encerrada

A Frente 2a (validação externa com Dias Pinto & da Rocha 2011, NCEP/NCAR R1) está **encerrada com o FAIL pré-registrado** (D-021, PR #58): 53 de 89 termos aceitáveis (59,6 %; limiar 60 %). Decisão do responsável em 08/11/2026, reconfirmada em 10/11/2026.

- `queue.md` ainda dizia "decisão de encerramento pendente"; o estado foi corrigido por nota datada, sem apagar o texto anterior do histórico.
- A subfrente 2a.1 (termos de fronteira) está concluída (PR #60, D-022).
- Fica no backlog, **sem abrir**: investigar CA e CE 19–25 % menores no Reg3.

## 2. Frente 2b — emendas e regras

- **Máquina (D-023, 09/11):** execução no servidor de processamento em vez do Mac Mini. O plano antigo não foi encontrado no repositório nem no servidor (busca com cobertura parcial: sem PDFs binários, sem `miniconda3/` e `MPAS-Research*`) e **não foi procurado no Mac Mini**.
- **Tetos de tráfego (D-024, 10/11):** 420 GB para as rotas exceto DestinE; 500 GB próprios da DestinE; máximo combinado de 920 GB. São tetos de tráfego, não de disco; a regra de disco (parar com menos de 1 TB livre) é separada. Consumo final: 200,9 GB e 133,4 GB, respectivamente.
- **Get a Subset (D-025a):** os defeitos do cliente do LECTK (HTTP 301 e caixas que cruzam Greenwich) são falhas do nosso cliente, não da fonte. Os resultados da campanha para essa rota não descrevem a fonte. Uma nova medição, depois da correção pelo agente, exige emenda ao pré-registro aprovada pelo responsável.
- **Conta CDS compartilhada (D-025b):** o CDS limita pedidos ativos por usuário; tempos de espera medidos com a conta em uso por outro trabalho podem estar inflados.

## 3. Ordem dos próximos passos (D-026, 13/11)

1. Registrar estas decisões.
2. Corrigir os 3 defeitos dos clientes NCAR (`src/utils/gdex_sources.py`) em frente própria, pelo agente.
3. Uma única rodada complementar — CDS com outra chave e as rotas NCAR corrigidas — com as mesmas caixas seladas, sob uma só emenda ao pré-registro aprovada antes de qualquer download.
4. Fase 2 (cálculos LEC e comparação de produtos).
5. Retenção (GOV-RET-01) dos dados da campanha só depois da fase 2, que depende deles.

Nenhum desses passos abre automaticamente uma frente.
