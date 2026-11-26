# Frente S — atualização do estado do projeto (26/11/2026)

**Natureza:** só documentação em `docs/evolution/`; nenhum código, teste, dependência ou CI mudou. Base: `develop` em `5c1fd96` (PR #90).

## 1. Problema

A fila e o panorama estavam desatualizados: a linha 2b da fila ainda começava por "Fase 1 concluída… pendentes", embora a 2b tenha sido concluída em 19/11/2026 (a correção autorizada nesse dia não foi executada); `README.md` e `overview.md` tinham corte em 07/11, sem a linha do tempo de 08/11 a 25/11.

## 2. Decisões do responsável registradas

1. **2b concluída em 19/11/2026** (PRs #61–#74, D-023 a D-028; fase 2 e retenção encerradas, retenção absorvida pela etapa 4). Registrada por nota na fila.
2. **"Operacionalizar o LEC"** é objetivo de outro projeto, fora do Projeto Evolução (26/11; D-045).
3. **Frente O** (controle remoto do agente pelo celular e arquivo de instruções do agente, aberta em 19/11) **descontinuada por enquanto**; sem registro anterior no repositório e nada integrado (26/11; D-046).
4. **Etapa 3 não aberta:** o bloco preparado em 19/11 nunca foi usado. Ao abrir, a escolha da fonte deve usar tempos de pedidos novos: o "~1 min" da WEkEO era cache; pedido novo ≈ 24 min ([4b·1](2026-11-23-stage-4b1-wekeo.md)). Registrada por nota na fila.
5. **Regra de modelo e esforço nos prompts para agentes** (19/11, já aplicada nas instruções do Project): primeiro execução correta, depois menor uso de créditos; faixas simples / intermediária / complexa; só opções disponíveis, sem fixar fornecedor. Resolve o item de D-019 sobre recomendar Codex e Claude Code (D-047).

## 3. Backlog acrescentado (pendências da etapa 4), não aberto

- `drop` → `drop_vars` em `src/utils/preprocessing.py` (já estava no backlog; mantido);
- robustez das tentativas de download das referências do kerchunk (código inalterado desde `f7569f2`);
- investigar o clone do repositório no servidor de processamento, que tem um arquivo não rastreado.

## 4. Estado resultante

1, 2a, 2a.1, 2b e 4 concluídas; 3 e 5 não abertas; 6 bloqueada até autorização específica. Histórico preservado: correções por nota datada, sem reescrita dos textos anteriores.

## 5. Fechamento (D-044)

Aprovado pelo responsável, PR #91.
