# Frente 2b — retenção dos dados, nota de correção e documentação do usuário (19/11/2026)

**Natureza:** registro documental das decisões do responsável de 19/11/2026 (D-028). **Nada executável muda:**
nenhum código, teste, script, dependência ou CI. Nada foi apagado. `main`, tags e releases continuam inalterados.

**Base:** `develop` em `b7e95a3` (PR #73).

## 1. Retenção dos dados das campanhas (decisão do responsável, 19/11/2026)

**Decisão:** todos os dados das campanhas da Frente 2b ficam guardados no servidor de processamento até o
encerramento da etapa 4. Isso inclui a fase 1, a rodada complementar e a fase 2.

- Nada é apagado antes disso.
- Depois do encerramento da etapa 4, os dados passam por revisão pela GOV-RET-01, com lista e aprovação
  específica.

**Motivo:** os manifestos e relatórios sustentam os números que irão para a documentação da etapa 4, e o espaço
ocupado é pequeno.

## 2. Espaço ocupado em disco

A medida foi feita com `du`, somente leitura, em 19/11/2026 às 20:51:55 UTC.

| Diretório da campanha | Espaço ocupado |
|---|---|
| Fase 1 (10–12/11) | 31 GB (33.164.367.524 bytes) |
| Rodada complementar (14–18/11) | 28 GB (29.804.646.344 bytes) |
| Fase 2 (19/11) | 175 MB (175.322.074 bytes) |
| **Total** | **≈ 63 GB** |

## 3. Correção por nota (D-020)

**Correção:** o [registro dos resultados da rodada complementar](2026-11-19-front-2b-round2-results.md) diz "cerca
de 237 GB no diretório da rodada", e isso está incorreto.

- **O que eram os 236,7 GB:** o **tráfego** da rodada, não o espaço em disco.
- **Espaço realmente ocupado:** o medido na §2, 28 GB para a rodada complementar.

**Texto original:** não é alterado.

## 4. Documentação do usuário no Read the Docs (decisão do responsável, 19/11/2026)

**Situação atual:** a documentação do usuário no Read the Docs (`docs/source/remote_data_acquisition.rst` e páginas
relacionadas) se baseia hoje na D-Bench de outubro e está parcialmente desatualizada.

**Decisão:**

- **Atualização:** essa documentação será atualizada na etapa 4, junto com a inclusão das fontes, usando os
  resultados da Frente 2b.
- **Sem aviso provisório:** nenhum aviso provisório será colocado antes disso.

## Fora deste registro

Merge, CI e releitura deste documento não podem constar dele.
