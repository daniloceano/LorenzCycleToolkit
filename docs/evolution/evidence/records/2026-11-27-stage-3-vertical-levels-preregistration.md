# Etapa 3 — Níveis verticais: pré-registro da medição simples (Etapa 0)

**Estado:** APROVADO em 27/11/2026; texto congelado antes de qualquer download ou cálculo. Resultados no [registro de resultados](2026-11-27-stage-3-vertical-levels-results.md).
**Base:** `develop` em `754ec2e` (PR #92); `main` em `48fe8e2`. Decisões: D-015, D-016, D-034, D-040, D-042, D-044, D-047.
**Natureza:** medição simples, amostra reduzida (10 casos, Atlântico Sul, 1979–1980). Não é a hipótese K.

## 1. Pergunta

Nos 10 casos da D-Bench, em caixa fixa, quanto muda o ciclo de energia de Lorenz quando o LECTK usa 13 níveis
(grade ARCO) em vez de todos os níveis ERA5 que ele aproveita? Antecedente (LZV, um caso, modo móvel): Ca −5,1 %
na média (29 → 13), Ce −1,8 %, Ck −0,2 %; extensão 32 → 29 < 0,05 %. Orienta, não conclui.

## 2. Caixa (fixada antes de rodar)

**As mesmas caixas seladas da Frente 2b**: retângulo que envolve a trajetória + 5° de margem, um por caso, já
calculadas e seladas por hash antes da 2b.

Por quê: (a) já são regra pré-registrada e independente do resultado; (b) permitem conferir os novos dados contra
os arquivos de 13 níveis retidos da 2b (§6); (c) mantêm a comparabilidade com a 2b. Limitação herdada: as bordas da
caixa coincidem com as dos dados (derivadas de um lado só nas paredes), igual nas duas configurações.

Horários: 00/06/12/18 UTC (6 h), os mesmos da 2b. Execução `-f -r`.

## 3. Conjuntos de níveis e referência

Um único download por caso com os **37 níveis** ERA5; os conjuntos são recortes do **mesmo** arquivo, escritos pelo
**mesmo** script (evita o ruído de caminho de leitura visto no C0 da LZV):

| Conjunto | Níveis | Papel |
|---|---|---|
| L32 | 1000–10 hPa (o LECTK descarta 1–7 hPa) | **Referência**: o que o LECTK faz por padrão com CDS/WEkEO |
| L29 | 1000–50 hPa, todos | Isola extensão × discretização |
| L13 | 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 50 | Grade ARCO / 2b |

Comparação principal: **L32 → L13**, decomposta em L32 → L29 (topo) e L29 → L13 (discretização na camada comum).
O agente confirma as contagens nos arquivos (37 no download; 32/29/13 efetivos no LEC, pelo log) antes do LEC.

## 4. Fonte

Precisa de ≥ 29 níveis: ARCO (13) não serve; DestinE (19, quantizada) excluída. Restam CDS, WEkEO e NCAR.

| Fonte | Pedido novo (registros) | Observação |
|---|---|---|
| **WEkEO** | ≈ 24 min (4b·1, 37 níveis; o "1 min" era cache) | Idêntica ao CDS na 2b (fase 2, grupo {CDS, WEkEO}); volume pequeno |
| NCAR kerchunk | 16 min com 13 níveis (2b); com 37, tráfego ≈ 2–3× maior (estimativa) | Difere do CDS só por arredondamento |
| CDS API | mediana 2,7 h (fila; 6,5 min–4,4 h) | Referência histórica, lenta |

**Recomendação: WEkEO** para os três conjuntos (mesma fonte, mesmo arquivo). Alternativa se a WEkEO falhar: kerchunk
para **todos** os casos (nunca misturar fontes entre conjuntos ou casos); a troca exige aviso ao responsável.

## 5. Grandezas e tolerâncias (fixadas antes de calcular)

Todos os termos do CSV do LEC: reservatórios (Az, Ae, Kz, Ke), conversões (Ca, Ce, Ck, Cz), geração/dissipação
estimadas e resíduos (Gz, Ge, Dz, De por resíduo), fronteiras (BAz, BAe, BKz, BKe, BΦZ, BΦE) e tendências.

Por caso e termo: média temporal da referência e de cada conjunto; Δ médio (média de X − X_ref) e Δ máximo
(maior |X − X_ref| num horário); correlação temporal r; inversão de sinal da média.

Escala de normalização: `S = max(|média_ref|, desvio-padrão temporal_ref)` (evita % infinitos em fronteiras e
resíduos próximos de zero). Diferença normalizada `d = |Δ médio| / S`.

| Classe | Critério (por termo, mediana dos 10 casos) |
|---|---|
| a — numérica | d ≤ 10⁻⁶ |
| b — pequena | 10⁻⁶ < d ≤ 2 % e sem inversão de sinal |
| c — moderada | 2 % < d ≤ 10 % |
| d — relevante | d > 10 %, ou inversão de sinal da média em ≥ 2 casos, ou r < 0,9 em ≥ 2 casos |

Relatar também a faixa (mín.–máx.) entre casos. Conclusão escrita como amostra reduzida; sem generalizar para
outras regiões, estações ou modo móvel.

## 6. Controles

- **C1 (repetição):** L32 e L13 de um caso processados duas vezes → CSV idêntico byte a byte. Se não, |Δ| relativo
  ≤ 10⁻¹² (ordem do ruído da LZV, 7 × 10⁻¹⁶). Falha de C1 = parar e informar antes de interpretar.
- **C2 (recorte):** L37 reescrito pelo script, sem remover níveis, vs. L32 → mesmo CSV (o recorte não altera o LEC).
- **C3 (fonte, conferência):** L13 novo vs. arquivo de 13 níveis da 2b retido (grupo {CDS, WEkEO}) → valores
  idênticos esperados. Divergência é registrada, mas não invalida a comparação de níveis (feita no mesmo arquivo).

## 7. Recursos e limites

- Download: 10 pedidos WEkEO, 37 níveis, caixas seladas, 6 h. Estimativa: ≈ 0,3 GB/caso (13 níveis na 2b: 0,11 GB),
  ≈ 3 GB no total. **Teto: 10 GB de tráfego; 20 GB de disco novo.** Tempo: ≈ 24 min por pedido, ≈ 4–6 h em
  sequência; **teto de 12 h** de relógio para os downloads.
- Cálculo: 10 casos × 3 conjuntos + controles ≈ 35 execuções LEC em arquivo local (segundos a poucos minutos cada).
- Reaproveitado (só leitura): caixas seladas, arquivos de 13 níveis e CSVs de referência da 2b. Nada retido é
  alterado, movido ou apagado (cópia única, até o Marco B).
- Resultados fora do Git, em nova pasta `/p1-swell/danilocs/lectk_stage3_levels_<data>/`, com sha256.

## 8. Portões

1. Aprovação deste pré-registro e dos limites.
2. Inspeção somente leitura na Swell (disponibilidade, `develop` do clone, env `lorenz-4d`, caixas seladas e
   arquivos de 13 níveis retidos, disco livre, processos em curso) → relatório ao responsável.
3. Confirmação do uso da credencial WEkEO pelo responsável → downloads.
4. LEC, controles e tabelas; C1 deve passar antes da interpretação.
5. Registro no repositório (branch única, diff mostrado antes do commit, PR → `develop`).

## 9. Fora do escopo

Hipótese K; 1000 hPa abaixo do relevo (D-015); eleger fonte; mudar o padrão de níveis do LECTK (se sugerido,
vira recomendação para decisão separada); alterar código.

Aprovado pelo responsável em 27/11/2026, no pedido de execução da etapa 3; uso da credencial WEkEO autorizado no mesmo dia, depois da inspeção somente leitura da Swell.
