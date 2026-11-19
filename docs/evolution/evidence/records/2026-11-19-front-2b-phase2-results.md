# Frente 2b — fase 2: resultados da comparação de produtos (19/11/2026)

**Natureza:** registro documental. Torna permanentes os resultados da fase 2 da Frente 2b, executada em 19/11/2026
conforme o [pré-registro](2026-11-19-front-2b-phase2-preregistration.md) (D-027). O relatório foi **aprovado pelo
responsável em 19/11/2026**. **Nada executável muda no repositório:** nenhum código, teste, script, dependência ou
CI. O harness, os manifestos, as tabelas e o relatório original ficam fora do Git, no servidor de processamento.
`main`, tags e releases continuam inalterados.

**Base:** `develop` em `fd30f83` (PR #72).

**Estado:** fase 2 CONCLUÍDA. A retenção dos dados das campanhas (D-026, passo 5) não foi aberta.

## Resumo em linguagem simples

**Pergunta:** os mesmos 10 casos de ERA5, baixados por 6 rotas, chegam com os mesmos números e dão o mesmo ciclo de
energia de Lorenz no LECTK?

**Método:** o LECTK rodou com versão e configuração fixas sobre o arquivo já guardado de cada rota, sem nenhum
download novo. Cada rota foi comparada com a CDS API, que é a referência por convenção ([D-015](../../decisions.md)).

**Resultado:** apareceram três grupos de produto.

- **CDS e WEkEO:** são idênticos.
- **ARCO, kerchunk e Get a Subset:** entregam o mesmo arquivo entre si e diferem do CDS só por arredondamento, com
  diferenças desprezíveis no ciclo.
- **DestinE:** os valores vêm arredondados em degraus. Isso passa da tolerância nos campos de entrada e muda entre
  1 % e 2,7 % alguns termos de fronteira, de geração e de resíduo. Os reservatórios e as conversões não mudam.

## 1. Execução

- **Rodadas:** 68, todas bem-sucedidas na primeira tentativa: 1 piloto, 10 de referência, 10 do controle C1 e 47 das
  rotas.
- **Código:** export de `develop` em `068c3a9`, rodado em caixa fixa com resíduos (`-f -r`), nas caixas seladas.
  Todas as rodadas usaram os 13 níveis, de 1000 a 50 hPa.
- **Seleção:** 57 arquivos, que são o primeiro download bem-sucedido de cada rota e caso. O Get a Subset entrou com 7
  casos e o NCSS ficou de fora.
- **Recursos:** 506 s de processamento no total e pico de 6 GB de memória. Não houve rede nem tarefa agendada.
- **Selos (fora do Git):**
  - antes de qualquer rodada: `9c510d78…`;
  - desvio 01: em `state/deviation01*`;
  - selo final: `8a2fee50…` (2.277 arquivos);
  - relatório: `a9651182…`.
- **Desvio 01 (autorizado pelo responsável em 19/11):**
  - **o que houve:** o autoteste do comparador (C3) falhou em 2 de 6 verificações, porque o leitor padrão de CSV do
    pandas não é exato no último bit;
  - **correção:** uma linha no comparador, para ler com `float_precision="round_trip"`;
  - **o que não mudou:** nenhum critério. Nenhum resultado de rota existia antes do desvio.
- **Interpretações registradas antes das rodadas:**
  - o teste de sinal só se aplica quando a média da referência passa do piso F;
  - o C3 usa o termo Az.

## 2. Controles

| Controle | Resultado |
|---|---|
| C3, autoteste do comparador (barreira) | PASS 6/6, depois do desvio 01 |
| C1, determinismo do LECTK (barreira) | PASS: 260/260 pares caso × coluna idênticos |
| C2, repetibilidade da entrega (informativo) | 0 diferenças: os 158 arquivos de downloads repetidos (56 pares rota × caso) são idênticos bit a bit ao primeiro da mesma rota |

## 3. Campos de entrada

| Rota | Casos | Classe por caso | Maior \|Δ\| em T | Maior \|Δ\| em z | u, v, ω |
|---|---|---|---|---|---|
| WEkEO | 10 | 10 idênticos | 0 | 0 | idênticos |
| ARCO | 10 | 10 pequena | 0,00024 K | 0,0625 m² s⁻² | idênticos |
| NCAR kerchunk | 10 | 10 pequena | 0,00024 K | 0,0625 m² s⁻² | idênticos |
| NCAR Get a Subset | 7 | 7 pequena | 0,00024 K | 0,0625 m² s⁻² | idênticos |
| DestinE | 10 | **10 relevante** | 0,125 K | 64 m² s⁻² | pequena |

**DestinE:**

- **T:** a diferença máxima é de 0,0625 K de 50 a 300 hPa e de 0,125 K de 400 a 1000 hPa.
- **z:** a diferença máxima é de 64 m² s⁻² de 50 a 150 hPa e cai até 2 m² s⁻² em 1000 hPa. Passa da tolerância de 50
  a 600 hPa.
- **viés:** é sempre pequeno (|viés| ≤ 0,0006 K e ≤ 0,49 m² s⁻²).

**Análise posterior** (feita depois de ver os resultados, não pré-registrada): pelos hashes, os arquivos de grade
comum formam três grupos idênticos em todos os casos: {CDS, WEkEO}, {ARCO, kerchunk, Get a Subset} e {DestinE}.

## 4. Termos do ciclo

Foram 24 termos com classe por caso; e_max é o erro máximo relativo à escala do termo na referência.

| Rota | Pares | Idêntico | Pequena | Relevante | Maior e_max | Sinal trocado |
|---|---|---|---|---|---|---|
| WEkEO | 240 | 240 | 0 | 0 | 0 | 0 |
| ARCO | 240 | 95 | 145 | 0 | 0,043 % | 0 |
| NCAR kerchunk | 240 | 95 | 145 | 0 | 0,043 % | 0 |
| NCAR Get a Subset | 168 | 57 | 111 | 0 | 0,043 % | 0 |
| DestinE | 240 | 0 | 222 | **18** | 2,67 % | 0 |

**Onde estão os 18 relevantes da DestinE:**

| Grupo de termos | Relevantes |
|---|---|
| fronteiras | 11 (BΦE ×6, BΦZ ×4, BAe ×1) |
| geração | 4 (Gz ×2, Ge ×2) |
| resíduos | 3 (RGz ×3) |
| reservatórios, conversões e tendências | nenhum (todos com diferença pequena) |

O e_max desses 18 pares fica entre 1,04 % e 2,67 %, e a diferença média no caso (e_med) fica em no máximo 0,37 %.

**[INTERPRETAÇÃO]** Os degraus da DestinE são compatíveis com valores guardados com precisão reduzida, o mecanismo não
foi verificado. Os termos afetados dependem de z nas paredes da caixa ou fecham o balanço por diferença.

## 5. Ressalvas

1. As bordas da caixa coincidem com as bordas dos dados, então as derivadas nas paredes são de um lado só, igualmente
   em todas as rotas. Isso pode amplificar o efeito nos termos de fronteira.
2. "Relevante" quer dizer acima do limiar pré-registrado (τ nas entradas, 1 % no ciclo), e não que uma conclusão
   meteorológica mude.
3. A amostra tem 10 casos do Atlântico Sul, de 1979–1980, com intervalo de 6 h e 13 níveis.
4. A referência é só uma convenção. Nenhum dos três grupos é declarado "o correto" (D-015).

## Fora deste registro

Fica de fora a retenção dos dados das campanhas (D-026, passo 5), e nada foi apagado. Também não constam merge, CI e
a releitura deste documento.
