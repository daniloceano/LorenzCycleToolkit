# Frente 2b — fase 2: pré-registro da comparação de produtos (aprovado em 19/11/2026)

**Natureza:** registro documental. Torna permanente o pré-registro da fase 2 da Frente 2b, autorizada pelo
responsável em 19/11/2026 ([D-026](../../decisions.md), passo 4; D-027) e **aprovado por ele na mesma data**, com a
CDS API como referência de comparação. **Nada foi calculado antes da aprovação:** nenhum cálculo do LECTK, nenhuma
comparação e nenhum download. **Nada executável muda no repositório:** nenhum código, teste, script, dependência ou CI. `main`,
tags e releases continuam inalterados.

**Base:** `develop` em `068c3a990114cd66b44460f16a19334dc6f5f2d8` (PR #71);
[pré-registro v1](2026-11-09-front-2b-stage0-preregistration.md);
[resultados da fase 1](2026-11-13-front-2b-phase1-results.md);
[resultados da rodada complementar](2026-11-19-front-2b-round2-results.md).

## Resumo em linguagem simples

As campanhas da Frente 2b baixaram os mesmos 10 casos de ERA5 por várias rotas. A fase 2 pergunta: **os dados
que chegam por rotas diferentes são os mesmos números, e dão o mesmo ciclo de energia de Lorenz no LECTK?**
Para responder, rodamos o LECTK, sempre com a mesma versão e a mesma configuração, sobre o arquivo de cada rota já
guardado no servidor. Depois medimos a diferença de cada rota para uma rota escolhida como referência de
comparação. A referência é só uma convenção: o resultado não diz qual fonte é "a correta" ([D-015](../../decisions.md)).
Antes de olhar qualquer resultado, ficam fixados aqui os arquivos usados, as grandezas comparadas, as tolerâncias
e um controle.

## 1. Pergunta

Para cada caso e cada rota com download bem-sucedido, qual a diferença em relação à rota de referência:

- **(A)** nos campos de entrada T, z, u, v e ω, por nível;
- **(B)** nos termos do ciclo de energia calculados pelo LECTK.

Não se elege a melhor fonte nem a "mais correta" (D-015). Uma diferença diz que as rotas entregam números
diferentes, e não qual delas está certa.

## 2. Dados que entram

**O que é usado:** o arquivo de grade comum que o harness gravou para cada pedido bem-sucedido. Esse arquivo é um
recorte do que a rota entregou, sem interpolação: as mesmas datas, os 13 níveis, os mesmos pontos de 0,25° da caixa
selada e as variáveis em float32. A mesma função gravou os arquivos nas duas campanhas. O hash sha256 de cada
arquivo está no manifesto do pedido.

**Regra de escolha:** para cada rota e cada caso, usa-se o **primeiro download bem-sucedido**, isto é, o de menor
horário de início, considerando a fase 1 e a rodada complementar juntas. Os outros downloads bem-sucedidos só
entram no controle C2 (§6).

**Seleção resultante:** P1 = fase 1, R2 = rodada complementar e r = repetição. Entre parênteses, o número de
downloads bem-sucedidos disponíveis.

| Caso | CDS API (referência) | ARCO | WEkEO | DestinE | kerchunk | Get a Subset |
|---|---|---|---|---|---|---|
| 19790166 | R2r2 (2) | P1r1 (3) | P1r1 (3) | P1r2 (2) | R2r1 (3) | R2r1 (2) |
| 19790205 | R2r2 (2) | P1r1 (3) | P1r1 (3) | P1r1 (3) | R2r1 (3) | — |
| 19790535 | R2r2 (2) | P1r1 (3) | P1r1 (3) | P1r1 (3) | P1r2 (4) | P1r3 (1) |
| 19790585 | P1r1 (3) | P1r1 (3) | P1r1 (3) | P1r1 (3) | R2r1 (3) | — |
| 19790644 | R2r2 (2) | P1r1 (3) | P1r1 (3) | P1r2 (2) | R2r1 (3) | R2r1 (2) |
| 19791008 | P1r1 (3) | P1r1 (3) | P1r1 (3) | P1r2 (2) | R2r1 (3) | R2r1 (2) |
| 19791150 | R2r2 (2) | P1r1 (3) | P1r1 (3) | P1r1 (3) | P1r2 (4) | R2r1 (2) |
| 19800256 | R2r2 (2) | P1r1 (3) | P1r1 (3) | P1r1 (3) | R2r1 (3) | R2r1 (3) |
| 19800296 | P1r1 (3) | P1r1 (3) | P1r1 (3) | P1r1 (3) | P1r2 (4) | P1r2 (3) |
| 19800318 | P1r1 (3) | P1r1 (3) | P1r1 (3) | P1r1 (3) | R2r1 (3) | — |

Ao todo, são 10 arquivos de referência e 47 arquivos a comparar.

**Rotas incompletas:**

- **Get a Subset (7 de 10 casos):** compara-se só nos 7 casos com dado. Os casos 19790205, 19790585 e 19800318
  ficam marcados "sem dado (falha no download)": nada é estimado nem preenchido. Todo total da rota informa o
  número de casos (n = 7).
- **NCSS (0 de 16 na rodada 2, 0 de 12 na fase 1):** fica fora da fase 2, porque não há nenhum download completo. As
  partes incompletas que ficaram no disco **não** são usadas.

## 3. Referência de comparação (convenção, não "fonte correta")

**Referência: CDS API.** É só a rota contra a qual as diferenças são medidas. Os motivos da escolha são
práticos:

- é a rota mais utilizada para obter o ERA5;
- é a fonte oficial de distribuição do Copernicus (Climate Data Store).

Os 10 arquivos de referência vêm de dois momentos: 4 casos da fase 1 (repetição 1) e 6 da rodada complementar
(repetição 2, com outra chave), sempre com o mesmo pedido selado (§2).

A escolha não declara o CDS mais correto (D-015). Com outra referência, mudariam o sinal das diferenças e a escala
usada nas diferenças relativas, mas não a diferença absoluta entre duas rotas.

## 4. Configuração fixa do LECTK

- **Código:** export (`git archive`) de `develop` em `068c3a990114cd66b44460f16a19334dc6f5f2d8`, com o sha256 do
  tar registrado. O export roda a partir do seu próprio diretório. O registro de cada rodada anota de onde o módulo
  `src` foi carregado, para provar que veio do export e não do checkout do servidor. O checkout do servidor não é
  tocado.
- **Ambiente:** o ambiente `lorenz` do servidor (Python 3.12.12, NumPy 2.0.0, xarray 2024.2.0, MetPy 1.6.2), com a
  lista completa de pacotes (`pip freeze`) registrada.
- **Modo:** caixa fixa com termos residuais (`-f -r`), sem figuras e sem download.
- **Caixa:** a caixa selada de cada caso (`boxes.json`, sha256 `791a06e5…`), igual para todas as rotas. Ela é o
  retângulo que envolve a trajetória, com 5° de margem.
- **Níveis:** os 13 níveis comuns, de 1000 a 50 hPa. Nenhum arquivo selecionado tem valor não finito (verificação
  da campanha), então se esperam os 13 níveis em todas as rodadas. Se o LECTK excluir algum nível, isso é
  registrado e a comparação daquele caso usa a mesma regra para todas as rotas.
- **Tempo:** os horários do pedido comum, 00/06/12/18 UTC.
- **Lista de variáveis do LECTK:** uma só, escrita para os arquivos de grade comum (t, z, u, v, w; longitude,
  latitude, time, level), com o sha256 registrado.
- **Ressalva:** as bordas da caixa coincidem com as bordas dos dados. Por isso, as derivadas nas bordas usam
  diferenças de um lado só, igualmente em todas as rotas. Isso não afeta a comparação, mas os valores absolutos
  desses termos não servem para interpretação física.

**Ajustes permitidos depois da aprovação:** só no piloto (§8, E1), que roda **apenas a referência** e, por isso, não
mostra nenhuma diferença entre rotas. Pode-se ajustar apenas o formato da lista de variáveis ou do arquivo de
caixa, e cada ajuste fica registrado. Se o piloto pedir uma mudança de código no LECTK, a execução para e o
responsável é consultado.

## 5. Grandezas comparadas

### 5.1 Campos de entrada (A)

As medidas são feitas por caso, variável (T, z, u, v, ω) e nível (13), sobre todos os pontos e horários da caixa,
com Δ = rota − referência:

- diferença máxima, max |Δ|;
- diferença média, média de Δ (viés);
- média de |Δ|;
- fração de valores exatamente iguais;
- contagem de valores não finitos, que se espera ser 0.

### 5.2 Termos do ciclo (B)

As medidas são feitas por caso, a partir do CSV de resultados do LECTK, com as 26 colunas do modo `-f -r`:

| Grupo | Termos |
|---|---|
| Reservatórios | Az, Ae, Kz, Ke |
| Conversões | Cz, Ca, Ck, Ce |
| Fronteiras | BAz, BAe, BKz, BKe, BΦZ, BΦE |
| Geração (cálculo direto) | Gz, Ge |
| Tendências | ∂Az/∂t, ∂Ae/∂t, ∂Kz/∂t, ∂Ke/∂t |
| Resíduos | RGz, RKz, RGe, RKe |
| Diagnósticos auxiliares | C_overturning, M (só reportados, sem classe) |

Para cada caso c, termo X e rota s, com Δ(t) = X_s(t) − X_ref(t) em cada horário t:

- **escala:** S = máximo entre a média no tempo de |X_ref(t)| e um piso F. O piso é F = 10³ J m⁻² para os
  reservatórios e F = 0,01 W m⁻² para os demais termos. Ele evita dividir por quase zero;
- **erro máximo:** e_max = max_t |Δ(t)| / S;
- **erro médio:** e_med = |média_t Δ(t)| / S;
- **sinal:** se o sinal da média no tempo concorda entre a rota e a referência.

## 6. Classes e tolerâncias (fixadas antes de qualquer resultado)

**Campos de entrada (A)**, por caso, variável e nível:

| Classe | Critério |
|---|---|
| Idêntico | todos os valores iguais bit a bit (fração de iguais = 100 %) |
| Diferença pequena | não idêntico, max \|Δ\| ≤ τ e \|viés\| ≤ τ/10 |
| Diferença relevante | qualquer outro caso, ou algum valor não finito |

Tolerâncias τ:

| Variável | τ |
|---|---|
| T | 0,1 K |
| z | 10 m² s⁻² (≈ 1 m) |
| u, v | 0,1 m s⁻¹ |
| ω | 0,01 Pa s⁻¹ |

Esses valores ficam acima do passo de empacotamento de 16 bits do GRIB do ERA5, que é calculado por campo. Pela
amplitude típica de cada campo, esse passo é da ordem de milésimos de K em T e de décimos de m² s⁻² em z
**[ESTIMATIVA, não medida]**. As tolerâncias ficam também bem abaixo das diferenças que mudariam a meteorologia de
um caso. A classe de um caso
e rota é a pior entre as 65 combinações de variável e nível; a tabela completa também é reportada.

**Termos do ciclo (B)**, por caso, termo e rota:

| Classe | Critério |
|---|---|
| Idêntico | Δ(t) = 0 em todos os horários (valores iguais no CSV) |
| Diferença pequena | não idêntico, e_max ≤ 1 % e o sinal da média concorda |
| Diferença relevante | e_max > 1 %, ou sinal da média diferente com média da referência acima do piso F |

O limiar de 1 % fica bem abaixo do efeito de escolhas de método já medidas no projeto. Por exemplo, reduzir de 17
para 13 níveis mudou Ca em cerca de 5 % na investigação LZV. Os valores de e_max e e_med são sempre reportados,
não só a classe.

**Agregação:** para cada rota, a contagem de pares caso × termo em cada classe, por grupo de termos e no total, sempre
com o n. As rotas não são ordenadas entre si (D-015).

**As tolerâncias acima ficaram fixadas na aprovação de 19/11/2026.** Não se mudam depois de ver resultados: um resultado que surpreenda é
reportado como está, com explicação.

## 7. Controles

- **C1, determinismo do LECTK (barreira):**
  - **O que se faz:** cada um dos 10 arquivos de referência é processado duas vezes, em diretórios separados;
  - **o que se espera:** classe "idêntico" em todos os termos e horários, nos 10 casos (diferença zero);
  - **se falhar:** para antes de qualquer comparação entre rotas, porque as diferenças seriam ruído do próprio
    cálculo.
- **C2, repetibilidade da entrega (informativo):**
  - **O que se faz:** para cada rota e caso com 2 ou mais downloads bem-sucedidos, compara-se o sha256 dos arquivos
    de grade comum entre as repetições, recalculado no disco e conferido com o manifesto;
  - **se os hashes diferirem:** aplicam-se as medidas da §5.1 entre as repetições, sem rodar o LECTK;
  - **papel:** mostra se a mesma rota devolve os mesmos números em dias diferentes. Não é barreira.
- **C3, autoteste do comparador (barreira):**
  - **referência contra ela mesma:** tudo deve dar zero e "idêntico";
  - **perturbação sintética:** T + 0,05 K num único ponto deve dar "diferença pequena" com max |Δ| = 0,05 K;
    T + 0,2 K deve dar "diferença relevante";
  - **termos do ciclo:** um CSV com um termo multiplicado por 1,005 deve dar "pequena", e por 1,02, "relevante".

## 8. Sequência de execução (só depois da aprovação)

1. **E0 — selagem:** gravar o sha256 deste pré-registro, da tabela de seleção (caminhos e sha256 dos 57 arquivos,
   conferidos no disco), do comparador, da lista de variáveis e do tar do export. A selagem vem antes de qualquer
   cálculo.
2. **E1 — piloto, só com a referência:** rodar o LECTK sobre o arquivo do CDS do caso 19800296, o menor. O piloto confere a
   configuração (13 níveis usados, caixa, colunas) e mede tempo e memória. Se o tempo projetado passar de 48 h, a
   execução para e o responsável é consultado.
3. **E2 — controles:** C3, depois as 10 rodadas de referência e C1.
4. **E3 — rodadas das rotas:** as 47 rodadas dos arquivos selecionados.
5. **E4 — comparação e relatório:** comparações A e B, C2 e relatório. O relatório é apresentado ao responsável,
   e nada é homologado sem a aprovação dele.

**Regras de parada e de falha:**

- **Falha de controle:** C1 ou C3 falhou → para.
- **Erro do LECTK num arquivo:** a rodada é registrada como "falha do LECTK", com o log. Não se muda o código nem
  a configuração e não há nova tentativa. A única exceção é uma interrupção externa, como queda do servidor, que
  permite uma nova tentativa.
- **Arquivo alterado:** se o sha256 de um arquivo selecionado não conferir com o selado, o arquivo fica de fora e
  isso é registrado. Não se troca o arquivo por outro download.
- **Disco ou tempo esgotado:** disco livre abaixo de 1 TB ou duração total acima de 72 h → para, com relatório
  parcial.
- **Sem rede:** a fase 2 não faz nenhum download nem acesso a fontes remotas.

## 9. Recursos no servidor

- **Diretório:** um diretório novo, `<diretório da fase 2 no servidor>`. Os dois diretórios das campanhas são
  **só lidos**. Os arquivos de entrada entram por links simbólicos, com o sha256 conferido antes e depois de cada
  rodada. Nada existente é alterado, movido ou apagado.
- **Disco:** até 5 GB de saída (CSVs, logs e tabelas; sem figuras).
- **Processos:** no máximo 4 rodadas do LECTK ao mesmo tempo, com `nice -n 10` e no máximo 2 threads numéricas cada
  uma, numa sessão tmux própria.
- **Agendamento:** nenhuma tarefa agendada (crontab) e nenhuma retomada automática. Se algo for interrompido, a
  retomada é feita à mão, a partir do estado registrado.
- **Memória:** o maior arquivo tem 0,71 GB (caso 19790205). Esperam-se menos de 20 GB por rodada, num servidor com
  cerca de 1 TiB.
- **Tempo:** estimado em menos de 1 dia para as 67 rodadas (10 de referência, 10 do C1 e 47 das rotas). O piloto
  confirma a estimativa, e o teto é de 72 h.

## 10. Proveniência e resultados (fora do Git)

Ficam no diretório da fase 2:

- os manifestos por rodada: comando exato, horários UTC, sha256 da entrada e de cada saída, origem do módulo
  `src`, versão dos pacotes e resultado;
- as tabelas de comparação (A, B e C2) em CSV;
- os logs;
- o relatório;
- um índice final de sha256 que sela tudo.

No Git entra só o registro documental dos resultados, depois da aprovação do relatório. Antes de qualquer commit, o
diff exato é mostrado ao responsável.

## 11. Fora da fase 2

- **Eleger fonte:** não se elege uma fonte nem se recomenda uma rota (D-015). Isso fica para a etapa 4 (inclusão das
  fontes), com os resultados das fases 1 e 2.
- **Níveis verticais:** ficam para a etapa 3 (D-016).
- **Retenção:** a retenção dos dados das campanhas vem só depois da fase 2 (D-026); nada é apagado.
- **Defeitos:** a correção dos defeitos da §12 não é feita na fase 2.

## 12. Pendências registradas para a etapa 4 (inclusão das fontes)

Por decisão do responsável de 19/11/2026, os três defeitos achados na rodada complementar ficam registrados como
pendências da etapa 4, **sem abrir frente agora**:

1. **Limite por arquivo do Get a Subset (cliente LECTK, `src/utils/gdex_sources.py`).** O limite por arquivo,
   calculado pelo próprio cliente, recusou os casos 19790205 e 19790585 nas três repetições.
2. **Proteção do CDS (`src/utils/cds_guard.py`).** Uma única falha de download, como o download de 0 B de 14/11,
   trava a proteção e, com ela, a rota inteira, até uma liberação manual.
3. **Harness das campanhas (fora do Git).** O harness purgava no NCAR só os pedidos de casos bem-sucedidos, e os
   pedidos de casos com falha ocupavam vagas. **Corrigir só se houver nova rodada.**

Os outros achados da §6 dos [resultados da rodada complementar](2026-11-19-front-2b-round2-results.md) continuam
como estão.

## Fora deste registro

Aprovação, merge, CI e releitura deste documento não podem constar dele.
