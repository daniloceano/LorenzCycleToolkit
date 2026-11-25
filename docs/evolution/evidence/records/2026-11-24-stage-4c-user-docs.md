# Etapa 4c — reescrita da documentação de usuário das rotas ERA5 (24/11/2026)

**Natureza:** registro documental. Só documentação de usuário (`docs/source`), registros e `CHANGELOG.md`. Nenhum
código, teste, dependência ou CI muda. Nenhum uso da Swell e nenhum download. `main`, tags e releases inalterados.

**Base:** `develop` em `dbd7fbc` (PR #85); `main` em `48fe8e2`. Etapa 3 "Não aberta" (conferido). Branch única da
etapa: `docs/stage-4c-user-docs` (D-042).

## Resumo em linguagem simples

A página `remote_data_acquisition.rst` foi reescrita para as 7 rotas ERA5. Para cada uma, o usuário encontra como
ativá-la, o que precisa (conta, token, pacote), o tempo e o volume que observamos, como os dados diferem dos do CDS,
as limitações e os melhores usos (D-012). Pela D-015, a página mostra as diferenças e não elege uma fonte correta. As
páginas ligadas foram atualizadas sem duplicar conteúdo, e as notas provisórias da etapa 4 foram removidas.

## 1. Decisão D-042 (responsável, 24/11/2026)

Com uma única frente ativa, a etapa usa uma branch só: cada parte é revisada pelo responsável e commitada na branch da
etapa; há um PR para `develop` no fim, com as notas de encerramento dentro dele. Com frentes em paralelo, volta o PR
por parte. Não é preciso esperar o CI de uma parte para começar a seguinte.

## 2. Tabela de fatos por rota (fonte de cada número da documentação)

"2b" = benchmark de 11/2026: 10 casos de 1979–1980, caixa = track + 5°, 13 níveis (1000–50 hPa), 6 h, 3 repetições;
mediana por caso ([pré-registro §2.3](2026-11-09-front-2b-stage0-preregistration.md)). "Teste da etapa 4" = caso
19790205 pelo comando principal, caixa de 15° (N−24 W−82 S−80 E76, 225 × 633), 34 horários de 6 h, níveis do produto,
uma execução.

| Rota | Opção / extra | Credencial | Níveis / passo | 2b: tempo e volume por caso | Teste da etapa 4 | Diferença em relação ao CDS | Limites | Estado |
|---|---|---|---|---|---|---|---|---|
| CDS API | `--cdsapi` / — | `~/.cdsapirc` (modo protegido: `CDSAPI_URL`+`CDSAPI_KEY` ou `CDSAPI_RC`) | 37 / 1 h | rep. 2: 2,7 h (6,5 min–4,4 h), 0,09 GB; rep. 3: 30 s (15–73 s), cache provável [NÃO CONFIRMADO] ([round2 §2–3](2026-11-19-front-2b-round2-results.md)) | não feito | referência por convenção (D-015) | fila; limites de custo da conta | oficial (D-029) |
| WEkEO | `--wekeo` / `[wekeo]` | `~/.hdarc` ou `HDA_USER`/`HDA_PASSWORD` | 37 / 1 h | 1,2 min (0,4–7,3), 0,11 GB, 30/30, cache possível (round2 §2) | fila 1 017 s, transferência 273 s, montagem 135 s, ≈ 24 min; 1,60 GiB; limite superior 3,54 GiB ([4b·1 §5.2](2026-11-23-stage-4b1-wekeo.md)) | 2b: idêntica (240/240, [fase 2 §4](2026-11-19-front-2b-phase2-results.md)); 4b·1: t ≤ 1,2·10⁻⁴ K, z ≤ 0,0029, ciclo ≤ 0,007 % (4b·1 §5.3–5.4) | 500 itens/h (doc. HDA); sem cancelamento; fila com prazo de 1 h (4b·1 §1.2, §2) | integrada (D-036), 1 caso |
| kerchunk | `--ncar` (padrão) / — | nenhuma; referências 14–28 MB | 37 / 1 h | 16 min (6,9–33), 7,5 GB, 30/30 (round2 §2) | aquisição 1 856 s, transferência 1 656 s; 9,85 GiB em 170 pedidos; dry-run 9,864 GiB ([4a §7.3](2026-11-19-stage-4a-routes-and-fixes.md)) | máx. 0,043 %; T 2,4·10⁻⁴ K; z 0,0625 (fase 2 §3–4) | globo inteiro por campo e hora; teto 20 GiB (D-035); conexão intermitente (4a §6) | oficial, padrão (D-029, D-030) |
| Get a Subset | `--ncar --ncar-method subset` / — | `GDEX_API_TOKEN` ou `~/.gdex_api_token` | 37 / todas as horas | 6,6 min (2,6–12,6), 1,1 GB, 13/20 rodados, sem a espera de vaga (20 esperas, 46 h; maiores 6–8 h) (round2 §2–3) | não feito | igual ao kerchunk (fase 2) | 8 abertos (LECTK para com 7 listados); dedup.; purga 7 dias | oficial; correção do limite por arquivo confirmada só num arquivo (4a §4–5) |
| NCSS | `--ncar --ncar-method ncss` / — | nenhuma | 37 / 1 h | 0/16 (round2 §2); 0/12 na fase 1 ([fase 1 §2](2026-11-13-front-2b-phase1-results.md)) | não feito | — | timeouts, 504 | experimental (D-029) |
| ARCO | `--arco` / `[arco]` | chave CDS (`CDSAPI_KEY`, `CDSAPI_RC`, `~/.cdsapirc`) | 13 / 6 h | 6,8 min (2,6–20), 2,5 GB, 30/30 (round2 §2) | 826 s total (aquisição 793,1 s); 4,60 GiB; 4 431 leituras; limite superior 12,16 GiB; 1ª execução 50 min ([4b·2 §5.1, §5.4](2026-11-24-stage-4b2-arco.md)) | máx. 0,043 % (fase 2); 0,007 % (4b·2 §5.3) | guia sem limites declarados; formatos podem mudar (4b·2 §1.2) | integrada (D-039), 1 caso |
| DestinE | `--destine` / `[destine]` | `~/.netrc` ou `DESTINE_EDH_TOKEN` | 19 (17 efetivos) / 1 h; store v2 congelado desde 06/07/2026 | 7,8 min (2,6–19), 3,1 GB, 27/30 (round2 §2) | 1 066 s total (aquisição 1 010,4 s); 1 528 requisições; 8,79 GiB; limite superior 45,23 GiB ([4b·3 §5.1](2026-11-24-stage-4b3-destine.md)) | 18/240 termos entre 1,04 % e 2,67 %; T ≤ 0,125 K; z ≤ 64 (fase 2 §3–4); 4b·3: só BΦE 1,21 % (§5.3) | 500 000 requisições/mês; `--max-requests` 10 000 (D-041) | integrada (D-041), 1 caso |

Outros números usados: efeito de 13 e 17 níveis (Ca −5,1 %, Ce −1,8/−1,9 %, Ck −0,2 %; figuras) da
[LZV §3.1](2026-11-02-lzv-era5-zarr-vertical-levels.md); Ck +12,5 % com os níveis da DestinE até 10 hPa contra 50 hPa,
[PRELIMINAR] (4b·3 §5.3); excesso do ARCO de 7,5× (LZV §3.2), 10,1× (fase 1 §3) e 3,9× (4b·2 §5.1); comparação NCAR ×
CDS com 37 níveis e tempos da D-Bench (texto anterior da página, mantido como histórico, com legendas datadas); datas de
cobertura: CDS conferido em 02/11, GDEX em 26/10 (texto anterior), ARCO 18/11 18 UTC lido em 24/11 (4b·2 §1.4).

## 3. O que estava desatualizado (levantamento da Parte 1)

- `remote_data_acquisition.rst`: aviso "under revision" e três notas provisórias; NCSS apresentado como o método do
  `--ncar`; kerchunk e Get a Subset "fora do comando principal" e "testados só offline"; módulo que "não baixa as
  referências"; ARCO "preparação externa" (com exemplo que lia `CDS_API_KEY`, nome que o código não usa); DestinE
  "não testada", com o formato de bloco do v3; WEkEO ausente; desempenho, falhas e precisão só com a D-Bench.
- `usage.rst`: download e lista de flags só com `--cdsapi` e `--ncar`, "experimental".
- `troubleshooting.rst`: nada sobre WEkEO, ARCO, DestinE, recusas por volume ou requisições, pacote ausente.
- `installation.rst`: sem os extras. `results.rst`: `provider` só `local`, `cds`, `ncar`.
- `examples_and_tutorials.rst`: "experimental `--ncar` route". `ncar_gdex_example.rst`: referências duplicadas.
- Teto de 8 GiB: nenhuma citação restante. Cópia do namelist: já corrigida em `configuration.rst` (D-038).

## 4. Mudanças por página

| Página | Mudança |
|---|---|
| `remote_data_acquisition.rst` | Reescrita: 1. Escolha rápida (duas tabelas: acesso e o que observamos; datas e passos; escolha por finalidade); 2. Uma seção por rota, com uma seção comum a todas (contrato do arquivo, pedido, níveis efetivos, tetos, reuso, proveniência, credenciais) e uma seção comum às três NCAR; 3. Diferenças nos dados (três grupos da fase 2; WEkEO = CDS só com pedido igual; aviso da DestinE em destaque; testes da etapa 4; níveis verticais da LZV; comparação NCAR × CDS de 37 níveis como histórico); 4. Ressalvas das medidas (medidas poucas, épocas, cache, tamanho do pedido; tabela da 2b; tabela dos testes da etapa 4; D-Bench com as figuras stageB como histórico); 5. Limitações conhecidas. Âncoras antigas mantidas; novas: `rda-quick-choice`, `rda-wekeo`, `rda-ncar`, `rda-destine-warning`, `rda-ncar-cds-37`, `rda-dbench-history`, `rda-limitations` |
| `usage.rst` | Exemplos das 5 opções de download; reuso do arquivo; flags `--wekeo`, `--ncar-method`, `--kerchunk-refs`, `--arco`, `--destine`, `--max-gib`, `--max-requests`; `--time-resolution` e `--delete-downloaded-input` para todas as rotas |
| `troubleshooting.rst` | Novas entradas: recusa por volume ou requisições; pacote opcional ausente; credencial ausente ou recusada (WEkEO, ARCO, DestinE, Get a Subset); DestinE HTTP 429 e store congelado; ARCO e `--time-resolution`; pedido WEkEO de estado desconhecido. Entradas existentes generalizadas para todas as rotas; FAQ de diferenças entre rotas atualizada |
| `installation.rst` | Seção dos extras `[wekeo]`, `[arco]`, `[destine]` |
| `configuration.rst` | Namelist com `--ncar` |
| `results.rst` | `provider` e fase `acquisition` para todas as rotas |
| `examples_and_tutorials.rst` | Referência às outras rotas |
| `ncar_gdex_example.rst`, `cdsapi_example.rst`, `era5_products_vertical_levels.rst`, `remote_data_acquisition_case_details.rst` | Páginas "movidas": referências cruzadas para as novas seções |
| `index.rst` | Sem mudança (a estrutura do índice continua válida) |

Ajustes à estrutura aprovada, com motivo: seção comum a todas as rotas (evita repetir 7 vezes o mesmo texto); seção
comum às três NCAR (idem); níveis verticais e comparação de 37 níveis dentro de "Diferenças nos dados" (conteúdo que
existia e se perderia); D-Bench e figuras stageB dentro de "Ressalvas das medidas", como histórico; a tabela da escolha
rápida virou duas tabelas (legibilidade). Rótulo de estado de WEkEO, ARCO e DestinE: "integrated; tested on one case".

Ajuste pedido pelo responsável após a prévia (25/11/2026): no início da escolha rápida, um bloco visual "Pros and cons
at a glance", com um cartão por rota (prós ✔ e contras ✖, estado e link para os detalhes). Cada item resume um fato já
presente na página e na tabela de fatos acima; nenhum número novo.

## 5. Build e links

- **Build (Sphinx 7.4.4, sphinx-rtd-theme 2.0.0, ambiente isolado fora do repositório, `-n`):** `develop` (antes): 0
  avisos; branch (depois): 0 avisos. Nenhum aviso novo e nenhum preexistente.
- **Linkcheck:** antes e depois, um link quebrado, preexistente: `../LICENSE` (link relativo de `license.rst`, fora da
  árvore do Sphinx). Redirecionamentos (funcionam): `doi.org/10.5065/BH6N-5N20`, a página do ARCO no Confluence, a página
  do dataset no CDS, um anexo do GitHub e, novos, `www.wekeo.eu/register` (→ `data.wekeo.copernicus.eu/register`) e a
  página de chaves da Earth Data Hub (→ login da DestinE). Nenhum endereço foi alterado por adivinhação.

## 6. Achados para decisão (não corrigidos: código)

- O `--help` de `lorenzcycletoolkit.py` ainda descreve `--cdsapi` como "(experimental)", o que contradiz a D-029.
- O `--help` de `--time-resolution` diz "for CDS API data download", mas a opção vale para todas as rotas.

## 7. Limitações

- Os números da documentação são os dos registros; nenhuma medida nova foi feita.
- A data mais recente servida pela WEkEO e a última data do store v2 da DestinE não foram conferidas; a página o diz.
- A revisão visual do responsável é sobre a prévia local; a publicação no Read the Docs não faz parte desta etapa.

## 8. Nota de fechamento (25/11/2026): aprovação, PR, CI, merge e releitura

- **Aprovação:** o responsável revisou a prévia renderizada e o diff, pediu um ajuste (cartões de prós e contras na
  escolha rápida, §4) e autorizou o commit em 25/11/2026.
- **Commit e PR:** commit `b37420d` na branch da etapa (`docs/stage-4c-user-docs`); PR único para a `develop`:
  [#86](https://github.com/daniloceano/LorenzCycleToolkit/pull/86).
- **CI:** verde: os dois `build` do workflow "Python application" e `ci/circleci: build_test`.
- **Merge:** feito pelo responsável em 25/11/2026, merge `30e6c2b` (`30e6c2bcfc7e136d2ed277d96672d3aac4844145`).
- **Releitura no SHA do merge:**
  - a árvore do merge é idêntica à do commit `b37420d` (diff vazio);
  - o build do Sphinx (`-n`, o mesmo ambiente isolado) a partir de um `git archive` de `30e6c2b` terminou com 0 avisos;
  - a página renderizada tem os 7 cartões de prós e contras;
  - em `docs/source` não há caminho da Swell nem credencial, só modelos como `<your-api-key>`.
- **Pendências da etapa 4 (sem mudança):** decisão sobre ler os stores Zarr v3 da DestinE, retenção dos dados na Swell
  (D-028) e resumo de transferência. Os dois achados do `--help` (§6) continuam abertos para decisão.
