# Investigação complementar pré-F — ERA5/Zarr, níveis verticais e correção do alinhamento temporal (02/11/2026)

**Estado: CONCLUÍDA TECNICAMENTE — AGUARDANDO INSPEÇÃO PRÉ-COMMIT E INTEGRAÇÃO.** Investigação **complementar**, fora da ordem principal: **não abre a Frente F**, não reabre D, D-Bench, E, G ou K e não altera `main`, tags, versões ou publicações. Este registro só deve ser lido como definitivo depois de integrado e relido em `develop`.

## 1. Natureza, escopo e perguntas

Investigação autônoma executada no servidor de processamento do grupo em 02/11/2026 (≈ 1 h 52 min de campanha), com quatro perguntas:

1. Os produtos ERA5 em Zarr (CDS ARCO, 13 níveis, 6 h; DestinE, 19 níveis, 1 h) trazem vantagem operacional para estudos de caso com o LECTK?
2. Quanto a redução do número e da distribuição dos níveis verticais altera Ca, Ce e Ck?
3. Quando cada configuração é cientificamente adequada?
4. Convém implementar uma rota Zarr ou priorizar outros métodos de aquisição?

Na mesma data, com autorização posterior do responsável, foi implementada a correção do defeito de alinhamento temporal encontrado pela investigação (§7). A página de usuário [ERA5 Products and Vertical Levels](../../../source/era5_products_vertical_levels.rst) resume os resultados para quem usa o toolkit.

## 2. Caso, protocolo e proveniência

- **Caso:** trajetória real `19791150` da D-Bench (29/12/1979–04/01/1980), modo `-t -r`, caixas de 15° × 15°, **27 horários de 6 h**. Referência: o arquivo CDS de 37 níveis e 3 h da D-Bench (sha256 `ad9de4037fcc5e807bdddad823e8a01a5a96c77dcb0ef0179cbc1cb5496b1a02`), reutilizado sem nova aquisição. Escolhido por ser o único caso completo com ≥ 4 dias cujo domínio não cruza Greenwich; a divergência NCAR × CDS desse caso **não** foi critério.
- **Protocolo:** pré-registro selado antes de medir (`PRE_REGISTRO_LZV_v1.md`, sha256 `3a47e732f18ac437b34e1c2c5d98ed576745bb7f4a5a8454210dc00935c66844`, 13:18:27 UTC) e **adendo v2** selado depois do desvio do controle positivo (`ADENDO_LZV_v2.md`, sha256 `1a08b51afbecb5183860f6c787a6bcc1640a793e2462a6f44e9031a292c7267a`); ver §4.
- **Código:** export imutável (`git archive`) de `develop` `835478e` (árvore `b959f5801b61f2f51f24a99f53a33becd53c2794`). Entre `835478e` e `b4a5f77` só mudaram documentos de `docs/evolution/`; `src/`, CLI e `inputs/` são idênticos.
- **Ambiente:** Python 3.12.12, NumPy 2.0.0, xarray 2024.2.0, pandas 2.2.2, MetPy 1.6.2, netCDF4 1.7.1.post1; leitura Zarr em ambiente isolado (zarr 2.18.4, numcodecs 0.15.1, fsspec 2024.6.1).
- **Orçamento:** ≈ 8,6 GiB de tráfego novo (teto 30 GiB), uma aquisição por vez, nenhum serviço pago, nenhuma credencial alternativa (as chaves exclusivas da D-Bench não foram usadas).

## 3. Resultados

### 3.1 Sensibilidade vertical controlada

Cinco conjuntos de níveis recortados do **mesmo** arquivo CDS, com valores idênticos bit a bit nos níveis mantidos: 32 (1000–10 hPa, referência), 17 (grade DestinE), 29 (1000–50 hPa), 14 (DestinE ≥ 50 hPa) e 13 (grade ARCO). Os 32 níveis efetivos da referência resultam do corte do LECTK abaixo de 10 hPa (1, 2, 3, 5 e 7 hPa descartados dos 37).

| Comparação | Natureza | Ca | Ce | Ck |
|---|---|---|---|---|
| 32 → 17 | discretização, 1000–10 hPa | média −5,1 % (7,692 → 7,293 W m⁻²); até 19,5 % num horário; r = 0,995 | −1,9 %; até 5,7 % | −0,2 %; até 12,9 % |
| 29 → 13 | discretização, camada comum 1000–50 hPa | −5,1 %; até 19,2 % | −1,8 %; até 5,9 % | −0,2 %; até 12,9 % |
| 32 → 29 | extensão (topo 10 → 50 hPa) | < 0,05 % | < 0,05 % | < 0,05 % |

- Nenhuma inversão de sinal nem deslocamento dos extremos; todas as comparações ficaram na classe pré-registrada **b** (acima da tolerância numérica estrita, sem marcador interpretativo). Uma reimplementação independente das métricas concordou até 5 × 10⁻¹².
- **Origem:** −0,396 dos −0,399 W m⁻² da queda média de Ca vêm da camada 1000–925 hPa, onde o espaçamento passa de 25 para 75 hPa e a estabilidade estática no nível inferior (diferença unilateral) muda. Em Ck, compensações de até ±0,11 W m⁻² entre camadas. Parte da caixa pode ter 1000 hPa abaixo da superfície (valores extrapolados): ressalva física não resolvida.
- **Hipóteses:** H1 parcialmente confirmada (Ca mais sensível; refutado que Ck seja mais sensível que Ce na média); H2 e H3 confirmadas; H4 refutada quanto aos bits e confirmada quanto à energética; H5 confirmada; H6 confirmada no sentido (7,5×, não 14×).

### 3.2 Produtos reais e precisão

- **CDS ARCO Zarr (store time-chunked):** funcionou com a chave CDS existente usada como token *bearer*; não se afirma acesso sem conta. u, v, w idênticos bit a bit ao NetCDF CDS; T com deslocamentos uniformes de +1 a +8 ULP float32 (até 2,44 × 10⁻⁴ K) em 92 de 351 fatias, só em 700–1000 hPa; z com diferenças pequenas. Energética × 13 níveis do CDS: 26/26 colunas dentro da tolerância estrita.
- **NCAR NCSS (37 níveis):** **idêntico bit a bit ao ARCO** nos cinco campos dos 13 níveis comuns. Com 32 níveis, os mesmos deslocamentos de T levaram Az, Ae, BAz, Gz, tendências e resíduos além da tolerância estrita (padrão da D-Bench). O NetCDF entregue pelo CDS é o divergente; **não foi demonstrado qual representação reproduz o GRIB original** (exigiria pedido GRIB, não feito).
- **DestinE:** BLOQUEADO (HTTP 401 inclusive em metadados; sem credencial). Estimativa não medida: ≈ 1425–2280 chunks (≈ 29–47 GB não comprimidos) para < 1 GB de valores neste caso.

### 3.3 Desempenho operacional (uma medição por configuração)

| Rota | Níveis | Tempo até o resultado | Recebido |
|---|---:|---:|---:|
| ARCO, `xarray` + fsspec, 4 fluxos | 13 | 6,6 min | 2,71 GB |
| ARCO, 1 fluxo, conexão persistente (repetição: 21 min) | 13 | 23,1 min | 2,71 GB |
| NCAR NCSS (`--ncar`), 1 fluxo | 32 | 19,7 min | 1,07 GB |
| NetCDF local já em disco | 32 / 13 | 11 s / 10 s | — |
| CDS, histórico D-Bench (30/10, 3 h; não controlado) | 32 | 65,6 min | 0,92 GB |

**Limites:** as linhas diferem em produto, cliente e concorrência; não constituem comparação controlada de formatos. O ARCO só foi mais rápido que a NCAR com paralelismo, transferiu 2,5× os bytes (faixas de latitude inteiras: 7,5× o volume útil) e sua leitura local foi ≈ 200× mais lenta que a do NetCDF convertido. O NetCDF CDS original (3 h, zlib em blocos de 19 horários) fez o LECTK levar 128 s em vez de ≈ 12 s; a descompressão de blocos inteiros é hipótese não isolada.

## 4. Desvio de protocolo e diagnóstico (preservado)

O controle positivo **C0** previa CSV idêntico byte a byte entre o arquivo original e o subconjunto de 32 níveis. **A previsão foi refutada** (diferenças de arredondamento: viés ≤ 7 × 10⁻¹⁶), e o driver **não respeitou a regra de parada** pré-registrada: C1–C3 já tinham rodado. O desvio foi registrado no adendo v2 **antes** dos diagnósticos: D1 (repetição de L32, idêntica), D2b (37 níveis, 6 h, sem compressão: idêntico a L32 e diferente de C0) e C0b (repetição de C0, idêntica). Conclusão do adendo: a diferença vem do caminho de leitura do arquivo original, e o instrumento é válido para C1–C3, cujas diferenças estão 10 a 13 ordens de grandeza acima desse ruído.

**Independência:** as métricas de C1–C3 foram calculadas antes do adendo, com critérios pré-registrados que não foram alterados; a validade do instrumento foi estabelecida depois, por diagnósticos definidos no adendo. Os diagnósticos não são independentes do mesmo código e dos mesmos dados. Nada disso constitui validação física: a correspondência entre execuções é verificação numérica. A primeira tentativa de D2 falhou por erro do instrumento e está registrada como inválida.

## 5. Avaliação científica do responsável

O responsável examinou os resultados e aprovou a interpretação: **neste estudo de caso, as diferenças quantitativas entre as configurações verticais examinadas não alteram de maneira substantiva a interpretação física da evolução energética** (mesma evolução temporal, sinais e extremos das três conversões). Isso não significa resultados numericamente idênticos nem adequação de níveis reduzidos a qualquer aplicação. A diferença de ≈ 5 % em Ca **não** fundamenta, por si só, recomendar sempre 37 níveis: a escolha depende do objetivo (evolução qualitativa, magnitudes absolutas, estrutura vertical, fechamentos ou reprodutibilidade), conforme a página de usuário.

## 6. Resultados negativos e o que permanece não demonstrado

- Vantagem geral de velocidade do Zarr: **não demonstrada**.
- DestinE e nova medição CDS: **não realizadas** (credencial ausente; cancelamento seguro de jobs CDS não demonstrável antes da submissão).
- Qual NetCDF (CDS ou ARCO/NCAR) reproduz o GRIB: **desconhecido**.
- Generalidade do viés de Ca (outros casos, estações, regiões, modo `--fixed`): **não examinada**.
- Peso físico da camada 1000–925 hPa onde 1000 hPa está abaixo da superfície: **não verificado**.
- Validação física independente do LECTK: **não demonstrada** (como já registrado em A).

## 7. Defeito temporal e correção (02/11/2026, pré-commit)

**Defeito:** em `process_data`, a validação do passo temporal (dados 6 h > track 3 h) e a do último horário rodavam **antes** da seleção do track aplicada a `--cdsapi`/`--ncar`; com arquivo local o caso nunca era tratado. Na investigação, `--ncar --time-resolution 6` com track de 3 h falhou **depois** de baixar 1,07 GB; reproduzido sem rede (D3).

**Correção (no checkout, não commitada):** novo módulo `src/utils/track_alignment.py`. Quando os dados são amostrados com menos frequência que o track, as posições originais nos horários dos dados são mantidas por correspondência exata de timestamp completo; nada é interpolado, extrapolado ou preenchido e o arquivo do track não é modificado. A mesma regra vale para `infile`, `--cdsapi` e `--ncar`. Continuam sendo erro: duplicatas ou desordem, intervalos irregulares (lacunas), intervalo dos dados que não seja múltiplo inteiro do do track, horários dos dados fora dos horários do track, horários na grade dos dados ausentes (lacuna ou truncamento) e menos de duas posições restantes. Um WARNING único, no terminal e no log, informa intervalos e contagens; `lec_metadata.json` ganha a chave aditiva `track_time_alignment`. Antes de qualquer download, `--cdsapi`/`--ncar` conferem o track contra os horários que serão pedidos. O framework móvel passou a manter as posições do track nos horários dos dados (antes, um filtro `hora % passo` só para `--cdsapi`). Nenhuma equação, referência energética ou tolerância foi alterada.

**Evidências:**
- Offline sintético, `develop` inalterado (`b4a5f77`): local, `--ncar` e `--cdsapi` falham com track de 3 h e dados de 6 h; `--ncar` só depois das 10 requisições.
- Caso real, servidor, sem rede (proxy inválido), arquivo NCAR verificado reaproveitado (sha256 `0571b3a3…`, inalterado): o comando D3 passa a terminar com sucesso, 55 → 27 posições (28 excluídas). As 26 saídas CSV e de track são **idênticas byte a byte** às da execução P2b feita na investigação com o track reamostrado à mão. O mesmo vale para o arquivo local com o track de 3 h.
- Suíte: 19 testes novos (`tests/test_track_alignment.py`), incluindo as três origens com os mesmos dados de 6 h e track de 3 h atravessando a meia-noite (resultados idênticos entre origens). Suíte completa com cache Natural Earth: **591 passed, 1 failed** na correção; **572 passed, 1 failed** em `b4a5f77`. A falha é a mesma nas duas (`test_packaging::test_installed_metadata_matches_pyproject`, metadados do pacote instalado no ambiente local) e preexiste à mudança. CI remota: **PENDENTE** (depende de push/PR).

## 8. Decisões e encaminhamento

- **Não implementar rota Zarr agora.** O ARCO fica documentado como preparação externa para quem aceita 13 níveis; o DestinE só faria sentido, com credencial, para séries longas em áreas pequenas.
- **Get a Subset** permanece o primeiro candidato recomendado para um eventual piloto da Frente F, condicionado à abertura expressa da frente ([E §6](../../fronts/front-E-era5-provider-architecture.md#6-encaminhamento-condicional-para-a-frente-f)). **F continua NÃO ABERTA.**
- Testes que mudariam as conclusões: mais casos; um pedido CDS em GRIB decodificado com ecCodes; DestinE com credencial para o uso específico acima.
- A medição de um caso **não** abre a hipótese K (sensibilidade/convergência versus níveis), que continua em backlog.

## 9. Retenção e evidências externas (GOV-RET-01)

Área externa: área privada da investigação no servidor de processamento (`lectk_zarr_vertical_20261102/`, fora do Git). Relatório `report/RELATORIO_LZV_2026-11-02.md` (sha256 `883d201419270921…`), índice `report/SHA256SUMS_entregaveis.txt` (sha256 `b67bce98aef15dd6…`).

| Rodada | Escopo | Resultado |
|---|---|---|
| 1 (02/11, campanha) | 1902 duplicatas byte a byte (cópia da referência, espelhos ARCO A2/mínimo/interrompido) | 3,44 GiB; 0 falhas; manifesto `RETENCAO_INVENTARIO_v2.json` (sha256 `e292af4b…`) |
| 3 (02/11, após a correção) | subconjuntos controlados e diagnóstico (regenerados bit a bit a partir do arquivo D-Bench antes da exclusão), export `835478e` (idêntico a `git archive 835478e`), arquivo de envio da correção, ambiente Zarr (recriável pelo *freeze*), bytecode | 2287 arquivos, 8,09 GB incluindo as cópias de prova; 0 falhas após uma parada registrada por permissão somente leitura do export, resolvida por `chmod` documentado; manifesto `RETENCAO_MANIFESTO_v3.json` (sha256 `7150f00d…`) |
| Notebook | export de `b4a5f77`, temporários do pytest e cache do Sphinx | 13 608 arquivos, 2,45 GB; 0 falhas |

A área do servidor passou de 7,6 GB para **3,7 GB**. **Retidos temporariamente** até a releitura desta integração em `develop` e, no máximo, até o Marco B: espelho ARCO A1 e NetCDF A1/A3 (bytes brutos e entrada do LECTK), arquivo NCAR com proveniência (entrada da verificação da correção) e export da correção (única cópia, até o commit). **Memória permanente externa:** protocolo, adendo, ledger HTTP, logs, instrumentos, resultados, figuras e execuções LECTK. Reproduzir C0–C3 exige o arquivo CDS retido em `dbench_private/`, cuja retenção pertence à D-Bench. Nada de `dbench_private/`, `LEC_Results/`, Frente D, D-Bench ou do job MPAS em execução foi tocado.

**Estado de retenção:** auditada e executada no escopo verificado; **não** declarar ENCERRADA — RETENÇÃO VALIDADA antes da integração e da releitura em `develop`.
