# Etapa 0 e etapa 4a — rotas ERA5, correções dos clientes NCAR e CDS e kerchunk como padrão do `--ncar`

**Data:** 19/11/2026. **Base:** `develop` em `d65e515` (PR #74); `main` em `48fe8e2`.
**Estado deste registro:** integrado em `develop` pelo PR #75 (merge `c0e3942`); diagnóstico na Swell **confirmado** (23/11/2026); verificação do kerchunk como padrão e encerramento na §6 (23/11/2026).

## 1. Decisões do responsável (19/11/2026)

| ID | Conteúdo |
|---|---|
| D-029 | CDS API, kerchunk e Get a Subset entram oficialmente; NCSS continua disponível como experimental; WEkEO, ARCO e DestinE terão integração nativa na 4b. |
| D-030 | O padrão do `--ncar` passa de NCSS para kerchunk, com o seletor `--ncar-method` da Frente E. Altera a parte "NCSS como padrão" da decisão da Frente E. As referências do kerchunk ficam em `inputs/kerchunk_refs` (ou `--kerchunk-refs`) e são baixadas uma vez quando faltam (escolha do responsável na Etapa 0 desta subetapa). |
| D-031 | ARCO Zarr terá integração nativa; substitui a orientação da LZV de "não implementar rota Zarr agora" (base: D-012). |
| D-032 | Etapa 4 dividida em 4a (correções e padrão), 4b (um cliente por PR, começando pela WEkEO) e 4c (reescrita completa da documentação, por último). |
| D-033 | Proteção do CDS: estado remoto conhecido e encerrado não bloqueia (no máximo uma nova tentativa de download, sem novo pedido); estado desconhecido ou pendente continua bloqueando. |
| D-034 | A Swell é o servidor de processamento preferido, sempre com consulta prévia ao responsável. |

## 2. Inventário das rotas após a 4a

| Rota | Como se usa | Situação |
|---|---|---|
| NetCDF local | `infile` | Fluxo principal (D-001) |
| CDS API | `--cdsapi` (modo protegido opcional `--cds-*`) | Oficial (D-029) |
| NCAR GDEX kerchunk | `--ncar` (padrão) ou `--ncar --ncar-method kerchunk` | Oficial e padrão do `--ncar` (D-029, D-030) |
| NCAR GDEX Get a Subset | `--ncar --ncar-method subset` | Oficial (D-029); exige token GDEX |
| NCAR GDEX NCSS | `--ncar --ncar-method ncss` | Experimental; 0/16 na rodada complementar da 2b |
| WEkEO, ARCO Zarr, DestinE | — | Integração nativa prevista na 4b (D-029, D-031) |

## 3. Mudanças da 4a

1. **Proteção do CDS (`src/utils/cds_guard.py`).** Novos estados no ledger: `download_retry` (nova tentativa de
   download do mesmo resultado) e `download_failed` (terminal, `blocks_new_submissions: false`). Os estados
   `failed`, `abandoned`, `interrupted` e pedidos sem estado final continuam bloqueando. Entradas antigas de ledgers
   já existentes não são reinterpretadas.
2. **Get a Subset (`src/utils/gdex_sources.py`).** Pedidos "Set for Purge" deixam de contar no limite de 7 pedidos
   listados; continuam impedindo um pedido idêntico (deduplicação do GDEX).
3. **Seletor `--ncar-method {kerchunk,subset,ncss}`**, válido só com `--ncar`; kerchunk é o padrão, também no
   comando separado `python -m src.utils.ncar_gdex`. Referências do kerchunk: procuradas na pasta; baixadas uma vez
   quando faltam, só as das variáveis e décadas do pedido; arquivo inválido ou parcial não é guardado; a
   proveniência registra a origem de cada referência.
4. **Limite por arquivo do Get a Subset (§4).** Um arquivo pode chegar com até 2× o tamanho listado
   (`SUBSET_SIZE_CEILING`) mais `HEADER_ALLOWANCE`, contando o enquadramento dos blocos. Maior que o listado:
   aviso no log e `listed_bytes` no estado do pedido, e o arquivo é aceito. Menor que o listado: arquivo
   incompleto, baixado de novo. Arquivo que não começa com assinatura NetCDF/HDF5: baixado de novo e, se
   persistir, recusado. O conteúdo (variável, tempos, níveis, latitudes e longitudes) continua conferido na
   leitura. Escolha do responsável: a diferença de tamanho deixa de interromper; fica só um teto alto.

## 4. Diagnóstico do limite por arquivo do Get a Subset (casos 19790205 e 19790585)

**Confirmado** em 23/11/2026, com um pedido real aprovado pelo responsável (só o caso 19790205, segmento
oeste).

**Proveniência.** Registros da 2b lidos sem alteração em
`/p1-swell/danilocs/lectk_bench2b_round2_20261113` (`logs/campaign_rep{1,2,3}_gas_19790205.log`,
`manifests/campaign/rep1/gas/19790205_case.json`, estados `subset_request*.json` das três repetições).
Pedido de teste: mesmo controle do segmento oeste da 2b (`d633000`, `197903181200/to/197903270000`,
T/Z/W/U/V, −34,75 a −70, −71,75 a −0,25), pedido GDEX 989514. Script e dados em
`/p1-swell/danilocs/lectk_stage4a_diag_20261122/` (fora do repositório; ≈ 101 MB); só o arquivo T de
19/03/1979 foi baixado, por um soquete TLS bruto com teto de 200 MB.

**Na 2b.** As três repetições falharam no mesmo arquivo
(`…128_130_t.ll025sc.1979031900_1979031923.nc`), quatro tentativas cada, com
`_ResponseTooLarge: more than 52853805 bytes read` (tamanho listado + 64 KiB); o arquivo da véspera (25 MiB)
passou.

**Medido.**

| Grandeza | Valor |
|---|---|
| Tamanho listado pela API | 52 788 269 B |
| Corpo recebido | 52 788 269 B (igual ao listado) |
| `Content-Length` | ausente; `Transfer-Encoding: chunked` |
| Blocos | 9 665 |
| Enquadramento dos blocos | 74 104 B |
| Linha de estado e cabeçalhos | 530 B |
| Total na conexão | 52 862 903 B (listado + 74 634 B) |
| Limite do LECTK na 2b | 52 853 805 B (listado + 65 536 B) |
| Assinatura do corpo | HDF5 (`89 48 44 46`); abre com netCDF4: T, 24 × 37 × 142 × 287 |

**Conclusão.** O arquivo entregue está certo. A falha vem da contagem: a conexão do LECTK conta o
enquadramento dos blocos, e na rota Get a Subset a folga era só `HEADER_ALLOWANCE`. O enquadramento cresce
com o tamanho do arquivo (≈ 0,14 %) e passa de 64 KiB perto de 45 MB, o que explica por que o arquivo de
25 MiB passou. As hipóteses "arquivo maior que o listado" e "cabeçalhos maiores que a folga" foram
descartadas para este arquivo. O caso 19790585 não foi testado; a mesma causa é **[INFERIDA]**.

## 5. Limitações

- O endereço do servidor de arquivos das referências do kerchunk
  (`https://tds.gdex.ucar.edu/thredds/fileServer/files/g/d633000/kerchunk/pressure_levels/`) foi derivado do
  catálogo público do GDEX e **ainda não foi conferido** contra o serviço real.
- Todos os testes novos são offline, com respostas simuladas. Só o arquivo T de 19/03/1979 do caso 19790205
  foi medido de novo; a recuperação completa do caso pela rota Get a Subset não foi refeita.
- O pedido de teste 989514 não pôde ser apagado (a conexão da Swell com `gdex.ucar.edu` deu timeout duas
  vezes no pedido de remoção); ficou listado como `Completed` e expira pelas regras do GDEX.
- Os testes de "antes" da correção do limite por arquivo rodaram com o `gdex_sources.py` da `develop`.
- Os testes de mapas não rodaram localmente (sem os arquivos Natural Earth); ficam para o CI.

## 6. 23/11/2026 — Verificação do kerchunk como padrão e encerramento da 4a

**Endereço das referências (corrige a §5).** O endereço derivado em 19/11
(`…/thredds/fileServer/files/g/d633000/kerchunk/pressure_levels/`) respondeu **HTTP 404** no serviço real. O catálogo
`catalog_d633000.xml` do GDEX lista as referências em `files/d633000/kerchunk/pressure_levels/`, sem o `g/`; o
endereço `https://tds.gdex.ucar.edu/thredds/fileServer/files/d633000/kerchunk/pressure_levels/` respondeu
**HTTP 200** (`temperature.1970-1979-remote.json`, 14 172 501 B), do Mac e da Swell. Com autorização do responsável,
`KERCHUNK_REFS_URL` foi corrigido, com um teste que falhava antes da correção e passa depois.

**Teste real na Swell (autorizado em 23/11/2026).** Caso 19790205, mesmo pedido selado da rodada 2 da 2b (caixa,
13 níveis, 35 horários), método padrão do `--ncar` (`DEFAULT_NCAR_METHOD`, kerchunk), referências numa pasta vazia,
teto de 12 GiB, um download por vez. Código: `c0e3942` mais a correção do endereço. Pasta:
`/p1-swell/danilocs/lectk_stage4a_kerchunk_20261123/` (fora do repositório; ≈ 1,7 GB).

| Grandeza | Valor |
|---|---|
| 1ª tentativa (11:55–12:09 UTC) | 3 de 5 referências baixadas; a de `u` falhou 4 vezes (1 resposta incompleta e 3 timeouts de conexão da Swell, `Errno 110`); parou antes de qualquer dado do caso |
| 2ª tentativa (12:22–12:52 UTC) | as 3 referências reaproveitadas; `u` e `v` baixadas do endereço corrigido; recuperação completa |
| Referências | 5 arquivos, 14,2 a 28,2 MB; sha256 iguais aos da rodada 2 da 2b |
| Tempo de recuperação (2ª tentativa) | 1 560,5 s (2b, rodada 2, rep. 1: 1 914,4 s) |
| Bytes recebidos (dados, contados pelo `ByteBudget`) | 10 889 525 245 B em 175 pedidos HTTP (igual à 2b) |
| Arquivo de grade comum | sha256 `567ccdce…`, **idêntico** ao da 2b |
| LEC (`-f -r`, mesma caixa e namelist da fase 2) | rodou até o fim; CSV de resultados com sha256 `ac5c7eec…`, **idêntico** ao da rota kerchunk da fase 2 (35 × 26 valores, diferença máxima 0) |

**Conclusão.** PASS: o kerchunk como padrão baixa as referências do endereço corrigido, monta o arquivo e
reproduz bit a bit a saída da 2b, dentro da tolerância do pré-registro da fase 2 (classe "idêntico").

**Achados para a 4b.**

- O `--ncar` do programa principal usa o teto fixo de 8 GiB (`DEFAULT_MAX_TOTAL_BYTES`) e não tem `--max-gib`;
  recusaria este caso (≈ 10,1 GiB). O teste usou o mesmo cliente com teto de 12 GiB.
- A conexão da Swell com `tds.gdex.ucar.edu` falhou de forma intermitente (respostas incompletas e timeouts); as 4
  tentativas por referência não bastaram na primeira vez. Nenhum arquivo parcial foi guardado.

**Notas do responsável (23/11/2026).**

- A folga de 2× do Get a Subset (`SUBSET_SIZE_CEILING`) fica mantida por escolha consciente; a causa medida (§4)
  pedia ≈ 0,14 %.
- O pedido de teste 989514 (§4, §5) foi apagado pelo responsável.
- O destino de `/p1-swell/danilocs/lectk_stage4a_diag_20261122/` (≈ 101 MB) e de
  `/p1-swell/danilocs/lectk_stage4a_kerchunk_20261123/` (≈ 1,7 GB) fica para o fim da etapa 4, junto com os dados
  da 2b (D-028).

## Fora deste registro

Aprovação, merge, CI e releitura deste documento não podem constar dele.
