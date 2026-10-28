# Frente D — Rotas alternativas de aquisição ERA5 e adaptador NCAR GDEX experimental

**Período:** 26–28/10/2026. **Estado:** rota NCAR **integrada em `develop` como funcionalidade EXPERIMENTAL**
([PR #38](https://github.com/daniloceano/LorenzCycleToolkit/pull/38), merge
`dec13fe6dcacc2e2cd7f5720160ea40d40e75403`, 28/10/2026). **Comparação de tempo com o CDS: não realizada**
(P2-B0 bloqueada antes de qualquer submissão).

**Encerramento formal (28/10/2026): PASS COM RESSALVAS ACEITAS.** Estado operacional: **ENCERRADA — RETENÇÃO
VALIDADA**, com limpeza parcial e retenções temporárias documentadas até a próxima versão/tag estável (ver
"Encerramento formal" e "Retenção remanescente"). Os benchmarks controlados previstos para D **não** foram
concluídos.

## Objeto

Avaliar rotas de aquisição de ERA5 que reduzam a dependência da fila de jobs do CDS, **sem alterar o núcleo
científico**, o arquivo local `infile` como fluxo principal (D-001) nem a integração opcional `--cdsapi` (D-002). A
fila aprovada (D-010) previa para D benchmarks controlados. O trabalho executado cobriu um levantamento de rotas,
um piloto, uma tentativa de comparação com o CDS e um adaptador opcional; **não houve campanha de benchmarks**.

## Cronologia e resultados

| Etapa | Data | Resultado |
|---|---|---|
| P1 — levantamento sem download | 26/10 | Rotas consideradas. ERA5 no Google (ARCO) exige baixar blocos globais (~4,25 GB para o caso E1) e usa unidade de nível que o LECTK não converte. OPeNDAP via `curl` exigiria um decodificador novo. Escolha: NCAR GDEX (d633000) pelo serviço THREDDS NCSS, que recorta a área no servidor. |
| P2-A — piloto NCSS autorizado | 26/10 | 78,80 MiB, custo zero. Os campos do caso C0 reobtidos ficaram idênticos bit a bit à amostra versionada; 49 regressões passaram. O caso E1 (2018, 13 e 37 níveis) rodou nos modos fixo e móvel. **PASS OPERACIONAL LIMITADO**, não PASS da frente. |
| P2-B0 — um caso NCAR × CDS | 26/10 | **BLOQUEADA, nenhum job criado** (ver "Resultado negativo"). |
| P3 — protótipo offline | 26/10 | Opção `--ncar` e comando `python -m src.utils.ncar_gdex`; a reprodução das respostas da P2-A gerou valores e CSV do LECTK idênticos. |
| P3-L — validação ao vivo, teto de 100 MiB | 27/10 | Correção do contador de bytes (cobra tentativas falhas e repetições). **PARCIAL 6/10**: o servidor parou de enviar aos 89 % de uma resposta, e o teto não comportava terminar. As 6 respostas recebidas ficaram idênticas bit a bit à P2-A. |
| P3-L2 — retomada, teto de 85 MiB | 27/10 | 4/4 respostas na primeira tentativa (39,4 MB, 69 s). Detalhes abaixo. |
| Ajustes e verificação pré-PR | 27/10 | Contagem dos bytes na conexão e trava contra duas recuperações do mesmo arquivo. Um teste adverso revelou e permitiu corrigir um **defeito real**: o prazo total da tentativa não era imposto dentro de uma leitura lenta (409 s contra 14 s permitidos, em simulação; 9,7 s contra 1 s, em loopback). |
| Revisão visual e integração | 27–28/10 | Diffs e documentação v4 aprovados pelo responsável; commit documental `348bb42` e abertura do PR #38 em 27/10; merge em 28/10. CI verde no PR (Actions push e PR: 417 testes; CircleCI) e em `develop` após o merge. |

**Resultados da P3-L2,** sobre o caso E1 de 09–10/08/2018, com 5 campos, 37 níveis e 129×129 pontos:

- arquivo com valores idênticos bit a bit à P2-A em 15/15 comparações no domínio comum de 41×41;
- 4 comandos congelados do LECTK com 24/24 CSV idênticos;
- reaproveitamento com `--ncar` sem nenhum pedido.

## O que foi integrado (PR #38)

| Arquivo | Mudança |
|---|---|
| `src/utils/ncar_gdex.py` | novo |
| `lorenzcycletoolkit.py` | +13: opção `--ncar`, exclusiva com `--cdsapi` |
| `src/utils/preprocessing.py` | +19/−4: escolha da fonte, mesmo tratamento de trajetória do `--cdsapi` |
| `tests/ncss_emulator.py` | novo: emulador sem rede |
| `tests/test_ncar_gdex.py` | novo: 115 casos |
| `docs/source/ncar_gdex_example.rst` | novo, mais a entrada no índice |

**Comportamento do adaptador:**

- pede a mesma área, níveis e horas do `--cdsapi`;
- confere cada resposta antes de usá-la;
- copia os valores sem alteração, em float32;
- só publica o arquivo completo e verificado, com registro de proveniência;
- retoma downloads interrompidos e reaproveita um arquivo existente apenas se a proveniência o confirmar;
- limita bytes, tempo por tentativa e repetições, e impede duas recuperações simultâneas do mesmo arquivo.

**Diretrizes aprovadas pelo responsável (27/10):**

- float32 nativo como padrão;
- float64 opcional apenas no comando independente;
- `--ncar` como opção experimental mais o comando independente;
- fluxo CDS, NetCDF local, núcleo científico e esquema de `lec_metadata.json` preservados.

## O que cada resultado demonstra, e o que não demonstra

| Conceito | Situação | Base |
|---|---|---|
| **Reprodutibilidade do LECTK** | **Demonstrada no caso E1.** Com os dados NCAR, os mesmos comandos reproduziram byte a byte os CSV da P2-A, e o reuso via `--ncar` reproduziu os da execução de base sobre o mesmo arquivo. | P2-A, P3-L2 |
| **Validação técnica do adaptador** | **Demonstrada no escopo testado:** requisição igual à do `--cdsapi`, verificação das respostas, integridade e proveniência, limites de bytes e de tempo, retomada, trava e reuso. Suíte e CI verdes. A operação ao vivo foi validada num único caso e numa única rede. | P3, P3-L, P3-L2, ajustes, PR #38 |
| **Equivalência direta com o CDS** | **Não demonstrada.** Os valores NCAR foram comparados só com outra cópia NCAR (P2-A), nunca com NetCDF do CDS. O tempo também não foi comparado (P2-B0 bloqueada). | P2-B0 |
| **Validação física independente** | **Não demonstrada.** Igualdade bit a bit, regressões e CI verificam computação e reprodução, não a física dos resultados. | [Visão geral](../overview.md) |
| **Conclusão técnica × encerramento** | **Conclusão técnica:** produto integrado em `develop` com CI verde. **Encerramento formal:** PASS COM RESSALVAS ACEITAS em 28/10/2026. **Encerramento operacional (GOV-RET-01):** retenção validada, com limpeza parcial e retenções temporárias até a próxima versão/tag. Benchmarks e comparação com o CDS **não** fazem parte do que foi concluído. | Este registro |

## Resultado negativo: comparação com o CDS (P2-B0)

A P2-B0 deveria medir, no caso E1, o tempo até um diagnóstico utilizável pelo NCAR e pelo `--cdsapi`. Ela foi
**bloqueada antes de qualquer submissão**, por dois motivos:

- **sem cancelamento seguro:** não havia meio autorizado de cancelar apenas os jobs do teste. O cliente CDS
  instalado não implementa cancelamento, e as alternativas estavam vedadas;
- **conta ocupada:** a vaga de requisição por usuário estava tomada por outra campanha.

**Consequências:**

- Não existe medição CDS do caso. A regra de decisão pré-aprovada **não permite** afirmar que o NCAR é mais rápido
  nem recomendar o adaptador com base em desempenho.
- Os números disponíveis são só contexto, não comparação. O lado NCAR do E1 levou cerca de 100 s, entre aquisição,
  preparo e LECTK. Os jobs CDS da mesma conta, de outro dataset, esperaram de 5h55 a 7h41 naquele dia.
- Destravar a comparação exige decisão do responsável (forma de cancelamento e janela sem a outra campanha) e nova
  autorização.

## Limitações e o que não foi demonstrado

- **Experimental:**
  - validação ao vivo em **um único caso** e uma única rede;
  - equivalência com NetCDF do CDS **não medida**;
  - nenhuma validação física. Testes verdes e igualdade com a cópia NCAR não são validação física (ver
    [visão geral](../overview.md)).
- **Datas:** o GDEX publica ERA5 meses depois do CDS (em 26/10/2026, cobertura até 31/07/2026). Datas recentes
  continuam exigindo `--cdsapi`.
- **Precisão:** o LECTK calcula no tipo da entrada. float32 contra float64 muda os resultados em até 2,7×10⁻⁷ da
  escala de cada saída, sem troca de sinal. As referências de regressão foram calculadas com entrada float64.
- **Contagem de bytes:**
  - cobre os bytes HTTP lidos depois do TLS;
  - os bytes que podem ter chegado sem ser lidos entram apenas como cota superior (buffer de recepção);
  - TLS, resposta de proxy e cabeçalhos TCP/IP não são observados.
- **Prazo:** a conexão inicial (TCP e TLS) fica limitada só pelo timeout de cada passo.
- **Trava:** sem `flock` (Windows), a trava é a criação exclusiva do arquivo de trava; não foi testada em Windows.
- **Volume:** o NCSS entrega NetCDF-3 sem compressão; em domínios grandes trafega mais bytes que o CDS comprimido
  (estimativa, não medida).
- **Pendência herdada:** NCEP-R2 continua sem referência numérica (Frente B); não foi necessária para a
  equivalência ERA5 desta frente.

## Evidências essenciais (fora do Git)

Os relatórios completos, ledgers, respostas brutas, pré-registros selados e prévias ficam na área de trabalho
externa da frente, sem caminhos pessoais aqui. Um inventário de retenção, segundo a
[GOV-RET-01](../retention-policy.md), foi preparado fora do Git em 28/10. Hashes SHA-256 dos documentos principais:

| Documento | SHA-256 |
|---|---|
| Relatório de preservação P2-A (GOV-RET-01) | `8123d86010089f47a0e90fd8257ec95a2c4b6b21333fca94c594c81057c1b90d` |
| Relatório P2-B0 (bloqueio) | `28d0c68e90be79248032dd587eaea3124233aa55cf10548006224f196d277949` |
| Relatório P3 (protótipo) | `27df8384bdef2ca40bd6d77ebb2630acf6db6312c065895261820f01b573cf5a` |
| Relatório P3-L + manifesto | `d03e2e68bf99c4b225d342ab672f69ed080f6c02c78461e3e7226e72450ef499` / `ba08d4649b39bb39da2094d9231479a1bc614f53903da40b3b0ecf321734352b` |
| Relatório P3-L2 + manifesto | `461abb21d1097cb543778c6d236ad0625d60e71316a1c79ddd7e8cded53581ec` / `e4ce4d466e5c8227344baf92455c0cd69cfe1ac50201f42728527e47aa38643e` |
| Relatório de ajustes finais + manifesto | `a9b8731f4de39dd706226e890504cc5cbdc969f6964b71558eef8b9ead1cef51` / `ba2991c4df49d4402796e7e9daef9ccfc004ae01cb06fb5648241ac072a43143` |
| Relatório de verificação pré-PR + manifesto | `27f10381bb2705d542c1c1fb20a542b7c56bed1a631f3e9d8261e78a60cc9d41` / `7668c2f014aa218e3991da626bc2f3c797eb7c91f94f469c931342eb8d77ab49` |
| PDF da documentação v4 aprovada | `16ac926fc0bb17a321e6c1abd01132b4779d01d41a3f7ea8e6a306894179099b` |

Os hashes provam a identidade dos documentos, não backup integral dos dados.

## Processo

Os 13 primeiros commits da branch precederam a revisão visual do responsável. A ocorrência está no
[registro de 27/10](../evidence/records/2026-10-27-front-D-precommit.md); o histórico foi preservado.

## Limpeza de 28/10 (GOV-RET-01, Marco A)

Executada com autorização expressa do responsável para os grupos G1, G3, G4 e G6. Submanifesto de 838 caminhos
(sha256 `b1ebb93e29ee8740176b9b5074fda17cbe91666910ab9e39c747bf420214789a`), executor com conferência de conteúdo
imediatamente antes de cada exclusão, remoção sem seguir links e registro durável (executor sha256
`2094a759f84e0661353a7152becba5c81debb32a43e92dbcba74ef0442c74c5e`).

| Grupo | Conteúdo | Resultado |
|---|---|---|
| G1 | temporários de pytest e sandbox | 35/35 caminhos, 49 956 arquivos removidos |
| G3 | builds de wheel e venvs recriáveis | 20/20 caminhos, 6 287 arquivos removidos |
| G4 | arquivos das cópias do LECTK idênticos ao blob do commit de origem | 782/782 caminhos, 4 904 arquivos removidos |
| G6 | exportações da revisão visual de 27/10 | **falha de permissão**: subpastas exportadas como somente leitura; 5 de 16 arquivos removidos, 11 permanecem |

- **Resultado:** PARCIAL (837/838), código de saída 3. A falha não foi contornada; a pasta G6 restante é recriável com
  `git show` de commits integrados e não contém evidência. Por decisão do responsável, a exclusão não foi repetida,
  as permissões não foram alteradas e não se gerou novo manifesto para recuperar cerca de 0,2 MiB: o resto de G6
  fica como **resíduo operacional recriável**, a reavaliar na próxima manutenção de retenção. O registro durável da
  execução e os hashes foram preservados.
- **Espaço:** estimativa de 10,98 GiB para os três grupos concluídos; o espaço livre medido no volume subiu cerca de
  11,0 GiB (medida global do disco, sujeita a outras atividades).
- **Preservados e conferidos:** G2 e G5 intactos; dados P2-A, P2-B0, P3-L e E1; produtos C0 e E-R9; `LEC_Results` e
  `inputs/namelist` das cópias; logs; ambiente de trabalho; os quatro `SHA256SUMS` das evidências sem divergência;
  worktrees, branches e stash.

## Retenção remanescente (GOV-RET-01)

Inspeção somente leitura de 28/10, com hash de conteúdo de cada conjunto registrado num inventário externo de
integridade. Marco de reavaliação de todos os itens: **próxima versão/tag estável** (Marco B).

| Conjunto | Tamanho | Classe | Finalidade, dependência e recriabilidade | Cópia independente |
|---|---:|---|---|---|
| G2 — 13 cópias de worktree idênticas a commits integrados | ≈ 0,59 GiB | DESCARTAR — CANDIDATO (próxima manutenção) | código igual, arquivo a arquivo, a commits em `develop`; resultados nos logs retidos | não exigida |
| G2 — 9 cópias com estados intermediários não commitados (versões de `ncar_gdex.py` e testes antes de correções, controle negativo, scripts de cobertura) | ≈ 0,42 GiB | RETER TEMPORARIAMENTE | únicos registros desses estados de código; os resultados estão nos logs retidos; não se apagam automaticamente | não exigida: os estados finais estão integrados e os resultados registrados |
| G5 — rodadas superadas do protótipo P3 | 0,79 GiB | RETER TEMPORARIAMENTE | superadas pelas rodadas `final2` (referência), mas **não idênticas** a elas arquivo a arquivo; reprocessáveis a partir das respostas P2-A | não exigida |
| P3-L — partes brutas da validação ao vivo | 56,4 MiB | RETER TEMPORARIAMENTE; DESCARTAR — CANDIDATO no Marco B | as 12 partes são idênticas byte a byte às partes da P3-L2; mantidas para que o manifesto da P3-L continue completo | a cópia da P3-L2 já existe |
| P3-L2 — `E1_P3L2_float64.nc` e `E1_13_P3L2.nc` | 58 MiB | RETER TEMPORARIAMENTE | análises citadas no manifesto da P3-L2; recriáveis a partir do E1 | não exigida |
| P2-A e E1 da P3-L2 — respostas brutas NCAR | ≈ 0,25 GiB | RETER TEMPORARIAMENTE | base das igualdades bit a bit; uma nova obtenção do serviço público reproduziu os mesmos bytes (P3-L 6/6, P3-L2 15/15), e os hashes estão nos relatórios | não exigida: recriabilidade demonstrada |
| Resto de G6 | 0,2 MiB | resíduo operacional recriável | exportação por `git show` de commits integrados | não exigida |

Nenhum item INVESTIGAR bloqueia o encerramento: a suspeita de duplicação P3-L × P3-L2 foi resolvida (duplicata
confirmada), e os NetCDF de análise foram classificados.

## Encerramento formal (28/10/2026) — PASS COM RESSALVAS ACEITAS

**O que se conclui:** o adaptador NCAR GDEX está integrado em `develop` como funcionalidade **experimental** e foi
testado no escopo descrito acima: suíte, emulador, controles negativos, validação ao vivo de um caso, reprodução
bit a bit dos CSV do LECTK no caso E1 e CI verde.

**Ressalvas aceitas:**

1. Os **benchmarks controlados** originalmente previstos para D (D-010) **não foram concluídos**.
2. A **comparação direta NCAR × CDS não foi realizada** (P2-B0 bloqueada): não há equivalência demonstrada com o
   NetCDF do CDS.
3. **Nenhuma superioridade de velocidade** foi demonstrada; não há recomendação baseada em desempenho.
4. **Nenhuma validação física independente.**
5. Permanecem as diferenças de precisão (float32 × float64, até 2,7×10⁻⁷ da escala), a defasagem de cobertura
   temporal do GDEX, a validação ao vivo de um único caso e uma única rede, e as limitações de contagem de bytes e
   de prazo descritas em "Limitações".
6. Limpeza de retenção parcial (G6) e retenções temporárias até a próxima versão/tag.

## Pendências após o encerramento

- Comparação NCAR × CDS e benchmarks controlados: exigem nova frente ou decisão do responsável, com forma segura de
  cancelamento de jobs e uma janela sem outra campanha na conta.
- Reavaliação das retenções no Marco B (próxima versão/tag estável).
- Remoção do worktree e dos branches local e remoto da Frente D: procedimento próprio, não executado.
- A Frente E (decisão de arquitetura de provedores) **não foi aberta**.
