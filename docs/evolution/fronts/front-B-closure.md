# Frente B — Encerramento por triagem de necessidade

**Decisão do responsável:** 24/10/2026. **Estado: ENCERRADA — TRIAGEM CONCLUÍDA, SEM EXECUÇÃO DE NOVAS REFERÊNCIAS.**

## Objeto e alcance da decisão

A Frente B foi prevista condicionalmente para tratar **lacunas de referências numéricas e testes de regressão do próprio LorenzCycleToolkit**, quando indispensáveis às comparações posteriores. A Frente A já havia sido encerrada com ressalvas aceitas. A avaliação desta frente foi **documental e restrita à necessidade atual**; não incluiu nova execução de cálculos ou testes.

O responsável decidiu **manter a falta de uma referência numérica NCEP-R2 como pendência documentada** e encerrar a Frente B caso não houvesse outra lacuna que exigisse tratamento imediato. A triagem na `develop` identificou os materiais a seguir.

## Evidências e lacunas

| Tema | Evidência na `develop` | Decisão |
|---|---|---|
| Referências de regressão ERA5 | **49 arquivos** em `tests/data/reference/`, documentados por `tests/data/reference/MANIFEST.tsv` e `PROVENANCE.md`: **24 para caixa fixa**, **25 para trajetória**; usados por `tests/test_ERA5_fixed.py` e `tests/test_ERA5_track.py` | **Suficientes para a continuidade documental da Frente C**, nos cenários efetivamente cobertos. Não equivalem a validação física |
| NCEP-R2 fixo e trajetória | `tests/test_R2_fixed.py` e `tests/test_R2_track.py` checam estrutura, dimensões, termos e finitude; `PROVENANCE.md` declara ausência de referência numérica aprovada para esses casos | **PENDENTE DOCUMENTADA**, não resolvida nem apresentada como PASS de regressão numérica; sem novos resultados nesta etapa |
| ERA5 37 níveis/3 h, ciclone real e validação física independente | Limites científicos já explicitados em [Frente A](front-A-closure.md) e [visão geral](../overview.md) | **Continuam não demonstrados**, mas não são automaticamente tarefas da Frente B; avaliar em escopos científicos futuros apropriados |

## Critério de fechamento

A triagem não encontrou **outra lacuna de referência que imponha produção de novos artefatos numéricos antes do levantamento documental de provedores ERA5 (Frente C)**. Isso **não** significa ausência universal de lacunas, suficiência para qualquer benchmark ou validação científica do NCEP-R2.

**Veredito administrativo:** Frente B **ENCERRADA POR DISPENSA DE EXECUÇÃO NESTE MOMENTO**, com a pendência NCEP-R2 explicitamente aceita pelo responsável. Não houve alteração de código, equações, resultados, testes ou datasets para encerrar B; não se aplica um resultado PASS de testes inexistentes.

## Gatilho de reavaliação e sucessão

- **Antes de definir benchmarks da Frente D**, decidir, de forma explícita, se as comparações planejadas exigem uma linha de base numérica NCEP-R2 ou novas referências; se sim, abrir trabalho específico com protocolo aprovado, sem reclassificar retroativamente o NCEP-R2 como validado.
- **Frente C é a próxima na fila**: investigação documental dos provedores ERA5 e requisitos científicos; este registro não autoriza automaticamente mudanças de arquitetura, novos provedores, processamento oneroso ou alteração do fluxo `infile`.
- Validação física independente e novos cenários científicos permanecem temas separados e não devem ser confundidos com a regressão computacional aqui tratada.

**Proveniência:** decisão expressa do responsável de 24/10/2026 e inspeção dos documentos e testes efetivamente versionados em `develop` no commit `61d37b5f3f860c89c932052979696fd42e80ba36` (estado anterior à formalização documental). `main` não foi alterada por esta triagem.
