# D-Bench — diagnóstico complementar pós-encerramento das divergências energéticas CDS × NCAR (01/11/2026)

**Estado científico: CONCLUSÕES ACEITAS pelo responsável em 01/11/2026; investigação complementar encerrada, sem correção de código.** Este registro é **posterior à campanha**, não modifica seu pré-registro nem seu [encerramento histórico](2026-11-01-dbench-closure.md). A Frente D e a D-Bench permanecem encerradas; **Frente E não aberta** e fila principal inalterada.

## Escopo, método e natureza da evidência

Investigação **exclusivamente offline**, com dados já preservados, dos quatro casos que continuavam sem causa energética demonstrada: `19791008`, `19800256`, `19790166` e `19791150`. O caso `19790644` foi usado como controle positivo do mecanismo de temperatura; casos equivalentes serviram de controles adicionais.

A documentação abaixo sintetiza o **relatório produzido pelo agente no servidor e aprovado pelo responsável**, e não afirma que a conversa de coordenação tenha reexecutado independentemente os artefatos externos. O agente relatou:
- 26/26 NetCDF comparados com os hashes SHA-256 da campanha, sem divergências; tracks e namelist conferidos.
- 14 reexecuções do cálculo LEC (`-t -r --time-resolution 3`), idênticas byte a byte às saídas originais, e 22 arquivos híbridos de substituição controlada; 36 execuções LEC no total, conforme inventário externo.
- Verificação de que os arquivos brutos Get a Subset do GDEX, nos casos verificados, eram idênticos aos arquivos montados pelo cliente experimental; os métodos NCAR haviam produzido campos idênticos na Etapa B.
- Reimplementação de comparação de colunas e reconstrução independente de estabilidade estática (σ) e energia potencial disponível zonal (Az), com concordância numérica relatada.

**Proveniência de código:** `develop` em `857b3ca84714ee03afaf960d67327cd470da7cc6`; o agente informou que o núcleo científico, a CLI e o preset de namelist da cópia congelada usada na campanha (`dc4ed748` com alterações locais) coincidiam byte a byte com a versão examinada. Ambiente relatado: Python 3.12.12; NumPy 2.0.0; xarray 2024.2.0; pandas 2.2.2; MetPy 1.6.2; netCDF4 1.7.1.post1.

## Resultados por caso

A razão abaixo é **máxima diferença energética / tolerância oficial da coluna**; **não** representa porcentagem da energia física. Todas as diferenças meteorológicas já haviam passado nas tolerâncias de campo pré-registradas.

| Caso | Razão máxima histórica | Componentes energéticas fora da tolerância | Conclusão aceita |
|---|---:|---|---|
| `19791008` | 1,24 × | ∂Az/∂t | A temperatura NCAR × CDS explica a divergência; sem defeito de processamento identificado. |
| `19800256` | 3,48 × | BAz, ∂Az/∂t, ∂Ae/∂t, RGz, RGe | Mesmo mecanismo de temperatura e propagação nas tendências/resíduos; sem defeito identificado. |
| `19790166` | 5,18 × | Az, BAz, ∂Az/∂t, ∂Ae/∂t, RGz | Mesmo mecanismo; sem defeito identificado. |
| `19791150` | 6,66 × | Az, Ae, Ca, BAz, Gz, ∂Az/∂t, ∂Ae/∂t, RGz, RGe | Mesmo mecanismo energético, mais disseminado; sem defeito identificado. |

**Experimentos causais relatados para os quatro casos:**
1. Substituir **somente a temperatura NCAR pela temperatura CDS** tornou a energética equivalente à referência CDS; a troca inversa (**CDS com temperatura NCAR**) reproduziu a divergência.
2. Trocar somente o geopotencial não resolveu o desacordo.
3. Substituir temperatura **e** geopotencial NCAR pelos campos CDS produziu saídas energéticas **idênticas byte a byte** à referência CDS; nenhuma contribuição relevante do geopotencial foi identificada nesse escopo.
4. O mecanismo já observado em `19790644` também foi reproduzido como controle. As diferenças dos ventos e da velocidade vertical não explicaram os desvios na amostra comparada.

## Interpretação científica

**Origem nos arquivos comparados (representação):** o relatório constatou diferenças muito pequenas na temperatura disponibilizada pelas duas rotas, até aproximadamente `2,4 × 10⁻⁴ K` no conjunto original, frequentemente próximas de deslocamentos uniformes por fatia horário–nível, e já presentes antes da montagem pelo cliente NCAR. O agente encontrou evidência numérica de deslocamentos subquantização e modelou as fatias divergentes em três casos (`19791008`, `19800256` e `19790166`) por um deslocamento escalar ε reproduzindo bit a bit suas temperaturas após arredondamento float32.

**Limite de atribuição:** o GRIB1 original e o algoritmo histórico exato de conversão NCAR/CISL não estavam disponíveis nesta investigação. Portanto, o modelo de representação é fortemente sustentado para os arquivos analisados, mas **não constitui prova direta de qual operação da conversão histórica introduziu o desvio** nem de que o CDS reproduz exatamente o GRIB1 original. Para `19791150`, o resultado causal energético foi relatado, mas o registro sintético não documenta validação bit a bit da modelagem ε para todas as fatias. Não tratar essas ressalvas como defeito confirmado.

**Amplificação energética:** pequenas diferenças na temperatura média em níveis vizinhos alteram o gradiente vertical e a estabilidade estática σ. Componentes que dependem de `1/σ` (incluindo Az) podem responder de forma relativamente maior. Variações temporais das pequenas diferenças de entrada se refletem nas derivadas por diferenças finitas; resíduos energéticos propagam essas diferenças algebricamente. O agente relatou reconstrução independente de σ e Az que explicou a discrepância de ∂Az/∂t com erro inferior a 0,1% no teste realizado. Este mecanismo é compatível com a matemática implementada, não caracteriza, por si só, defeito.

**Magnitude física e limites:** o agente relatou discrepâncias inferiores a 0,3% da magnitude média das componentes e a 0,23% de sua variabilidade temporal nas métricas examinadas. São métricas **diagnósticas posteriores** à campanha, não seus critérios oficiais. Isso sugere impacto quantitativo pequeno nos casos avaliados, mas **não comprova irrelevância física em toda aplicação**, especialmente ao analisar tendências, sinais, compensações e fechamentos de balanço. A divergência de coluna originalmente registrada continua válida.

**Ressalva secundária não bloqueante:** o relatório do agente apresentou `0,44 ×` para o controle `19800296`, enquanto a documentação histórica da Etapa B registra `0,66 ×`. A discrepância auxiliar não foi reconciliada; o **valor oficial de 0,66 × permanece inalterado** e a classificação do caso permanece equivalente. Não usar o número auxiliar para alterar o encerramento da campanha.

## Decisão aprovada e consequências

1. **Causa das quatro divergências energéticas:** diferenças de temperatura entre os arquivos CDS e NCAR e sua amplificação numericamente consistente no cálculo da estabilidade estática, componentes energéticas, tendências e resíduos. A hipótese inicialmente remanescente do caso `19790644` foi estendida por experimentos causais aos outros quatro casos.
2. **Defeito do LECTK ou dos clientes:** **não identificado no escopo testado**. Não afirmar ausência universal de defeitos, nem afirmar que a origem da diferença é exclusivamente uma operação histórica de conversão externa.
3. **Correção:** **dispensada por ausência de defeito demonstrado**. Não modificar código, equações, piso da estabilidade estática, referências ou tolerâncias.
4. **Intercambialidade:** arquivos CDS e NCAR **não são automaticamente intercambiáveis** para análises energéticas que exigem a tolerância original de `10⁻⁴` da escala por coluna, ainda que os campos meteorológicos passem nas tolerâncias originais. Evitar misturar fontes sem avaliar a sensibilidade pretendida; a conclusão abrange os casos históricos examinados e não estabelece uma propriedade universal dos serviços.
5. **Sem novas aquisições ou investigação adicional obrigatória:** o responsável aceitou as ressalvas acima e decidiu encerrar este diagnóstico sem novos testes. A conclusão não reclassifica os resultados de desempenho ou de cobertura da D-Bench original; a Frente E não abre automaticamente.

## Preservação, cronologia e referências

- Evidências completas, arquivos híbridos, scripts diagnósticos, resultados e hashes foram relatados pelo agente como preservados **fora do Git**, na área privada de diagnóstico pós-encerramento, com `MANIFEST.sha256` (112 entradas) e `RELATORIO.md`. A acessibilidade atual desses bytes é responsabilidade do ambiente de investigação; este registro não substitui a retenção externa.
- Permanecem preservados `dbench_private/`, `LEC_Results/` e os registros históricos de [pré-registro](2026-10-29-dbench-preregistration.md), [Etapa B](../../../../source/remote_data_acquisition.rst) e [encerramento de 01/11](2026-11-01-dbench-closure.md).
- O presente documento é um **adendo posterior**, aprovado cientificamente pelo responsável, não uma edição retroativa dos resultados da campanha. Não altera a fila principal, `main`, tags, releases ou publicações.
