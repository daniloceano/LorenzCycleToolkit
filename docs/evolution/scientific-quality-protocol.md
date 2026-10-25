# Protocolo de rigor científico e verificação computacional

**Estado: APROVADO em 25/10/2026; aplicação prospectiva.** Não é aplicado retroativamente a frentes homologadas. Subordinado às decisões científicas explícitas, ao escopo autorizado, ao protocolo de coordenação e à política de desenvolvimento vigente.

## 1. Planejamento e portões verificáveis

1. Formular problema, hipótese, previsão e variável observável **antes da medição**, quando aplicável. Indicar população/cenário, instrumento, configuração, commit, versão dos dados, tolerâncias e tratamento de resultados ausentes.
2. Definir critérios de PASS/FAIL ou de aceitação **alcançáveis** para a população-alvo. Verificar por raciocínio e, quando necessário, aritmética, se o teste pode de fato reprovar. Não criar critérios que sempre passem, que exijam erros deliberados do produto ou que dependam exclusivamente de resposta de terceiros.
3. Separar **investigação/Etapa 0** e **implementação** quando a hipótese de causa não estiver suficientemente demonstrada. Resultado negativo/inconclusivo deve ser mantido; não reescrever previsão após ver os números.
4. Registrar claramente se uma checagem é **prospectiva** ou **retrospectiva**. Uma referência escolhida depois de observar o resultado não equivale a um teste independente pré-registrado.

## 2. Baseline, regressão e testes informativos

1. Sempre que comparar desempenhos, colocar ao lado uma **linha de base adequada** (comportamento anterior, método trivial ou referência previamente aprovada) e identificar as mesmas condições de medição. Não pressupor que um número histórico copiado represente o código atual.
2. Em mudanças que prometem comportamento padrão inalterado, comparar antes/depois **no mesmo ambiente, configuração, amostra, gerador e tolerâncias**, quando tecnicamente viável. Distinguir igualdade bit a bit de equivalência dentro de tolerâncias.
3. Conferir conteúdo, valores, unidades, dimensões, fronteiras temporais/espaciais, alinhamento e integridade — não apenas a existência de arquivos e ausência de erro.
4. Para testes críticos, quando seguro e proporcional, incluir **controle de sensibilidade/falha** em ambiente descartável: demonstrar que uma violação real da regra testada é detectada. O agente de código implementa toda mutação controlada; não tocar em dados/branches canônicos nem alterar teste para produzir um FAIL artificial.
5. Um CI verde é evidência de passagem dos jobs executados, não prova automática da validade da ciência. Verificar a receita real da CI (pacote instalado versus import local, ambiente limpo, versões importantes, jobs e skips), sempre que isso afetar a alegação.

## 3. Independência e dados reservados

1. **Regressão numérica**, **reprodução computacional**, **aprovação da formulação** e **validação física independente** são afirmações diferentes e devem ser rotuladas separadamente.
2. Quando houver alegação de validação independente, não usar apenas o mesmo script, resultado ou hipótese do agente como verificador. Empregar implementação, fonte, derivação, dados ou protocolo suficientemente independente, conforme escopo e risco; documentar dependências compartilhadas e limites restantes.
3. Havendo dados reservados de teste/validação, identificá-los e mantê-los isolados. A exposição antecipada de rótulos, resultados ou cortes invalida o uso posterior como prova independente **sem ressalva explícita**. Não presumir a existência de conjuntos reservados onde não foram definidos.

## 4. Proveniência e registro de erros

1. Para quantidades novas e relevantes, guardar método de cálculo ou script de geração, entradas, unidades, versões, SHA do código, amostra e ambiente. Regenerar números quando a mudança atual exigir; não aceitar transcrição de memória como verificação suficiente.
2. Se dois instrumentos discordarem, investigar também o instrumento, a conversão de unidades, as definições e a população antes de concluir que mudou o sistema.
3. Documentar hipóteses refutadas, falhas de testes, correções e erros da conversa-filha, do agente ou de terceiros de maneira atribuída à ação e verificável, sem ocultar achados negativos ou dados de baixa confiança.
4. Nenhuma alteração técnica ou científica pode ser apresentada como implementada ou homologada sem evidência adequada à alegação. O desenho é aprovado separadamente quando a governança exigir.

## 5. Responsabilidades e proporcionalidade

- **Filha:** especifica critérios, consulta código, executa testes existentes de maneira segura, interpreta resultados e documentação; **não edita código, testes-fonte, scripts, notebooks executáveis ou CI**.
- **Agente de código:** implementa mudanças executáveis após autorização, executa experimentos e controles previstos, registra evidências e verifica escopo e CI.
- **Mãe:** confere a coerência do fechamento e a evidência documental, sem declarar validação física independente apenas por ter relido o relatório do agente.
- Aplicar critérios de maneira **proporcional ao risco e ao objeto**. Revisão meramente documental não precisa de experimento científico artificial. Casos sem medição quantitativa devem usar critérios qualitativos observáveis.

**Este protocolo não abre novas frentes, não altera equações ou convenções e não autoriza publicação, exclusão, alterações em `main` ou bypass da governança Git.**
