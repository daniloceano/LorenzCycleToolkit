# Frente A — auditoria, critérios e cronologia

**Objeto auditado:** LECTK 2.0.0, commit `d38cda7e37d8e8a3a937a5919640a94bef19e34a`. Auditoria original concluída **19/10/2026 sem PASS global limpo**. A matriz usada para **fechamento retrospectivo** em 23/10 foi a §20.B do relatório consolidado; não foi comprovada igualdade literal com lista externa original.

| # | Critério | Revisão retrospectiva 23/10 | Observação |
|---|---|---|---|
| 1 | Inventário | PASS | Código, wheel, dependências e Natural Earth; pendências de packaging separadas |
| 2 | Suíte integral | PASS | 66 testes originais; CI posterior 302 develop, 291 main em SHAs diferentes |
| 3 | Referência fixa | PASS | 24 arquivos `ERA5_fixed` |
| 4 | Referência trajetória | **PASS com ressalva** | 25 arquivos `ERA5_track`; trajetórias estacionária/sintética, sem ciclone real |
| 5 | Reprodução congelada | PASS | Bit a bit no cenário 2D; entre plataformas tolerâncias |
| 6 | Unidades e dimensões | PASS | **Parcial em 19/10**, completado pela Onda 2 (23 termos) |
| 7 | Domínio geográfico | PASS | **Parcial para track em 19/10**, Onda 3 passou a registrar área efetiva |
| 8 | Rastreabilidade E/S | **PASS com ressalva** | Versão de resultados ainda pode constar `2.0.0`, igual ao PyPI antigo |
| 9 | Falhas classificadas | PASS com correção temporal | Onda 4 **posteriormente confirmou** o defeito E0-21b, corrigido |
| 10 | Preservação do ambiente original | **PASS com ressalva** | R1: `mtime` interno `.git` por autofetch, sem conteúdo alterado segundo evidências |

**Resultado retrospectivo: 10/10 PASS, três com ressalvas**, sem reexecutar ciência na Etapa 0 documental.

**Defeito descoberto depois:** Onda 4 identificou volume vertical variável entre instantes, com tendências artificiais no track. Correção de conjunto vertical consistente homologada lote 8 e integrada PR #25. **Onda 6:** quatro FAIL pré-registrados (V1/V2/V3/V5) devem permanecer. Referência 2D: 154 arquivos, 148 científicos; 49 referências numéricas versionadas relatadas preservadas.

**Limites:** cálculos e CI não implicam validação física independente; dados ERA5 37 níveis/3 h, NCEP-R2 e ciclone real rastreado permanecem não demonstrados. [Encerramento](front-A-closure.md), [proveniência](../evidence/primary-artifacts.md).
