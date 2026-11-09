# Subfrente 2a.1 — conferência dos termos de fronteira e do fechamento do resíduo (09/11/2026)

**Natureza:** testes sintéticos pré-registrados; nenhuma mudança de comportamento em `src/`.
**Base:** `develop` em `4840ae55c61ff62622728f2f56dccacb7bb05d27`. Ligada à [Frente 2a](2026-11-08-frente-2a-dias-pinto-2011.md) (D-021).

## 1. Objetivo

Verificar, com campos analíticos e uma referência contínua independente do código (quadratura de
Gauss-Legendre), se BAZ, BAE, BKZ, BKE, BΦZ e BΦE de `src/analysis/boundary_terms.py` reproduzem as
expressões de área limitada de Brennan & Vincent (1980, Apêndice) e Muench (1965, Eqs. 6–11), e se os
resíduos fecham como RGZ = GZ, RGE = GE, RKZ = BΦZ + DZ, RKE = BΦE + DE. O código de produção é usado
como está (`BoxData` → `BoundaryTerms`; `calc_budget_diff` → `calc_residuals`).

## 2. Conclusão

BAZ, BAE, BKZ, BKE, BΦE e os quatro resíduos conferem com Brennan & Vincent (1980) e Muench (1965).
BΦZ segue a forma com Φ* em todas as faces ([D-022](../../decisions.md)), invariante a uma constante
aditiva em Φ; a forma de B&V com Φ completo, montada só no teste, muda pela quantidade prevista
analiticamente. Nenhum erro de fórmula foi encontrado.

## 3. Resultados (critérios pré-registrados)

| Teste | Critério | Resultado |
|---|---|---|
| T1 | Faixa global de 360° com coluna 0° = 360° repetida: parte L/O < 1e-10 × magnitude de uma parede | PASS |
| T2 | u, v, ω constantes; T′ e Φ′ polinomiais: erro relativo ≤ 1 % em 0,25° (termos nulos: ≤ 1e-10 × escala) | PASS |
| T3c | Cada uma das 6 paredes isolada × 6 termos (36 combinações); refinamento conjunto 1°/2000 Pa → 0,5°/1000 Pa → 0,25°/500 Pa; ordens observadas em [1,8; 2,2] e erro do valor extrapolado por Richardson (p = 2) ≤ 0,1 % | PASS 36/36; ordens 1,968–2,008; erro de Richardson máx. 8,10e-5 |
| T4 | BΦZ(Φ + C) × BΦZ(Φ): variação relativa < 1e-10 | PASS |
| T5 | RGZ, RGE, RKZ, RKE reproduzem GZ, GE, BΦZ + DZ, BΦE + DE com erro ≤ 1e-10 | PASS |
| Controles | Mutantes em memória: M1 (sinal da parede oeste), M2 (cos φ → \|sen φ\|) devem reprovar T2 e T3c | Ambos detectados |

## 4. Resultados negativos registrados

- **T3 — FAIL (4/36).** Critério original: erro relativo ≤ 1 % por termo em 0,25°. Topo e base de BAE e
  BKE ficaram em 1,11 %, com convergência de 2ª ordem (4,4 % → 1,11 % → 0,28 %): tolerância apertada
  demais para a grade, não erro de fórmula.
- **T3b — FAIL (4/36).** Mesmo critério de ordem e Richardson do T3c, mas refinando só a horizontal
  (Δp fixo em 1000 Pa). Falharam L/O de BAE e N/S de BAZ pela ordem observada. Causa: patamar de
  discretização vertical de 2ª ordem (≈ 1,1e-4, proporcional a Δp²) nos termos ponderados por 1/σ, que
  o T3b não refinava. O T3c, com refinamento conjunto, foi a tentativa final pré-registrada.

## 5. Implicação

A FAIL da Frente 2a ([D-021](../../decisions.md)) não decorre de erro na implementação dos termos de
fronteira nem do fechamento dos resíduos. A questão de CA e CE menores no Reg3 segue no backlog da
[fila](../../queue.md).

## 6. Alterações integradas

- `tests/test_boundary_terms_synthetic.py`: T1, T2, T3c, T4, T5 e controles M1/M2.
- `src/utils/thermodynamics.py`: só docstring de `StaticStability` (SIGMA_FLOOR também afeta BAZ e BAE).
- `docs/evolution/decisions.md` (D-022) e `docs/evolution/queue.md` (estado da 2a.1).

## 7. Proveniência

- Base: `develop` `4840ae55c61ff62622728f2f56dccacb7bb05d27`.
- sha256 de `tests/test_boundary_terms_synthetic.py`: `90056185735ff4445b645e6589f1020edad0206e216618cf6e7dbf3875b3956d`.
- Ambiente: Python 3.12.15 (arm64), numpy 2.0.0, xarray 2024.2.0, pandas 2.2.2, metpy 1.6.2, pint 0.24.3,
  scipy 1.14.0, pytest 8.2.2.
- Suíte completa: única falha em `test_packaging::test_installed_metadata_matches_pyproject`, preexistente
  e também presente antes das mudanças.
- Evidências completas (JSON por combinação, logs, snapshots e diffs de T3, T3b e T3c, diagnóstico do
  patamar vertical) mantidas fora do Git, na área de trabalho da frente, sob
  [GOV-RET-01](../../retention-policy.md).
