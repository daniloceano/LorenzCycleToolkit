# Frente 2a — Passo 5: regras operacionais (congeladas antes de abrir as saídas do LECTK)

Escrito em 08/11/2026, antes de qualquer leitura de valores dos CSV de resultado do LECTK. Este documento não
altera o pré-registro aprovado; ele só torna explícito como cada item foi operacionalizado onde o texto admite
mais de uma leitura. Onde houve escolha, ela é marcada com **[INTERPRETAÇÃO]** para revisão do responsável.

## R1. Mapeamento dos 16 termos (colunas do CSV `-r -f`)

| Termo do artigo | Coluna LECTK | Sinal |
|---|---|---|
| ∂AZ/∂t, ∂AE/∂t, ∂KZ/∂t, ∂KE/∂t | `∂Az/∂t (finite diff.)` … `∂Ke/∂t (finite diff.)` | igual |
| CZ, CA, CE, CK | `Cz`, `Ca`, `Ce`, `Ck` | igual |
| BAZ, BAE, BKZ, BKE | `BAz`, `BAe`, `BKz`, `BKe` | igual |
| RGZ, RGE, RKZ, RKE | `RGz`, `RGe`, `RKz`, `RKe` | igual |

Justificativa (código, não saídas): `src/utils/calc_budget_and_residual.py::calc_residuals` define
RGz = ∂Az/∂t + Cz + Ca − BAz; RKz = ∂Kz/∂t − Cz − Ck − BKz; RGe = ∂Ae/∂t − Ca + Ce − BAe;
RKe = ∂Ke/∂t − Ce + Ck − BKe. Isso equivale às eqs. (4)–(7) do artigo com CZ: AZ→KZ, CA: AZ→AE, CE: AE→KE,
CK: KE→KZ e B positivo entrando. RK do LECTK não subtrai BΦ, logo RK = BΦ − D + ε, como nas eqs. (10)–(11)
do artigo. Nenhuma troca de sinal é aplicada.

## R2. Médias por fase

Média simples das linhas de 6 h com início ≤ t ≤ fim (UTC), fases de `inputs/periods_RegN.csv`
(= legendas das Figs. 6, 10, 14). Contagens esperadas: Reg1 8/8/9; Reg2 6/4/4; Reg3 6/5/4.

## R3. Execução de referência e envelope de topo

L = execução principal (`-r -f`, 12 níveis, 1000–100 hPa). S = sensibilidade (omega = 0 de 70 a 10 hPa,
17 níveis). E = |L − S| por termo e fase (pré-registro, item 4).

## R4. Derivada alternativa (item 5, “troca da derivada centrada por diferença avançada”)

O LECTK usa `np.gradient` (centrada no interior). Alternativa: diferença avançada
X'(t) = [X(t+Δt) − X(t)]/Δt sobre as séries de energia da própria execução (colunas `Az, Ae, Kz, Ke`);
os resíduos são recalculados com as fórmulas de R1 trocando só a tendência. Conversões e termos B não mudam.
Todas as fases têm o horário seguinte disponível na série.

## R5. Classificação de cada termo (A = artigo, L = LECTK, T = max(0,25·|A|; 0,5))

Ordem de aplicação (leitura literal do item 5) **[INTERPRETAÇÃO]**:
1. |A| < 0,5 → **informativo** (só o sinal; não entra na contagem).
2. sinal(L) = sinal(A) e |L − A| ≤ T + E → **concordância**.
3. sinal oposto e |A| ≥ 1 → **falha** (sem possibilidade de explicação).
4. **discordância explicada** se (a) coberta por E: existe x em [L − E, L + E] com sinal(x) = sinal(A) e
   |x − A| ≤ T; ou (b) coberta pela derivada: L_fwd com sinal de A e |L_fwd − A| ≤ T + E.
5. caso contrário → **falha**.

Como E já está somado na faixa de concordância, (a) só tem efeito para sinal oposto com 0,5 ≤ |A| < 1.

## R6. Percentual

Denominador: termos com |A| ≥ 1. Numerador: concordância + discordância explicada entre eles.
Calculado por caso (48 termos) e global (144 termos).

## R7. Mecanismos dominantes (avaliados nos valores L; também conferidos no artigo) **[INTERPRETAÇÃO]**

Conversões consideradas: CZ, CA, CE, CK.
- Reg1 formação, “CK < 0 e relevante”: **igual** se CK < 0 e |CK| ≥ 1 W m⁻² (mesmo limiar de relevância do
  item 6); **invertido** se CK > 0 e |CK| ≥ 1; senão, não igual (sem inversão).
- Reg2 maturidade, “CE dominante”: **igual** se CE > 0 e |CE| é o maior dos quatro |C|; **invertido** se CE < 0
  com |CE| ≥ 1; senão, não igual.
- Reg3, “CA e CE dominantes” (sem fase indicada → aplicado às três fases): **igual** se, em cada fase, CA > 0 e
  CE > 0 são os dois maiores |C|; **invertido** se em alguma fase CA ou CE é negativo com módulo ≥ 1; senão,
  não igual.

## R8. Vereditos

- Por caso: **PASS** se percentual ≥ 80 % e mecanismo igual; **FAIL** se percentual < 60 % ou mecanismo
  invertido; senão **INCONCLUSIVO**.
- Global numérico: o mesmo, com o percentual dos 144 termos e os três mecanismos (todos iguais para PASS;
  qualquer invertido para FAIL).
- Controles (item 7): se algum controle não for majoritariamente falha, o veredito global passa a
  **INCONCLUSIVO**; o valor anterior à regra é registrado ao lado **[INTERPRETAÇÃO: a regra vale também sobre
  FAIL, por leitura literal]**.
- Item 8: o veredito é relatado como calculado; a aprovação final depende da inspeção visual do responsável.

## R9. Controles de falha

- `ctrl_Reg1_shift15E` comparado com os valores do artigo para o Reg1; `ctrl_Reg3_2007` com os do Reg3, nas
  mesmas datas de calendário (ano 2007).
- Não há sensibilidade de topo para os controles: usa-se o E do caso correspondente **[INTERPRETAÇÃO:
  escolha conservadora, que torna mais difícil o controle falhar]**. Derivada alternativa: da própria série.
- “Majoritariamente falha”: mais de 50 % dos termos com |A| ≥ 1 classificados como falha.

## R10. Transcrição

Duas transcrições independentes (A: sessão principal; B: subagente isolado, sem acesso às saídas do LECTK, à
transcrição A nem às digitalizações antigas da auditoria de 2026-08). Fechamento por caixa:
AZ: −CZ − CA + BAZ + RGZ; KZ: CZ + CK + BKZ + RKZ; AE: CA − CE + BAE + RGE; KE: CE − CK + BKE + RKE;
|tendência − soma| ≤ 0,1 W m⁻². Discrepâncias A×B resolvidas por nova leitura ampliada antes de abrir as saídas.
