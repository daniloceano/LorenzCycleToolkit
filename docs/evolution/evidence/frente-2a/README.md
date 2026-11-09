# Frente 2a — evidências versionadas (Dias Pinto & da Rocha, 2011)

Arquivos mínimos para reproduzir a comparação registrada em
[2026-11-08-frente-2a-dias-pinto-2011.md](../records/2026-11-08-frente-2a-dias-pinto-2011.md). Dados R1,
saídas brutas, logs, figuras e recortes do artigo ficam fora do Git, sob GOV-RET-01.

| Arquivo | Conteúdo |
|---|---|
| `REGRAS_OPERACIONAIS_P5.md` | Regras congeladas antes de abrir as saídas do LECTK (cópia idêntica, SHA-256 `248d57ad…`) |
| `p5_tabela_144.csv` | 144 linhas: caso, fase, termo, artigo, LECTK, derivada avançada, sensibilidade, E, faixa, classe (SHA-256 `a055ded1…`) |
| `periods_Reg{1,2,3}.csv` | Fases das legendas das Figs. 6, 10 e 14, no formato do `periods.csv` do LECTK |
| `box_limits_Reg{1,2,3}` | Caixas fixas das legendas da Fig. 2 |

## Como reproduzir

1. Obter NCEP/NCAR R1 (2,5°, 6 h; `air`, `hgt`, `uwnd`, `vwnd`, `omega`) para a caixa e o período de cada
   caso e montar um NetCDF por caso. O ômega do R1 só existe em 12 níveis (1000–100 hPa); a execução
   principal usa esses 12 níveis.
2. Usar `inputs/namelist_NCEP-R1` como `inputs/namelist`.
3. Rodar o LECTK em caixa fixa com resíduos:
   `python lorenzcycletoolkit.py <arquivo_RegN.nc> -r -f --box_limits <box_limits_RegN> --outname RegN`.
4. Sensibilidade de topo: repetir com 17 níveis (até 10 hPa), impondo ômega = 0 de 70 a 10 hPa.
5. Calcular as médias por fase com `periods_RegN.csv` (limites inclusivos, média simples dos horários de 6 h)
   e classificar cada termo conforme `REGRAS_OPERACIONAIS_P5.md`. O resultado deve reproduzir
   `p5_tabela_144.csv`.
6. Controles: caixa do Reg1 deslocada 15° para leste; caixa e datas do Reg3 com dados de 2007.

Os scripts de análise (`p5a`–`p5d`) não são versionados; seus hashes estão no registro.
