# Replanejamento do projeto de modernização — decisões do responsável (07/11/2026)

**Natureza:** registro documental. As decisões abaixo foram tomadas **pelo responsável** em conversa de coordenação de 06–07/11/2026 e existiam apenas no chat; este registro as torna permanentes. **Nada executável muda:** nenhum código, teste, script, dependência ou CI. Nenhuma frente é aberta por este registro, e `main`, tags e releases continuam inalterados.

**Base:** `develop` em `96b1da8e100e8df0bfcfc46c369b2343270c7c94` (PR #56); `main` em `48fe8e2e75d0c1527214cedda145924ab13ae9ea`.

## 1. Objetivo e entregável (D-012)

O projeto de modernização tem como entregável **incluir no LECTK novas fontes de download de ERA5 e documentar cada uma**: limitações, tempos de espera estimados e melhores candidatas para usos distintos.

## 2. Versão única ao fim do projeto (D-013)

Todo o trabalho em `develop` e nas branches da modernização sai numa **única nova versão oficial, ao fim do projeto**. Esse marco passa a ser o gatilho de reavaliação das retenções temporárias de GOV-RET-01 hoje vinculadas à “próxima versão/tag” (Marco B). **A tag continua sendo gatilho de avaliação, não autorização de exclusão**, e a criação da versão exige autorização específica.

## 3. Validação científica externa (D-014)

Reproduzir os três ciclones (Reg1, Reg2, Reg3) de Dias Pinto & da Rocha (2011), *J. Geophys. Res.*, 116, D14112, doi:10.1029/2011JD016217, com a reanálise NCEP/NCAR R1, as caixas fixas e as fases definidas no artigo.

- As figuras 6, 10 e 14 do artigo trazem valores médios por fase; as séries das figuras 5, 9 e 13 exigiriam digitalização.
- **As implementações são independentes:** o algoritmo do artigo vem do grupo de estudos climáticos da USP; o LECTK vem do doutorado do responsável, em outro grupo. Por isso a comparação conta como verificação externa, e não como regressão.
- O repositório já contém `inputs/namelist_NCEP-R1`, `inputs/box_limits_Reg1` e `inputs/track_Reg1-Representative` (presença conferida em `96b1da8`). **A correspondência desses arquivos com o artigo ainda não foi verificada.**
- Critérios de comparação devem ser fixados antes de rodar, conforme o [protocolo de qualidade científica](../../scientific-quality-protocol.md).

## 4. Fora do escopo (D-015)

- Determinar qual reanálise ou fonte é “mais correta”: apenas documentar que os produtos diferem; a escolha é do usuário.
- Investigação aprofundada de níveis verticais: a hipótese K continua fora.
- O efeito do nível de 1000 hPa abaixo do relevo.

## 5. Níveis verticais: medição simples (D-016)

Medição simples e documentada com os 10 casos da D-Bench em **modo de caixa fixa**, usando a fonte mais eficiente apontada pelo benchmark de fontes, relatando as diferenças médias como **amostra reduzida**. A regra de definição da caixa de cada caso deve ser fixada antes de rodar. Não equivale à hipótese K.

## 6. Premissas do benchmark unificado de fontes (D-017)

Estado: **Etapa 0, nada executado.**

- Cada fonte é medida como ela é, incluindo o volume que obriga a baixar.
- Medir a degradação de desempenho em pedidos seguidos (pedidos em sequência, com horário registrado).
- Condições comuns: a menor resolução temporal entre as fontes, os níveis da fonte mais restrita e topo em 50 hPa.
- **Somente fontes gratuitas.**
- **Somente credenciais que um usuário não desenvolvedor configure seguindo tutorial simples**; fonte que não cumprir sai do benchmark.
- Modo protegido CDS (`src/utils/cds_guard.py`) **obrigatório**, e regras de parada no código para todas as fontes.
- Execução prevista no Mac Mini do responsável, sem supervisão, por vários dias.
- O plano detalhado anterior pode estar só no servidor de processamento ou no Mac Mini; localizá-lo é tarefa de frente futura.

## 7. Nova ordem de trabalho (D-018)

Substitui a fila antiga a partir da Frente F; G–J seguem como backlog.

1. Este registro.
2. Em paralelo: (a) validação Dias Pinto & da Rocha; (b) benchmark de fontes.
3. Níveis verticais (depende de 2b).
4. Inclusão das fontes escolhidas e documentação de cada uma (depende de 2b; a Frente F é absorvida aqui).
5. Arrumação final.
6. Nova versão (exige autorização específica).

Fechar uma etapa não abre automaticamente a seguinte.

## 8. Pendências de arrumação (D-019) — registrar, não executar

| Pendência | Situação verificada em `96b1da8` |
|---|---|
| Teste `test_packaging::test_installed_metadata_matches_pyproject` falhando há tempos | Falha preexistente também relatada na releitura da LZV (03/11) |
| Mensagens de erro para domínio inválido | Pendente |
| Aviso no CHANGELOG da rejeição de tracks irregulares (commit `7e170ac`) | **Já presente** no `CHANGELOG.md` (entrada sobre alinhamento track × dados cita tempos irregulares e menos de duas posições); manter apenas a conferência na arrumação final |
| Harmonizar a regra que exige recomendar Codex e Claude Code com a orientação de neutralidade de fornecedor | Pendente; a regra de 25/10 continua vigente até decisão |
| Limpeza opcional de branches antigas | Operação destrutiva: só com autorização específica |

## 9. Correções factuais (D-020)

Correções por nota; os registros originais não foram reescritos.

- O limite de tempo/tentativas do cliente CDS, descrito em 29/10 como “mudança de código ainda a autorizar e implementar”, **foi implementado como modo protegido CDS no PR #47** (commit `7d72890`, 01/11/2026; `src/utils/cds_guard.py`).
- O responsável **revisou visualmente** a reorganização da documentação do usuário feita no PR #55 (declaração do responsável; revisão humana não verificável no Git).

## Fora deste registro

Merge, CI e releitura deste documento não podem constar dele; ficam no relatório de encerramento da Frente R e no handoff à conversa-mãe.
