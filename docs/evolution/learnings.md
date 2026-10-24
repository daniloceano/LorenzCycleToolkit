# Aprendizados, negativos e limites

1. Em **19/10** a auditoria original não tinha PASS global limpo. Unidades por nível e observação de domínio real no modo trajetória eram parciais; Ondas 2 e 3 completaram essas verificações posteriormente.
2. A frase da auditoria original “nenhum defeito científico confirmado” era datada: a **Onda 4** confirmou **E0-21b**, variação do volume de controle vertical entre instantes no modo track que gerava tendências artificiais. A correção foi homologada no lote 8 e integrada via PR #25. Não apagar a descoberta nem atribuí-la retroativamente à auditoria.
3. A **Onda 6** teve **V1, V2, V3 e V5 FAIL** nos critérios pré-registrados; **V4 e V6 PASS**. Uma verificação do instrumento foi corrigida *post hoc* com explicitação; quatro termos Kz, Ke, Cz, Ce foram relatados reproduzidos com desvio até 5,5×10⁻¹⁶ no caso delimitado. Isso **não valida fisicamente** os resultados.
4. A linha de base 2D continha **154 arquivos, 148 científicos**. Entre plataformas, reprodução por **tolerâncias**, não necessariamente igualdade bit a bit. **49 referências** em `tests/data/reference/` foram relatadas intactas no corte de 23/10.
5. A CI registrada do `develop` `81cc000` teve **302 testes** em Actions/CircleCI; `main` `48fe8e2` teve **291**. CI verde não substitui comparação física independente.
6. Lacunas: ciclone real rastreado, validação física externa, ERA5 produtivo de 37 níveis/3 h, NCEP-R2, Dz/De; falha de mensagem para domínio inválido (quatro de oito casos com erro interno), possíveis problemas de empacotamento e rótulo `lectk_version` ainda `2.0.0`.
7. Hipótese K de convergência com menos níveis é **não medida**: distinguir economia computacional, desvio da referência e validade científica.

**Método:** definir perguntas, entradas, SHA, ambiente, unidades, controles e tolerâncias antes de medir; guardar FAILs e mudanças pós-teste; testar o pacote instalado efetivamente contra o checkout; diferenciar números relatados por documentos de execução independente. Consulte [auditoria A](fronts/front-A-audit.md) e [artefatos](evidence/primary-artifacts.md).
