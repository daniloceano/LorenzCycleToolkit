# Evidência primária e conferência

| Conjunto | Resultado |
|---|---|
| Sources do Project | 19 inventariadas; duplicação S10 = S11 identificada |
| Anexos formais lote 15 | Quatro recebidos e conferidos |
| Primeiro ZIP documental | 20 das 33 evidências indexadas; 94 registros de coleta |
| ZIP complementar | Dez documentos de revisão consolidada/lotes 5, 8 e 10–14 |
| Sources já existentes | Três evidências indexadas |
| **Total** | **33/33 MATCH, 0 MISMATCH, 0 ausentes** |

### Identificadores de integridade

- ZIP inicial `LECTK_RECUPERACAO_DOCUMENTAL_2026-10-24.zip`: SHA-256 `2f179e1f52270fb55c3d12aa9d03c71bff74124543ebcd9026b54ed53b8aa9ee`.
- ZIP complementar `LECTK_COMPLEMENTO_DOCUMENTAL_2026-10-24.zip`: SHA-256 `3bd2dc9c83ff5687015971a57073cc2cc11cf815937d7d02a72f7b7b8d5be6a2`.
- Índice original `EVIDENCIAS_SHA256.txt`: SHA-256 `facba873cc9bdecbd2caf15dc2f1c363558c6b3556ef32b2819b0ef90a0fc6a9`.
- Registro técnico lote 15 `RELATORIO_FECHAMENTO_FRENTE_A.md`: SHA-256 `bf6e33fb7766d2f390ffbe6814f79c25c94b5ab569e5fbd0f01657a2d20a74ba`.

[Índice integral de 33 SHA-256](sha256-index.tsv) com caminhos relativos saneados. Os hashes identificam bytes privados, não implicam publicação dos artefatos correspondentes no Git.

### Estado Git atestado no corte

PR #25 `main` `48fe8e2e75d0c1527214cedda145924ab13ae9ea`; PR #26 `develop` `911122dcbee49a4779413acdae64cc9cd513459d`; PR #27 `develop` `65820a341ed974d7212e32e9d60030defc79d488`; PR #28 `develop` `81cc00091e546a136ab2ceb0e2c372a444076d97`. Esse estado é datado, não o SHA de futuras versões.

**Limite:** a reconciliação não reexecutou testes atmosféricos, não recuperou todos os bytes dos diretórios de dezenas de GB e não realizou validação física independente. Documentos brutos, datasets, logs e arquivos não rastreados requerem retenção/backup separado. Migração da memória não autoriza apagar diretórios dos computadores.
