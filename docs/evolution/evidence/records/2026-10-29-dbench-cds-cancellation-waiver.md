# D-Bench — dispensa do cancelamento seletivo CDS (29/10/2026)

**Decisão:** opção 2, escolhida pelo responsável em 29/10/2026, para o portão 2 do [plano D-Bench](../../fronts/front-D-bench-cds-ncar.md#pré-registro-e-portões).

## Fatos declarados pelo responsável (não verificados pelo agente)

- A D-Bench usará uma credencial CDS que o responsável designou, armazenada fora do repositório. Nenhum conteúdo, identificador ou caminho dela é registrado aqui.
- A credencial é válida e de uso exclusivo da D-Bench.
- Não há meio de cancelar seletivamente as requisições CDS.

O agente verificou apenas que existe um arquivo local, que ele tem permissões restritas ao dono e que a entrada designada tem formato de token. Validade, licenças, cota e lista de jobs da conta não foram consultadas, porque nenhum contato com o CDS é autorizado nesta etapa.

## Mitigações que substituem o cancelamento seletivo (todas obrigatórias)

1. **Um pedido CDS por vez**, sem concorrência. A comparação CDS sob concorrência continua bloqueada.
2. **Ledger privado criado antes de cada submissão**, fora do Git, com caso, repetição, ordem e horário UTC; anexar identificador remoto do job **assim que recebido** (ele não pode ser conhecido antes da submissão). Os identificadores ficam somente nesse ledger.
3. **Limites no cliente:** um teto de tempo por caso e um número pequeno de tentativas, em vez do atual `timeout=600, retry_max=500`. É uma mudança de código que **ainda depende de aprovação de escopo e implementação**, com testes, antes do piloto.
4. **Regra de parada:** depois de um timeout ou de um erro de cota, autenticação ou política, nenhuma nova submissão naquela sessão. O caso é registrado como falha ou censurado, nunca reaproveitado nem imputado.
5. **Risco residual explicitamente aceito:** um job remoto poderá continuar após timeout ou encerramento local e consumir cota até terminar ou expirar; limitar tentativas locais **não garante** cancelamento remoto. Antes de cada nova submissão, confirmar por meio seguro, no portal ou na API, que não há job pendente da campanha; se isso não puder ser confirmado, **não submeter novo pedido**.

## Uso isolado da credencial

A credencial deverá ser entregue ao cliente **somente em configuração isolada e compatível com a versão instalada** (avaliar `CDSAPI_RC` apenas após confirmar suporte efetivo, ou parâmetros privados do cliente). O `~/.cdsapirc` e as credenciais globais não são alterados, e a credencial nunca entra em logs, `run_status.json`, commits ou relatórios.

## O que este registro não autoriza

Nenhuma submissão CDS, download ou início do piloto. Esta dispensa de uma salvaguarda **não demonstra** validade da credencial, isolamento de cota, ausência de jobs pendentes ou disponibilidade do serviço. Os outros portões do plano continuam pendentes: orçamento, integridade e tolerâncias, verificação efetiva da compatibilidade dos namelists (divergência **não comprovada**) e proteção contra reuso CDS sem proveniência. Não há autorização para abrir solicitações por tentativa e erro.
