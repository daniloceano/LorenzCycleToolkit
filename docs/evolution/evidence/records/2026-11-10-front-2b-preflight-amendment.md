# Frente 2b — pré-voo da fase 1, emenda 1 ao pré-registro e início da campanha (10/11/2026)

**Natureza:** registro documental. Torna permanentes o resultado do pré-voo da fase 1, a emenda 1 ao
[pré-registro v1](2026-11-09-front-2b-stage0-preregistration.md) e as decisões do responsável de 10/11/2026
que autorizaram a campanha. **Nada executável muda no repositório:** nenhum código, teste, script,
dependência ou CI. O harness, os manifestos e os dados ficam fora do Git, em `<diretório de trabalho no
servidor>`. `main`, tags e releases continuam inalterados.

**Base:** `develop` em `cac34e96a84fbe693b89a7d300797b8b41e2b3e0` (PR #61), exportado com `git archive`
para o diretório de trabalho; o clone do servidor não foi tocado.

O pré-registro v1 não é reescrito: as mudanças abaixo são emendas e decisões datadas, tomadas **antes**
de qualquer medida da campanha.

## 1. Pré-voo (10/11/2026, 01:35–01:51 UTC)

- **Caixas seladas** antes do 1º download: `boxes.json` sha256
  `791a06e5d9ddaa85d55848e2366c5d67860df49d19807ef726a47bf28efbe89c`; trajetórias conferidas com os sha256
  da D-Bench; 2,32 GB úteis (float32) por repetição; menor caixa: caso 19800296.
- **As 7 rotas passaram** no teste de acesso sem dados e no pedido mínimo (1 horário, 500 hPa, caixa
  2°×3°, 5 variáveis), com conteúdo verificado. Nenhuma rota foi excluída.
- **WEkEO:** a ordem do `bbox` é [W, S, E, N], confirmada por um pedido mínimo que devolveu exatamente a
  grade pedida.
- **DestinE (opção A, §2.8 do pré-registro):** o caso 19800296 recebeu 1,28 GB para 0,03 GB úteis, em 269
  requisições HTTP e 161 s. Como ficou ≤ 10 GB/caso, vale o **ramo completo** (10 casos, teto próprio de
  500 GB, parada em 150.000 requisições no mês).
- **Tráfego do pré-voo:** 1,62 GB (limite: 60 GB). Nenhuma repetição da campanha foi feita no pré-voo.
- **Estimativa revista** (volumes a partir de metadados; tempos inferidos): ≈ 125 GB por repetição sem a
  DestinE, dos quais ≈ 71 GB do kerchunk; ≈ 376 GB nas 3 repetições, acima do teto de 300 GB do §2.7.
  DestinE ≈ 37 GB por repetição. Duração inferida: 2 a 4 dias, mais se a fila da CDS for longa.

### Achados do pré-voo

- **Defeito do LECTK — cliente Get a Subset (`src/utils/gdex_sources.py`):** a API do GDEX passou a
  responder `301 Moved Permanently` em `status/<id>` e `get_req_files/<id>`, redirecionando para o mesmo
  caminho com barra final; o cliente recusa redirecionamentos e falha. **Registrado, não corrigido nesta
  frente.** O harness acrescenta a barra final sem mexer em `src/` (decisão (b), §3).
- **Conta CDS compartilhada:** havia pedidos de outro trabalho na fila e pedidos rejeitados que não são
  deste benchmark; nenhum foi tocado. Condição de medida: o limite por usuário da CDS pode alongar a fila.
- **Observação, não comparação de produtos:** no pedido mínimo, CDS, WEkEO, ARCO e as três rotas NCAR
  deram valores idênticos bit a bit (CDS e WEkEO entregaram o mesmo arquivo); a DestinE difere em até
  0,062 K em T e 15,6 m² s⁻² em z, compatível com armazenamento quantizado. A comparação de produtos
  continua reservada à fase 2.

## 2. Emenda 1 ao pré-registro (10/11/2026)

| Item (§2.7) | Pré-registro v1 | Emenda 1 |
|---|---|---|
| Tráfego total da campanha, exceto DestinE | > 300 GB → a campanha para | **> 420 GB → a campanha para** |

**Motivo:** a estimativa revista (≈ 376 GB) excederia o teto, e a parada global truncaria a repetição 3.
**Decidida pelo responsável antes do início da campanha**, sem nenhum dado de medida. O teto próprio da
DestinE (500 GB) e as demais regras não mudam.

## 3. Decisões do responsável (10/11/2026)

- **(b) Get a Subset — aceito:** o adaptador de barra final fica no harness. Se o pedido da repetição
  anterior ainda estiver listado no serviço e o cliente parar por deduplicação, isso conta como falha da
  fonte, pela regra pré-registrada.
- **(c) Retomada automática — autorizada:** um vigia agendado no servidor (na reinicialização e a cada
  15 min) reinicia o controlador da campanha quando ele não está rodando, salvo se a campanha terminou ou
  parou por regra. Implementa o §2.9.
- **(d) Clientes e desenho dos pedidos — aceitos:**
  - Clientes kerchunk e Get a Subset de `src/utils/gdex_sources.py` do export de `develop`, e não as
    cópias congeladas da D-Bench (que dependem de caminhos e grades da D-Bench e ficam guardadas para
    proveniência). **Desvio da instrução original do Portão 2, aceito.**
  - CDS: 1 pedido por caso (intervalo de datas), e não o fatiamento diário do LECTK.
  - WEkEO: 1 pedido por mês do caso.
  - NCAR: NCSS e kerchunk entregam 37 níveis, e Get a Subset, horas cheias; o excesso entra na métrica.
  - Zarr (ARCO, DestinE): 4 conexões HTTP em paralelo dentro de cada pedido; NCSS com 1 conexão.
  - "Um pedido por vez no servidor inteiro" vale no nível do pedido.

## 4. Início da campanha

Com aprovação expressa do responsável, a campanha começou em **10/11/2026, 10:30 UTC**, com:

- harness v5 sha256 `cc1538a090f3e117b05b2777174a15336abdecd340ab38cc956c5436e214db42` (fora do Git);
- código do LECTK: export de `develop` em `cac34e9`;
- caixas `791a06e5…` (§1);
- semente da ordem sorteada: 20261109.

Ordem das rotas sorteada:

| Repetição | Ordem |
|---|---|
| 1 | destine → wekeo → gas → ncss → kerchunk → cds → arco |
| 2 | gas → ncss → kerchunk → wekeo → arco → cds → destine |
| 3 | gas → cds → ncss → arco → kerchunk → wekeo → destine |

## Fora deste registro

Os resultados da campanha, o relatório da fase 1 e a correção do defeito do cliente Get a Subset ficam
para registros posteriores. Merge, CI e releitura deste documento não podem constar dele.
