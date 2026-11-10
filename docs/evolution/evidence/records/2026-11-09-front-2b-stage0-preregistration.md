# Frente 2b — benchmark unificado de fontes ERA5: Etapa 0 e pré-registro v1 (09/11/2026)

**Natureza:** registro documental. Torna permanentes as decisões do responsável de 09/11/2026 sobre a
Frente 2b e o pré-registro v1 da fase 1, aprovado na mesma data. **Nada executável muda:** nenhum código,
teste, script, dependência ou CI. Nenhum download foi feito para este registro, e `main`, tags e releases
continuam inalterados.

**Base:** `develop` em `498ecec73b176ed56a9167ed8d5d762b0e15609e` (PR #60); `main` em
`48fe8e2e75d0c1527214cedda145924ab13ae9ea`. Premissas de origem: [D-017](../../decisions.md)
([registro de 07/11](2026-11-07-replanning-decisions.md), §6); ordem de trabalho: D-018.

## 1. Decisões do responsável (09/11/2026)

- **(a) Emenda a D-017 — local de execução.** A campanha é executada no servidor de processamento, e não
  no Mac Mini do responsável, que está indisponível. As demais premissas de D-017 permanecem.
- **(b) Busca do plano antigo encerrada.** O plano detalhado anterior não foi encontrado no repositório
  nem no servidor de processamento. A busca teve cobertura parcial: não incluiu PDFs binários nem os
  diretórios `miniconda3/` e `MPAS-Research*`; o Mac Mini não foi verificado. O plano foi refeito a
  partir de D-017 (pré-registro v1, §2).
- **(c) Dados no disco local do servidor**, em `<diretório de trabalho no servidor>`.
- **(d) Retomada automática após queda** do servidor ou da conexão (regras em §2.9).

Com isso, a Etapa 0 da 2b está concluída; a fase 1 entra em pré-voo.

## 2. Pré-registro v1 (aprovado pelo responsável em 09/11/2026)

### 2.1 Pergunta

Para cada rota gratuita de acesso ao ERA5, sob pedido comum: tempo de fila, tempo de download, volume
baixado × útil, taxa de falhas e degradação em pedidos seguidos. Não se elege a melhor fonte nem a
"mais correta" ([D-015](../../decisions.md)).

### 2.2 Rotas

| Rota | Situação |
|---|---|
| CDS API | Incluída (sempre em modo protegido, `src/utils/cds_guard.py`) |
| CDS ARCO Zarr | Incluída |
| NCAR Get a Subset | Incluída |
| NCAR kerchunk | Incluída |
| NCAR NCSS | Incluída |
| WEkEO HDA (dataset `EO:ECMWF:DAT:REANALYSIS_ERA5_PRESSURE_LEVELS`) | Incluída |
| DestinE Earth Data Hub | Incluída, com regra própria (§2.8) |
| CDS timeseries | **Excluída:** séries pontuais, sem grade 3D |
| Google ARCO | **Excluída:** sem geopotencial em níveis de pressão |
| Earth Engine | **Excluída:** sem geopotencial em níveis de pressão |

**Resultados do teste de conta:**

- **WEkEO** — a credencial é o usuário/senha da conta, em `~/.hdarc`; o caminho da API não é evidente no
  portal (achado de usabilidade a documentar).
- **DestinE** — chave obtida na hora; cota de 500.000 requisições/mês, com renovação no dia 1º.

### 2.3 Pedido comum

- 10 casos da D-Bench.
- Caixa fixa = retângulo que envolve a trajetória + 5° de margem; as caixas são calculadas e seladas com
  sha256 antes do 1º download.
- Horários 00/06/12/18 UTC.
- 13 níveis, de 1000 a 50 hPa (grade ARCO).
- Resolução horizontal de 0,25°.
- Variáveis T, z, u, v, ω.
- Cada rota baixa como permite; o excesso entra na métrica.

### 2.4 Métricas por pedido

- Tempo de fila (quando a rota informa).
- Tempo de download e tempo total (relógio externo).
- Bytes baixados × bytes úteis (grade comum em float32).
- Resultado: sucesso / falha da fonte / falha de rede / interrupção externa.
- Verificação de conteúdo: variáveis, níveis, horários e caixa completos.

### 2.5 Degradação e concorrência

- Por repetição, os 10 casos de cada rota em sequência, em ordem sorteada com semente registrada.
- Tendência do tempo total ao longo da sequência.
- Repetições separadas por ≥ 12 h; 3 repetições.
- **Concorrência:** um pedido por vez no servidor inteiro.

### 2.6 Registro

Manifesto JSON por pedido (UTC de cada evento, parâmetros exatos, versão/commit do harness e dos
pacotes, tamanho e sha256 de cada arquivo, resultado, erros) e índice de hashes por repetição; **fora do
Git**.

### 2.7 Regras de parada

- 2 falhas consecutivas da fonte numa rota → a rota para naquela repetição.
- Teto por pedido de 12 h (CDS) e 2 h (demais) → falha.
- Tráfego total > 300 GB (exceto DestinE) → a campanha para.
- Disco livre < 1 TB ou duração > 7 dias → a campanha para.
- CDS sempre em modo protegido (`cds_guard`).

### 2.8 DestinE (opção A)

- O pré-voo baixa 1 caso e mede requisições e bytes reais.
- Se ≤ 10 GB/caso → campanha completa, com teto próprio de 500 GB.
- Se > 10 GB/caso → 3 casos × 3 repetições, marcado como **amostra parcial**.
- Em qualquer caso, parada ao atingir 150.000 requisições no mês.
- Se a campanha cruzar a virada do mês, registrar contadores separados.

### 2.9 Quedas

- Retomada automática.
- Pedidos com hash verificado não são repetidos.
- Pedido interrompido = "interrupção externa": não conta para a regra de falhas e é refeito.
- A repetição afetada fica marcada, e a degradação é reportada com e sem ela.
- Sem retentativas além das do próprio cliente.

### 2.10 Validade

Uma rota é documentável se ≥ 2 de 3 repetições terminarem sem parada global. Resultados negativos e
exclusões são conclusões registradas.

### 2.11 Para depois

- Fase 2: cálculos LEC e comparação de produtos, só após a 2a.1.
- Níveis verticais: etapa 3 (D-016).
- Inclusão no LECTK: etapa 4.

## 3. Sequência da fase 1

1. **Portão 1 (este registro):** documentação apenas.
2. **Portão 2:** harness fora do Git e pré-voo no servidor (caixas seladas, testes de acesso sem dados,
   um pedido mínimo por rota, caso único do DestinE pela opção A e estimativa revista de volume e tempo),
   com limite de ≤ 60 GB de tráfego e nenhuma repetição da campanha. **A campanha só começa com
   aprovação expressa do responsável após o relatório do pré-voo.**

Fora do escopo da fase 1: cálculos LEC, comparação de produtos e qualquer mudança em `src/`, testes, CI
ou dependências do repositório.

## Fora deste registro

Merge, CI e releitura deste documento não podem constar dele; ficam no relatório do Portão 1 e no
handoff à conversa-mãe.
