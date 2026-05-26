# Rchimerism — Workflow

Este repositório contém o pacote R **Rchimerism** para análise de quimerismo pós-transplante de células-tronco hematopoéticas, e um **conversor Python** para preparar exportações do GeneMapper para o formato esperado pelo Rchimerism.

---

## Visão geral do fluxo

```
GeneMapper (.txt)
       │
       ▼
genemapper-converter     ←── markers.csv
       │
       ▼
  ddata.txt / rdata.txt / sdata.txt
       │
       ▼
   Rchimerism()
       │
       ▼
  Resultados de quimerismo
```

---

## 1. Pré-processamento — Conversor GeneMapper

### Exportação do GeneMapper

Exportar como **"Sample/Peak report"** (Table Setting com as colunas: `Dye/Sample Peak`, `Sample File Name`, `Marker`, `Allele`, `Size`, `Height`, `Area`, `Data Point`). O arquivo exportado já está no formato longo correto — o conversor apenas divide as amostras em arquivos separados.

### Instalação

Requer Python ≥ 3.12 e [uv](https://github.com/astral-sh/uv).

```bash
cd converter
uv sync
```

### Comandos

#### Listar amostras (com detecção automática de papéis)

```bash
uv run python -m genemapper_converter list caminho/para/arquivo.txt
```

Exemplo de saída:
```
Found 3 sample(s):
  '0028MD01-071020 01D SG DOADOR_B02.fsa'  → auto-detected as: donor
  '0028MR-071020 01D SG PRE_A02.fsa'       → auto-detected as: recipient
  'CP_QUI353_F03.fsa'                       → auto-detected as: sample
```

#### Converter com detecção automática

Detecta papéis automaticamente por padrões no `Sample File Name`:

| Padrão no nome | Papel detectado | Arquivo gerado |
|---|---|---|
| `DOADOR`, `DONOR`, `DON` | Doador | `ddata.txt` |
| `PRE`, `RECEP`, `RECEPTOR` | Receptor | `rdata.txt` |
| `QUI`, `CHIM`, `POS`, `POST` | Amostra quimera | `sdata.txt` |

```bash
uv run python -m genemapper_converter convert arquivo.txt \
    --output-dir ./saida \
    --auto-detect
```

#### Converter com nomes explícitos

Use quando o `Sample File Name` não segue os padrões acima:

```bash
# Doador único
uv run python -m genemapper_converter convert arquivo.txt \
    --output-dir ./saida \
    --donor   "NOME_EXATO_B02.fsa" \
    --recipient "NOME_EXATO_A02.fsa" \
    --sample  "NOME_EXATO_F03.fsa"

# Duplo doador (cordão umbilical)
uv run python -m genemapper_converter convert arquivo.txt \
    --output-dir ./saida \
    --donor1  "NOME_DOADOR1_B02.fsa" \
    --donor2  "NOME_DOADOR2_C02.fsa" \
    --recipient "NOME_RECEPTOR_A02.fsa" \
    --sample  "NOME_AMOSTRA_F03.fsa"
```

### O que o conversor faz

- Lê o export do GeneMapper (já no formato longo — 1 linha por pico)
- Divide o arquivo em arquivos separados por papel (doador, receptor, amostra)
- Preserva todos os picos, incluindo picos sem marcador atribuído

### Estrutura dos arquivos gerados

```
saida/
├── ddata.txt    ← dados pré-transplante do doador
├── rdata.txt    ← dados pré-transplante do receptor
└── sdata.txt    ← dados pós-transplante (amostra quimera)
```

Para **duplo doador**, são gerados `d1data.txt` e `d2data.txt` no lugar de `ddata.txt`.

---

## 2. Arquivo de marcadores (`markers.csv`)

O `markers.csv` define quais marcadores serão usados na análise de quimerismo. Deve estar no formato CSV de uma única linha, com os nomes **exatamente iguais** aos que aparecem nos arquivos de dados convertidos.

### Painel GlobalFiler (21 marcadores autossômicos)

```
D3S1358,vWA,D16S539,CSF1PO,TPOX,D8S1179,D21S11,D18S51,D2S441,D19S433,TH01,FGA,D22S1045,D5S818,D13S317,D7S820,SE33,D10S1248,D1S1656,D12S391,D2S1338
```

O arquivo pronto está em: [`example/convertedExample/markers.csv`](example/convertedExample/markers.csv)

> **Atenção — três marcadores ligados ao sexo devem ser excluídos:**
> - **AMELOGENINA**: alelos `X` e `Y` são alfabéticos; o Rchimerism os filtra e gera erro `'AMEL' from markers not found in input data`.
> - **DYS391**: marcador do cromossomo Y; ausente em amostras femininas, causa o mesmo erro.
> - **Yindel**: marcador do cromossomo Y; mesmo problema.
>
> Esses três marcadores são pré-carregados como excluídos no padrão GlobalFiler bundled no pacote.

> **Nomenclatura:**
> Os nomes devem corresponder exatamente ao que o GeneMapper exporta:
> - `vWA` (não `VWA`)

---

## 3. Análise no Rchimerism

### Instalação do Rchimerism

```r
install.packages("devtools")
install.packages("shinyFiles")
devtools::install_github("joaohccampos/Rchimerism", subdir = "source")
library(Rchimerism)
```

### Execução

```r
Rchimerism()
```

Na interface Shiny:

1. Carregar `ddata.txt` (ou `d1data.txt` + `d2data.txt` para duplo doador)
2. Carregar `rdata.txt`
3. Carregar `sdata.txt`
4. Clicar em **"Read input files"**

> O painel GlobalFiler (21 marcadores autossômicos) já está pré-carregado. O campo "Marker File" pode ser ignorado, ou substituído por outro `markers.csv` se necessário.

---

## Exemplo completo (dados reais)

Os arquivos de exemplo convertidos estão em [`example/convertedExample/`](example/convertedExample/):

| Arquivo | Amostra original | Papel |
|---|---|---|
| `ddata.txt` | `0028MD01-071020 01D SG DOADOR` | Doador |
| `rdata.txt` | `0028MR-071020 01D SG PRE` | Receptor (pré-transplante) |
| `sdata.txt` | `CP_QUI353` | Amostra quimera (pós-transplante) |
| `markers.csv` | — | Painel GlobalFiler (21 marcadores autossômicos) |

```bash
# Reproduzir a conversão do exemplo
cd converter
uv run python -m genemapper_converter convert ../genemapper_export.txt \
    --output-dir ../example/convertedExample \
    --auto-detect
```
