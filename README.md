# Rchimerism

Pacote R para análise de quimerismo pós-transplante de células-tronco hematopoéticas, com interface Shiny interativa.

---

## Instalação

```r
install.packages("devtools")
devtools::install_github("joaohccampos/Rchimerism", subdir = "source")
```

## Execução

```r
library(Rchimerism)
Rchimerism()
```

---

## Modos de entrada

O app oferece três modos de entrada, selecionáveis na parte superior da interface:

### Normal Mode

Upload manual dos arquivos separados por papel:

| Campo | Arquivo | Conteúdo |
|---|---|---|
| Donor Data | `ddata.txt` | Perfil pré-transplante do doador |
| Recipient Data | `rdata.txt` | Perfil pré-transplante do receptor |
| Sample Data | `sdata.txt` | Amostra pós-transplante (quimera) |

Para **duplo doador**, aparece um campo adicional para `d2data.txt`.

### Directory Mode

Aponta para um diretório com a estrutura esperada pelo Rchimerism. Os arquivos `ddata.txt`, `rdata.txt` e `sdata.txt` são localizados automaticamente pelos caminhos relativos padrão.

### GeneMapper Mode

Importação direta do export bruto do GeneMapper — sem etapa de pré-processamento externo.

1. Exportar do GeneMapper como **"Sample/Peak report"** com as colunas: `Dye/Sample Peak`, `Sample File Name`, `Marker`, `Allele`, `Size`, `Height`, `Area`, `Data Point`
2. Fazer upload do arquivo `.txt` exportado
3. O app lista todas as amostras encontradas e tenta detectar o papel de cada uma automaticamente pelos padrões abaixo:

| Padrão no `Sample File Name` | Papel detectado |
|---|---|
| `DOADOR`, `DONOR`, `DON` | Doador |
| `PRE`, `RECEP`, `RECEPTOR`, `RECIPIENT` | Receptor |
| `QUI`, `CHIM`, `POS`, `POST` | Amostra quimera |

4. Confirmar ou corrigir o papel nos menus e clicar em **"Read input files"**

---

## Arquivo de marcadores

O campo "Marker File" é opcional. O painel **GlobalFiler (21 marcadores autossômicos)** já está pré-carregado.

Para usar um painel customizado, forneça um `.csv` de uma única linha com os nomes dos marcadores separados por vírgula, **exatamente** como aparecem no export do GeneMapper (ex: `vWA`, não `VWA`).

> **Três marcadores ligados ao sexo devem ser excluídos do painel:**
> - **AMELOGENINA** — alelos `X`/`Y` são alfabéticos e filtrados pelo Rchimerism
> - **DYS391** — ausente em amostras femininas
> - **Yindel** — ausente em amostras femininas
>
> O painel GlobalFiler bundled já exclui esses três.

---

## Resultados

Após clicar em **"Read input files"**, o app exibe:

- Tabela de resultados por marcador com percentual de quimerismo
- Estatísticas agregadas: Donor% Mean, SD e CV (outliers excluídos automaticamente)
- Matrizes de alelos (doador, receptor, amostra)
- Botão **"Download Results Excel File"** — exporta resultados e estatísticas em `.xlsx`
- Botão **"Download Check File"** — exporta matrizes completas em `.txt` para auditoria

Use **"New Analysis"** para limpar todos os campos e iniciar uma nova análise sem reiniciar o R.
