# Proyecto LINCS_TFG_USP19

Este proyecto está orientado al análisis de perturbaciones utilizando la línea celular MDA-MB-231, en la que se ha anulado el gen **USP19**. Se utilizan diferentes scripts para buscar y analizar firmas genéticas conectadas con la perturbación de USP19, utilizando bases de datos de **LINCS** y el motor de busqueda **SigCom**. El objetivo final es identificar compuestos químicos o perturbaciones que puedan revertir o mimetizar el estado patológico de las células.

## Índice

1. [Descripción General](#descripción-general)
2. [Requisitos](#requisitos)
3. [Estructura del Proyecto](#estructura-del-proyecto)
   - [perturbation_Search.py](#perturbation_searchpy)
   - [applied_Filters](#applied_filters)
   - [metadata_Search](#metadata_search)
   - [z_scores](#z_scores)
4. [Uso del Proyecto](#uso-del-proyecto)
5. [Contacto](#contacto)

## Descripción General

El proyecto emplea un análisis computacional para investigar las respuestas de las células a distintas perturbaciones, tomando en cuenta la expresión génica diferencial en la línea celular **MDA-MB-231** con la anulación del gen **USP19**. Mediante el uso de bases de datos **SigCom** y **LINCS**, se busca identificar mimickers (compuestos que mimetizan la perturbación) y reversers (compuestos que revierten el estado patológico).

El flujo de trabajo abarca desde la búsqueda de perturbaciones hasta el análisis de la metadata de estas perturbaciones, utilizando un conjunto de filtros y métricas que ayudan a priorizar los resultados.

## Requisitos

- Python 3.x
- R
- Librerías necesarias:
  - `numpy`
  - `pandas`
  - `requests`
  - `openxlsx` (R)
  - `dplyr` (R)
  - `ggplot2` (R)
  
## Estructura del Proyecto

### perturbation_Search.py

Este script ejecuta la búsqueda de perturbaciones relacionadas con la anulación de **USP19**. Utiliza los valores de *log2FoldChange* para filtrar los genes significativamente regulados de forma positiva y negativa, y luego ejecuta una búsqueda en bases de datos de **SigCom LINCS**. Los resultados son almacenados en un archivo **Excel** llamado `all_scenarios_limit5.xlsx`.

- **Input:** `USP19_res_table.txt` (archivo con expresión génica)
- **Output:** `all_scenarios_limit5.xlsx` (archivo con los ID de las perturbaciones y los valores Z asociados)
- **Librerías:** `numpy`, `pandas`, `requests`
  
### applied_Filters

Este archivo de R filtra y calcula métricas adicionales basadas en los resultados de la búsqueda de perturbaciones de **USP19**.

- **Input:** `all_scenarios_limit5.xlsx`
- **Output:** `filtered_data5.xlsx` (archivo con las perturbaciones filtradas y métricas adicionales)
- **Métricas Calculadas:**
  - **quantity_values:** Número total de umbrales pasados por cada perturbación.
  - **last_col_present:** Si la perturbación pasa o no el último umbral.
  - **last_5_cols:** Número de umbrales pasados dentro de los 5 más estrictos.
  - **last_10_cols:** Número de umbrales pasados dentro de los 10 más estrictos.
  - **max_consecutive:** Cantidad máxima de umbrales consecutivos pasados.
  - **weighted_sum:** Suma ponderada según el umbral pasado y los valores Z asociados.
  - **absolute_sum:** Suma absoluta de los valores Z.
  - **mimicker_or_reverser:** Clasificación de la perturbación como mimicker o reverser.

### metadata_Search

Este archivo de R extrae la metadata relacionada con las perturbaciones filtradas en el paso anterior. Utiliza las tablas **cellinfo_beta** y **siginfo_beta** descargadas de **CLUE.io** para enriquecer los datos con información adicional sobre las líneas celulares y los compuestos.

- **Input:** `filtered_data5.xlsx`, `cellinfo_beta`, `siginfo_beta`
- **Output:** `metadata_top_filtered.xlsx`
- **Metadata Extraída:**
  - **cell_iname:** Nombre de la línea celular.
  - **cell_type:** Tipo de célula (ej. epitelial, mesenquimal, etc.).
  - **cell_lineage:** Linaje celular (origen de la célula).
  - **donor_ethnicity:** Etnicidad del donante de la línea celular.
  - **subtype:** Subtipo de la célula, si aplica.
  - **primary_disease:** Enfermedad primaria asociada a la célula.
  - **pert_id:** ID del compuesto/perturbación.
  - **pert_dose:** Dosis del compuesto utilizado.
  - **pert_dose_unit:** Unidad de medida de la dosis (ej. micromolar).
  - **pert_time:** Tiempo de exposición a la perturbación.
  - **pert_time_unit:** Unidad de tiempo (ej. horas).
  - **pert_type:** Tipo de perturbación (químico, ligando o anticuerpo).
  - **MoR:** Indicador de si la perturbación es un mimicker o un reverser.

### z_scores

Este archivo de R se encarga de calcular y visualizar la distribución de los valores Z obtenidos en los análisis anteriores.

- **Input:** `all_scenarios_limit5.xlsx`
- **Output:** `perturbations_Zdistribution_limit5.csv`, `Zplot_limit5.png`
- **Métricas Calculadas:**
  - **Maximum:** Valor Z máximo de mimickers y reversers.
  - **Minimum:** Valor Z mínimo de mimickers y reversers.
  - **Mean:** Valor medio de los Z-scores de mimickers y reversers.
  - **Standard Deviation:** Desviación estándar de los Z-scores de mimickers y reversers.

El script también genera una visualización de la distribución de los valores Z para mimickers y reversers en un gráfico de histograma.

## Uso del Proyecto

1. Ejecutar el script `perturbation_Search.py` para buscar las perturbaciones utilizando los datos de expresión génica de USP19.
2. Filtrar y calcular métricas adicionales usando el archivo `applied_Filters`.
3. Calcular y visualizar las distribuciones de valores Z usando el script `z_scores`.
4. Extraer la metadata relacionada a las perturbaciones filtradas con el script `metadata_Search`

## Contacto

Para más información sobre este proyecto, contacta a:

**Nombre:** Agustina Verschoor, Camila Contestabile. 
**Email:** verschooragustina@gmail.com, camilacontestabile@gmail.com

