# LINCS\_TFG\_USP19 Project

This project is focused on the analysis of perturbations using the **MDA-MB-231** cell line, in which the **USP19** gene has been knocked out. Different scripts are used to search for and analyze gene signatures connected with the perturbation of USP19, using **LINCS** databases and the **SigCom** search engine. The ultimate goal is to identify chemical compounds or perturbations that could reverse or mimic the pathological state of the cells.

## Index

1. [General Description](#general-description)
2. [Requirements](#requirements)
3. [Project Structure](#project-structure)

   * [perturbation\_Search.py](#perturbation_searchpy)
   * [applied\_Filters](#applied_filters)
   * [metadata\_Search](#metadata_search)
   * [z\_scores](#z_scores)
4. [Project Usage](#project-usage)

## General Description

The project employs a computational analysis to investigate cellular responses to various perturbations, considering differential gene expression in the **MDA-MB-231** cell line with the **USP19** gene knockout. Using **SigCom** and **LINCS** databases, it aims to identify **mimickers** (compounds that mimic the perturbation) and **reversers** (compounds that reverse the pathological state).

The workflow ranges from the search for perturbations to the analysis of their metadata, using a set of filters and metrics that help prioritize results.

## Requirements

* Python 3.x
* R
* Required libraries:

  * `numpy`
  * `pandas`
  * `requests`
  * `openxlsx` (R)
  * `dplyr` (R)
  * `ggplot2` (R)

## Project Structure

### perturbation\_Search.py

This script performs the search for perturbations related to **USP19** knockout. It uses *log2FoldChange* values to filter significantly upregulated and downregulated genes, then runs a search on **SigCom LINCS** databases. The results are stored in an **Excel** file named `all_scenarios_limit5.xlsx`.

* **Input:** `USP19_res_table.txt` (gene expression file)
* **Output:** `all_scenarios_limit5.xlsx` (file with perturbation IDs and associated Z values)
* **Libraries:** `numpy`, `pandas`, `requests`

### applied\_Filters

This R script filters and calculates additional metrics based on the results of the **USP19** perturbation search.

* **Input:** `all_scenarios_limit5.xlsx`
* **Output:** `filtered_data5.xlsx` (file with filtered perturbations and additional metrics)
* **Calculated Metrics:**

  * **quantity\_values:** Total number of thresholds passed by each perturbation.
  * **last\_col\_present:** Whether the perturbation passes the last threshold.
  * **last\_5\_cols:** Number of thresholds passed within the 5 most stringent.
  * **last\_10\_cols:** Number of thresholds passed within the 10 most stringent.
  * **max\_consecutive:** Maximum number of consecutive thresholds passed.
  * **weighted\_sum:** Weighted sum according to the threshold passed and the associated Z values.
  * **absolute\_sum:** Absolute sum of Z values.
  * **mimicker\_or\_reverser:** Classification of the perturbation as mimicker or reverser.

### metadata\_Search

This R script extracts metadata related to the perturbations filtered in the previous step. It uses **cellinfo\_beta** and **siginfo\_beta** tables downloaded from **CLUE.io** to enrich the data with additional information about cell lines and compounds.

* **Input:** `filtered_data5.xlsx`, `cellinfo_beta`, `siginfo_beta`
* **Output:** `metadata_top_filtered.xlsx`
* **Extracted Metadata:**

  * **cell\_iname:** Name of the cell line.
  * **cell\_type:** Type of cell (e.g., epithelial, mesenchymal, etc.).
  * **cell\_lineage:** Cell lineage (origin of the cell).
  * **donor\_ethnicity:** Ethnicity of the cell line donor.
  * **subtype:** Subtype of the cell, if applicable.
  * **primary\_disease:** Primary disease associated with the cell.
  * **pert\_id:** Compound/perturbation ID.
  * **pert\_dose:** Dose of the compound used.
  * **pert\_dose\_unit:** Unit of dose measurement (e.g., micromolar).
  * **pert\_time:** Exposure time to the perturbation.
  * **pert\_time\_unit:** Time unit (e.g., hours).
  * **pert\_type:** Type of perturbation (chemical, ligand, or antibody).
  * **MoR:** Indicator of whether the perturbation is a mimicker or a reverser.

### z\_scores

This R script calculates and visualizes the distribution of Z values obtained in previous analyses.

* **Input:** `all_scenarios_limit5.xlsx`
* **Output:** `perturbations_Zdistribution_limit5.csv`, `Zplot_limit5.png`
* **Calculated Metrics:**

  * **Maximum:** Maximum Z value for mimickers and reversers.
  * **Minimum:** Minimum Z value for mimickers and reversers.
  * **Mean:** Average Z-score for mimickers and reversers.
  * **Standard Deviation:** Standard deviation of Z-scores for mimickers and reversers.

The script also generates a histogram showing the distribution of Z values for mimickers and reversers.

## Project Usage

1. Run the `perturbation_Search.py` script to search for perturbations using USP19 gene expression data.
2. Filter and calculate additional metrics using the `applied_Filters` file.
3. Calculate and visualize Z-score distributions using the `z_scores` script.
4. Extract metadata related to the filtered perturbations with the `metadata_Search` script.
