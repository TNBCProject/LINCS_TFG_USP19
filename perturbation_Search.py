import pandas as pd
import numpy as np
import requests

#Lectura del archivo. Contiene la expresión génica de MDA-MB-231 con la anulación de USP19
USP19_df = pd.read_csv('USP19_res_table.txt', delimiter='\t')

#Filtrado de valores estadísticamente significativos
USP19_df = USP19_df.loc[USP19_df["pvalue"] < .05]

#Selección de la columna conteniendo los valores de fold change
log2foldChange = USP19_df["log2FoldChange"]

#Separación de los valores de fold change entre genes regulados positiva y negativamente 
foldChange_UP = [x for x in log2foldChange if x > 0]
foldChange_DOWN = [x for x in log2foldChange if x < 0]

#Calculo de la mediana para ambos grupos de genes
median_UP = np.median(foldChange_UP)
median_DOWN = np.median(foldChange_DOWN)

# Rango de valores del umbral
values_UP = np.arange(median_UP, median_UP + 0.6, 0.01)
values_DOWN = np.arange(median_DOWN, median_DOWN - 0.3, -0.005)


all_scenarios = []
ftime = True
temp = 60

# Definición de las APIs a utilizar
METADATA_API = "https://maayanlab.cloud/sigcom-lincs/metadata-api/"
DATA_API = "https://maayanlab.cloud/sigcom-lincs/data-api/api/v1/"

# Definición de las bases de datos
all_databases = ["human_GEO", "gtex_age_sigs", "LINCS chemical perturbagen signatures",
                     "l1000_shRNA", "l1000_siRNA", "LINCS gene overexpression signatures",
                     "l1000_lig", "LINCS consensus gene (CGS) knockdown signatures", "l1000_mean_xpr",
                     "l1000_mean_cp", "l1000_cp", "l1000_aby", "l1000_xpr", "l1000_oe"]

# Bucle para laa búsqueda y guardado de perturbaciones. 
# Plantea 60 escenarios desde la mediana hasta que solo un gen pasa el umbral.
for i in range(temp):

    # Umbral perteneciente al escenario correspondiente
    threshold_UP = values_UP[i]
    threshold_DOWN = values_DOWN[i]

    # Selección de los genes
    up_genes = USP19_df.loc[USP19_df["log2FoldChange"] > threshold_UP]["gene"].tolist()
    down_genes = USP19_df.loc[USP19_df["log2FoldChange"] < threshold_DOWN]["gene"].tolist()

    input_gene_set = {
        "up_genes": up_genes,
        "down_genes": down_genes
    }

    all_genes = input_gene_set["up_genes"] + input_gene_set["down_genes"]
    

    first_time = True
    all_dfs = []

    # Bucle para buscar en las bases de datos de SigCom
    for mydb in all_databases:

        #Conversión de los nombres de los genes a UUIDs
        payload = {
            "filter": {
                "where": {
                    "meta.symbol": {
                        "inq": all_genes
                    }
                },
                "fields": ["id", "meta.symbol"]
            }
        }

        res = requests.post(METADATA_API + "entities/find", json=payload)

        entities = res.json()

        for_enrichment = {
            "up_entities": [],
            "down_entities": []
        }

        for e in entities:
            symbol = e["meta"]["symbol"]
            if symbol in input_gene_set["up_genes"]:
                for_enrichment["up_entities"].append(e["id"])
            elif symbol in input_gene_set["down_genes"]:
                for_enrichment["down_entities"].append(e["id"])

        # Búsqueda de perturbaciones usando ranktwosided
        query = {
            **for_enrichment,
            "limit": 5,
            "database": mydb
        }

        res = requests.post(DATA_API + "enrich/ranktwosided", json=query)
        results = res.json()

        # Cambia los valores Z de los reversers para que sean negativos
        for ii in results["results"]:
            ii["z-down"] = -ii["z-down"]
            ii["direction-down"] = -ii["direction-down"]


        sigids = {i["uuid"]: i for i in results["results"]}

        # Recuperación de la metadata
        payload = {
            "filter": {
                "where": {
                    "id": {
                        "inq": list(sigids.keys())
                    }
                }
            }
        }

        res = requests.post(METADATA_API + "signatures/find", json=payload)
        signatures = res.json()

        # Conectar los valores Z con la metadata
        for sig in signatures:
            uid = sig["id"]
            scores = sigids[uid]
            scores.pop("uuid")
            sig["scores"] = scores

        # Guarda las perturbaciones y sus valores Z en una lista.
        # Contiene los resultados de todas las bases de datos
        all_dfs.append(pd.json_normalize(signatures))

    # Selección de las 3 bases de datos a estudiar
    ''' 
    Bases de datos utilizadas: 
        "LINCS L1000 Antibody Perturbations (2021)": "l1000_aby"
        "LINCS L1000 Ligand Perturbations (2021)": "l1000_lig"
        "LINCS L1000 Chemical Perturbations (2021)": "l1000_cp"
    '''
    df_merged_1 = pd.merge(all_dfs[6], all_dfs[10], how='outer')
    df_merged = pd.merge(df_merged_1, all_dfs[11], how='outer')

    # Guardado de los datos encontrados en un archivo excel. Contiene el ID de las perturbaciones, su valor Z y el escenario asociado.

    selected_columns = ["meta.tissue",
                        "meta.cmap_id",
                        "meta.disease",
                        "meta.cell_line",
                        "meta.pert_dose",
                        "meta.pert_name",
                        "meta.pert_time",
                        "meta.pert_type",
                        "scores.z-sum"]

    df_merged = df_merged[selected_columns]

    if ftime:
        all_scenarios = df_merged[["meta.cmap_id", "scores.z-sum"]]
        all_scenarios = all_scenarios.rename(columns={'scores.z-sum': f"[{round(values_UP[i], 3)};({round(values_DOWN[i], 3)})]"})
        ftime = False
    else:
        df_merged = df_merged[["meta.cmap_id", "scores.z-sum"]]
        df_merged = df_merged.rename(columns={'scores.z-sum': f"[{round(values_UP[i], 3)};({round(values_DOWN[i], 3)})]"})

        all_scenarios = pd.merge(all_scenarios, df_merged, how="outer")

    all_scenarios.to_excel("all_scenarios_limit5.xlsx", index=False)
