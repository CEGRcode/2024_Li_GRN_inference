import numpy as np
import pandas as pd
import math
import os

# Read the Excel file and specify the sheet name
df = pd.read_excel("./data/41586_2021_3314_MOESM3_ESM.xlsx", sheet_name="Supplementary Data 1", header=0)

# Define the allowed values
allowed_values = ["01_RP", "02_STM", "03_TFO", "04_UNB", "05_NoPIC", "06_tRNAprox", "07_ChExMix_extreme", "08_Hyper-variable"]

# Filter the dataframe
df_filtered = df[df["Feature class Level 1"].isin(allowed_values)]

# columns to output
columns_annotation = [
    "Chrom",
    "Strand",
    "Experiment_Left",
    "Experiment_Right",
    "Systematic ID",
    "Common Name",
    "NFR/NDR_Left",
    "NFR/NDR_Right"
]

columns_plus_one_nuc = ["PlusOne_ID", "Median Occupancy", "Median Variance"]
columns_minus_one_nuc = ["MinusOne_ID", "Median Occupancy.1", "Median Variance.1"]

renamed_columns_minus_one_nuc = {
    "MinusOne_ID": "MinusOne_ID",
    "Median Occupancy.1": "Median Occupancy",
    "Median Variance.1": "Median Variance"
}

df_genome_annotation = df_filtered[columns_annotation]
df_plus_one_nuc = df[columns_plus_one_nuc]
df_minus_one_nuc = df[columns_minus_one_nuc]
df_minus_one_nuc = df_minus_one_nuc.rename(columns=renamed_columns_minus_one_nuc)

# Save to tab-delimited text file
df_genome_annotation.to_csv("./result/Sc_genome_annotations.txt", sep="\t", index=False)
df_plus_one_nuc.to_csv("./result/PlusOneNuc.txt", sep="\t", index=False)
df_minus_one_nuc.to_csv("./result/MinusOneNuc.txt", sep="\t", index=False)
