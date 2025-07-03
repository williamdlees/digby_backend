import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, Polygon
import re
import tempfile
from itertools import product
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly.io as pio
from collections import defaultdict

GENE_loc = {
    'IGH': ['IGHV7-81', 'IGHV5-78', 'IGHV3-74', 'IGHV3-73', 'IGHV3-72', 'IGHV3-71', 'IGHV2-70', 'IGHV1-69D', 'IGHV1-f', 'IGHV2-70D', 'IGHV1-69', 'IGHV1-68', 'IGHV1-69-2', 'IGHV3-69-1', 'IGHV3-66', 'IGHV3-64', 'IGHV3-63', 'IGHV3-62', 'IGHV4-61', 'IGHV4-59', 'IGHV1-58', 'IGHV4-55', 'IGHV3-54', 'IGHV3-53', 'IGHV3-52', 'IGHV5-51', 'IGHV3-49', 'IGHV3-48', 'IGHV3-47', 'IGHV1-46', 'IGHV1-45', 'IGHV3-43', 'IGHV3-43D', 'IGHV7-40', 'IGHV4-39', 'IGHV1-38-4', 'IGHV3-38-3', 'IGHV4-38-2', 'IGHV3-38', 'IGHV3-35', 'IGHV7-34-1', 'IGHV4-34', 'IGHV3-33-2', 'IGHV3-33', 'IGHV3-32', 'IGHV4-31', 'IGHV3-30-52', 'IGHV3-30-5', 'IGHV3-30-42', 'IGHV4-30-4', 'IGHV3-30-33', 'IGHV3-30-3', 'IGHV3-30-22', 'IGHV4-30-2', 'IGHV4-30-1', 'IGHV3-30-2', 'IGHV3-30', 'IGHV3-29', 'IGHV4-28', 'IGHV2-26', 'IGHV3-25', 'IGHV1-24', 'IGHV3-23', 'IGHV3-23D', 'IGHV3-22', 'IGHV3-21', 'IGHV3-20', 'IGHV3-19', 'IGHV1-18', 'IGHV3-16', 'IGHV3-15', 'IGHV3-13', 'IGHV3-11', 'IGHV2-10', 'IGHV3-9', 'IGHV5-10-1', 'IGHV1-8', 'IGHV3-64D', 'IGHV3-7', 'IGHV2-5', 'IGHV7-4-1', 'IGHV4-4', 'IGHV1-3', 'IGHV1-2', 'IGHV6-1', 'IGHD1-1', 'IGHD2-2', 'IGHD3-3', 'IGHD6-6', 'IGHD1-7', 'IGHD2-8', 'IGHD3-9', 'IGHD3-10', 'IGHD4-11', 'IGHD5-12', 'IGHD6-13', 'IGHD1-14', 'IGHD2-15', 'IGHD3-16', 'IGHD4-17', 'IGHD5-18', 'IGHD6-19', 'IGHD1-20', 'IGHD2-21', 'IGHD3-22', 'IGHD4-23', 'IGHD5-24', 'IGHD6-25', 'IGHD1-26', 'IGHD7-27', 'IGHJ1', 'IGHJ2', 'IGHJ3', 'IGHJ4', 'IGHJ5', 'IGHJ6'],
    'IGK': ['IGKJ5', 'IGKJ4', 'IGKJ3', 'IGKJ2', 'IGKJ1', 'IGKV4-1', 'IGKV5-2', 'IGKV7-3', 'IGKV2-4', 'IGKV1-5', 'IGKV1-6', 'IGKV3-7', 'IGKV1-8', 'IGKV1-9', 'IGKV2-10', 'IGKV3-11', 'IGKV1-12', 'IGKV1-13', 'IGKV2-14', 'IGKV3-15', 'IGKV1-16', 'IGKV1-17', 'IGKV2-18', 'IGKV2-19', 'IGKV3-20', 'IGKV6-21', 'IGKV1-22', 'IGKV2-23', 'IGKV2-24', 'IGKV3-25', 'IGKV2-26', 'IGKV1-27', 'IGKV2-28', 'IGKV2-29', 'IGKV2-30', 'IGKV3-31', 'IGKV1-32', 'IGKV1-33', 'IGKV3-34', 'IGKV1-35', 'IGKV2-36', 'IGKV1-37', 'IGKV2-38', 'IGKV1-39', 'IGKV2-40', 'IGKV2D-40', 'IGKV1D-39', 'IGKV2D-38', 'IGKV1D-37', 'IGKV2D-36', 'IGKV1D-35', 'IGKV3D-34', 'IGKV1D-33', 'IGKV1D-32', 'IGKV3D-31', 'IGKV2D-30', 'IGKV2D-29', 'IGKV2D-28', 'IGKV1D-27', 'IGKV2D-26', 'IGKV3D-25', 'IGKV2D-24', 'IGKV2D-23', 'IGKV1D-22', 'IGKV6D-21', 'IGKV3D-20', 'IGKV2D-19', 'IGKV2D-18', 'IGKV6D-41', 'IGKV1D-17', 'IGKV1D-16', 'IGKV3D-15', 'IGKV2D-14', 'IGKV1D-13', 'IGKV1D-12', 'IGKV3D-11', 'IGKV2D-10', 'IGKV1D-42', 'IGKV1D-43', 'IGKV1D-8', 'IGKV3D-7', 'IGKV1-NL1'],
    'IGL': ['IGLV4-69', 'IGLV8-61', 'IGLV4-60', 'IGLV6-57', 'IGLV10-54', 'IGLV5-52', 'IGLV1-51', 'IGLV9-49', 'IGLV1-47', 'IGLV7-46', 'IGLV5-45', 'IGLV1-44', 'IGLV7-43', 'IGLV1-40', 'IGLV5-37', 'IGLV1-36', 'IGLV3-27', 'IGLV3-25', 'IGLV2-23', 'IGLV3-22', 'IGLV3-21', 'IGLV3-19', 'IGLV2-18', 'IGLV3-16', 'IGLV2-14', 'IGLV3-12', 'IGLV2-11', 'IGLV3-10', 'IGLV3-9', 'IGLV2-8', 'IGLV4-3', 'IGLV3-1', 'IGLJ1', 'IGLJ2', 'IGLJ3', 'IGLJ6', 'IGLJ7']
}

PSEUDO = {
    'IGH': ["IGHV2-10", "IGHV3-52", "IGHV3-47", "IGHV3-71", "IGHV3-22", "IGHV4-55", "IGHV1-68", "IGHV2-10", "IGHV5-78", "IGHV3-32", "IGHV3-33-2", "IGHV3-38-3", "IGHV3-25", "IGHV3-19", "IGHV7-40", "IGHV3-63", "IGHV3-62", "IGHV3-29", "IGHV3-54", "IGHV1-38-4", "IGHV7-34-1", "IGHV1-38-4", "IGHV3-30-2", "IGHV3-69-1", "IGHV3-30-22", "IGHV1-f", "IGHV3-30-33", "IGHV3-38", "IGHV7-81", "IGHV3-35", "IGHV3-16", "IGHV3-30-52", "IGHV1-69D", "IGHD1-14", "IGHV3-30-42"],
    'IGK': ["IGKV7-3", "IGKV6D-41", "IGKV3D-34", "IGKV3D-31", "IGKV3D-25", "IGKV3/OR22-2", "IGKV3/OR2-5", "IGKV3/OR2-268", "IGKV3-34", "IGKV3-31", "IGKV3-25", "IGKV2D-38", "IGKV2D-36", "IGKV2D-24", "IGKV2D-23", "IGKV2D-19", "IGKV2D-18", "IGKV2D-14", "IGKV2D-10", "IGKV2/OR22-4", "IGKV2/OR22-3", "IGKV2/OR2-8", "IGKV2/OR2-7D", "IGKV2/OR2-7", "IGKV2/OR2-4", "IGKV2/OR2-2", "IGKV2/OR2-10", "IGKV2/OR2-1", "IGKV2-4", "IGKV2-38", "IGKV2-36", "IGKV2-26", "IGKV2-23", "IGKV2-19", "IGKV2-18", "IGKV2-14", "IGKV2-10", "IGKV1D-42", "IGKV1D-37", "IGKV1D-35", "IGKV1D-32", "IGKV1D-27", "IGKV1D-22", "IGKV1/ORY-1", "IGKV1/OR9-2", "IGKV1/OR9-1", "IGKV1/OR22-5", "IGKV1/OR22-1", "IGKV1/OR2-9", "IGKV1/OR2-6", "IGKV1/OR2-3", "IGKV1/OR2-2", "IGKV1/OR2-118", "IGKV1/OR2-11", "IGKV1/OR2-108", "IGKV1/OR2-1", "IGKV1/OR2-0", "IGKV1/OR15-118", "IGKV1/OR10-1", "IGKV1/OR1-1", "IGKV1/OR-4", "IGKV1/OR-3", "IGKV1/OR-2", "IGKV1-37", "IGKV1-35", "IGKV1-32", "IGKV1-22"],
    'IGL': ["IGLV7-35", "IGLV3-7", "IGLV3-6", "IGLV3-4", "IGLV3-32", "IGLV3-31", "IGLV3-30", "IGLV3-29", "IGLV3-26", "IGLV3-24", "IGLV3-2", "IGLV3-17", "IGLV3-15", "IGLV3-13", "IGLV2-NL1", "IGLV2-5", "IGLV2-34", "IGLV2-33", "IGLV2-28", "IGLV11-55", "IGLV10-67", "IGLV1-62", "IGLV1-50", "IGLV(VII)-41-1", "IGLV(VI)-25-1", "IGLV(VI)-22-1", "IGLV(V)-66", "IGLV(V)-58", "IGLV(IV)/OR22-2", "IGLV(IV)/OR22-1", "IGLV(IV)-66-1", "IGLV(IV)-65", "IGLV(IV)-64", "IGLV(IV)-59", "IGLV(IV)-53", "IGLV(I)-70", "IGLV(I)-68", "IGLV(I)-63", "IGLV(I)-56", "IGLV(I)-42", "IGLV(I)-38", "IGLV(I)-20", "IGLL4", "IGLL2", "IGLJ5", "IGLJ4", "IGLC5", "IGLC4", "IGLC/OR22-2", "IGLC/OR22-1"]
}

# Allele palette colors
ALLELE_PALETTE = {
    "01": "#f5bc6e", "02": "#9d69f4", "03": "#598200", "04": "#e375fd", 
    "05": "#01a452", "06": "#ff6ed5", "07": "#60de8a", "08": "#e3006d", 
    "09": "#0edec9", "10": "#c43e00", "11": "#0197f3", "12": "#ffa629", 
    "13": "#7c82ff", "14": "#8dd979", "15": "#9a0f79", "16": "#008844", 
    "17": "#ff64aa", "18": "#ff5d68", "19": "#ff88ce", "20": "#00633c", 
    "21": "#355d06", "22": "#009f8b", "23": "#ff7d4e", "24": "#02a8ec", 
    "25": "#7d4408", "26": "#89c8ff", "i01": "#ff897a", "i02": "#58d8e0", 
    "i03": "#ff8793", "i04": "#82d89e", "NA": "#803a6a"
}

def allele_palette(hap_alleles, nra=True):
    """
    Create the allele color palette for haplotype graphical output.
    
    Parameters:
    -----------
    hap_alleles : list
        List of haplotype alleles
    nra : bool
        Include NRA (Non-Reliable Allele)
        
    Returns:
    --------
    dict
        Color palette for alleles
    """
    # Extract alleles with digits
    alleles = [a for a in set(hap_alleles) if re.search(r'[012]', a)]
    
    # Get unique allele base names
    allele_col_tmp = sorted(set([a.split('_')[0] for a in alleles]))
    
    # Assign colors to base alleles
    tmp_col = {allele: ALLELE_PALETTE.get(allele, "#000000") for allele in allele_col_tmp}
    
    # Handle novel alleles (with underscore)
    novels = [a for a in alleles if '_' in a]
    allele_col = tmp_col.copy()
    
    if novels:
        # Assign colors to novel alleles
        novel_colors = {novel: ALLELE_PALETTE.get(novel.split('_')[0], "#000000") for novel in novels}
        allele_col.update(novel_colors)
        
    # Add special colors
    allele_col.update({
        "Unk": "#dedede", 
        "Del": "#6d6d6d", 
        "NR": "#000000", 
        "NRA": "#fbf7f5"
    })
    
    # Remove alleles not in the input
    for allele in ["NR", "Del", "Unk", "NRA"]:
        if allele not in hap_alleles:
            allele_col.pop(allele, None)
    
    # Calculate transparency for alleles
    transper = {}
    for allele in allele_col:
        if '_' in allele:
            mom_allele = allele.split('_')[0]
            all_novel = [a for a in allele_col if a.startswith(f"{mom_allele}_")]
            if len(all_novel) == 1:
                transper[allele] = 0.5
            elif len(all_novel) == 2:
                m = all_novel.index(allele)
                transper[allele] = 0.6 if m == 0 else 0.3
            elif len(all_novel) == 3:
                m = all_novel.index(allele)
                if m == 0:
                    transper[allele] = 0.6
                else:
                    transper[allele] = 0.4 if m == 0 else 0.2
            elif len(all_novel) > 9:
                m = all_novel.index(allele)
                if m == 0:
                    transper[allele] = 1
                else:
                    transper[allele] = 1 - m/20
            elif len(all_novel) > 3:
                m = all_novel.index(allele)
                if m == 0:
                    transper[allele] = 0.85
                else:
                    transper[allele] = 0.85 - m/10
        else:
            transper[allele] = 1
            
    # Filter alleles
    special = [a for a in ["Unk", "Del", "NR", "NRA"] if a in allele_col]
    allele_col_filtered = {a: allele_col[a] for a in allele_col if a in alleles + special}
    transper_filtered = {a: transper[a] for a in transper if a in allele_col_filtered}
    
    return {
        "transper": transper_filtered,
        "AlleleCol": allele_col_filtered
    }


def sort_df_by_gene(data, chain="IGH", method="position", remove_chain=False, 
                    geno=False, pseudo_remove=False, orf_remove=False):
    """
    Sort the dataframe by gene names or position.
    
    Parameters:
    -----------
    data : DataFrame
        Data frame to sort
    chain : str
        The Ig chain: IGH, IGK, IGL
    method : str
        Sorting method: 'name' or 'position'
    removeIGH : bool
        If True, remove IGH/IGK/IGL prefix from gene names
    geno : bool
        For genotype data
    pseudo_remove : bool
        Remove pseudo genes when False (matches R function behavior)
    orf_remove : bool
        Remove ORF genes when False (matches R function behavior)
    
    Returns:
    --------
    DataFrame
        Sorted data frame
    """
    # Get gene locations for the specified chain
    GENE_loc_tmp = GENE_loc[chain].copy()
    
    vs = [g for g in GENE_loc_tmp if "V" in g]
    ps = [g for g in PSEUDO[chain] if "V" in g]
    ps = [g for g in ps if g not in vs]
    orf = [g for g in data["gene"].unique() if isinstance(g, str) and re.search(r"OR|NL", g)]
    non_v_genes = [g for g in GENE_loc_tmp if not re.search(r"V", g)]
    GENE_loc_tmp = non_v_genes + vs + ps + orf
    
    # Note: Logic reversed to match R function
    if pseudo_remove:
        data = data[~data["gene"].isin(PSEUDO[chain])]
        GENE_loc_tmp = [g for g in GENE_loc_tmp if g not in ps]
    
    if orf_remove:
        data = data[~data["gene"].isin(orf)]
        GENE_loc_tmp = [g for g in GENE_loc_tmp if g not in orf]
    
    if method == "name":
        # Sort by gene name
        unique_genes = sorted(data["gene"].unique())
        data["gene"] = pd.Categorical(data["gene"], categories=unique_genes)
        data = data.sort_values("gene")
        if remove_chain:
            data["gene"] = data["gene"].str.replace(r"IG[HKL]", "", regex=True)
            unique_genes = sorted(data["gene"].unique())
            data["gene"] = pd.Categorical(data["gene"], categories=unique_genes)
            data = data.sort_values("gene")
            if not geno and "hapBy" in data.columns:
                data["hapBy"] = data["hapBy"].str.replace(r"IG[HKL]", "", regex=True)
    else:
        # Sort by position
        data["gene"] = pd.Categorical(data["gene"], categories=GENE_loc_tmp)
        data = data.sort_values("gene")
        if remove_chain:
            # Remove IGH/IGK/IGL prefix from gene locations
            GENE_loc_tmp = [re.sub(r"IG[HKL]", "", g) for g in GENE_loc_tmp]
            data["gene"] = data["gene"].str.replace(r"IG[HKL]", "", regex=True)
            data["gene"] = pd.Categorical(data["gene"], categories=GENE_loc_tmp)
            data = data.sort_values("gene")
            if not geno and "hapBy" in data.columns:
                data["hapBy"] = data["hapBy"].str.replace(r"IG[HKL]", "", regex=True)
    return data

def explode_hap_table(hap_table,hapBy_cols):
    records = []
    for _, row in hap_table.iterrows():
        subject = row["subject"]
        gene = row["gene"]
        alleles = str(row["alleles"]).split(",") if pd.notna(row["alleles"]) else []
        col1_alleles = str(row[hapBy_cols[0]]).split(",") if pd.notna(row[hapBy_cols[0]]) else []
        col2_alleles = str(row[hapBy_cols[1]]).split(",") if pd.notna(row[hapBy_cols[1]]) else []
        for i, allele in enumerate(alleles):
            count_col = f"counts{i+1}"
            k_col = f"k{i+1}"
            count_val = row.get(count_col, np.nan)
            k_val = row.get(k_col, np.nan)
            # Split count into parts (if string like "2,0")
            if pd.notna(count_val) and isinstance(count_val, str) and "," in count_val:
                count_val = count_val.split(",")
            else:
                count_val = [count_val, count_val]  # duplicate single value for both sides
            allele_in_1 = allele in col1_alleles
            allele_in_2 = allele in col2_alleles
            if allele_in_1:
                records.append({
                    "subject": subject,
                    "gene": gene,
                    "hap_col": hapBy_cols[0],
                    "allele": allele,
                    "count": count_val[0],
                    "k": k_val
                })
            if allele_in_2:
                records.append({
                    "subject": subject,
                    "gene": gene,
                    "hap_col": hapBy_cols[1],
                    "allele": allele,
                    "count": count_val[1],
                    "k": k_val
                })
        # Also include "Unk" / "Del" that are in hap cols but NOT in alleles
        for hap_col, hap_alleles in [(hapBy_cols[0], col1_alleles), (hapBy_cols[1], col2_alleles)]:
            for hap_allele in hap_alleles:
                if hap_allele not in alleles:
                    records.append({
                        "subject": subject,
                        "gene": gene,
                        "hap_col": hap_col,
                        "allele": hap_allele,
                        "count": 0,
                        "k": 0
                    })
    return pd.DataFrame(records)


def generate_haplotype_plot(sample_name,hap_table,chain = "IGH", remove_chain = True,genes_order="position", pseudo_genes=False, ORF_genes=False, html=False, n_line=4, file= None):
    print(sample_name)
    if chain not in ["IGH", "IGK", "IGL", "TRB", "TRA"]:
        raise ValueError("Invalid chain input")

    if not hap_table["gene"].str.startswith(chain).all():
        raise ValueError("The chain input does not match the genes in the dataset")

    if genes_order is None:
        if GENE_loc and PSEUDO:
            genes_order = [gene for gene in GENE_loc.get(chain, []) if gene not in PSEUDO.get(chain, [])]
        else:
            genes_order = sorted(hap_table["gene"].unique())
    else:
        if not all(gene.startswith(chain) for gene in genes_order):
            raise ValueError("The chain input does not match the genes_order")

    id_GENE_col = list(hap_table.columns).index("gene")
    hapBy_col_id = [id_GENE_col + 1, id_GENE_col + 2]
    hapBy_cols = hap_table.columns[hapBy_col_id]
    hapBy_alleles = [col.replace("_", "*") for col in hapBy_cols]
    samples = hap_table["subject"].unique()
    
    # Ensuring all samples have all genes, filling missing values with "Unk"
    all_subjects = hap_table['subject'].unique()
    all_genes = hap_table['gene'].unique()
    index = pd.MultiIndex.from_product([all_subjects, all_genes], names=["subject", "gene"])
    index = pd.MultiIndex.from_product([all_subjects, all_genes], names=["subject", "gene"])
    complete_df = pd.DataFrame(index=index).reset_index()
    hap_table = pd.merge(complete_df, hap_table, on=["subject", "gene"], how="left")
    hap_table[hapBy_cols] = hap_table[hapBy_cols].replace({"": "Unk", None: "Unk"}).fillna("Unk")
    hap_table["alleles"] = hap_table["alleles"].replace({"": "Unk", None: "Unk"}).fillna("Unk")
    hap_table[["k1", "k2", "k3", "k4"]] = hap_table[["k1", "k2", "k3", "k4"]].fillna(0)
    # Sorting the data
    hap_table["order"] = hap_table["gene"].map({g: i for i, g in enumerate(genes_order)})
    hap_table = hap_table.dropna(subset=["order"]).sort_values("order").reset_index(drop=True)
    
    if remove_chain:
        genes_order = [gene.replace(chain, "") for gene in genes_order]
        hap_table["gene"] = hap_table["gene"].str.replace(chain, "", regex=False)
        hap_table["gene"] = pd.Categorical(hap_table["gene"], categories=[g.replace(chain, "") for g in genes_order], ordered=True)
    else:
        hap_table["gene"] = pd.Categorical(hap_table["gene"], categories=genes_order, ordered=True)
    
    # Fix NA in k columns
    hap_table[["k1", "k2", "k3", "k4"]] = hap_table[["k1", "k2", "k3", "k4"]].fillna(0)
    hap_table[["counts1", "counts2", "counts3", "counts4"]] = hap_table[["counts1", "counts2", "counts3", "counts4"]].fillna("0,0")
    
    # Melt haplotype columns to one
    hap_table = explode_hap_table(hap_table,hapBy_cols)
    unique_genes_sorted = sorted(hap_table["gene"].unique(), key=lambda x: list(genes_order).index(x))
    gene_loc = {gene: i+1 for i, gene in enumerate(unique_genes_sorted)}
    hap_table["GENE_LOC"] = hap_table["gene"].map(gene_loc)
    hap_table.rename(columns={"allele": "ALLELES","hap_col":"hapBy"}, inplace=True)
    # Copy alleles for grouping
    hap_table["ALLELES_G"] = hap_table["ALLELES"]
    hap_table["text"] = ""
    hap_table["text_bottom"] = hap_table["ALLELES"]
    # Handle numeric-range alleles (NRA)
    id_nra = hap_table["ALLELES"].str.contains(r"i?\d{2}_i?\d{2}", regex=True)
    nra = False
    if any(id_nra):
        hap_table.loc[id_nra, "ALLELES"] = "NRA"
        nra = True
    
    # Generate palette for alleles
    allele_palette_data = allele_palette(hap_table["ALLELES"].tolist())
    # Convert alleles to categorical with defined order
    allele_cats = pd.Categorical(hap_table["ALLELES"], categories=list(allele_palette_data["AlleleCol"].keys()))
    hap_table["ALLELES"] = allele_cats
    # Get unique samples and genes
    samples = hap_table["subject"].unique()
    samples_n = len(samples)
    genes = hap_table["gene"].unique()
    genes_n = len(genes)
    # Sort by subject and gene location
    hap_table = hap_table.sort_values(["subject", "GENE_LOC"])
    # Calculate line height
    hap_table["n"] = hap_table.groupby(["subject", "gene","hapBy"])["gene"].transform("count")
    hap_table["line"] = 12 / hap_table["n"]
    
    # Create allele code mapping
    clean_allele_list = []
    #Extract and clean alleles
    for allele in allele_palette_data["AlleleCol"].keys():
        code = allele.split("_")[0] if "_" in allele else allele
        clean_allele_list.append(code)
    
    # Keep only numeric values
    clean_allele_list = [item for item in clean_allele_list if item.isdigit()]
    # Assign a unique numeric code from 1 to N
    unique_clean_alleles = sorted(set(clean_allele_list))  # Sort for consistent ordering
    clean_allele_to_code = {allele: str(i + 1) for i, allele in enumerate(unique_clean_alleles)}
    
    non_numeric_keys = [key for key, value in allele_palette_data["AlleleCol"].items() if (not str(key).isdigit() and "_" not in key)]
    if non_numeric_keys:
        numeric_values = [int(value) for value in clean_allele_to_code.values() if value]
        last_numeric = max(numeric_values, default=0) + 1  # Default to 1 if no numbers exist
        for i, key in enumerate(non_numeric_keys):
            clean_allele_to_code[key] = str(last_numeric + i)
    
    # Map cleaned alleles to their unique numeric code  
    allele_code = {}
    for allele in allele_palette_data["AlleleCol"].keys():
        code = allele.split("_")[0] if "_" in allele else allele
        allele_code[allele] = clean_allele_to_code.get(code, None)  # Use None instead of "UNK"
    
    # Apply allele code mapping efficiently
    def get_allele_code(x):
        key = x.split("_")[0] if "_" in x else x
        return allele_code.get(key, x)
    
    hap_table["A_CODE"] = hap_table["ALLELES"].apply(get_allele_code)
    
    # Handle special NRA cases
    nra_pattern = hap_table["ALLELES"].astype(str).str.contains(r"i?\d{2}_i?\d{2}", regex=True)
    if "NRA" in allele_code:
        hap_table.loc[nra_pattern, "A_CODE"] = allele_code["NRA"]
    # Sort by subject, gene location, and allele code
    hap_table = hap_table.sort_values(["subject", "GENE_LOC", "A_CODE"])
    # Assign ID within each subject-gene group
    hap_table["id"] = hap_table.groupby(["subject", "gene","hapBy"]).cumcount() + 1
    
    hap_1 = hap_table[hap_table["hapBy"] == hapBy_cols[0]]
    hap_2 = hap_table[hap_table["hapBy"] == hapBy_cols[1]]
    
    hap_table_f = []
    for subject, gene, gene_loc, alleles_g, a_code, text_bottom, n_lines, id_, hapBy,n, k,count in zip(
        hap_table["subject"], 
        hap_table["gene"], 
        hap_table["GENE_LOC"], 
        hap_table["ALLELES_G"], 
        hap_table["A_CODE"], 
        hap_table["text_bottom"], 
        hap_table["line"],
        hap_table["id"],
        hap_table["hapBy"],
        hap_table["n"],
        hap_table["k"],
        hap_table["count"]
    ):
        n_lines = int(n_lines)
        n = int(n)
        for i in range(1, int(n_lines) + 1):
            row = {
                "subject": subject,
                "gene": gene,
                "GENE_LOC": gene_loc,
                "ALLELES_G": alleles_g,
                "A_CODE": a_code,
                "text_bottom": text_bottom,
                "n_line": i,
                "id": id_,
                "hapBy": hapBy,
                "k": k,
                "count": count
            }
            hap_table_f.append(row)
            if n == 5 and ((id_ == 1 and i == 1) or (id_ == 5 and i == 1)):
                hap_table_f.append(row.copy())
    
    hap_table_f = pd.DataFrame(hap_table_f)
    # Deviding the data into two dataframes
    hap_table_1 = hap_table_f[hap_table_f["hapBy"] == hapBy_cols[0]]
    hap_table_2 = hap_table_f[hap_table_f["hapBy"] == hapBy_cols[1]]
    
    hap_table_1 = hap_table_1.reset_index(drop=True)
    m_hap_1 = np.full((genes_n, 12 * len(samples)), "", dtype=object)
    for i in range(genes_n):
        for j in range(12 * len(samples)):
            a_code = str(hap_table_1.loc[(i * 12 * len(samples)) + j, "A_CODE"])
            m_hap_1[i, j] = a_code
    
    hap_table_2 = hap_table_2.reset_index(drop=True)
    m_hap_2 = np.full((genes_n, 12 * len(samples)), "", dtype=object)
    for i in range(genes_n):
        for j in range(12 * len(samples)):
            a_code = str(hap_table_2.loc[(i * 12 * len(samples)) + j, "A_CODE"])
            m_hap_2[i, j] = a_code
    
    # Combine alleles and NRA text_bottom entries into a single legend
    bottom_annot_mask = hap_table["text_bottom"].str.contains(r"i?\d+_i?\d+", regex=True)
    bottom_annot = hap_table.loc[bottom_annot_mask, "text_bottom"].unique().tolist()
    # Assign color and code for each NRA bottom annotation if not already in allele_code
    nra_color = allele_palette_data["AlleleCol"].get("NRA", "#fbf7f5")
    
    for i, text in enumerate(bottom_annot):
        if text not in allele_code:
            allele_code[text] = allele_code["NRA"]
            allele_palette_data["AlleleCol"][text] = nra_color
    
    # Get full list of alleles including NRA pseudo-alleles
    alleles = list(allele_palette_data["AlleleCol"].keys())
    n_alleles = len(alleles)
    
    # Map alleles to numeric code values (for color block rendering)
    allele_code_map = {}
    for allele in alleles:
        if allele in allele_code and allele_code[allele] and allele_code[allele].isdigit():
            allele_code_map[allele] = int(allele_code[allele])  # +1 to avoid 0
    
    # Sort alleles by code value for consistent order
    sorted_alleles = sorted(allele_code_map.items(), key=lambda x: x[1])
    # Layout settings
    longest_allele = max(len(a) for a in alleles) * 3
    legend_width = genes_n * 12        # Total available width
    num_columns = 1
    # Recalculate layout based on this dynamic number of columns
    items_per_column = n_alleles
    total_rows = items_per_column
    # Create new legend matrix
    m2 = np.zeros((total_rows, legend_width))
    for i, (allele, code_val) in enumerate(sorted_alleles):
        column = i // items_per_column
        row = i % items_per_column
        col_start = column * (legend_width // num_columns)
        col_width = longest_allele  # Width of color block
        m2[row, col_start:legend_width] = code_val
    
    longest_allele = max(len(allele) for allele in allele_palette_data["AlleleCol"].keys()) * 3 + 40
    legend_alleles = {allele: color for allele, color in allele_palette_data["AlleleCol"].items() if "_" not in allele}
    # Convert m to numeric matrix using allele_to_index
    m_hap_1_numeric = m_hap_1.astype(float)
    m_hap_2_numeric = m_hap_2.astype(float)
    m2_numeric = m2.astype(float)
    samples_ordered = hap_table_f["subject"].unique()
    
    colors_k = ['#f7fbff', '#deebf7', '#c6dbef', '#9ecae1', '#6baed6','#4292c6', '#2171b5', '#08519c', '#08306b']
    bins = [0, 1, 2, 3, 4, 5, 10, 20, 50, np.inf]
    labels_k = ["0-1", "1-2", "2-3", "3-4", "4-5", "5-10", "10-20", "20-50", "50+"]
    cmap_k = ListedColormap(colors_k)
    bin_lower_bounds = bins[:-1]  # lower bound of each bin
    
    m_hap_1_k = np.full((genes_n, 12 * len(samples)), "", dtype=object)
    # Ensure numeric
    hap_table_1["k"] = pd.to_numeric(hap_table_1["k"], errors="coerce")
    # Reset index to ensure sequential access
    hap_table_1 = hap_table_1.reset_index(drop=True)
    # Fill the matrix
    for i in range(genes_n):
        for j in range(12 * len(samples)):
            row_index = i * 12 * len(samples) + j
            k_val = hap_table_1.loc[row_index, "k"]
            if pd.isna(k_val):
                m_hap_1_k[i, j] = ""
            else:
                bin_index = np.digitize(k_val, bins) - 1  # 0-based index
                m_hap_1_k[i, j] = bin_lower_bounds[bin_index]
    
    m_hap_2_k = np.full((genes_n, 12 * len(samples)), "", dtype=object)
    # Ensure numeric
    hap_table_2["k"] = pd.to_numeric(hap_table_2["k"], errors="coerce")
    # Reset index to ensure sequential access
    hap_table_2 = hap_table_2.reset_index(drop=True)
    # Fill the matrix
    for i in range(genes_n):
        for j in range(12 * len(samples)):
            row_index = i * 12 * len(samples) + j
            k_val = hap_table_2.loc[row_index, "k"]
            if pd.isna(k_val):
                m_hap_2_k[i, j] = ""
            else:
                bin_index = np.digitize(k_val, bins) - 1  # 0-based index
                m_hap_2_k[i, j] = bin_lower_bounds[bin_index]
    
    m_hap_1_k_numeric = m_hap_1_k.astype(float)
    m_hap_2_k_numeric = m_hap_2_k.astype(float)
    # Handle novel alleles
    novel_count = 1
    nra_count = 1
    novel_map = {}
    nra_map = {}
    hap_table = hap_table.sort_values(by=["GENE_LOC", "subject"]).reset_index(drop=True)
    for idx, row in hap_table.iterrows():
        allele_g = str(row["ALLELES_G"])
        if "_" in allele_g:
            if re.match(r"i?\d+_i?\d+", allele_g):  # NRA-style
                if allele_g not in nra_map:
                    nra_map[allele_g] = nra_count
                    nra_count += 1
                hap_table.at[idx, "text_bottom"] = f"*{nra_map[allele_g]}"
            else:  # Novel allele
                if allele_g not in novel_map:
                    novel_map[allele_g] = novel_count
                    novel_count += 1
                hap_table.at[idx, "text_bottom"] = f"^{novel_map[allele_g]}"
    
    allele_to_label = hap_table.dropna(subset=["text_bottom"])[["ALLELES_G", "text_bottom"]] \
        .drop_duplicates() \
        .set_index("ALLELES_G")["text_bottom"].to_dict()
    
    # Ensure 'count' is numeric
    hap_table_f['count'] = pd.to_numeric(hap_table_f['count'], errors='coerce').fillna(0)
    # Group by gene and A_CODE
    gene_allele_count_df = (
        hap_table_f.groupby(['gene', 'A_CODE',"hapBy"])['count']
        .sum()
        .reset_index()
    )
    
    gene_allele_count_df = hap_table_f
    # Map A_CODE to color
    gene_allele_count_df['color'] = gene_allele_count_df['ALLELES_G'].map(allele_palette_data["AlleleCol"])
    # Add log10 count for plotting
    gene_allele_count_df['log_count'] = np.log10(gene_allele_count_df['count'] + 1)
    gene_allele_count_df.loc[
        gene_allele_count_df['hapBy'] == hapBy_cols[0],
        'log_count'
    ] *= -1
    unique_genes = hap_table["gene"].unique()
    
    if html:
        # Create subplots layout
        fig = make_subplots(
            rows=1, cols=8,
            column_widths=[0.4,0.015,0.2,0.2,0.015,0.2,0.2,0.3],
            specs=[[{"type": "bar"},None,{"type": "heatmap"},{"type": "heatmap"},None,{"type": "heatmap"}, {"type": "heatmap"},{"type": "table"}]],
            horizontal_spacing=0.01,
            subplot_titles=["", hapBy_cols[0], hapBy_cols[1], hapBy_cols[0], hapBy_cols[1],""]  # Add the titles in correct col positions
        )
        # Expanded list of gene labels so each bar has a unique y-value
        filtered_df = gene_allele_count_df.dropna(subset=['log_count', 'color'])
        x_vals = filtered_df['log_count'].tolist()
        # Generate a unique y-axis label for each bar to force stacking
        y_vals = [
            f"{gene} ({i})" for i, gene in enumerate(filtered_df['gene'])
        ]
        colors = filtered_df['color'].tolist()
        hover_texts = [
            f"Individual: {subject}<br />Gene: {gene}<br />Allele: {ALLELES_G}<br />Kdiff: {k}<br>Count: {count:.0f}<br>Log10 count: {log_count:.0f}"
            for subject,gene, ALLELES_G, k, count, log_count in zip(filtered_df['subject'], filtered_df['gene'],filtered_df['ALLELES_G'], filtered_df['k'], filtered_df['count'], filtered_df['log_count'])
        ]
        # Add single bar trace with unique y values
        fig.add_trace(go.Bar(
            x=x_vals,
            y=y_vals,
            orientation='h',
            marker=dict(
                color=colors,
                line=dict(color='black', width=0)
            ),
            showlegend=False,
            hovertext=hover_texts,
            hoverinfo="text"
        ), row=1, col=1)
        
        # Unique genes in order
        y_labels = hap_table['gene'].unique().tolist()
        # One label every 24 positions, centered at position +12
        tickvals = [i * 24 + 12 for i in range(len(y_labels))]
        ticktext = y_labels
        # Apply to plot
        fig.update_yaxes(
            row=1, col=1,
            tickvals=tickvals,
            ticktext=ticktext,
            #autorange='reversed',
            title="Gene"
        )
        fig.update_xaxes(
            row=1, col=1,
            title='log₁₀(Count+1)',
            zeroline=True,
            zerolinewidth=1,
            zerolinecolor='black'
        )
        # Main Allele Heatmap with Discrete Color Mapping
        allele_labels = list(legend_alleles.keys())
        allele_colors = list(legend_alleles.values())
        n = len(allele_colors)
        discrete_colorscale = []
        for i in range(n):
            start = i / n
            end = (i + 1) / n
            color = allele_colors[i]
            discrete_colorscale.append([start, color])
            discrete_colorscale.append([end, color])
        
        tickvals = [(i + 0.5)/n * (n+1) for i in range(n)]
        ticktext = allele_labels
        # Prepare hover text
        hover_texts = [
            f"Individual: {subject}<br />Gene: {gene}<br />Allele: {ALLELES_G}<br />Kdiff: {k}<br>Count: {count}"
            for subject,gene, ALLELES_G, k, count in zip(hap_table_1['subject'], hap_table_1['gene'],hap_table_1['ALLELES_G'], hap_table_1['k'], hap_table_1['count'])
        ]
        hover_texts_2d = np.array(hover_texts).reshape(m_hap_1_numeric.shape)
        fig.add_trace(go.Heatmap(
            z=m_hap_1_numeric,
            x=[],
            y=genes,
            hovertext=hover_texts_2d,
            hoverinfo="text",
            zmin=0,
            zmax=n+1,
            colorscale=discrete_colorscale,
            showscale=False,
            colorbar=dict(
                title="",
                tickmode="array",
                len=0.8,
                tickvals=tickvals,
                ticktext=ticktext
            )
        ), row=1, col=3)
        
        # Create subject and gene index mappings
        subject_to_index = {subject: idx for idx, subject in enumerate(samples_ordered)}
        gene_to_index = {gene: idx for idx, gene in enumerate(unique_genes)}
        # Add annotations for nra_text (e.g., 'i12_i3')
        nra_text = hap_1[hap_1["ALLELES_G"].str.match(r"i?\d+_i?\d+")]
        for _, row in nra_text.iterrows():
            if row["subject"] in subject_to_index and row["gene"] in gene_to_index:
                y = gene_to_index[row["gene"]]
                x = (row["id"] - 0.5) * (12 / row["n"])
                label = allele_to_label.get(row["ALLELES_G"], "")
                fig.add_annotation(
                    x=x, y=row["gene"], text=label,
                    showarrow=False, font=dict(size=8),
                    xanchor="center", yanchor="middle",
                    row=1, col=3
                )
        
        # Add annotations for novel alleles (e.g., ones with "_")
        novel_text = hap_1[hap_1["ALLELES_G"].str.contains("_")]
        for _, row in novel_text.iterrows():
            if row["subject"] in subject_to_index and row["gene"] in gene_to_index:
                y = gene_to_index[row["gene"]]
                x = (row["id"] - 0.5) * (12 / row["n"])
                label = allele_to_label.get(row["ALLELES_G"], "")
                fig.add_annotation(
                    x=x, y=row["gene"], text=label,
                    showarrow=False, font=dict(size=8),
                    xanchor="center", yanchor="middle",
                    row=1, col=3
                )
        
        # hap_2 heatmap
        allele_labels = list(legend_alleles.keys())
        allele_colors = list(legend_alleles.values())
        n = len(allele_colors)
        discrete_colorscale = []
        for i in range(n):
            start = i / n
            end = (i + 1) / n
            color = allele_colors[i]
            discrete_colorscale.append([start, color])
            discrete_colorscale.append([end, color])
        
        tickvals = [(i + 0.5)/n * (n+1) for i in range(n)]
        ticktext = allele_labels
        # Prepare hover text
        hover_texts = [
            f"Individual: {subject}<br />Gene: {gene}<br />Allele: {ALLELES_G}<br />Kdiff: {k}<br>Count: {count}"
            for subject,gene, ALLELES_G, k, count in zip(hap_table_2['subject'], hap_table_2['gene'],hap_table_2['ALLELES_G'], hap_table_2['k'], hap_table_2['count'])
        ]
        hover_texts_2d = np.array(hover_texts).reshape(m_hap_2_numeric.shape)
        fig.add_trace(go.Heatmap(
            z=m_hap_2_numeric,
            x=[],
            y=genes,
            hovertext=hover_texts_2d,
            hoverinfo="text",
            zmin=0,
            zmax=n+1,
            colorscale=discrete_colorscale,
            showscale=False,
            colorbar=dict(
                title="",
                tickmode="array",
                len=0.8,
                tickvals=tickvals,
                ticktext=ticktext
            )
        ), row=1, col=4)
        
        # Create subject and gene index mappings
        subject_to_index = {subject: idx for idx, subject in enumerate(samples_ordered)}
        gene_to_index = {gene: idx for idx, gene in enumerate(unique_genes)}
        # Add annotations for nra_text (e.g., 'i12_i3')
        nra_text = hap_2[hap_2["ALLELES_G"].str.match(r"i?\d+_i?\d+")]
        for _, row in nra_text.iterrows():
            if row["subject"] in subject_to_index and row["gene"] in gene_to_index:
                y = gene_to_index[row["gene"]]
                x = (row["id"] - 0.5) * (12 / row["n"])
                label = allele_to_label.get(row["ALLELES_G"], "")
                fig.add_annotation(
                    x=x, y=row["gene"], text=label,
                    showarrow=False, font=dict(size=8),
                    xanchor="center", yanchor="middle",
                    row=1, col=4
                )

        # Add annotations for novel alleles (e.g., ones with "_")
        novel_text = hap_2[hap_2["ALLELES_G"].str.contains("_")]
        for _, row in novel_text.iterrows():
            if row["subject"] in subject_to_index and row["gene"] in gene_to_index:
                y = gene_to_index[row["gene"]]
                x = (row["id"] - 0.5) * (12 / row["n"])
                label = allele_to_label.get(row["ALLELES_G"], "")
                fig.add_annotation(
                    x=x, y=row["gene"], text=label,
                    showarrow=False, font=dict(size=8),
                    xanchor="center", yanchor="middle",
                    row=1, col=4
                )
        
        # K_diff Heatmap
        cmap_k = {i: c for i, c in enumerate(colors_k)}
        k = len(colors_k)
        k_discrete_colorscale = []
        for i in range(k):
            start = i / k
            end = (i + 1) / k
            color = colors_k[i]
            k_discrete_colorscale.append([start, color])
            k_discrete_colorscale.append([end, color])
        
        k_tickvals = [(i + 0.5)/k * (k+1) for i in range(k)]
        hover_texts = [
            f"Individual: {subject}<br />Gene: {gene}<br />Allele: {ALLELES_G}<br />Kdiff: {k}<br>Count: {count}"
            for subject,gene, ALLELES_G, k, count in zip(hap_table_1['subject'], hap_table_1['gene'],hap_table_1['ALLELES_G'], hap_table_1['k'], hap_table_1['count'])
        ]
        hover_texts_2d = np.array(hover_texts).reshape(m_hap_1_numeric.shape)
        fig.add_trace(go.Heatmap(
            z=m_hap_1_k_numeric,
            x=[],
            y=genes,
            colorscale=k_discrete_colorscale,
            zmin=0,
            zmax=k+1,
            hovertext=hover_texts_2d,
            hoverinfo="text",
            showscale=False,
            colorbar=dict(
                title="K_diff",
                tickmode="array",
                tickvals=k_tickvals,
                ticktext=labels_k,
                len=0.8
            )
        ), row=1, col=6)
        
        cmap_k = {i: c for i, c in enumerate(colors_k)}
        k = len(colors_k)
        k_discrete_colorscale = []
        for i in range(k):
            start = i / k
            end = (i + 1) / k
            color = colors_k[i]
            k_discrete_colorscale.append([start, color])
            k_discrete_colorscale.append([end, color])
        
        k_tickvals = [(i + 0.5)/k * (k+1) for i in range(k)]
        hover_texts = [
            f"Individual: {subject}<br />Gene: {gene}<br />Allele: {ALLELES_G}<br />Kdiff: {k}<br>Count: {count}"
            for subject,gene, ALLELES_G, k, count in zip(hap_table_2['subject'], hap_table_2['gene'],hap_table_2['ALLELES_G'], hap_table_2['k'], hap_table_2['count'])
        ]
        hover_texts_2d = np.array(hover_texts).reshape(m_hap_2_k_numeric.shape)
        fig.add_trace(go.Heatmap(
            z=m_hap_2_k_numeric,
            x=[],
            y=genes,
            colorscale=k_discrete_colorscale,
            zmin=0,
            zmax=k+1,
            hovertext=hover_texts_2d,
            hoverinfo="text",
            showscale=False,
            colorbar=dict(
                title="K_diff",
                tickmode="array",
                tickvals=k_tickvals,
                ticktext=labels_k,
                len=0.8
            )
        ), row=1, col=7)
        legend_colors = [allele_palette_data["AlleleCol"].get(k[0], "#ffffff") for k in sorted_alleles]
        def wrap_text(text, width=26):
            return '<br>'.join([text[i:i+width] for i in range(0, len(text), width)])
        fig.add_trace(go.Table(
            header=dict(
                values=["<b style='font-size:16px'>Alleles</b>"],
                fill_color='white',
                align='center'
            ),
            cells=dict(
                values=[[  # legend values
                    f"{allele_to_label.get(allele, '')} - {wrap_text(allele)}"
                    if allele_to_label.get(allele, '') != allele
                    else wrap_text(allele)
                    for allele, _ in sorted_alleles
                ]],
                fill_color=[legend_colors],
                align='center'
            ),
            domain=dict(x=[0.83, 1], y=[0.5, 1])
        ))
        fig.add_trace(go.Table(
            header=dict(
                values=["<b style='font-size:16px'>K_diff</b>"],
                fill_color='white',
                align='center'
            ),
            cells=dict(
                values=[labels_k],
                fill_color=[colors_k],
                align='center'
            ),
            domain=dict(x=[0.83, 1], y=[0, 0.5])
        ))
        
        fig.update_layout(
            height=1500,
            width=1500,
            title=dict(
                text=sample_name,
                x=0.5,               # Center the title horizontally
                xanchor="center"     # Anchor the title at the center
            ),
            xaxis2=dict(showticklabels=False),  # Hide x-axis of first heatmap
            yaxis2=dict(showticklabels=False),   # Hide y-axis of second heatmap
            xaxis3=dict(showticklabels=False),   # Hide x-axis of second heatmap
            yaxis3=dict(showticklabels=False),   # Hide y-axis of second heatmap
            xaxis4=dict(showticklabels=False),  # Hide x-axis of first heatmap
            yaxis4=dict(showticklabels=False),   # Hide y-axis of second heatmap
            xaxis5=dict(showticklabels=False),   # Hide x-axis of second heatmap
            yaxis5=dict(showticklabels=False)   # Hide y-axis of second heatmap
        )
        # Save to HTML
        pio.write_html(fig, file=file, auto_open=False)
    else:
        # ----- PDF GENERATION -----
        size = 20
        with PdfPages(file) as pdf:  
            fig = plt.figure(figsize=(size, size))
            plt.suptitle(sample_name, fontsize=16)
            gs = plt.GridSpec(3, 6, width_ratios=[3,2, 2, 1,1,3], height_ratios=[1, 0.1, 1])
            # --- Main Heatmap ---
            colors =list(legend_alleles.values())
            cmap = ListedColormap(colors)
            # ax7 = barplot colored by A_CODE per row in gene_allele_count_df
            ax7 = fig.add_subplot(gs[:, 0])
            # Assign a y-position to each gene
            gene_to_y = {gene: i for i, gene in enumerate(unique_genes)}
            for _, row in gene_allele_count_df.iterrows():
                gene = row['gene']
                y = gene_to_y.get(gene)
                if y is None or pd.isna(row['log_count']):
                    continue
                ax7.barh(
                    y, 
                    row['log_count'], 
                    color=row['color'], 
                    edgecolor='black',
                    height=0.8
                )
            # Add grid lines
            # Axis styling
            ax7.set_yticks(range(len(unique_genes)))
            ax7.set_yticklabels(unique_genes)
            ax7.set_xlabel('log$_{10}$(Count+1)')
            ax7.invert_yaxis()
            ax7.margins(x=0)
            ax7.set_ylim(-0.5, len(unique_genes) - 0.5)
            ax7.invert_yaxis()
            ax1 = fig.add_subplot(gs[:, 1])
            im = ax1.imshow(m_hap_1_numeric, cmap=cmap, aspect='auto', interpolation='none')
            # Vertical grid lines
            for g in range(1, genes_n):
                ax1.axhline(y=g-0.5, color='white', linestyle='-', linewidth=1)
            # X and Y labels
            ax1.set_yticks([])
            ax1.set_yticklabels([])
            ax1.set_xticks([])
            ax1.set_xticklabels([])
            ax1.set_title(hapBy_cols[0])
            nra_text = hap_1[hap_1["ALLELES_G"].str.match(r"i?\d+_i?\d+")]
            subject_to_index = {subject: idx for idx, subject in enumerate(samples_ordered)}
            gene_to_index = {gene: idx for idx, gene in enumerate(unique_genes)}
            for _, row in nra_text.iterrows():
                if row["subject"] in subject_to_index:
                    y = gene_to_index[row["gene"]]
                    x = (row["id"] - 0.5) * (12 / row["n"])
                    label = allele_to_label.get(row["ALLELES_G"], "")
                    ax1.text(x,y,label, fontsize=5, ha='center', va='center')
            novel_text = hap_1[hap_1["ALLELES_G"].str.contains("_")]
            subject_to_index = {subject: idx for idx, subject in enumerate(samples_ordered)}
            gene_to_index = {gene: idx for idx, gene in enumerate(unique_genes)}
            for _, row in novel_text.iterrows():
                if row["subject"] in subject_to_index:
                    y = gene_to_index[row["gene"]]
                    x = (row["id"] - 0.5) * (12 / row["n"])
                    label = allele_to_label.get(row["ALLELES_G"], "")
                    ax1.text(x,y,label, fontsize=5, ha='center', va='center')
            
            ax2 = fig.add_subplot(gs[:, 2])
            im = ax2.imshow(m_hap_2_numeric, cmap=cmap, aspect='auto', interpolation='none')
            # Vertical grid lines
            for g in range(1, genes_n):
                ax2.axhline(y=g-0.5, color='white', linestyle='-', linewidth=1)
            # X and Y labels
            ax2.set_yticks([])
            ax2.set_yticklabels([])
            ax2.set_xticks([])
            ax2.set_xticklabels([])
            ax2.set_title(hapBy_cols[1])
            nra_text = hap_2[hap_2["ALLELES_G"].str.match(r"i?\d+_i?\d+")]
            subject_to_index = {subject: idx for idx, subject in enumerate(samples_ordered)}
            gene_to_index = {gene: idx for idx, gene in enumerate(unique_genes)}
            for _, row in nra_text.iterrows():
                if row["subject"] in subject_to_index:
                    y = gene_to_index[row["gene"]]
                    x = (row["id"] - 0.5) * (12 / row["n"])
                    label = allele_to_label.get(row["ALLELES_G"], "")
                    ax2.text(x,y,label, fontsize=5, ha='center', va='center')
            novel_text = hap_2[hap_2["ALLELES_G"].str.contains("_")]
            subject_to_index = {subject: idx for idx, subject in enumerate(samples_ordered)}
            gene_to_index = {gene: idx for idx, gene in enumerate(unique_genes)}
            for _, row in novel_text.iterrows():
                if row["subject"] in subject_to_index:
                    y = gene_to_index[row["gene"]]
                    x = (row["id"] - 0.5) * (12 / row["n"])
                    label = allele_to_label.get(row["ALLELES_G"], "")
                    ax2.text(x,y,label, fontsize=5, ha='center', va='center') 
            # --- Legend ---
            # Modify the legend plotting code
            ax5 = fig.add_subplot(gs[0, 5]) 
            #colors = ["white"] +list(legend_alleles.values())
            colors = list(legend_alleles.values())
            cmap = ListedColormap(colors)
            ax5.imshow(m2_numeric, cmap=cmap, aspect='auto', interpolation='none')
            ax5.set_title("Alleles", fontsize=10)
            # Add horizontal grid lines
            for r in range(1, total_rows):
                ax5.axhline(y=r-0.5, color='black', linestyle='--', linewidth=0.5)
            # Add vertical grid lines to separate columns
            for c in range(1, num_columns):
                col_pos = c * (genes_n * 12 // num_columns) - 0.5
                ax5.axvline(x=col_pos, color='black', linestyle='-', linewidth=1)
            # Add text labels for each allele
            for i, (allele, code_val) in enumerate(sorted_alleles):
                column = i // items_per_column
                row = i % items_per_column
                # Calculate text position (to right of colored area)
                col_start = column * (genes_n * 12 // num_columns)
                # Add the text label
                    # Default label
                full_label = allele
                # Check if allele is in ALLELES_G with underscore
                label = allele_to_label.get(allele, "")
                if allele == label:
                    label = ""
                #full_label = f"{text_bottom} {allele}"
                full_label = f"{label} - {allele}" if label else allele
                #ax5.text(legend_width/2, row, full_label, fontsize=8, ha='center', va='center')
                ax5.text(
                    5,
                    row,
                    full_label,
                    fontsize=6,               # slightly smaller font helps
                    ha='left',                # align to left to avoid center-squeeze
                    va='center',
                    wrap=False,               # make sure wrap is off
                    fontdict={'family': 'monospace'},  # monospace font keeps it tidy
                    usetex=False              # avoid TeX interpretation of underscores
                )
            ax5.set_xticks([])
            ax5.set_yticks([])
            
            # --- K Heatmap ---
            ax3 = fig.add_subplot(gs[:, 3])
            im = ax3.imshow(m_hap_1_k_numeric, cmap=cmap_k, aspect='auto')
            # Tick labels
            #ax3.set_xticks(np.arange(len(m_hap_1_k.columns)))
            ax3.set_xticks([])
            ax3.set_xticklabels([])
            ax3.set_yticks([])
            ax3.set_yticklabels([])
            ax3.set_title(hapBy_cols[0])
            
            ax4 = fig.add_subplot(gs[:, 4])
            im = ax4.imshow(m_hap_2_k_numeric, cmap=cmap_k, aspect='auto')
            # Tick labels
            ax4.set_xticks([])
            ax4.set_xticklabels([])
            ax4.set_yticks([])
            ax4.set_yticklabels([])
            ax4.set_title(hapBy_cols[1])
            
            # ax6 = legend_k
            ax6 = fig.add_subplot(gs[2, 5])
            for i, (color, label) in enumerate(zip(colors_k, labels_k)):
                ax6.add_patch(plt.Rectangle((0, i), 1, 1, color=color))
                ax6.text(0.5, i + 0.5, label, ha='center', va='center', fontsize=8, color='black')
            ax6.set_ylim(0, len(colors_k))
            ax6.set_xlim(0, 1)
            ax6.set_xticks([])
            ax6.set_yticks([])
            ax6.set_title("K_diff", fontsize=10)
            plt.subplots_adjust(
                top=0.95,      # Move subplots closer to the top (0-1 scale)
                bottom=0.05,   # Move subplots closer to the bottom
                left=0.05,     # Move subplots closer to the left
                right=0.95,    # Move subplots closer to the right
                hspace=0     # Control vertical spacing between subplots
            )
            # Save the figure
            pdf.savefig(fig, bbox_inches='tight')  # 'tight' removes excess whitespace
            plt.close()
    
    return file

