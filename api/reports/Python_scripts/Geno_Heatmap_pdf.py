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
from collections import defaultdict

GENE_LOCATIONS = {
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
    """
    alleles = hap_alleles
    base_alleles = sorted(set(a.split('_')[0] for a in alleles))
    
    # Base color assignment
    allele_col = {a: ALLELE_PALETTE.get(a, "#000000") for a in base_alleles}

    # Add novel allele colors (based on base allele color)
    novels = [a for a in alleles if '_' in a]
    for novel in novels:
        base = novel.split('_')[0]
        allele_col[novel] = ALLELE_PALETTE.get(base, "#000000")
        
        if base not in allele_col:
            allele_col[base] = ALLELE_PALETTE.get(base, "#000000")

    # Add special alleles
    special_alleles = {"Unk": "#dedede", "Del": "#6d6d6d", "NR": "#000000", "NRA": "#fbf7f5"}
    allele_col.update({k: v for k, v in special_alleles.items() if k in hap_alleles})

    # Transparency calculation
    novel_groups = defaultdict(list)
    for allele in allele_col:
        if '_' in allele:
            novel_groups[allele.split('_')[0]].append(allele)

    transper = {}
    for allele in allele_col:
        if '_' in allele:
            base = allele.split('_')[0]
            group = novel_groups[base]
            m = group.index(allele)
            n = len(group)
            if n == 1:
                transper[allele] = 0.5
            elif n == 2:
                transper[allele] = 0.6 if m == 0 else 0.3
            elif n == 3:
                transper[allele] = 0.6 if m == 0 else 0.4 if m == 1 else 0.2
            elif n > 9:
                transper[allele] = 1 if m == 0 else 1 - m / 20
            else:
                transper[allele] = 0.85 if m == 0 else 0.85 - m / 10
        else:
            transper[allele] = 1

    # Filter output
    included = set(alleles) | set(special_alleles.keys()) if nra else set(alleles)
    final_col = {a: c for a, c in allele_col.items() if a in included}
    final_trans = {a: t for a, t in transper.items() if a in final_col}

    return {
        "transper": final_trans,
        "AlleleCol": final_col
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
    GENE_loc_tmp = GENE_LOCATIONS[chain].copy()
    
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


def generate_heatmap_pdf(geno_table,chain = "IGH", remove_chain = True,gene_sort="position", pseudo_genes=False, ORF_genes=False, html=False, file= None):
    if file is None:
        file = os.path.join(tempfile.gettempdir(), "genotype_heatmap.html")
    
    # Subset and rename columns
    geno_db = geno_table[["subject", "gene", "GENOTYPED_ALLELES", "k_diff"]].copy()
    if not all(gene.startswith(chain) for gene in geno_db["gene"]):
        raise ValueError(f"The chain input {chain} does not match the genes")
    
    geno_db.rename(columns={"GENOTYPED_ALLELES": "ALLELES", "k_diff": "K"}, inplace=True)
    # Replace "Deletion" with "Del"
    geno_db["ALLELES"] = geno_db["ALLELES"].str.replace("Deletion", "Del")
    # Create complete dataframe with all subject-gene combinations
    all_subjects = geno_db["subject"].unique()
    all_genes = geno_db["gene"].unique()
    index = pd.MultiIndex.from_product([all_subjects, all_genes], names=["subject", "gene"])
    complete_df = pd.DataFrame(index=index).reset_index()
    geno_db = pd.merge(complete_df, geno_db, on=["subject", "gene"], how="left")
    # Fill missing values
    geno_db.loc[geno_db["ALLELES"].isna(), ["ALLELES", "K"]] = ["Unk", np.nan]
    geno_db.loc[geno_db["ALLELES"].str.contains("Del"), "K"] = np.nan
    # Split comma-separated alleles into separate rows
    geno_db = geno_db.assign(ALLELES=geno_db['ALLELES'].str.split(',')).explode('ALLELES')
    geno_db = sort_df_by_gene(geno_db, chain=chain, method=gene_sort, remove_chain=remove_chain, 
                                geno=True, pseudo_remove=pseudo_genes, orf_remove=ORF_genes)
    geno_db = geno_db.dropna(subset=["gene"])
    # Get unique samples and genes
    samples = geno_table["subject"].unique()
    samples_n = len(samples)
    genes = geno_db["gene"].unique()
    genes_n = len(genes)
    # Create gene location mapping
    gene_loc = {gene: i+1 for i, gene in enumerate(genes)}
    geno_db["GENE_LOC"] = geno_db["gene"].map(gene_loc)
    # Count alleles per subject-gene combination
    geno_db["n"] = geno_db.groupby(["subject", "gene"])["gene"].transform("count")
    # Copy alleles for grouping
    geno_db["ALLELES_G"] = geno_db["ALLELES"]
    geno_db["text"] = ""
    geno_db["text_bottom"] = geno_db["ALLELES"]
    # Handle numeric-range alleles (NRA)
    id_nra = geno_db["ALLELES"].str.contains(r"i?\d{2}_i?\d{2}", regex=True)
    if any(id_nra):
        geno_db.loc[id_nra, "ALLELES"] = "NRA"
    
    # Generate palette for alleles
    allele_palette_data = allele_palette(geno_db["ALLELES"].tolist())
    # Convert alleles to categorical with defined order
    allele_cats = pd.Categorical(geno_db["ALLELES"], categories=list(allele_palette_data["AlleleCol"].keys()))
    geno_db["ALLELES"] = allele_cats
    
    # Sort by subject and gene location
    geno_db = geno_db.sort_values(["subject", "GENE_LOC"])
    # Calculate line height
    geno_db["line"] = 12 / geno_db["n"]
    
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
    
    geno_db["A_CODE"] = geno_db["ALLELES"].apply(get_allele_code)

    # Handle special NRA cases
    nra_pattern = geno_db["ALLELES"].astype(str).str.contains(r"i?\d{2}_i?\d{2}", regex=True)
    if "NRA" in allele_code:
        geno_db.loc[nra_pattern, "A_CODE"] = allele_code["NRA"]
    # Sort by subject, gene location, and allele code
    geno_db = geno_db.sort_values(["subject", "GENE_LOC", "A_CODE"])
    # Assign ID within each subject-gene group
    geno_db["id"] = geno_db.groupby(["subject", "gene"]).cumcount() + 1
    geno_db_f = []
    for subject, gene, gene_loc, alleles_g, a_code, text_bottom, k, n_lines, id_, n in zip(
        geno_db["subject"], 
        geno_db["gene"], 
        geno_db["GENE_LOC"], 
        geno_db["ALLELES_G"], 
        geno_db["A_CODE"], 
        geno_db["text_bottom"], 
        geno_db["K"], 
        geno_db["line"],
        geno_db["id"],
        geno_db["n"]
    ):
        n_lines = int(n_lines)
        n = int(n)
        for i in range(1, n_lines + 1):
            row = {
                "subject": subject,
                "gene": gene,
                "GENE_LOC": gene_loc,
                "ALLELES_G": alleles_g,
                "A_CODE": a_code,
                "text_bottom": text_bottom,
                "K": k,
                "n_line": i,
                "id": id_
            }
            geno_db_f.append(row)
            if n == 5 and ((id_ == 1 and i == 1) or (id_ == 5 and i == 1)):
                geno_db_f.append(row.copy())
    
    geno_db_f = pd.DataFrame(geno_db_f)
    
    n_subjects = len(samples)
    n_cols = 12 * genes_n
    # Extract the A_CODE column directly, convert to string array and reshape
    a_codes = geno_db_f["A_CODE"].astype(str).to_numpy().reshape(n_subjects, n_cols)
    # Create the matrix directly
    m = np.array(a_codes, dtype=object)
    
    # Combine alleles and NRA text_bottom entries into a single legend
    bottom_annot_mask = geno_db["text_bottom"].str.contains(r"i?\d{2}_i?\d{2}", regex=True)
    bottom_annot = geno_db.loc[bottom_annot_mask, "text_bottom"].unique().tolist()
    
    # Assign color and code for each NRA bottom annotation if not already in allele_code
    nra_color = allele_palette_data["AlleleCol"].get("NRA", "#fbf7f5")
    for i, text in enumerate(bottom_annot):
        if text not in allele_code:
            allele_code[text] = allele_code["NRA"]
            allele_palette_data["AlleleCol"][text] = nra_color
    
    # Get full list of alleles including NRA pseudo-alleles
    alleles = list(allele_palette_data["AlleleCol"].keys())
    n_alleles = len(alleles)
    sorted_alleles = sorted(allele_code.items(), key=lambda x: x[1])
    longest_allele = int(max(len(a) for a in alleles)) + 3
    col_width_pixels = longest_allele  # Estimated width needed per allele
    heatmap_width = int(genes_n * 12 *0.5)        # Total available width
    num_columns = max(1, heatmap_width // col_width_pixels)
    # Recalculate layout based on this dynamic number of columns
    items_per_column = (n_alleles + (num_columns-1)) // num_columns
    # Create new legend matrix
    legend_width = int(genes_n * 12)
    total_rows = (n_alleles + (num_columns - 1)) // num_columns
    col_width = 5  # Width of color block
    column_width = legend_width // num_columns
    # Initialize output matrix
    m2 = np.zeros((int(total_rows), int(legend_width)))
    # Create arrays of indices
    allele_indices = np.arange(len(sorted_alleles))
    rows = allele_indices // num_columns
    cols = allele_indices % num_columns
    col_starts = cols * column_width
    # Fill matrix
    for idx, (row, col_start) in enumerate(zip(rows, col_starts)):
        code_val = sorted_alleles[idx][1]
        m2[row, col_start:col_start + col_width] = code_val
    
    # Set the height and width of plot
    height = samples_n * 0.5 + items_per_column * 0.1
    width = genes_n *12 * 0.05

    # ----- SETUP -----
    legend_alleles = {allele: color for allele, color in allele_palette_data["AlleleCol"].items() if "_" not in allele}
    print(legend_alleles)
    m_numeric = m.astype(float)
    m2_numeric = m2.astype(float)
    samples_ordered = geno_db_f["subject"].unique()
    # ----- PDF GENERATION -----
    with PdfPages(file) as pdf:  
        fig = plt.figure(figsize=(width, height))
        x = True
        gs = plt.GridSpec(2, 1, height_ratios=[3, 1], hspace=5/round(samples_n, -1))
        # --- Main Heatmap ---
        ax1 = plt.subplot(gs[0])
        colors =list(legend_alleles.values())
        cmap = ListedColormap(colors)
        ax1.imshow(m_numeric, cmap=cmap, aspect='auto', interpolation='none')
        # Vertical grid lines
        for g in range(1, genes_n):
            ax1.axvline(x=g*12-0.5, color='white', linestyle='-', linewidth=1)
        # X and Y labels
        ax1.set_xticks([(g*12+6-0.5) for g in range(genes_n)])
        ax1.set_xticklabels(list({gene: i+1 for i, gene in enumerate(genes)}), rotation=90)
        ax1.set_yticks(range(len(samples)))
        ax1.set_yticklabels(samples_ordered, fontsize=8)
        subject_to_index = {subject: idx for idx, subject in enumerate(samples_ordered)}
        # Helper function to add symbols
        def plot_symbols(df, symbol):
            df = df[df["subject"].isin(subject_to_index)]
            if df.empty:
                return
            y_vals = df["subject"].map(subject_to_index).to_numpy()
            x_vals = ((df["GENE_LOC"] - 1) * 12 + (df["id"] - 0.5) * (12 / df["n"])).to_numpy()
            texts = [symbol] * len(df)
            for x, y, text in zip(x_vals, y_vals, texts):
                ax1.text(x, y, text, fontsize=5, ha='center', va='center')
        
        # Plot NRA and novel markers
        plot_symbols(geno_db[geno_db["text_bottom"].str.match(r"i?\d{2}_i?\d{2}")], "*")
        plot_symbols(geno_db[geno_db["ALLELES"].str.contains("_")], "^")
        
        # --- Legend ---
        # Modify the legend plotting code
        ax2 = plt.subplot(gs[1])
        colors = ["#FFFFFF"] + list(legend_alleles.values()) if "#FFFFFF" not in legend_alleles.values() else list(legend_alleles.values())
        cmap = ListedColormap(colors)
        ax2.imshow(m2_numeric, cmap=cmap, aspect='auto', interpolation='none')
        # Add horizontal grid lines
        for r in range(1, total_rows):
            ax2.axhline(y=r-0.5, color='black', linestyle='--', linewidth=0.5)
        # Add vertical grid lines to separate columns
        for c in range(1, num_columns):
            col_pos = c * (genes_n * 12 // num_columns) - 0.5
            ax2.axvline(x=col_pos, color='black', linestyle='-', linewidth=1)
        # Add text labels for each allele
        for i, (allele, code_val) in enumerate(sorted_alleles):
            row = i // num_columns
            column = i % num_columns
            # Calculate text position (to right of colored area)
            col_start = column * (genes_n * 12 // num_columns)
            text_x = col_start + 6  # Position text after colored area
            # Add the text label
            ax2.text(text_x, row, allele, fontsize=8, ha='left', va='center')
        ax2.set_xticks([])
        ax2.set_yticks([])
        plt.subplots_adjust(
            top=0.95,      # Move subplots closer to the top (0-1 scale)
            bottom=0.05,   # Move subplots closer to the bottom
            left=0.05,     # Move subplots closer to the left
            right=0.95,    # Move subplots closer to the right
            hspace=0.05     # Control vertical spacing between subplots
        )
        # Save the figure
        pdf.savefig(fig, bbox_inches='tight')  # 'tight' removes excess whitespace
        plt.close()
        
    return plt
