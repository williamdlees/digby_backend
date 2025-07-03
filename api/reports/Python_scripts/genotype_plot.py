import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, Polygon
import re
from itertools import product
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly.io as pio

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



def generate_genotype_plot(geno_table,chain = "IGH", remove_chain = True, lk_cutoff= 1, mark_low_lk = True,gene_sort="position", pseudo_genes=False, ORF_genes=False, html=False, n_line=4, line_length=60, file= None):
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
    # Handle pseudo and ORF genes
    color_pes_orf = []
    if pseudo_genes:
        color_pes_orf.extend([g for g in PSEUDO[chain] if 'V' in g])
    if ORF_genes:
        color_pes_orf.extend(geno_db[geno_db['gene'].str.contains('OR|NL')]['gene'].unique())
    # Sort the data
    geno_db = sort_df_by_gene(geno_db, chain=chain, method=gene_sort, remove_chain=remove_chain, 
                                geno=True, pseudo_remove=pseudo_genes, orf_remove=ORF_genes)
    geno_db = geno_db.dropna(subset=["gene"])
    # Remove chain prefix if requested
    if remove_chain:
        geno_db["gene"] = geno_db["gene"].str.replace(chain, "")
    
    # Create gene location mapping
    unique_genes = geno_db["gene"].unique()
    gene_loc = {gene: i+1 for i, gene in enumerate(unique_genes)}
    geno_db["GENE_LOC"] = geno_db["gene"].map(gene_loc)
    # Count alleles per subject-gene combination
    geno_db["n"] = geno_db.groupby(["subject", "gene"])["gene"].transform("count")
    # Copy alleles for grouping
    geno_db["ALLELES_G"] = geno_db["ALLELES"]
    geno_db["text"] = ""
    geno_db["text_bottom"] = geno_db["ALLELES"]
    # Handle numeric-range alleles (NRA)
    id_nra = geno_db["ALLELES"].str.contains(r"i?\d+_i?\d+", regex=True)
    if any(id_nra):
        geno_db.loc[id_nra, "ALLELES"] = "NRA"
    
    # Generate palette for alleles
    allele_palette_data = allele_palette(geno_db["ALLELES"].tolist())
    # Convert alleles to categorical with defined order
    allele_cats = pd.Categorical(geno_db["ALLELES"], categories=list(allele_palette_data["AlleleCol"].keys()))
    geno_db["ALLELES"] = allele_cats
    # Get unique samples and genes
    samples = geno_table["subject"].unique()
    samples_n = len(samples)
    genes = geno_db["gene"].unique()
    genes_n = len(genes)
    # Sort by subject and gene location
    geno_db = geno_db.sort_values(["subject", "GENE_LOC"])
    # Calculate line height
    geno_db["line"] = 12 / geno_db["n"]
    # Create allele code mapping
    clean_allele_list = []
    # Extract and clean alleles
    for allele in allele_palette_data["AlleleCol"].keys():
        clean_allele = re.sub(r"\^[0-9]+-", "", allele)  
        code = clean_allele.split("_")[0] if "_" in clean_allele else clean_allele
        clean_allele_list.append(code)
    
    # Keep only numeric values
    clean_allele_list = [item for item in clean_allele_list if item.isdigit()]
    # Assign a unique numeric code from 1 to N
    unique_clean_alleles = sorted(set(clean_allele_list))  
    clean_allele_to_code = {allele: str(i + 1) for i, allele in enumerate(unique_clean_alleles)}
    # Map cleaned alleles to their unique numeric code
    allele_code = {}
    for allele in allele_palette_data["AlleleCol"].keys():
        clean_allele = re.sub(r"\^[0-9]+-", "", allele)
        code = clean_allele.split("_")[0] if "_" in clean_allele else clean_allele
        allele_code[clean_allele] = clean_allele_to_code.get(code, None)
    
    # Handle non-numeric codes
    non_numeric_keys = [key for key, value in allele_code.items() if not (value and value.isdigit())]
    if non_numeric_keys:
        numeric_values = [int(value) for value in allele_code.values() if value and value.isdigit()]
        last_numeric = max(numeric_values, default=0) + 1  
        for i, key in enumerate(non_numeric_keys):
            allele_code[key] = str(last_numeric + i)
    
    geno_db["A_CODE"] = geno_db["ALLELES"].apply(
        lambda x: allele_code[re.sub(r"\^[0-9]+[-]", "", str(x))] 
        if pd.notna(x) else x
    )
    
    # Handle special NRA cases
    nra_pattern = geno_db["ALLELES"].astype(str).str.contains(r"[0-9]_[0-9]")
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
    
    m = np.full((genes_n, 12 * len(samples)), "", dtype=object)
    for i in range(genes_n):
        for j in range(12 * len(samples)):
            a_code = str(geno_db_f.loc[(i * 12 * len(samples)) + j, "A_CODE"])
            m[i, j] = a_code
    
    # Combine alleles and NRA text_bottom entries into a single legend
    bottom_annot_mask = geno_db["text_bottom"].str.contains(r"i?\d+_i?\d+", regex=True)
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
    # Fill legend matrix with color codes
    for i, (allele, code_val) in enumerate(sorted_alleles):
        column = i // items_per_column
        row = i % items_per_column
        col_start = column * (legend_width // num_columns)
        col_width = longest_allele  # Width of color block
        m2[row, col_start:legend_width] = code_val
    
    longest_allele = max(len(allele) for allele in allele_palette_data["AlleleCol"].keys()) * 3 + 40

    legend_alleles = {allele: color for allele, color in allele_palette_data["AlleleCol"].items() if "_" not in allele}
    # Convert m to numeric matrix using allele_to_index
    m_numeric = m.astype(float)
    m2_numeric = m2.astype(float)
    samples_ordered = geno_db_f["subject"].unique()
    
    geno_db["K"] = pd.to_numeric(geno_db["K"], errors='coerce')
    # Step 1: Aggregate unique K per (gene, subject)
    df_k = geno_db.groupby(["gene", "subject"])["K"].unique().reset_index()
    # Step 2: Pivot to matrix form
    k_matrix_df = df_k.pivot(index="gene", columns="subject", values="K").fillna(0)
    # Step 3: Bin the values manually
    colors_k = ['#f7fbff', '#deebf7', '#c6dbef', '#9ecae1', '#6baed6','#4292c6', '#2171b5', '#08519c', '#08306b']
    bins = [0, 1, 2, 3, 4, 5, 10, 20, 50, np.inf]
    labels_k = ["0-1", "1-2", "2-3", "3-4", "4-5", "5-10", "10-20", "20-50", "50+"]
    matrix_binned = pd.DataFrame(
        np.digitize(k_matrix_df, bins) - 1,  # bin indices from 0
        index=k_matrix_df.index,
        columns=k_matrix_df.columns
    )
    cmap_k = ListedColormap(colors_k)
    # Handle novel alleles
    novel_count = 1
    nra_count = 1
    novel_map = {}
    nra_map = {}
    geno_db = geno_db.sort_values(by=["GENE_LOC", "subject"]).reset_index(drop=True)
    for idx, row in geno_db.iterrows():
        allele_g = str(row["ALLELES_G"])
        if "_" in allele_g:
            if re.match(r"i?\d+_i?\d+", allele_g):  # NRA-style
                if allele_g not in nra_map:
                    nra_map[allele_g] = nra_count
                    nra_count += 1
                geno_db.at[idx, "text_bottom"] = f"*{nra_map[allele_g]}"
            else:  # Novel allele
                if allele_g not in novel_map:
                    novel_map[allele_g] = novel_count
                    novel_count += 1
                geno_db.at[idx, "text_bottom"] = f"^{novel_map[allele_g]}"
    
    allele_to_label = geno_db.dropna(subset=["text_bottom"])[["ALLELES_G", "text_bottom"]] \
        .drop_duplicates() \
        .set_index("ALLELES_G")["text_bottom"].to_dict()
    
    if html:
        # Prepare matrix_binned
        matrix_binned_values = matrix_binned.loc[genes].values.astype(int)
        # Create subplots layout
        fig = make_subplots(
            rows=1, cols=3,
            column_widths=[0.6, 0.2,0.2],
            specs=[[{"type": "heatmap"}, {"type": "heatmap"},{"type": "table"}]],
            horizontal_spacing=0.05
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
        text_column = []
        for subject,gene,text_bottom,k in zip(
                geno_db_f['subject'],geno_db_f['gene'],geno_db_f['text_bottom'],geno_db_f['K']
        ):
            text_column.append(f"Individual: {subject}<br />Gene: {gene}<br />Allele: {text_bottom}<br />Kdiff: {k}")
        
        hover_texts = [
            f"Individual: {subject}<br />Gene: {gene}<br />Allele: <span style='font-family:monospace'>{ALLELES}</span><br />Kdiff: {k}"
            for subject, gene, ALLELES, k in zip(
                geno_db_f['subject'], geno_db_f['gene'], geno_db_f['text_bottom'], geno_db_f['K']
            )
        ]
        hover_texts_2d = np.array(hover_texts).reshape(m.shape)
        geno_db_f['text'] = text_column
        # Fill hover matrix with corresponding text
        geno_db_dict = {}
        for subject, n_line, GENE_LOC, text in zip(
            geno_db_f["subject"], geno_db_f["n_line"], geno_db_f["GENE_LOC"], geno_db_f["text"]):
            geno_db_dict[subject, n_line, GENE_LOC] = text
        
        samples_ordered = geno_db_f["subject"].unique()
        conditions_text = np.full((genes_n, 12 * len(samples)), "", dtype=object)
        for i in range(genes_n):  # gene index
            for j in range(12 * len(samples)):  # allele slot across all subjects
                subject_index = j // 12
                n_line = (j % 12) + 1  # 1-based
                gene_loc = i + 1       # assuming gene locations are 1-based and ordered
                subject = samples[subject_index]
                conditions_text[i, j] = geno_db_dict.get((subject, n_line, gene_loc), "")
        
        fig.add_trace(go.Heatmap(
            z=m_numeric,
            x=[],
            y=genes,
            text=hover_texts_2d,
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
        ), row=1, col=1)
        
        # Create subject and gene index mappings
        subject_to_index = {subject: idx for idx, subject in enumerate(samples_ordered)}
        gene_to_index = {gene: idx for idx, gene in enumerate(unique_genes)}
        # Add annotations for nra_text (e.g., 'i12_i3')
        nra_text = geno_db[geno_db["ALLELES_G"].str.match(r"i?\d+_i?\d+")]
        for _, row in nra_text.iterrows():
            if row["subject"] in subject_to_index and row["gene"] in gene_to_index:
                y = gene_to_index[row["gene"]]
                x = (row["id"] - 0.5) * (12 / row["n"])
                label = allele_to_label.get(row["ALLELES_G"], "")
                fig.add_annotation(
                    x=x, y=row["gene"], text=label,
                    showarrow=False, font=dict(size=8),
                    xanchor="center", yanchor="middle",
                    row=1, col=1
                )
        
        # Add annotations for novel alleles (e.g., ones with "_")
        novel_text = geno_db[geno_db["ALLELES_G"].str.contains("_")]
        for _, row in novel_text.iterrows():
            if row["subject"] in subject_to_index and row["gene"] in gene_to_index:
                y = gene_to_index[row["gene"]]
                x = (row["id"] - 0.5) * (12 / row["n"])
                label = allele_to_label.get(row["ALLELES_G"], "")
                fig.add_annotation(
                    x=x, y=row["gene"], text=label,
                    showarrow=False, font=dict(size=8),
                    xanchor="center", yanchor="middle",
                    row=1, col=1
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
        hover_text = k_matrix_df.applymap(lambda x: f"Individual: {subject}<br />Gene: {gene}<br />Kdiff: {x}")
        fig.add_trace(go.Heatmap(
            z=matrix_binned_values.astype(float),
            x=matrix_binned.columns.astype(str),
            y=genes,
            colorscale=k_discrete_colorscale,
            zmin=0,
            zmax=k+1,
            text=hover_text.values,
            hoverinfo="text",
            showscale=False,
            colorbar=dict(
                title="K_diff",
                tickmode="array",
                tickvals=k_tickvals,
                ticktext=labels_k,
                len=0.8
            )
        ), row=1, col=2)
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
            domain=dict(x=[0.82, 1], y=[0.5, 1])
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
            domain=dict(x=[0.82, 1], y=[0, 0.5])
        ))
        fig.update_layout(
            height=1500,
            width=1500,
            title=dict(
                text=", ".join(samples_ordered),
                x=0.5,               # Center the title horizontally
                xanchor="center"     # Anchor the title at the center
            ),
            xaxis1=dict(showticklabels=False),  # Hide x-axis of first heatmap
            xaxis2=dict(showticklabels=False),   # Hide x-axis of second heatmap
            yaxis2=dict(showticklabels=False)   # Hide x-axis of second heatmap
        )
        # Save to HTML
        pio.write_html(fig, file=file, auto_open=False)
    
    else:
        # ----- PDF GENERATION -----
        size = 12
        with PdfPages(file) as pdf:  
            fig = plt.figure(figsize=(size, size))
            plt.suptitle(", ".join(samples_ordered), fontsize=16)
            gs = plt.GridSpec(3, 3, width_ratios=[3, 1, 1], height_ratios=[1, 0.1, 1])
            # --- Main Heatmap ---
            ax1 = fig.add_subplot(gs[:, 0])
            colors =list(legend_alleles.values())
            cmap = ListedColormap(colors)
            im = ax1.imshow(m_numeric, cmap=cmap, aspect='auto', interpolation='none')
            # Vertical grid lines
            for g in range(1, genes_n):
                ax1.axhline(y=g-0.5, color='white', linestyle='-', linewidth=1)
            # X and Y labels
            ax1.set_yticks([(g) for g in range(genes_n)])
            ax1.set_yticklabels(list({gene: i+1 for i, gene in enumerate(unique_genes)}))
            ax1.set_xticks([])
            ax1.set_xticklabels([])
            nra_text = geno_db[geno_db["ALLELES_G"].str.match(r"i?\d+_i?\d+")]
            subject_to_index = {subject: idx for idx, subject in enumerate(samples_ordered)}
            gene_to_index = {gene: idx for idx, gene in enumerate(unique_genes)}
            for _, row in nra_text.iterrows():
                if row["subject"] in subject_to_index:
                    y = gene_to_index[row["gene"]]
                    x = (row["id"] - 0.5) * (12 / row["n"])
                    label = allele_to_label.get(row["ALLELES_G"], "")
                    ax1.text(x,y,label, fontsize=5, ha='center', va='center')
            novel_text = geno_db[geno_db["ALLELES_G"].str.contains("_")]
            subject_to_index = {subject: idx for idx, subject in enumerate(samples_ordered)}
            gene_to_index = {gene: idx for idx, gene in enumerate(unique_genes)}
            for _, row in novel_text.iterrows():
                if row["subject"] in subject_to_index:
                    y = gene_to_index[row["gene"]]
                    x = (row["id"] - 0.5) * (12 / row["n"])
                    label = allele_to_label.get(row["ALLELES_G"], "")
                    ax1.text(x,y,label, fontsize=5, ha='center', va='center') 
            # --- Legend ---
            # Modify the legend plotting code
            ax2 = fig.add_subplot(gs[0, 2]) 
            #colors = ["white"] +list(legend_alleles.values())
            colors = list(legend_alleles.values())
            cmap = ListedColormap(colors)
            ax2.imshow(m2_numeric, cmap=cmap, aspect='auto', interpolation='none')
            ax2.set_title("Alleles", fontsize=10)
            # Add horizontal grid lines
            for r in range(1, total_rows):
                ax2.axhline(y=r-0.5, color='black', linestyle='--', linewidth=0.5)
            # Add vertical grid lines to separate columns
            for c in range(1, num_columns):
                col_pos = c * (genes_n * 12 // num_columns) - 0.5
                ax2.axvline(x=col_pos, color='black', linestyle='-', linewidth=1)
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
                ax2.text(
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
            ax2.set_xticks([])
            ax2.set_yticks([])
            # --- K Heatmap ---
            ax3 = fig.add_subplot(gs[:, 1])
            im = ax3.imshow(matrix_binned.values, cmap=cmap_k, aspect='auto')
            # Tick labels
            ax3.set_xticks(np.arange(len(matrix_binned.columns)))
            ax3.set_xticklabels([])
            ax3.set_yticks([])
            ax3.set_yticklabels([])
            # ax4 = legend_k
            ax4 = fig.add_subplot(gs[2, 2])
            for i, (color, label) in enumerate(zip(colors_k, labels_k)):
                ax4.add_patch(plt.Rectangle((0, i), 1, 1, color=color))
                ax4.text(0.5, i + 0.5, label, ha='center', va='center', fontsize=8, color='black')
            ax4.set_ylim(0, len(colors_k))
            ax4.set_xlim(0, 1)
            ax4.set_xticks([])
            ax4.set_yticks([])
            #ax4.set_yticklabels(labels, fontsize=8)
            ax4.set_title("K_diff", fontsize=10)
            #ax4.tick_params(left=False)
            plt.subplots_adjust(
                top=0.95,      # Move subplots closer to the top (0-1 scale)
                bottom=0.05,   # Move subplots closer to the bottom
                left=0.05,     # Move subplots closer to the left
                right=0.95,    # Move subplots closer to the right
                hspace=0     # Control vertical spacing between subplots
            )
            # Save the figure
            plt.rcParams['text.usetex'] = False  
            pdf.savefig(fig, bbox_inches='tight')  # 'tight' removes excess whitespace
            plt.close()
            
    return file