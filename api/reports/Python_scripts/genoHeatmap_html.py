
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.offline import plot
import os
import re
from typing import List, Dict, Union, Optional, Tuple
import tempfile
from dataclasses import dataclass
from typing import Optional
from werkzeug.exceptions import BadRequest
import time
import threading
import traceback
from dataclasses import dataclass
from typing import Optional
from werkzeug.exceptions import BadRequest
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

def allele_palette_old(hap_alleles, nra=True):
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
    #alleles = [a for a in set(hap_alleles) if re.search(r'[012]', a)]
    
    alleles = hap_alleles
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

import time
@dataclass
class ProcessTimer:
    process_name: str
    test_mode: bool
    process_id: int
    start_time: float
    last_step_time: float
    total_pages: Optional[int] = None
    
    @classmethod
    def start(cls, name: str, test_mode: bool = False):
        start = time.time()
        return cls(
            process_name=name,
            test_mode=test_mode,
            process_id=threading.get_ident(),
            start_time=start,
            last_step_time=start
        )
    
    def log_step(self, step_name: str):
        current = time.time()
        elapsed = current - self.last_step_time
        total = current - self.start_time
        print(f"{step_name}: {elapsed:.2f}s (Total: {total:.2f}s)")
        self.last_step_time = current
    
    def set_pages(self, total: int):
        self.total_pages = total
        print(f"/nGenerating {total} pages...")
    
    def log_page(self, page_num: int):
        if not self.total_pages:
            return
        print(f"Page {page_num}/{self.total_pages} completed in {time.time() - self.last_step_time:.2f}s")
        self.last_step_time = time.time()
    
    def finish(self):
        print(f"\n{'-'*20}")
        mode = 'TEST' if self.test_mode else 'Production'
        total = time.time() - self.start_time
        print(f"{mode} {self.process_name} completed in {total:.2f}s")
        print(f"{'-'*20}\n")

def generate_heatmap_html(geno_table,chain = "IGH", remove_chain = True, lk_cutoff= 1, mark_low_lk = True,gene_sort="position", pseudo_genes=False, ORF_genes=False, file= None):
    """
    Create an interactive heatmap of genotype data
    
    Parameters:
    -----------
    geno_table: pandas DataFrame
        DataFrame containing genotype data with columns: subject, gene, GENOTYPED_ALLELES, k_diff, Freq_by_Clone
    chain: str
        Immunoglobulin chain to filter by (default: "IGH")
    remove_chain: bool
        Whether to remove the chain prefix from gene names
    lk_cutoff: float
        K-difference cutoff for marking low likelihood genotypes
    mark_low_lk: bool
        Whether to mark low likelihood genotypes
    file: str, optional
        Output file path (default is a temporary directory)
    n_line: int
        Number of lines per gene
    
    Returns:
    --------
    fig: plotly.graph_objects.Figure
        Interactive heatmap figure
    """
    timer = ProcessTimer.start('Total Report Generation', test_mode=False)
    timer.log_step('Start')
    # Set default file path if not provided
    if file is None:
        file = os.path.join(tempfile.gettempdir(), "genotype_heatmap.html")
        
    # Subset and rename columns | we query more columns than we need for the plot, edit this in the query call
    geno_db = geno_table[["subject", "gene", "GENOTYPED_ALLELES", "k_diff"]].copy()
    
    if not all(gene.startswith(chain) for gene in geno_db["gene"]):
        raise ValueError(f"The chain input {chain} does not match the genes")
    
    geno_db.rename(columns={"GENOTYPED_ALLELES": "ALLELES", "k_diff": "K"}, inplace=True)
    geno_db["ALLELES"] = geno_db["ALLELES"].str.replace("Deletion", "Del")
    all_subjects = geno_db["subject"].unique()
    all_genes = geno_db["gene"].unique()
    index = pd.MultiIndex.from_product([all_subjects, all_genes], names=["subject", "gene"])
    complete_df = pd.DataFrame(index=index).reset_index()
    geno_db = pd.merge(complete_df, geno_db, on=["subject", "gene"], how="left")
    geno_db.loc[geno_db["ALLELES"].isna(), ["ALLELES", "K"]] = ["Unk", np.nan]
    geno_db.loc[geno_db["ALLELES"].str.contains("Del"), "K"] = np.nan
    geno_db = geno_db.assign(ALLELES=geno_db['ALLELES'].str.split(',')).explode('ALLELES')
    timer.log_step('Data Preparation complete')
    # Sort the data
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
        #clean_allele = re.sub(r"\^[0-9]+-", "", allele)  # Fixed regex pattern
        #code = clean_allele.split("_")[0] if "_" in clean_allele else clean_allele
        code = allele.split("_")[0] if "_" in allele else allele
        clean_allele_list.append(code)
    # Keep only numeric values
    clean_allele_list = [item for item in clean_allele_list if item.isdigit()]
    # Assign a unique numeric code from 1 to N
    unique_clean_alleles = sorted(set(clean_allele_list))  # Sort for consistent ordering
    clean_allele_to_code = {allele: str(i + 1) for i, allele in enumerate(unique_clean_alleles)}

    # Map cleaned alleles to their unique numeric code
    allele_code = {}
    for allele in allele_palette_data["AlleleCol"].keys():
        #clean_allele = re.sub(r"\^[0-9]+-", "", allele)
        #code = clean_allele.split("_")[0] if "_" in clean_allele else clean_allele        
        #allele_code[clean_allele] = clean_allele_to_code.get(code, None)  # Use None instead of "UNK"
        code = allele.split("_")[0] if "_" in allele else allele
        allele_code[allele] = clean_allele_to_code.get(code, None)  # Use None instead of "UNK"

    # Handle non-numeric codes
    non_numeric_keys = [key for key, value in allele_code.items() if not (value and value.isdigit())]

    if non_numeric_keys:
        numeric_values = [int(value) for value in allele_code.values() if value and value.isdigit()]
        last_numeric = max(numeric_values, default=0) + 1  # Default to 1 if no numbers exist
        for i, key in enumerate(non_numeric_keys):
            allele_code[key] = str(last_numeric + i)
    
    #geno_db["A_CODE"] = geno_db["ALLELES"].apply(
    #lambda x: allele_code[allele.split("_")[0] if "_" in x else x] 
    #if pd.notna(x) else x
    #)
    
    # Apply allele code mapping efficiently
    def get_allele_code(x):
        if pd.isna(x):
            return x
        key = x.split("_")[0] if "_" in x else x
        return allele_code.get(key, x)
    
    geno_db["A_CODE"] = geno_db["ALLELES"].apply(get_allele_code)
    
    
    
    nra_pattern_compiled = re.compile(r"i?\d+_i?\d+")
    # Handle special NRA cases
    if "NRA" in allele_code:
        geno_db.loc[geno_db["ALLELES"].astype(str).str.contains(nra_pattern_compiled), "A_CODE"]  = allele_code["NRA"]
    
    # Sort by subject, gene location, and allele code
    geno_db = geno_db.sort_values(["subject", "GENE_LOC", "A_CODE"])
    # Assign ID within each subject-gene group
    geno_db["id"] = geno_db.groupby(["subject", "gene"]).cumcount() + 1
    timer.log_step('Allele Code Mapping complete')
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
     
    timer.log_step('DataFrame f Creation complete')
    # Define vertical line function for plotly
    def vline(x=0, color="white"):
        return {
            "type": "line",
            "y0": 0, 
            "y1": 1, 
            "yref": "paper", 
            "x0": x, 
            "x1": x, 
            "line": {"color": color}
        }
    
    n_subjects = len(samples)
    n_cols = 12 * genes_n

    # Extract the A_CODE column directly, convert to string array and reshape
    a_codes = geno_db_f["A_CODE"].astype(str).to_numpy().reshape(n_subjects, n_cols)

    # Create the matrix directly
    m = np.array(a_codes, dtype=object)
    
    hover_texts = [
        f"Individual: {subject}<br />Gene: {gene}<br />Allele: {ALLELES}<br />Kdiff: {k}"
        for subject,gene, ALLELES, k in zip(geno_db_f['subject'], geno_db_f['gene'],geno_db_f['text_bottom'], geno_db_f['K'])
    ]
    hover_texts_2d = np.array(hover_texts).reshape(m.shape)
    samples_ordered = geno_db_f["subject"].unique() 
    # Fill the matrix    
    timer.log_step('Matrix Creation complete')
    # Create grid lines
    gridlines = [vline(x=x+0.5) for x in range(0, genes_n * 12, 12)]
    
    # Set plot dimensions
    plot_height = 500 + 12 * (m.shape[0])
    plot_width = 100 + 2 * (m.shape[1])
    
    legend_alleles = {allele: color for allele, color in allele_palette_data["AlleleCol"].items() if "_" not in allele}
    
    allele_labels = list(legend_alleles.keys())
    allele_colors = list(legend_alleles.values())
    # Number of unique values
    n = len(allele_colors)
    # Build a discrete colorscale with normalized ranges
    discrete_colorscale = []
    for i in range(n):
        start = i /n
        end = (i + 1)/n
        color = allele_colors[i]
        discrete_colorscale.append([start, color])
        discrete_colorscale.append([end, color])
    # Map tick values to middle of each color band
    tickvals = [(i + 0.5)/n *(n+1) for i in range(n)]
    ticktext = allele_labels
    timer.log_step('plot starting')
    # Create the heatmap
    fig = go.Figure(data=go.Heatmap(
        z=m,
        zmin=0,              # new!
        zmax=n+1,              # new! must match number of colors
        opacity=0.9,
        colorscale=discrete_colorscale,
        showscale=True,
        colorbar=dict(
            title="",
            tickmode="array",
            len=0.8,
            #tickvals=list(sorted(all_values, key=int)),
            #ticktext=list(legend_alleles.keys()),
            tickvals=tickvals,
            ticktext=ticktext,
        ),
        hoverinfo="text",
         hovertext=hover_texts_2d
    ))
    
    special_alleles = ["Del","Unk","NR","NRA"]
    # Add text annotations for special alleles
    #ids_text = geno_db[geno_db["text_bottom"].str.match(r"[0-9][0-9]_[0-9][0-9]")]
    ids_text = geno_db[geno_db["text_bottom"].str.match(r"^i?\d{2}_i?\d{2}$")]
    if not ids_text.empty:
        subject_to_index = {subject: idx for idx, subject in enumerate(samples_ordered)}
        ids_text["y"] = ids_text["subject"].map(subject_to_index)
        ids_text = ids_text.dropna(subset=["y"])
        gene_loc = ids_text["GENE_LOC"].to_numpy()
        gene_id = ids_text["id"].to_numpy()
        gene_n = ids_text["n"].to_numpy()
        y = ids_text["y"].to_numpy().astype(int)
        x = (gene_loc - 1) * 12 + (gene_id - 0.5) * ((12 / gene_n) - 1.5)
        
        # Pre-build common annotation parts
        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode="text",
            text=["*"] * len(x),
            textfont=dict(color="black", size=10),
            showlegend=False,
            hoverinfo="skip",
            cliponaxis=True  # Prevents plot from expanding for text
        ))
    
    novel_text = geno_db[geno_db["ALLELES"].str.contains("_", na=False)].copy()
    if not novel_text.empty:
        subject_to_index = {subject: idx for idx, subject in enumerate(samples_ordered)}
        novel_text["y"] = novel_text["subject"].map(subject_to_index)
        novel_text = novel_text.dropna(subset=["y"])

        gene_loc = novel_text["GENE_LOC"].to_numpy()
        gene_id = novel_text["id"].to_numpy()
        gene_n = novel_text["n"].to_numpy()
        y = novel_text["y"].to_numpy().astype(int)
        x = (gene_loc - 1) * 12 + (gene_id - 0.5) * ((12 / gene_n) - 1.5)
        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode="text",
            text=["^"] * len(x),
            textfont=dict(color="black", size=10),
            showlegend=False,
            hoverinfo="skip",
            cliponaxis=True  # Prevents plot from expanding for text
        ))
    
    # Add range slider
    # Configure layout
    fig.update_layout(
        width=plot_width,
        height=plot_height,
        yaxis=dict(
            dtick=1,
            ticktext=samples_ordered,
            tickmode="array",
            tickvals=list(range(len(samples)))
        ),
        xaxis=dict(
            dtick=1,
            ticktext=[gene for gene in genes for _ in range(1)],
            tickmode="array",
            tickvals=list(range(5, 12 * genes_n, 12)),
            rangeslider=dict(visible=True),
            type="linear"
        )
    )
    fig.update_yaxes(range=[-0.5, len(samples_ordered) - 0.5])
    fig.update_xaxes(range=[-0.5, 12 * genes_n - 0.5])
    timer.log_step('Heatmap Creation complete')
    # Save or return figure
    if file:
        plot(fig, filename=file, auto_open=False)
        return fig
    else:
        return fig
    
    timer.finish()