
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

def allele_palette(hap_alleles, nra=True):
    """
    Create the allele color palette for haplotype graphical output.
    """
    alleles = hap_alleles
    base_alleles = sorted(set(a.split('_')[0] for a in alleles))
    # Add base alleles from novel alleles if not in hap_alleles
    novel_bases = {a.split('_')[0] for a in alleles if '_' in a}
    missing_bases = novel_bases - set(alleles)
    alleles += list(missing_bases)
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

GENE_loc =  {
"IGH": ["IGHV7-81", "IGHV4-80", "IGHV3-79", "IGHV5-78", "IGHV7-77", "IGHV3-76", "IGHV3-75", "IGHV3-74", "IGHV3-73", "IGHV3-72", "IGHV3-71", "IGHV2-70", "IGHV1-69D", "IGHV1-69-2", "IGHV3-69-1", "IGHV2-70D", "IGHV1-69", "IGHV1-68", "IGHV1-67", "IGHV3-66", "IGHV3-65", "IGHV3-64", "IGHV3-63", "IGHV3-62", "IGHV4-61", "IGHV3-60", "IGHV4-59", "IGHV1-58", "IGHV3-57", "IGHV7-56", "IGHV4-55", "IGHV3-54", "IGHV3-53", "IGHV3-52", "IGHV3-51-1", "IGHV5-51", "IGHV3-50", "IGHV3-49", "IGHV3-48", "IGHV3-47", "IGHV1-46", "IGHV1-45", "IGHV3-43", "IGHV3-42", "IGHV3-41", "IGHV7-40", "IGHV4-39", "IGHV1-38-4", "IGHV3-38-3", "IGHV3-43D", "IGHV3-42D", "IGHV7-40D", "IGHV4-38-2", "IGHV3-38", "IGHV3-37", "IGHV3-36", "IGHV3-35", "IGHV7-34-1", "IGHV4-34", "IGHV3-33-2", "IGHV3-33", "IGHV3-32", "IGHV4-31", "IGHV3-30-52", "IGHV3-30-5", "IGHV3-30-42", "IGHV4-30-4", "IGHV3-30-33", "IGHV3-30-3", "IGHV3-30-22", "IGHV4-30-2", "IGHV4-30-1", "IGHV3-30-2", "IGHV3-30", "IGHV3-29", "IGHV4-28", "IGHV7-27", "IGHV2-26", "IGHV3-25", "IGHV1-24", "IGHV3-23D", "IGHV3-23", "IGHV3-22", "IGHV3-21", "IGHV3-20", "IGHV3-19", "IGHV1-18", "IGHV1-17", "IGHV3-16", "IGHV3-15", "IGHV1-14", "IGHV3-13", "IGHV1-12", "IGHV3-11", "IGHV2-10", "IGHV3-9", "IGHV1-8", "IGHV5-10-1", "IGHV3-64D", "IGHV3-7", "IGHV3-6", "IGHV2-5", "IGHV7-4-1", "IGHV4-4", "IGHV1-3", "IGHV1-2", "IGHV6-1", "IGHD1-1", "IGHD2-2", "IGHD3-3", "IGHD4-4", "IGHD5-5", "IGHD6-6", "IGHD1-7", "IGHD2-8", "IGHD3-9", "IGHD3-10", "IGHD4-11", "IGHD5-12", "IGHD6-13", "IGHD1-14", "IGHD2-15", "IGHD3-16", "IGHD4-17", "IGHD5-18", "IGHD6-19", "IGHD1-20", "IGHD2-21", "IGHD3-22", "IGHD4-23", "IGHD5-24", "IGHD6-25", "IGHD1-26", "IGHD7-27", "IGHJ1", "IGHJ2", "IGHJ3", "IGHJ4", "IGHJ5", "IGHJ6"],
"IGKO": ["IGKJ5", "IGKJ4", "IGKJ3", "IGKJ2", "IGKJ1", "IGKV4-1", "IGKV5-2", "IGKV7-3", "IGKV2-4", "IGKV1-5", "IGKV1-6", "IGKV3-7", "IGKV1-8", "IGKV1-9", "IGKV2-10", "IGKV3-11", "IGKV1-12", "IGKV1-13", "IGKV2-14", "IGKV3-15", "IGKV1-16", "IGKV1-17", "IGKV2-18", "IGKV2-19", "IGKV3-20", "IGKV6-21", "IGKV1-22", "IGKV2-23", "IGKV2-24", "IGKV3-25", "IGKV2-26", "IGKV1-27", "IGKV2-28", "IGKV2-29", "IGKV2-30", "IGKV3-31", "IGKV1-32", "IGKV1-33", "IGKV3-34", "IGKV1-35", "IGKV2-36", "IGKV1-37", "IGKV2-38", "IGKV1-39", "IGKV2-40", "IGKV2D-40", "IGKV1D-39", "IGKV2D-38", "IGKV1D-37", "IGKV2D-36", "IGKV1D-35", "IGKV3D-34", "IGKV1D-33", "IGKV1D-32", "IGKV3D-31", "IGKV2D-30", "IGKV2D-29", "IGKV2D-28", "IGKV1D-27", "IGKV2D-26", "IGKV3D-25", "IGKV2D-24", "IGKV2D-23", "IGKV1D-22", "IGKV6D-21", "IGKV3D-20", "IGKV2D-19", "IGKV2D-18", "IGKV6D-41", "IGKV1D-17", "IGKV1D-16", "IGKV3D-15", "IGKV2D-14", "IGKV1D-13", "IGKV1D-12", "IGKV3D-11", "IGKV2D-10", "IGKV1D-42", "IGKV1D-43", "IGKV1D-8", "IGKV3D-7", "IGKV1-NL1"],
"IGK": ["IGKJ5", "IGKJ4", "IGKJ3", "IGKJ2", "IGKJ1", "IGKV4-1", "IGKV5-2", "IGKV7-3", "IGKV2-4", "IGKV1-5", "IGKV1-6", "IGKV3-7", "IGKV1-8", "IGKV1-9", "IGKV2-10", "IGKV3-11", "IGKV1E-12", "IGKV1E-13", "IGKV2-14", "IGKV3-15", "IGKV1-16", "IGKV1-17", "IGKV2-18", "IGKV2-19", "IGKV3-20", "IGKV6E-21", "IGKV1-22", "IGKV2-23", "IGKV2-24", "IGKV3-25", "IGKV2-26", "IGKV1-27", "IGKV2E-28", "IGKV2-29", "IGKV2-30", "IGKV3-31", "IGKV1-32", "IGKV1E-33", "IGKV3-34", "IGKV1-35", "IGKV2-36", "IGKV1E-37", "IGKV2-38", "IGKV1E-39", "IGKV2E-40", "IGKV2D-38", "IGKV2D-36", "IGKV1D-35", "IGKV3D-34", "IGKV1D-32", "IGKV3D-31", "IGKV2D-30", "IGKV2D-29", "IGKV1D-27", "IGKV2D-26", "IGKV3D-25", "IGKV2D-24", "IGKV2D-23", "IGKV1D-22", "IGKV3D-20", "IGKV2D-19", "IGKV2D-18", "IGKV6D-41", "IGKV1D-17", "IGKV1D-16", "IGKV3D-15", "IGKV2D-14", "IGKV3D-11", "IGKV2D-10", "IGKV1D-42", "IGKV1D-43", "IGKV1D-8", "IGKV3D-7", "IGKV1-NL1"],
"IGL": ["IGLV(I)-70", "IGLV4-69", "IGLV(I)-68", "IGLV10-67", "IGLV(IV)-66-1", "IGLV(V)-66", "IGLV(IV)-65", "IGLV(IV)-64", "IGLV(I)-63", "IGLV1-62", "IGLV8-61", "IGLV4-60", "IGLV(IV)-59", "IGLV(V)-58", "IGLV6-57", "IGLV(I)-56", "IGLV11-55", "IGLV10-54", "IGLV(IV)-53", "IGLV5-52", "IGLV1-51", "IGLV1-50", "IGLV9-49", "IGLV5-48", "IGLV1-47", "IGLV7-46", "IGLV5-45", "IGLV1-44", "IGLV7-43", "IGLV(I)-42", "IGLV(VII)-41-1", "IGLV1-41", "IGLV1-40", "IGLV5-39", "IGLV(I)-38", "IGLV5-37", "IGLV1-36", "IGLV7-35", "IGLV2-34", "IGLV2-33", "IGLV3-32", "IGLV3-31", "IGLV3-30", "IGLV3-29", "IGLV2-28", "IGLV3-27", "IGLV3-26", "IGLV(VI)-25-1", "IGLV3-25", "IGLV3-24", "IGLV2-23", "IGLV(VI)-22-1", "IGLV3-22", "IGLV3-21", "IGLV(I)-20", "IGLV3-19", "IGLV2-18", "IGLV3-17", "IGLV3-16", "IGLV3-15", "IGLV2-14", "IGLV3-13", "IGLV3-12", "IGLV2-11", "IGLV3-10", "IGLV3-9", "IGLV2-8", "IGLV3-7", "IGLV3-6", "IGLV2-5", "IGLV3-4", "IGLV4-3", "IGLV3-2", "IGLV3-1", "IGLJ1", "IGLJ2", "IGLJ23", "IGLJ3", "IGLJ4", "IGLJ5", "IGLJ6", "IGLJ7"],
"TRB": ["TRBV1", "TRBV2", "TRBV3-1", "TRBV3-12", "TRBV4-1", "TRBV5-1", "TRBV6-1", "TRBV7-1", "TRBV4-2", "TRBV6-23", "TRBV3-2", "TRBV4-3", "TRBV7-2", "TRBV8-1", "TRBV5-2", "TRBV6-4", "TRBV7-3", "TRBV8-2", "TRBV5-3", "TRBV9", "TRBV10-1", "TRBV11-1", "TRBV12-1", "TRBV10-2", "TRBV11-2", "TRBV12-2", "TRBV6-5", "TRBV6-56", "TRBV7-4", "TRBV5-4", "TRBV6-6", "TRBV7-5", "TRBV5-5", "TRBV6-7", "TRBV7-6", "TRBV5-6", "TRBV6-8", "TRBV7-7", "TRBV5-7", "TRBV6-9", "TRBV7-8", "TRBV5-8", "TRBV7-9", "TRBV13", "TRBV10-3", "TRBV11-3", "TRBV12-3", "TRBV12-4", "TRBV12-34", "TRBV12-5", "TRBV14", "TRBV15", "TRBV16", "TRBV17", "TRBV18", "TRBV19", "TRBV20-1", "TRBV21-1", "TRBV22-1", "TRBV23-1", "TRBV24-1", "TRBV25-1", "TRBV26", "TRBV27", "TRBV28", "TRBV29-1", "TRBD1", "TRBJ1-1", "TRBJ1-2", "TRBJ1-3", "TRBJ1-4", "TRBJ1-5", "TRBJ1-6", "TRBD2", "TRBJ2-1", "TRBJ2-2", "TRBJ2-2P", "TRBJ2-3", "TRBJ2-4", "TRBJ2-5", "TRBJ2-6", "TRBJ2-7", "TRBV30"],
"TRA": ["TRAV1-1", "TRAV1-2", "TRAV2", "TRAV3", "TRAV4", "TRAV5", "TRAV6", "TRAV7", "TRAVA", "TRAV8-1", "TRAV9-1", "TRAV10", "TRAV11", "TRAV12-1", "TRAV8-2", "TRAV8-3", "TRAVB", "TRAV13-1", "TRAV14-1", "TRAV11-1", "TRAV12-2", "TRAV8-4", "TRAV8-5", "TRAV13-2", "TRAV14/DV4", "TRAV9-2", "TRAV15", "TRAV12-3", "TRAV8-6", "TRAV16", "TRAV17", "TRAV18", "TRAV19", "TRAVC", "TRAV20", "TRAV21", "TRAV8-6-1", "TRAV22", "TRAV23/DV6", "TRDV1", "TRAV24", "TRAV25", "TRAV26-1", "TRAV8-7", "TRAV27", "TRAV28", "TRAV29/DV5", "TRAV30", "TRAV31", "TRAV32", "TRAV33", "TRAV26-2", "TRAV34", "TRAV35", "TRAV36/DV7", "TRAV37", "TRAV38-1", "TRAV38-2/DV8", "TRAV39", "TRAV40", "TRAV41", "TRAV46", "TRDV2", "TRDD1", "TRDD2", "TRDD3", "TRDJ1", "TRDJ4", "TRDJ2", "TRDJ3", "TRDC", "TRDV3", "TRAJ61", "TRAJ60", "TRAJ59", "TRAJ58", "TRAJ57", "TRAJ56", "TRAJ55", "TRAJ54", "TRAJ53", "TRAJ52", "TRAJ51", "TRAJ50", "TRAJ49", "TRAJ48", "TRAJ47", "TRAJ46", "TRAJ45", "TRAJ44", "TRAJ43", "TRAJ42", "TRAJ41", "TRAJ40", "TRAJ39", "TRAJ38", "TRAJ37", "TRAJ36", "TRAJ35", "TRAJ34", "TRAJ33", "TRAJ32", "TRAJ31", "TRAJ30", "TRAJ29", "TRAJ28", "TRAJ27", "TRAJ26", "TRAJ25", "TRAJ24", "TRAJ23", "TRAJ22", "TRAJ21", "TRAJ20", "TRAJ19", "TRAJ18", "TRAJ17", "TRAJ16", "TRAJ15", "TRAJ14", "TRAJ13", "TRAJ12", "TRAJ11", "TRAJ10", "TRAJ9", "TRAJ8", "TRAJ7", "TRAJ6", "TRAJ5", "TRAJ4", "TRAJ3", "TRAJ2", "TRAJ1"]
}

PSEUDO = {
"IGH": ["IGHV2-10", "IGHV3-52", "IGHV3-47", "IGHV3-71", "IGHV3-22", "IGHV4-55", "IGHV1-68", "IGHV2-10", "IGHV5-78", "IGHV3-32", "IGHV3-33-2", "IGHV3-38-3", "IGHV3-25", "IGHV3-19", "IGHV7-40", "IGHV3-63", "IGHV3-62", "IGHV3-29", "IGHV3-54", "IGHV1-38-4", "IGHV7-34-1", "IGHV1-38-4", "IGHV3-30-2", "IGHV3-69-1", "IGHV3-30-22", "IGHV1-f", "IGHV3-30-33", "IGHV3-38", "IGHV7-81", "IGHV3-35", "IGHV3-16", "IGHV3-30-52", "IGHV1-69D", "IGHD1-14", "IGHV3-30-42"],
"IGK": ["IGKV7-3", "IGKV6D-41", "IGKV3D-34", "IGKV3D-31", "IGKV3D-25", "IGKV3/OR22-2", "IGKV3/OR2-5", "IGKV3/OR2-268", "IGKV3-34", "IGKV3-31", "IGKV3-25", "IGKV2D-38", "IGKV2D-36", "IGKV2D-24", "IGKV2D-23", "IGKV2D-19", "IGKV2D-18", "IGKV2D-14", "IGKV2D-10", "IGKV2/OR22-4", "IGKV2/OR22-3", "IGKV2/OR2-8", "IGKV2/OR2-7D", "IGKV2/OR2-7", "IGKV2/OR2-4", "IGKV2/OR2-2", "IGKV2/OR2-10", "IGKV2/OR2-1", "IGKV2-4", "IGKV2-38", "IGKV2-36", "IGKV2-26", "IGKV2-23", "IGKV2-19", "IGKV2-18", "IGKV2-14", "IGKV2-10", "IGKV1D-42", "IGKV1D-37", "IGKV1D-35", "IGKV1D-32", "IGKV1D-27", "IGKV1D-22", "IGKV1/ORY-1", "IGKV1/OR9-2", "IGKV1/OR9-1", "IGKV1/OR22-5", "IGKV1/OR22-1", "IGKV1/OR2-9", "IGKV1/OR2-6", "IGKV1/OR2-3", "IGKV1/OR2-2", "IGKV1/OR2-118", "IGKV1/OR2-11", "IGKV1/OR2-108", "IGKV1/OR2-1", "IGKV1/OR2-0", "IGKV1/OR15-118", "IGKV1/OR10-1", "IGKV1/OR1-1", "IGKV1/OR-4", "IGKV1/OR-3", "IGKV1/OR-2", "IGKV1-37", "IGKV1-35", "IGKV1-32", "IGKV1-22"],
"IGL": ["IGLV7-35", "IGLV3-7", "IGLV3-6", "IGLV3-4", "IGLV3-32", "IGLV3-31", "IGLV3-30", "IGLV3-29", "IGLV3-26", "IGLV3-24", "IGLV3-2", "IGLV3-17", "IGLV3-15", "IGLV3-13", "IGLV2-NL1", "IGLV2-5", "IGLV2-34", "IGLV2-33", "IGLV2-28", "IGLV11-55", "IGLV10-67", "IGLV1-62", "IGLV1-50", "IGLV(VII)-41-1", "IGLV(VI)-25-1", "IGLV(VI)-22-1", "IGLV(V)-66", "IGLV(V)-58", "IGLV(IV)/OR22-2", "IGLV(IV)/OR22-1", "IGLV(IV)-66-1", "IGLV(IV)-65", "IGLV(IV)-64", "IGLV(IV)-59", "IGLV(IV)-53", "IGLV(I)-70", "IGLV(I)-68", "IGLV(I)-63", "IGLV(I)-56", "IGLV(I)-42", "IGLV(I)-38", "IGLV(I)-20", "IGLL4", "IGLL2", "IGLJ5", "IGLJ4", "IGLC5", "IGLC4", "IGLC/OR22-2", "IGLC/OR22-1"],
"TRB": ["TRBV1", "TRBV7-1", "TRBV3-2", "TRBV8-1", "TRBV5-2", "TRBV8-2", "TRBV5-3", "TRBV12-1", "TRBV12-2", "TRBV7-5", "TRBV5-7", "TRBV6-7", "TRBV5-7", "TRBV17", "TRBV21-1", "TRBV22-1", "TRBV23-1", "TRBV26", "TRBJ2-2P"],
"TRA": ["TRAVA", "TRAV11", "TRAVB", "TRAV14-1", "TRAV11-1", "TRAV8-5", "TRAV15", "TRAVC", "TRAV8-6-1", "TRAV8-7", "TRAV28", "TRAV31", "TRAV32", "TRAV33", "TRAV37", "TRAV46"]
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


def hapHeatmap(hap_table,
               chain="IGH",
               genes_order=None,
               removeIGH=True,
               lk_cutoff=1,
               mark_low_lk=True,
               size_annot=1.5,
               color_y=None,
               order_subject=None,
               file=None,
               size_text=None,
               ylabel_size=1):
    """
    Generates a graphical heatmap output of alleles per gene in multiple samples.
    """
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


    # Selecting relevant columns
    hap_table = hap_table[["subject", "gene", *hapBy_cols, "alleles", "k1", "k2", "k3", "k4"]]

    # Ensuring all samples have all genes, filling missing values with "Unk"
    all_subjects = hap_table['subject'].unique()
    all_genes = hap_table['gene'].unique()
    index = pd.MultiIndex.from_product([all_subjects, all_genes], names=["subject", "gene"])

    index = pd.MultiIndex.from_product([all_subjects, all_genes], names=["subject", "gene"])
    complete_df = pd.DataFrame(index=index).reset_index()

    hap_table = pd.merge(complete_df, hap_table, on=["subject", "gene"], how="left")

    #hap_table = hap_table.set_index(["subject", "gene"]).reindex(complete_index).reset_index()
    hap_table[hapBy_cols] = hap_table[hapBy_cols].replace({"": "Unk", None: "Unk"}).fillna("Unk")
    hap_table["alleles"] = hap_table["alleles"].replace({"": "Unk", None: "Unk"}).fillna("Unk")

    # Sorting the data
    hap_table["order"] = hap_table["gene"].map({g: i for i, g in enumerate(genes_order)})
    hap_table = hap_table.dropna(subset=["order"]).sort_values("order").reset_index(drop=True)

    if removeIGH:
        genes_order = [gene.replace(chain, "") for gene in genes_order]
        hap_table["gene"] = hap_table["gene"].str.replace(chain, "", regex=False)
        hap_table["gene"] = pd.Categorical(hap_table["gene"], categories=[g.replace(chain, "") for g in genes_order], ordered=True)
    else:
        hap_table["gene"] = pd.Categorical(hap_table["gene"], categories=genes_order, ordered=True)

    # Rename genes to numbers
    unique_genes_sorted = sorted(hap_table["gene"].unique(), key=lambda x: list(genes_order).index(x))
    gene_loc = {gene: i+1 for i, gene in enumerate(unique_genes_sorted)}
    hap_table["GENE_LOC"] = hap_table["gene"].map(gene_loc)

    # Fix NA in k columns
    hap_table.fillna(np.inf, inplace=True)

    # Melt haplotype columns to one
    hap_table = hap_table.melt(
        id_vars=["subject", "gene", "GENE_LOC"],
        value_vars=hapBy_cols,
        var_name="hapBy",
        value_name="ALLELES"
    )

    hap_table = hap_table.assign(ALLELES=hap_table['ALLELES'].str.split(',')).explode('ALLELES')

    distinct_n = (
        hap_table.groupby(["subject", "gene", "hapBy"])
        .size()
        .reset_index(name="n")
    )

    distinct_n["n"] = distinct_n["n"]

    hap_table = hap_table.merge(distinct_n, on=["subject", "gene", "hapBy"], how="left")

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
    print(non_numeric_keys)
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
    
    print(allele_code)
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
    for subject, gene, gene_loc, alleles_g, a_code, text_bottom, n_lines, id_, hapBy,n in zip(
        hap_table["subject"], 
        hap_table["gene"], 
        hap_table["GENE_LOC"], 
        hap_table["ALLELES_G"], 
        hap_table["A_CODE"], 
        hap_table["text_bottom"], 
        hap_table["line"],
        hap_table["id"],
        hap_table["hapBy"],
        hap_table["n"]
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
                "hapBy": hapBy
            }
            hap_table_f.append(row)
            if n == 5 and ((id_ == 1 and i == 1) or (id_ == 5 and i == 1)):
                hap_table_f.append(row.copy())
                    
    hap_table_f = pd.DataFrame(hap_table_f)
    
    # Deviding the data into two dataframes
    hap_table_1 = hap_table_f[hap_table_f["hapBy"] == hapBy_cols[0]]
    hap_table_2 = hap_table_f[hap_table_f["hapBy"] == hapBy_cols[1]]
    
    #m1 = np.full((len(samples), 12 * genes_n), "", dtype=object)
    #hap_table_1 = hap_table_1.reset_index(drop=True)
    #for i, subject in enumerate(samples):
    #    for j in range(12 * genes_n):
    #        a_code = str(hap_table_1.loc[(i*genes_n*12)+j, "A_CODE"])
    #        m1[i, j] = str(a_code)
    #
    
    n_subjects = len(samples)
    n_cols = 12 * genes_n
    
    # Extract the A_CODE column directly, convert to string array and reshape
    a_codes = hap_table_1["A_CODE"].astype(str).to_numpy().reshape(n_subjects, n_cols)
    # Create the matrix directly
    m1 = np.array(a_codes, dtype=object)
    
    # Extract the A_CODE column directly, convert to string array and reshape
    a_codes = hap_table_2["A_CODE"].astype(str).to_numpy().reshape(n_subjects, n_cols)
    # Create the matrix directly
    m2 = np.array(a_codes, dtype=object)

    #m2 = np.full((len(samples), 12 * genes_n), "", dtype=object)
    #hap_table_2 = hap_table_2.reset_index(drop=True)
    #for i, subject in enumerate(samples):
    #    for j in range(12 * genes_n):
    #        a_code = str(hap_table_2.loc[int(i*genes_n*12)+j, "A_CODE"])
    #        m2[i, j] = str(a_code)

    # Combine alleles and NRA text_bottom entries into a single legend
    bottom_annot_mask = hap_table["text_bottom"].str.contains(r"i?\d{2}_i?\d{2}", regex=True)
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

    # Sort alleles by code value for consistent order
    sorted_alleles = sorted(allele_code.items(), key=lambda x: x[1])
    # Layout settings
    items_per_column = (n_alleles + 2) // 3  # 3 columns
    total_rows = items_per_column
    legend_width = genes_n * 12
    m3 = np.zeros((total_rows, legend_width))
    # Fill the legend matrix with color codes
    for i, (allele, code_val) in enumerate(sorted_alleles):
        column = i // items_per_column
        row = i % items_per_column
        col_start = column * (legend_width // 3)
        col_width = 5  # Width of color block
        m3[row, col_start:col_start + col_width] = code_val

    longest_allele = max(len(allele) for allele in allele_palette_data["AlleleCol"].keys()) * 3 + 40
    # Set the height and width of plot
    height = (samples_n * 0.1 + 2 + m3.shape[0] * 0.5) *2
    width = genes_n * 0.3 + 1.5
    size_text = m1.shape[0] / (height * width)
    size_text_leg = m3.shape[1] / (width * longest_allele) + 2
    # ----- SETUP -----
    legend_alleles = {allele: color for allele, color in allele_palette_data["AlleleCol"].items() if "_" not in allele}
    # Convert m to numeric matrix using allele_to_index
    m1_numeric = m1.astype(float)
    m2_numeric = m2.astype(float)
    m3_numeric = m3.astype(float)
    unique_genes = hap_table_f["gene"].unique()
    samples_ordered = hap_table_f["subject"].unique()

    # ----- PDF GENERATION -----
    with PdfPages(file) as pdf:  
        fig = plt.figure(figsize=(width, height))
        x = True
        gs = plt.GridSpec(3, 1, height_ratios=[3, 3, 1], hspace=0.05)
        # --- Main Heatmap ---
        ax1 = plt.subplot(gs[0])
        colors =list(legend_alleles.values())
        cmap = ListedColormap(colors)
        im = ax1.imshow(m1_numeric, cmap=cmap, aspect='auto', interpolation='none')
        # Vertical grid lines
        for g in range(1, genes_n):
            ax1.axvline(x=g*12-0.5, color='white', linestyle='-', linewidth=1)
        # X and Y labels
        ax1.set_xticks([(g*12+6-0.5) for g in range(genes_n)])
        ax1.set_xticklabels(list({gene: i+1 for i, gene in enumerate(unique_genes)}), rotation=90)
        ax1.set_yticks(range(len(samples)))
        ax1.set_yticklabels(samples_ordered, fontsize=8)
        ax1.title.set_text(hapBy_cols[0])
        # Optional: draw 
        nra_text = hap_1[hap_1["text_bottom"].str.match(r"[0-9][0-9]_[0-9][0-9]")]
        subject_to_index = {subject: idx for idx, subject in enumerate(samples_ordered)}
        for _, row in nra_text.iterrows():
            if row["subject"] in subject_to_index:
                y = subject_to_index[row["subject"]]
                x = (row["GENE_LOC"] - 1) * 12 + (row["id"] - 0.5) * (12 / row["n"])
                ax1.text(x,y,"*", fontsize=5, ha='center', va='center')
        novel_text = hap_1[hap_1["ALLELES"].str.contains("_")]
        subject_to_index = {subject: idx for idx, subject in enumerate(samples_ordered)}
        for _, row in novel_text.iterrows():
            if row["subject"] in subject_to_index:
                y = subject_to_index[row["subject"]]
                x = (row["GENE_LOC"] - 1) * 12 + (row["id"] - 0.5) * (12 / row["n"])
                ax1.text(x,y,"^", fontsize=5, ha='center', va='center')
                
        ax2 = plt.subplot(gs[1])
        colors =list(legend_alleles.values())
        cmap = ListedColormap(colors)
        im = ax2.imshow(m2_numeric, cmap=cmap, aspect='auto', interpolation='none')
        # Vertical grid lines
        for g in range(1, genes_n):
            ax2.axvline(x=g*12-0.5, color='white', linestyle='-', linewidth=1)
        # X and Y labels
        ax2.set_xticks([(g*12+6-0.5) for g in range(genes_n)])
        ax2.set_xticklabels(list({gene: i+1 for i, gene in enumerate(unique_genes)}), rotation=90)
        ax2.set_yticks(range(len(samples)))
        ax2.set_yticklabels(samples_ordered, fontsize=8)
        ax2.title.set_text(hapBy_cols[1])
        # Optional: draw 
        nra_text = hap_2[hap_2["text_bottom"].str.match(r"[0-9][0-9]_[0-9][0-9]")]
        subject_to_index = {subject: idx for idx, subject in enumerate(samples_ordered)}
        for _, row in nra_text.iterrows():
            if row["subject"] in subject_to_index:
                y = subject_to_index[row["subject"]]
                x = (row["GENE_LOC"] - 1) * 12 + (row["id"] - 0.5) * (1 / row["n"])
                ax2.text(x,y,"*", fontsize=5, ha='center', va='center')
        novel_text = hap_2[hap_2["ALLELES"].str.contains("_")]
        subject_to_index = {subject: idx for idx, subject in enumerate(samples_ordered)}
        for _, row in novel_text.iterrows():
            if row["subject"] in subject_to_index:
                y = subject_to_index[row["subject"]]
                x = (row["GENE_LOC"] - 1) * 12 + (row["id"] - 0.5) * (12 / row["n"])
                ax2.text(x,y,"^", fontsize=5, ha='center', va='center') 
        # --- Legend ---
        # Modify the legend plotting code
        ax3 = plt.subplot(gs[2])
        colors = ["white"] +list(legend_alleles.values())
        cmap = ListedColormap(colors)
        ax3.imshow(m3_numeric, cmap=cmap, aspect='auto', interpolation='none')
        # Add horizontal grid lines
        for r in range(1, total_rows):
            ax3.axhline(y=r-0.5, color='black', linestyle='--', linewidth=0.5)
        # Add vertical grid lines to separate columns
        for c in range(1, 3):
            col_pos = c * (genes_n * 12 // 3) - 0.5
            ax3.axvline(x=col_pos, color='black', linestyle='-', linewidth=1)
        # Add text labels for each allele
        for i, (allele, code_val) in enumerate(sorted_alleles):
            column = i // items_per_column
            row = i % items_per_column
            # Calculate text position (to right of colored area)
            col_start = column * (genes_n * 12 // 3)
            text_x = col_start + 6  # Position text after colored area
            # Add the text label
            ax3.text(text_x, row, allele, fontsize=8, ha='left', va='center')
        ax3.set_xticks([])
        ax3.set_yticks([])
        plt.subplots_adjust(
            top=0.95,      # Move subplots closer to the top (0-1 scale)
            bottom=0.05,   # Move subplots closer to the bottom
            left=0.05,     # Move subplots closer to the left
            right=0.95,    # Move subplots closer to the right
            hspace=0.05     # Control vertical spacing between subplots
        )
        # Save the figure
        #pdf.savefig(fig, bbox_inches='tight')  # 'tight' removes excess whitespace
        pdf.savefig(fig, bbox_inches='tight')  # 'tight' removes excess whitespace
        plt.close()

    return plt