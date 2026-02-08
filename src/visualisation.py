"""
This module provides visualization tool for metagenomic analyses
Creates several plots:
1) Contig Length Distribution
2) Contig Length Distribution according to Circularity
3) Contig Length Boxplot
4) Scatterplot of Relationship between Contig Length and Coverage
5) Comparison of quality of large circular contigs
6) Comparison of large circular contigs identified for each Phylum
"""

from pathlib import Path
import sys
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from analysis import compile_data

# create absolute path to root folder in which this code is
HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent
RESULTS_DIR = PROJECT_ROOT / "results"

# append folder into sys.path
sys.path.append(str(HERE))

# file path dictionary for every file
file_paths = {
    'headers_mylo': RESULTS_DIR / 'myloasm_assembly_headers.txt',
    'headers_meta': RESULTS_DIR / 'metamdbg_assembly_headers.txt',
    'checkm2_mylo': RESULTS_DIR / 'checkm2/myloasm/quality_report.tsv',
    'checkm2_meta': RESULTS_DIR / 'checkm2/metamdbg/quality_report.tsv',
    'gtdb_mylo_ar': RESULTS_DIR / 'gtdbtk/myloasm/classify/gtdbtk.ar53.summary.tsv',
    'gtdb_mylo_bac': RESULTS_DIR / 'gtdbtk/myloasm/classify/gtdbtk.bac120.summary.tsv',
    'gtdb_meta_ar': RESULTS_DIR / 'gtdbtk/metamdbg/classify/gtdbtk.ar53.summary.tsv',
    'gtdb_meta_bac': RESULTS_DIR / 'gtdbtk/metamdbg/classify/gtdbtk.bac120.summary.tsv',}

# load data
df = compile_data(file_paths)

# set visualization theme so that graphs look nice
sns.set_theme(style="whitegrid", context="talk")

# use different colours for assemblers
PALETTE = {"myloasm": "red", "metaMDBG": "blue"}

# set a definitive order for assemblers
ASSEMBLER_ORDER = ["myloasm", "metaMDBG"]

# contig length distrubution
print("Plotting Contig Length Distribution...")
plt.figure(figsize=(12, 6))

sns.kdeplot(
    data=df,
    x="length",
    hue="assembler",    # different colours for assemblers
    hue_order=ASSEMBLER_ORDER,
    palette=PALETTE,    # red and blue
    log_scale=True,      # apply logarithmic scale
    common_norm=False,   # each assembler is normalised separately
    fill=False,          # remove filling
    linewidth=4)          # make the line thicker so its better visible

plt.title("Contig Length Distribution")
plt.xlabel("Contig Length (bp, log)")
plt.ylabel("Density")
plt.tight_layout()

# distribution according to circularity
print("Plotting Contig Length Distribution According to Circularity...")
plt.figure(figsize=(12, 7))

# circular (full line)
sns.kdeplot(
    data=df[df['is_circular'] == True], # extract only circular contigs
    x="length",
    hue="assembler",    # colour based on assembler
    hue_order=ASSEMBLER_ORDER,
    palette=PALETTE,
    log_scale=True, # apply logarithmic scale
    common_norm=False,  # normalise each assembler separately
    linewidth=3,
    linestyle="-")  # full line

# non-circular (dashed line)
sns.kdeplot(
    data=df[df['is_circular'] == False],    # extract only non-circular
    x="length",
    hue="assembler",    # colour based on assembler
    hue_order=ASSEMBLER_ORDER,
    palette=PALETTE,
    log_scale=True, # apply logarithmic scale
    common_norm=False,  # normalise each assembler separately
    linewidth=2,
    linestyle="--") # dashed line

# define elements which will be displayed in the legend
legend_elements = [
    Line2D([0], [0], color=PALETTE['myloasm'], lw=3, label='myloasm'),
    Line2D([0], [0], color=PALETTE['metaMDBG'], lw=3, label='metaMDBG'),
    Line2D([0], [0], color='black', lw=2, linestyle='-', label='Circular'),
    Line2D([0], [0], color='black', lw=2, linestyle='--', label='Non-circular')]

plt.legend(handles=legend_elements, loc='upper right', title="Legend")
plt.title("Contig Length Distribution")
plt.xlabel("Contig Length (bp, log)")
plt.ylabel("Density")
plt.tight_layout()

# boxplot comparing different assemblers
print("Plotting Boxplot of Contig Lengths...")

# make a copy of our data
df_box = df.copy()
# logarithm is applied to Length column
df_box['log_length'] = np.log10(df_box['length'])

plt.figure(figsize=(10, 6))

sns.boxplot(
    data=df_box,
    x="assembler",
    y="log_length",
    palette=PALETTE,
    hue="assembler", # colour based boxplots
    legend=False,
    whis=75)    # set whiskers length

plt.title("Boxplot Comparing Assemblers Based on Contig Length")
plt.xlabel("Assembler")
plt.ylabel("Contig Length (bp, log)")
plt.tight_layout()
plt.show()

# correlation between contig length and coverage
print("Plotting Correlation between Contig Length and Coverage...")

# create a template for 2 graphs to be placed next to each other
# sharey=True ensures both graphs have the same Y-axis
fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharey=True)

# scatterplot for myloasm
sns.scatterplot(
    data=df[df['assembler'] == 'myloasm'],  # take only data for myloasm
    x="length",
    y="coverage",
    ax=axes[0], # plot in the first "column"
    color=PALETTE['myloasm'],
    alpha=0.5,  # make it half-transparent
    s=20,
    edgecolor=None) # remove edges

axes[0].set_title("myloasm: Length vs Coverage")
axes[0].set_xscale("log")
axes[0].set_yscale("log")
axes[0].set_xlabel("Contig Length (bp, log)")
axes[0].set_ylabel("Coverage (log)")

# scatterplot for metaMDBG
sns.scatterplot(
    data=df[df['assembler'] == 'metaMDBG'], # take only data for metaMDBG
    x="length",
    y="coverage",
    ax=axes[1], # plot in the second "column"
    color=PALETTE['metaMDBG'],
    alpha=0.5,  # make it half-transparent
    s=20,
    edgecolor=None) # remove edges

axes[1].set_title("metaMDBG: Length vs Coverage")
axes[1].set_xscale("log")
axes[1].set_yscale("log")
axes[1].set_xlabel("Contig Length (bp, log)")
axes[1].set_ylabel("")

plt.suptitle("Comparison of Relationship between Contig Length and Coverage", size=25)
plt.tight_layout()

# extract only contigs that are large and circular
large_circ = df[df['is_large_circular'] == True].copy()

if len(large_circ) > 0:
    # comparison of quality of large circular contigs
    print("Plotting quality of large circular contigs...")
    plt.figure(figsize=(12, 6))

    # barchart showing counts in each category
    sns.countplot(
        data=large_circ,
        x="quality_category",
        hue="assembler",    # colour based on assembler
        order=["High", "Medium", "Low"],
        hue_order=ASSEMBLER_ORDER,
        palette=PALETTE)

    plt.title("Quality of Large Circular Contigs (>500kb)")
    plt.xlabel("Quality Category")
    plt.ylabel("Number of large circular contigs")
    plt.tight_layout()

    # comparison of large circular contigs identified for each Phylum
    print("Plotting number of large circular contigs per Phylum...")
    plt.figure(figsize=(12, 10))

    # set order of Phylums from the most abundant to the least abundant
    order = large_circ['Phylum'].value_counts().index

    # barchart showing number of large circular contigs for each Phylum
    sns.countplot(
        data=large_circ,    # take into account only large circular contigs
        y="Phylum",
        hue="assembler",    # colour based on assembler
        order=order,
        hue_order=ASSEMBLER_ORDER,
        palette=PALETTE)

    plt.title("Number of Large Circular Contigs per Phylum", size=20)
    plt.xlabel("Number of large circular contigs")
    plt.ylabel("Phylum")
    plt.tight_layout()
    plt.subplots_adjust(left=0.18, right=0.92)
else:
    print("No large circular contigs found.")

plt.show()
