#!/usr/bin/env python3
"""
CTC Cluster vs Single CTC - Neutrophil-related Chemokine Heatmap Analysis
Based on GSE51827 (Aceto et al., 2014, Cell)
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Set plot style
plt.style.use('seaborn-whitegrid')
plt.rcParams['font.size'] = 12
plt.rcParams['figure.dpi'] = 150

# ============================================
# 1. Load and Process Data
# ============================================
print("=" * 60)
print("Loading GSE51827 data...")
print("=" * 60)

# Read the expression data
data_path = "/root/tumor_cell/data/GSE51827_readCounts.xls"
df = pd.read_excel(data_path, index_col=0)
print(f"Expression data shape: {df.shape}")

# Read the platform annotation
platform_path = "/root/tumor_cell/data/GSE51827_platform.xls"
platform = pd.read_excel(platform_path)
print(f"Platform annotation shape: {platform.shape}")

# Create ID to symbol mapping
id_to_symbol = dict(zip(platform['ID'], platform['symbol']))

# Map gene IDs to symbols
df.index = df.index.map(lambda x: id_to_symbol.get(x, x))
df = df[~df.index.isna()]  # Remove rows with no symbol mapping

# Aggregate duplicate gene symbols by sum
df = df.groupby(df.index).sum()

print(f"Data shape after symbol mapping: {df.shape}")
print(f"\nSample gene names: {df.index[:20].tolist()}")

# ============================================
# 2. Identify Sample Groups
# ============================================
print("\n" + "=" * 60)
print("Identifying sample groups...")
print("=" * 60)

# Based on naming convention: CL = Cluster, SC = Single Cell
cluster_samples = [col for col in df.columns if '_CL' in col or 'CL#' in col]
single_samples = [col for col in df.columns if '_SC' in col or 'SC#' in col]

print(f"\nCTC Clusters samples ({len(cluster_samples)}):")
print(cluster_samples)
print(f"\nSingle CTC samples ({len(single_samples)}):")
print(single_samples)

# ============================================
# 3. Define Neutrophil-related Chemokines
# ============================================
print("\n" + "=" * 60)
print("Defining neutrophil-related chemokines...")
print("=" * 60)

# Comprehensive list of neutrophil-related chemokines and their receptors
neutrophil_chemokines = [
    # CXC chemokines (primary neutrophil attractants)
    'CXCL1', 'CXCL2', 'CXCL3', 'CXCL5', 'CXCL6', 'CXCL7', 'CXCL8',  # IL8
    'CXCL9', 'CXCL10', 'CXCL11', 'CXCL12', 'CXCL13', 'CXCL14', 'CXCL16',
    # CC chemokines
    'CCL2', 'CCL3', 'CCL4', 'CCL5', 'CCL7', 'CCL8', 'CCL11', 'CCL20',
    # Chemokine receptors
    'CXCR1', 'CXCR2', 'CXCR3', 'CXCR4', 'CCR1', 'CCR2', 'CCR5', 'CCR7',
    # IL8 alternative names
    'IL8',
    # Additional inflammatory factors
    'IL6', 'IL1B', 'IL1A', 'TNF', 'TNFAIP3', 'CSF3', 'CSF2',
    'S100A8', 'S100A9', 'S100A12',
    # Adhesion molecules
    'ICAM1', 'VCAM1', 'SELE', 'SELP', 'ITGB2', 'ITGAM',
    # Other neutrophil-related genes
    'FPR1', 'FPR2', 'C5AR1', 'PTGS2', 'MMP9', 'MMP2',
    # Plakoglobin (key gene in this study)
    'JUP'
]

# Find which chemokines are in the dataset
available_chemokines = [gene for gene in neutrophil_chemokines if gene in df.index]
missing_chemokines = [gene for gene in neutrophil_chemokines if gene not in df.index]

print(f"\nAvailable chemokines/genes in dataset ({len(available_chemokines)}):")
print(available_chemokines)
print(f"\nMissing genes ({len(missing_chemokines)}):")
print(missing_chemokines)

# ============================================
# 4. Extract Chemokine Expression Data
# ============================================
print("\n" + "=" * 60)
print("Extracting chemokine expression data...")
print("=" * 60)

if len(available_chemokines) > 0:
    # Extract data for available chemokines
    chemokine_df = df.loc[available_chemokines].copy()
    
    # Calculate log2(count + 1) for better visualization
    chemokine_log = np.log2(chemokine_df + 1)
    
    print(f"\nChemokine expression matrix shape: {chemokine_log.shape}")
    print(chemokine_log.head())
    
    # ============================================
    # 5. Statistical Analysis: Cluster vs Single
    # ============================================
    print("\n" + "=" * 60)
    print("Performing differential expression analysis...")
    print("=" * 60)
    
    results = []
    for gene in available_chemokines:
        cluster_expr = chemokine_log.loc[gene, cluster_samples].values
        single_expr = chemokine_log.loc[gene, single_samples].values
        
        # Calculate statistics
        cluster_mean = np.mean(cluster_expr)
        single_mean = np.mean(single_expr)
        fold_change = cluster_mean - single_mean  # log2 fold change
        
        # T-test
        t_stat, p_value = stats.ttest_ind(cluster_expr, single_expr)
        
        results.append({
            'Gene': gene,
            'Cluster_Mean': round(cluster_mean, 3),
            'Single_Mean': round(single_mean, 3),
            'Log2FC': round(fold_change, 3),
            'P_value': round(p_value, 4) if not np.isnan(p_value) else 1.0
        })
    
    results_df = pd.DataFrame(results)
    results_df['Significant'] = results_df['P_value'] < 0.05
    results_df = results_df.sort_values('Log2FC', ascending=False)
    
    print("\nDifferential Expression Results (sorted by Log2FC):")
    print(results_df.to_string(index=False))
    
    # Save results
    results_df.to_csv('/root/tumor_cell/analysis/chemokine_DE_results.csv', index=False)
    print("\n[Saved] Differential expression results to: /root/tumor_cell/analysis/chemokine_DE_results.csv")
    
    # ============================================
    # 6. Create Heatmap
    # ============================================
    print("\n" + "=" * 60)
    print("Creating heatmap...")
    print("=" * 60)
    
    # Prepare data for heatmap
    # Reorder columns: Clusters first, then Singles
    ordered_samples = cluster_samples + single_samples
    heatmap_data = chemokine_log[ordered_samples]
    
    # Z-score normalization by row (gene)
    row_means = heatmap_data.mean(axis=1)
    row_stds = heatmap_data.std(axis=1)
    # Avoid division by zero
    row_stds = row_stds.replace(0, 1)
    heatmap_zscore = heatmap_data.subtract(row_means, axis=0).divide(row_stds, axis=0)
    
    # Sort genes by fold change
    gene_order = results_df['Gene'].tolist()
    heatmap_zscore = heatmap_zscore.loc[gene_order]
    
    # Create figure
    fig_height = max(8, len(available_chemokines) * 0.4)
    fig, ax = plt.subplots(figsize=(16, fig_height))
    
    # Create heatmap
    sns.heatmap(
        heatmap_zscore,
        cmap='RdBu_r',
        center=0,
        vmin=-2, vmax=2,
        xticklabels=True,
        yticklabels=True,
        ax=ax,
        cbar_kws={'label': 'Z-score', 'shrink': 0.5}
    )
    
    # Add vertical line to separate groups
    ax.axvline(x=len(cluster_samples), color='black', linewidth=2)
    
    # Styling
    ax.set_title('Neutrophil-related Chemokines: CTC Clusters vs Single CTCs\n(GSE51827, Aceto et al. 2014, Cell)', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel('Samples', fontsize=12)
    ax.set_ylabel('Genes', fontsize=12)
    
    # Rotate x labels
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.yticks(fontsize=10)
    
    # Add text annotations for groups
    ax.text(len(cluster_samples)/2, -0.5, f'CTC Clusters\n(n={len(cluster_samples)})', 
            ha='center', va='bottom', fontsize=10, fontweight='bold', color='#E74C3C')
    ax.text(len(cluster_samples) + len(single_samples)/2, -0.5, f'Single CTCs\n(n={len(single_samples)})', 
            ha='center', va='bottom', fontsize=10, fontweight='bold', color='#3498DB')
    
    plt.tight_layout()
    
    # Save heatmap
    output_path = '/root/tumor_cell/analysis/chemokine_heatmap.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\n[Saved] Heatmap to: {output_path}")
    
    # Also save as PDF for publication quality
    plt.savefig('/root/tumor_cell/analysis/chemokine_heatmap.pdf', bbox_inches='tight', facecolor='white')
    print("[Saved] PDF version to: /root/tumor_cell/analysis/chemokine_heatmap.pdf")
    
    plt.close()
    
    # ============================================
    # 7. Create Summary Bar Plot
    # ============================================
    print("\n" + "=" * 60)
    print("Creating fold change bar plot...")
    print("=" * 60)
    
    fig, ax = plt.subplots(figsize=(10, max(6, len(available_chemokines) * 0.35)))
    
    # Color by significance and direction
    colors = []
    for _, row in results_df.iterrows():
        if row['Significant'] and row['Log2FC'] > 0:
            colors.append('#E74C3C')  # Red - upregulated in clusters
        elif row['Significant'] and row['Log2FC'] < 0:
            colors.append('#3498DB')  # Blue - downregulated in clusters
        else:
            colors.append('#95A5A6')  # Gray - not significant
    
    y_pos = range(len(results_df))
    bars = ax.barh(y_pos, results_df['Log2FC'], color=colors, edgecolor='black', linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(results_df['Gene'])
    ax.axvline(x=0, color='black', linewidth=1)
    ax.set_xlabel('Log2 Fold Change (Cluster vs Single)', fontsize=12)
    ax.set_ylabel('Genes', fontsize=12)
    ax.set_title('Neutrophil-related Genes: Differential Expression\n(CTC Clusters vs Single CTCs)', 
                 fontsize=14, fontweight='bold')
    
    # Add significance markers
    for i, (_, row) in enumerate(results_df.iterrows()):
        if row['Significant']:
            x_pos = row['Log2FC'] + 0.1 if row['Log2FC'] > 0 else row['Log2FC'] - 0.3
            ax.text(x_pos, i, '*', fontsize=12, va='center', fontweight='bold')
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#E74C3C', edgecolor='black', label='Upregulated in Clusters (p<0.05)'),
        Patch(facecolor='#3498DB', edgecolor='black', label='Downregulated in Clusters (p<0.05)'),
        Patch(facecolor='#95A5A6', edgecolor='black', label='Not Significant (p≥0.05)')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
    
    # Set x-axis limits symmetrically
    max_fc = max(abs(results_df['Log2FC'].max()), abs(results_df['Log2FC'].min())) + 1
    ax.set_xlim(-max_fc, max_fc)
    
    plt.tight_layout()
    plt.savefig('/root/tumor_cell/analysis/chemokine_foldchange.png', dpi=300, bbox_inches='tight', facecolor='white')
    print("[Saved] Fold change plot to: /root/tumor_cell/analysis/chemokine_foldchange.png")
    plt.close()
    
    # ============================================
    # 8. Print Summary
    # ============================================
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    sig_up = results_df[(results_df['Significant']) & (results_df['Log2FC'] > 0)]
    sig_down = results_df[(results_df['Significant']) & (results_df['Log2FC'] < 0)]
    
    print(f"\nTotal genes analyzed: {len(available_chemokines)}")
    print(f"\n🔺 Significantly UPREGULATED in CTC clusters: {len(sig_up)}")
    if len(sig_up) > 0:
        for _, row in sig_up.iterrows():
            print(f"   • {row['Gene']}: Log2FC = {row['Log2FC']}, p = {row['P_value']}")
    
    print(f"\n🔻 Significantly DOWNREGULATED in CTC clusters: {len(sig_down)}")
    if len(sig_down) > 0:
        for _, row in sig_down.iterrows():
            print(f"   • {row['Gene']}: Log2FC = {row['Log2FC']}, p = {row['P_value']}")
    
    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print("Output files:")
    print("  1. /root/tumor_cell/analysis/chemokine_heatmap.png")
    print("  2. /root/tumor_cell/analysis/chemokine_heatmap.pdf")
    print("  3. /root/tumor_cell/analysis/chemokine_foldchange.png")
    print("  4. /root/tumor_cell/analysis/chemokine_DE_results.csv")
    print("=" * 60)

else:
    print("\nERROR: No chemokines found in the dataset!")
    print("Checking available gene names that contain 'CXC' or 'CCL'...")
    matching = [g for g in df.index if 'CXC' in str(g).upper() or 'CCL' in str(g).upper()]
    print(matching[:20])
