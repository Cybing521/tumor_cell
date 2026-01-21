#!/usr/bin/env python3
"""
CTC Cluster vs Single CTC - Neutrophil-related Chemokine Heatmap Analysis
Based on GSE111065 (Gkountela et al., 2019, Cell)
"Circulating Tumor Cell Clustering Shapes DNA Methylation to Enable Metastasis Seeding"
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
# 1. Load Expression Data
# ============================================
print("=" * 60)
print("Loading GSE111065 data (Gkountela et al., 2019)...")
print("=" * 60)

# Read the expression matrix
data_path = "/root/tumor_cell/data/GSE111065_matrix.txt"
df = pd.read_csv(data_path, sep='\t', index_col=0)
print(f"Expression data shape: {df.shape}")

# ============================================
# 2. Load Sample Metadata and Identify Groups
# ============================================
print("\n" + "=" * 60)
print("Loading sample metadata...")
print("=" * 60)

# Read series matrix to get sample types
# Parse the sample type information manually from the data we saw
# Line 7 contains: sample type: CTC-single, CTC-cluster, etc.

# Let's create sample classification based on the matrix file
# We'll read the full series matrix to extract sample types
series_matrix_path = "/root/tumor_cell/data/GSE111065_series_matrix.txt"

# Read the series matrix to find sample types
with open(series_matrix_path, 'r') as f:
    lines = f.readlines()

# Find the sample type line
sample_titles = []
sample_types = []
for line in lines:
    if line.startswith('!Sample_title'):
        parts = line.strip().split('\t')
        sample_titles = [p.strip('"') for p in parts[1:]]
    if 'sample type:' in line:
        parts = line.strip().split('\t')
        sample_types = [p.split(': ')[1].strip('"') if ': ' in p else 'NA' for p in parts[1:]]
        break

# Create sample info dataframe
sample_info = pd.DataFrame({
    'sample_title': sample_titles,
    'sample_type': sample_types
})

# Map sample titles to column names in expression matrix
# The expression matrix columns should match sample titles
print(f"\nTotal samples in metadata: {len(sample_info)}")
print(f"Sample types found:")
print(sample_info['sample_type'].value_counts())

# ============================================
# 3. Filter for CTC samples only (exclude cell lines)
# ============================================
print("\n" + "=" * 60)
print("Filtering for CTC samples...")
print("=" * 60)

# Identify CTC clusters and single CTCs
cluster_mask = sample_info['sample_type'].str.contains('CTC-cluster', na=False)
single_mask = sample_info['sample_type'] == 'CTC-single'

cluster_samples = sample_info[cluster_mask]['sample_title'].tolist()
single_samples = sample_info[single_mask]['sample_title'].tolist()

# Match to expression matrix columns
# Some column names might have slight differences
expr_columns = df.columns.tolist()

# Find matching columns
def find_matching_columns(sample_list, expr_cols):
    matched = []
    for sample in sample_list:
        # Try exact match first
        if sample in expr_cols:
            matched.append(sample)
        else:
            # Try partial match
            for col in expr_cols:
                if sample.replace(' #2', '') in col or col in sample:
                    matched.append(col)
                    break
    return list(set(matched))

cluster_cols = find_matching_columns(cluster_samples, expr_columns)
single_cols = find_matching_columns(single_samples, expr_columns)

print(f"\nCTC Clusters samples found in expression data: {len(cluster_cols)}")
print(f"Single CTC samples found in expression data: {len(single_cols)}")

if len(cluster_cols) == 0 or len(single_cols) == 0:
    print("\nWarning: Sample matching issue. Trying alternative matching...")
    # Alternative: use column names directly
    cluster_cols = [col for col in expr_columns if any(s.replace(' #2', '') == col for s in cluster_samples)]
    single_cols = [col for col in expr_columns if any(s.replace(' #2', '') == col for s in single_samples)]
    
    if len(cluster_cols) == 0:
        # Direct match from sample type
        cluster_cols = []
        single_cols = []
        for i, (title, stype) in enumerate(zip(sample_titles, sample_types)):
            clean_title = title.replace(' #2', '')
            if clean_title in expr_columns:
                if 'CTC-cluster' in stype:
                    cluster_cols.append(clean_title)
                elif stype == 'CTC-single':
                    single_cols.append(clean_title)
    
    print(f"After alternative matching:")
    print(f"  CTC Clusters: {len(cluster_cols)}")
    print(f"  Single CTCs: {len(single_cols)}")

# ============================================
# 4. Define Neutrophil-related Chemokines
# ============================================
print("\n" + "=" * 60)
print("Defining neutrophil-related chemokines...")
print("=" * 60)

neutrophil_chemokines = [
    # CXC chemokines (primary neutrophil attractants)
    'CXCL1', 'CXCL2', 'CXCL3', 'CXCL5', 'CXCL6', 'CXCL7', 'CXCL8',
    'CXCL9', 'CXCL10', 'CXCL11', 'CXCL12', 'CXCL13', 'CXCL14', 'CXCL16',
    # CC chemokines
    'CCL2', 'CCL3', 'CCL4', 'CCL5', 'CCL7', 'CCL8', 'CCL11', 'CCL20',
    # Chemokine receptors
    'CXCR1', 'CXCR2', 'CXCR3', 'CXCR4', 'CCR1', 'CCR2', 'CCR5', 'CCR7',
    # IL8 alternative
    'IL8',
    # Inflammatory factors
    'IL6', 'IL1B', 'IL1A', 'TNF', 'TNFAIP3', 'CSF3', 'CSF2',
    'S100A8', 'S100A9', 'S100A12',
    # Adhesion molecules
    'ICAM1', 'VCAM1', 'SELE', 'SELP', 'ITGB2', 'ITGAM',
    # Other neutrophil-related
    'FPR1', 'FPR2', 'C5AR1', 'PTGS2', 'MMP9', 'MMP2',
    # Key genes from Gkountela study
    'JUP', 'NANOG', 'OCT4', 'SOX2'
]

# Find available chemokines
available_chemokines = [gene for gene in neutrophil_chemokines if gene in df.index]
missing = [gene for gene in neutrophil_chemokines if gene not in df.index]

print(f"\nAvailable genes in dataset: {len(available_chemokines)}")
print(available_chemokines)
print(f"\nMissing genes: {len(missing)}")

# ============================================
# 5. Extract and Analyze Data
# ============================================
if len(available_chemokines) > 0 and len(cluster_cols) > 0 and len(single_cols) > 0:
    print("\n" + "=" * 60)
    print("Performing differential expression analysis...")
    print("=" * 60)
    
    # Extract chemokine data
    chemokine_df = df.loc[available_chemokines, cluster_cols + single_cols].copy()
    
    # Data is already log-normalized based on the values we saw
    print(f"\nChemokine expression matrix shape: {chemokine_df.shape}")
    
    # Calculate differential expression
    results = []
    for gene in available_chemokines:
        cluster_expr = chemokine_df.loc[gene, cluster_cols].values.astype(float)
        single_expr = chemokine_df.loc[gene, single_cols].values.astype(float)
        
        # Remove NaN values
        cluster_expr = cluster_expr[~np.isnan(cluster_expr)]
        single_expr = single_expr[~np.isnan(single_expr)]
        
        if len(cluster_expr) > 0 and len(single_expr) > 0:
            cluster_mean = np.mean(cluster_expr)
            single_mean = np.mean(single_expr)
            fold_change = cluster_mean - single_mean  # Already log-scale
            
            # T-test
            if len(cluster_expr) > 1 and len(single_expr) > 1:
                t_stat, p_value = stats.ttest_ind(cluster_expr, single_expr, nan_policy='omit')
            else:
                p_value = 1.0
            
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
    results_df.to_csv('/root/tumor_cell/analysis/gkountela_chemokine_DE_results.csv', index=False)
    print("\n[Saved] Results to: /root/tumor_cell/analysis/gkountela_chemokine_DE_results.csv")
    
    # ============================================
    # 6. Create Heatmap
    # ============================================
    print("\n" + "=" * 60)
    print("Creating heatmap...")
    print("=" * 60)
    
    # Prepare data
    ordered_samples = cluster_cols + single_cols
    heatmap_data = chemokine_df[ordered_samples].astype(float)
    
    # Z-score normalization
    row_means = heatmap_data.mean(axis=1)
    row_stds = heatmap_data.std(axis=1)
    row_stds = row_stds.replace(0, 1)
    heatmap_zscore = heatmap_data.subtract(row_means, axis=0).divide(row_stds, axis=0)
    
    # Sort genes by fold change
    gene_order = results_df['Gene'].tolist()
    heatmap_zscore = heatmap_zscore.loc[gene_order]
    
    # Create figure
    fig_height = max(10, len(available_chemokines) * 0.4)
    fig, ax = plt.subplots(figsize=(20, fig_height))
    
    # Heatmap
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
    
    # Add vertical line
    ax.axvline(x=len(cluster_cols), color='black', linewidth=2)
    
    # Title and labels
    ax.set_title('Neutrophil-related Chemokines: CTC Clusters vs Single CTCs\n(GSE111065, Gkountela et al. 2019, Cell)', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel('Samples', fontsize=12)
    ax.set_ylabel('Genes', fontsize=12)
    
    plt.xticks(rotation=90, fontsize=6)
    plt.yticks(fontsize=10)
    
    # Group labels
    ax.text(len(cluster_cols)/2, -0.5, f'CTC Clusters\n(n={len(cluster_cols)})', 
            ha='center', va='bottom', fontsize=10, fontweight='bold', color='#E74C3C')
    ax.text(len(cluster_cols) + len(single_cols)/2, -0.5, f'Single CTCs\n(n={len(single_cols)})', 
            ha='center', va='bottom', fontsize=10, fontweight='bold', color='#3498DB')
    
    plt.tight_layout()
    
    # Save
    plt.savefig('/root/tumor_cell/analysis/gkountela_chemokine_heatmap.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('/root/tumor_cell/analysis/gkountela_chemokine_heatmap.pdf', bbox_inches='tight', facecolor='white')
    print("[Saved] Heatmap to: /root/tumor_cell/analysis/gkountela_chemokine_heatmap.png")
    plt.close()
    
    # ============================================
    # 7. Create Fold Change Plot
    # ============================================
    fig, ax = plt.subplots(figsize=(10, max(6, len(available_chemokines) * 0.35)))
    
    colors = []
    for _, row in results_df.iterrows():
        if row['Significant'] and row['Log2FC'] > 0:
            colors.append('#E74C3C')
        elif row['Significant'] and row['Log2FC'] < 0:
            colors.append('#3498DB')
        else:
            colors.append('#95A5A6')
    
    y_pos = range(len(results_df))
    ax.barh(y_pos, results_df['Log2FC'], color=colors, edgecolor='black', linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(results_df['Gene'])
    ax.axvline(x=0, color='black', linewidth=1)
    ax.set_xlabel('Log2 Fold Change (Cluster vs Single)', fontsize=12)
    ax.set_ylabel('Genes', fontsize=12)
    ax.set_title('Neutrophil-related Genes: Differential Expression\n(Gkountela et al. 2019)', fontsize=14, fontweight='bold')
    
    # Significance markers
    for i, (_, row) in enumerate(results_df.iterrows()):
        if row['Significant']:
            x_pos = row['Log2FC'] + 0.1 if row['Log2FC'] > 0 else row['Log2FC'] - 0.3
            ax.text(x_pos, i, '*', fontsize=12, va='center', fontweight='bold')
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#E74C3C', edgecolor='black', label='Upregulated in Clusters (p<0.05)'),
        Patch(facecolor='#3498DB', edgecolor='black', label='Downregulated in Clusters (p<0.05)'),
        Patch(facecolor='#95A5A6', edgecolor='black', label='Not Significant')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('/root/tumor_cell/analysis/gkountela_chemokine_foldchange.png', dpi=300, bbox_inches='tight', facecolor='white')
    print("[Saved] Fold change plot to: /root/tumor_cell/analysis/gkountela_chemokine_foldchange.png")
    plt.close()
    
    # ============================================
    # 8. Summary
    # ============================================
    print("\n" + "=" * 60)
    print("SUMMARY - Gkountela 2019 (GSE111065)")
    print("=" * 60)
    
    sig_up = results_df[(results_df['Significant']) & (results_df['Log2FC'] > 0)]
    sig_down = results_df[(results_df['Significant']) & (results_df['Log2FC'] < 0)]
    
    print(f"\nTotal genes analyzed: {len(available_chemokines)}")
    print(f"CTC Clusters samples: {len(cluster_cols)}")
    print(f"Single CTC samples: {len(single_cols)}")
    
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
    print("=" * 60)

else:
    print("\nERROR: Could not find enough data for analysis!")
    print(f"Available chemokines: {len(available_chemokines)}")
    print(f"Cluster samples: {len(cluster_cols)}")
    print(f"Single samples: {len(single_cols)}")
