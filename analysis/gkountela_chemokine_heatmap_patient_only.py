#!/usr/bin/env python3
"""
CTC Cluster vs Single CTC - Neutrophil-related Chemokine Heatmap Analysis
Based on GSE111065 (Gkountela et al., 2019, Cell)
*** PATIENT SAMPLES ONLY - 排除动物模型和细胞系 ***
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
print("*** PATIENT SAMPLES ONLY ***")
print("=" * 60)

# Read the expression matrix
data_path = "/root/tumor_cell/data/GSE111065_matrix.txt"
df = pd.read_csv(data_path, sep='\t', index_col=0)
print(f"Expression data shape: {df.shape}")

# ============================================
# 2. Load Sample Metadata and Filter for PATIENT ONLY
# ============================================
print("\n" + "=" * 60)
print("Loading sample metadata and filtering for PATIENT samples...")
print("=" * 60)

series_matrix_path = "/root/tumor_cell/data/GSE111065_series_matrix.txt"

# Read the series matrix to find sample types and origins
with open(series_matrix_path, 'r') as f:
    lines = f.readlines()

# Find the sample info
sample_titles = []
sample_types = []
sample_origins = []

for line in lines:
    if line.startswith('!Sample_title'):
        parts = line.strip().split('\t')
        sample_titles = [p.strip('"') for p in parts[1:]]
    if line.startswith('!Sample_characteristics_ch1') and 'origin:' in line:
        parts = line.strip().split('\t')
        sample_origins = [p.split(': ')[1].strip('"') if ': ' in p else 'NA' for p in parts[1:]]
    if line.startswith('!Sample_characteristics_ch1') and 'sample type:' in line:
        parts = line.strip().split('\t')
        sample_types = [p.split(': ')[1].strip('"') if ': ' in p else 'NA' for p in parts[1:]]

# Create sample info dataframe
sample_info = pd.DataFrame({
    'sample_title': sample_titles,
    'sample_type': sample_types,
    'origin': sample_origins
})

print(f"\nTotal samples in metadata: {len(sample_info)}")
print("\nSample origins distribution:")
print(sample_info['origin'].value_counts())
print("\nSample types distribution:")
print(sample_info['sample_type'].value_counts())

# ============================================
# 3. Filter for PATIENT CTCs only (exclude cell lines and xenografts)
# ============================================
print("\n" + "=" * 60)
print("Filtering for PATIENT CTC samples ONLY...")
print("=" * 60)

# Only patient origin
patient_samples = sample_info[sample_info['origin'] == 'patient'].copy()
print(f"\nPatient-derived samples: {len(patient_samples)}")
print(patient_samples['sample_type'].value_counts())

# Identify CTC clusters and single CTCs from PATIENT samples only
cluster_mask = patient_samples['sample_type'].str.contains('CTC-cluster', na=False)
single_mask = patient_samples['sample_type'] == 'CTC-single'

cluster_samples = patient_samples[cluster_mask]['sample_title'].tolist()
single_samples = patient_samples[single_mask]['sample_title'].tolist()

print(f"\nPatient CTC Clusters: {len(cluster_samples)}")
print(f"Patient Single CTCs: {len(single_samples)}")

# Match to expression matrix columns
expr_columns = df.columns.tolist()

def find_matching_columns(sample_list, expr_cols):
    matched = []
    for sample in sample_list:
        # Try exact match first
        if sample in expr_cols:
            matched.append(sample)
        else:
            # Try partial match (remove #2 suffix)
            clean_sample = sample.replace(' #2', '')
            if clean_sample in expr_cols:
                matched.append(clean_sample)
            else:
                for col in expr_cols:
                    if clean_sample in col or col in clean_sample:
                        matched.append(col)
                        break
    return list(set(matched))

cluster_cols = find_matching_columns(cluster_samples, expr_columns)
single_cols = find_matching_columns(single_samples, expr_columns)

print(f"\nMatched in expression data:")
print(f"  CTC Clusters (Patient): {len(cluster_cols)}")
print(f"  Single CTCs (Patient): {len(single_cols)}")

# Show matched samples
print(f"\nCluster samples: {cluster_cols[:5]}..." if len(cluster_cols) > 5 else f"\nCluster samples: {cluster_cols}")
print(f"Single samples: {single_cols[:5]}..." if len(single_cols) > 5 else f"Single samples: {single_cols}")

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
print(f"Missing genes: {len(missing)}")

# ============================================
# 5. Extract and Analyze Data
# ============================================
if len(available_chemokines) > 0 and len(cluster_cols) > 0 and len(single_cols) > 0:
    print("\n" + "=" * 60)
    print("Performing differential expression analysis (PATIENT ONLY)...")
    print("=" * 60)
    
    # Extract chemokine data
    chemokine_df = df.loc[available_chemokines, cluster_cols + single_cols].copy()
    
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
    results_df.to_csv('/root/tumor_cell/analysis/gkountela_chemokine_DE_results_patient_only.csv', index=False)
    print("\n[Saved] Results to: gkountela_chemokine_DE_results_patient_only.csv")
    
    # ============================================
    # 6. Create Heatmap
    # ============================================
    print("\n" + "=" * 60)
    print("Creating heatmap (PATIENT ONLY)...")
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
    fig, ax = plt.subplots(figsize=(16, fig_height))
    
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
    
    # Mark significant genes with asterisks
    for i, gene in enumerate(gene_order):
        if results_df[results_df['Gene'] == gene]['Significant'].values[0]:
            ax.text(-0.5, i + 0.5, '*', fontsize=14, ha='right', va='center', fontweight='bold', color='red')
    
    # Title and labels
    ax.set_title('Neutrophil-related Chemokines: CTC Clusters vs Single CTCs\n(GSE111065, Gkountela et al. 2019 - PATIENT SAMPLES ONLY)', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel('Samples (Patient-derived CTCs)', fontsize=12)
    ax.set_ylabel('Genes (* = p<0.05)', fontsize=12)
    
    plt.xticks(rotation=90, fontsize=7)
    plt.yticks(fontsize=10)
    
    # Group labels
    ax.text(len(cluster_cols)/2, -0.5, f'CTC Clusters\n(n={len(cluster_cols)})', 
            ha='center', va='bottom', fontsize=10, fontweight='bold', color='#E74C3C')
    ax.text(len(cluster_cols) + len(single_cols)/2, -0.5, f'Single CTCs\n(n={len(single_cols)})', 
            ha='center', va='bottom', fontsize=10, fontweight='bold', color='#3498DB')
    
    plt.tight_layout()
    
    # Save
    plt.savefig('/root/tumor_cell/analysis/gkountela_chemokine_heatmap_patient_only.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('/root/tumor_cell/analysis/gkountela_chemokine_heatmap_patient_only.pdf', bbox_inches='tight', facecolor='white')
    print("[Saved] Heatmap to: gkountela_chemokine_heatmap_patient_only.png/pdf")
    plt.close()
    
    # ============================================
    # 7. Create Fold Change Plot
    # ============================================
    fig, ax = plt.subplots(figsize=(10, max(6, len(available_chemokines) * 0.35)))
    
    colors = []
    for _, row in results_df.iterrows():
        if row['Significant'] and row['Log2FC'] > 0:
            colors.append('#E74C3C')  # Red for significant up
        elif row['Significant'] and row['Log2FC'] < 0:
            colors.append('#3498DB')  # Blue for significant down
        else:
            colors.append('#95A5A6')  # Gray for not significant
    
    y_pos = range(len(results_df))
    ax.barh(y_pos, results_df['Log2FC'], color=colors, edgecolor='black', linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(results_df['Gene'])
    ax.axvline(x=0, color='black', linewidth=1)
    ax.set_xlabel('Log2 Fold Change (Cluster vs Single)', fontsize=12)
    ax.set_ylabel('Genes', fontsize=12)
    ax.set_title('Neutrophil-related Genes: Differential Expression\n(Gkountela 2019 - PATIENT SAMPLES ONLY)', fontsize=14, fontweight='bold')
    
    # Significance markers
    for i, (_, row) in enumerate(results_df.iterrows()):
        if row['Significant']:
            x_pos = row['Log2FC'] + 0.1 if row['Log2FC'] > 0 else row['Log2FC'] - 0.3
            ax.text(x_pos, i, '*', fontsize=12, va='center', fontweight='bold')
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#E74C3C', edgecolor='black', label='Upregulated (p<0.05)'),
        Patch(facecolor='#3498DB', edgecolor='black', label='Downregulated (p<0.05)'),
        Patch(facecolor='#95A5A6', edgecolor='black', label='Not Significant')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('/root/tumor_cell/analysis/gkountela_chemokine_foldchange_patient_only.png', dpi=300, bbox_inches='tight', facecolor='white')
    print("[Saved] Fold change plot to: gkountela_chemokine_foldchange_patient_only.png")
    plt.close()
    
    # ============================================
    # 8. Summary
    # ============================================
    print("\n" + "=" * 60)
    print("SUMMARY - Gkountela 2019 (GSE111065) - PATIENT ONLY")
    print("=" * 60)
    
    sig_up = results_df[(results_df['Significant']) & (results_df['Log2FC'] > 0)]
    sig_down = results_df[(results_df['Significant']) & (results_df['Log2FC'] < 0)]
    
    print(f"\nTotal genes analyzed: {len(available_chemokines)}")
    print(f"Patient CTC Clusters samples: {len(cluster_cols)}")
    print(f"Patient Single CTC samples: {len(single_cols)}")
    
    print(f"\n🔺 Significantly UPREGULATED in CTC clusters (p<0.05): {len(sig_up)}")
    if len(sig_up) > 0:
        for _, row in sig_up.iterrows():
            print(f"   • {row['Gene']}: Log2FC = {row['Log2FC']}, p = {row['P_value']}")
    
    print(f"\n🔻 Significantly DOWNREGULATED in CTC clusters (p<0.05): {len(sig_down)}")
    if len(sig_down) > 0:
        for _, row in sig_down.iterrows():
            print(f"   • {row['Gene']}: Log2FC = {row['Log2FC']}, p = {row['P_value']}")
    
    print("\n" + "=" * 60)
    print("Analysis Complete! (Patient samples only)")
    print("=" * 60)

else:
    print("\nERROR: Could not find enough data for analysis!")
    print(f"Available chemokines: {len(available_chemokines)}")
    print(f"Cluster samples: {len(cluster_cols)}")
    print(f"Single samples: {len(single_cols)}")
