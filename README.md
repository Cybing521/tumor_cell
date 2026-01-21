# CTC Clusters vs Single CTCs: Neutrophil-related Chemokine Analysis

## 项目简介

本项目分析循环肿瘤细胞簇（CTC Clusters）与单个循环肿瘤细胞（Single CTCs）之间中性粒细胞相关趋化因子的差异表达，探索中性粒细胞在CTC簇形成和肿瘤转移中的潜在作用。

## 研究背景

循环肿瘤细胞簇（CTC clusters）相比单个循环肿瘤细胞具有更高的转移潜能（约23-50倍）。近年研究发现，中性粒细胞可能通过与CTC形成复合体来促进肿瘤转移。本项目旨在通过RNA-seq数据分析，鉴定CTC簇中上调的中性粒细胞相关趋化因子。

## 数据来源

### GSE111065 (Gkountela et al., 2019, Cell)
- **论文**: "Circulating Tumor Cell Clustering Shapes DNA Methylation to Enable Metastasis Seeding"
- **样本**: 乳腺癌患者外周血
- **CTC Clusters**: 45个样本
- **Single CTCs**: 56个样本
- **数据类型**: 单细胞RNA测序

## 主要发现

### 在CTC簇中显著上调的中性粒细胞相关基因 (P < 0.05)

| 基因 | Log2FC | P值 | 功能 |
|------|--------|-----|------|
| **CXCL16** | 1.886 | 0.0013 | 趋化因子，招募中性粒细胞 |
| **TNF** | 1.400 | 0.0056 | 促炎因子，激活中性粒细胞 |
| **JUP** | 1.332 | 0.0277 | 连接蛋白，CTC簇形成关键分子 |
| **ITGB2** | 1.141 | 0.0096 | 整合素β2(CD18)，中性粒细胞关键黏附分子 |
| **SOX2** | 0.782 | 0.0449 | 干细胞转录因子 |

### 生物学意义

1. **CXCL16** 和 **ITGB2** 的显著上调支持"中性粒细胞护送CTC"的假说
2. **TNF** 等炎症因子的上调提示CTC簇具有更强的促炎微环境
3. **JUP** 的上调与CTC簇细胞间黏附机制一致

## 项目结构

```
tumor_cell/
├── README.md                    # 项目说明文档
├── analysis/                    # 分析脚本和结果
│   ├── gkountela_chemokine_heatmap.py    # 主分析脚本
│   ├── gkountela_chemokine_heatmap.png   # 热图PNG
│   ├── gkountela_chemokine_heatmap.pdf   # 热图PDF
│   ├── gkountela_chemokine_foldchange.png # 差异倍数图
│   └── gkountela_chemokine_DE_results.csv # 差异表达结果
├── data/                        # 原始数据
│   ├── GSE111065_matrix.txt     # 表达矩阵
│   └── GSE111065_series_matrix.txt # 样本信息
├── deliverables/                # 交付物
│   └── CTC_neutrophil_chemokine_heatmap/
│       ├── CTC_cluster_vs_single_chemokine_heatmap.pdf
│       ├── CTC_cluster_vs_single_foldchange.png
│       ├── differential_expression_results.csv
│       ├── source_data_GSE111065_expression_matrix.txt
│       └── README_热图说明.txt
└── material/                    # 参考文献PDF
    ├── PIIS0092867414009271.pdf # Aceto 2014
    └── gkountela2019.pdf        # Gkountela 2019
```

## 分析方法

1. **数据预处理**: 使用GEO数据库提供的标准化表达矩阵
2. **差异表达分析**: 
   - 计算Log2 Fold Change
   - 独立样本t检验
   - 显著性阈值: P < 0.05
3. **可视化**: 
   - Z-score标准化热图
   - Log2FC条形图

## 依赖环境

```python
pandas>=1.3.0
numpy>=1.21.0
scipy>=1.7.0
seaborn>=0.11.0
matplotlib>=3.4.0
```

## 使用方法

```bash
# 运行分析脚本
cd analysis
python3 gkountela_chemokine_heatmap.py
```

## 参考文献

1. Gkountela S, et al. (2019) Circulating Tumor Cell Clustering Shapes DNA Methylation to Enable Metastasis Seeding. *Cell*, 176(1-2):98-112.

2. Aceto N, et al. (2014) Circulating Tumor Cell Clusters Are Oligoclonal Precursors of Breast Cancer Metastasis. *Cell*, 158(5):1110-1122.

3. Szczerba BM, et al. (2019) Neutrophils Escort Circulating Tumour Cells to Enable Cell Cycle Progression. *Nature*, 566(7745):553-557.

## 作者

Tumor Cell Analysis Project

## 许可证

MIT License
