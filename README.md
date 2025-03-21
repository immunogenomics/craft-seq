# craft-seq
This repo provides code for the analysis of CRAFT-seq data corresponding to [Baglaenko, et al, biorxiv, 2024](https://www.biorxiv.org/content/10.1101/2024.03.28.587175v1).

Contact: Yuriy Baglaenko yuriy.baglaenko@cchmc.org, Michelle Curtis curtism@broadinstitute.org

Contained in this repo are figures/analysis notebooks or scripts for each of the following experiments.

| Experiment  |  Figure # | Notebooks/Scripts | Notes |
|--------|----------|----------|----------|
| PTEN  | Figure 1 | [PTEN Notebook.ipynb](./figures/PTEN/PTEN_FBXO11_DQB1_Github.ipynb) | Joint Notebook with DQB1/FBXO11. No CRISPR editing |
| FBXO11 | Figure 2 | [FBXO11_Run 1 Notebook.ipynb](./figures/FBXO11/PTEN_FBXO11_DQB1_Github.ipynb) <br> [FBXO11_Run 2 - Joint Analysis.ipynb](./figures/FBXO11/FBXO11_Run2_Joint_Github.ipynb)| Run 1 and Run 2 merged and analyzed together. Alternative splicing to induce a knockout with CRISPR |
| DQB1 | Figure 2 | [DQB1_Notebook.ipynb](./figures/DQB1/PTEN_FBXO11_DQB1_Github.ipynb) <br> [DQB1_Heatmap.ipynb](./figures/DQB1/DQB1_Heatmap_Github.ipynb)| CRISPR-Cas cutting and HDR in HLA-DQB1. |
| PTPRC | Figure 3 | [PTPRC_Notebook.ipynb](./figures/PTPRC/PTPRC_Github.ipynb) <br> [GSEA_PTPRC_IL2RA.ipynb](./figures/Misc/GSEA_PTPRC_IL2RA.ipynb) <br> [Bootstrapping CD45 Data.ipynb](./figures/PTPRC/PTPRC_Analysis_Bootstrap_Laters.ipynb)| Induction of ESC with CRISPR base editors|
| RPL8 | Figure 4 | [RPL8_Notebook.ipynb](./figures/RPL8/RPL8_Github.ipynb)| RPL8 eQTL CRISPR base editing | 
| IL2RA | Figure 4 | [IL2RA_Run 1.ipynb](./figures/IL2RA/IL2RA_Run1_Github.ipynb) <br> [IL2RA_Run 2.ipynb](./figures/IL2RA/IL2RA_Run2_Github.ipynb) <br> [IL2RA_Joint Analysis Notebook 2.ipynb](./figures/IL2RA/IL2RA_joint_Github.ipynb) <br> [GSEA_PTPRC_IL2RA.ipynb](./figures/Misc/GSEA_PTPRC_IL2RA.ipynb)| CRISPR editing of IL2RA variant in Th1 and Treg polarized naive CD4 T cells
| PAX5 | Figure 5 | [PAX5_genotype_analysis.R](./figures/PAX5/PAX5_genotype_analysis.R) <br> [PAX5_RNA_analysis.R](./figures/PAX5/PAX5_RNA_analysis.R) <br> [PAX5_RNA_glmNB.R](./figures/PAX5/PAX5_RNA_glmNB.R) <br> [PAX5_ADT_glmNB.R](./figures/PAX5/PAX5_ADT_glmNB.R) <br>  [PAX5_genotyping.ipynb](./figures/PAX5/PAX5_AlleleCalling_Github.ipynb)| PAX5 multiplexed editing with CRIPSR base editors |
| Misc |  | [Mixscape_Comparison_Latest.ipynb](./figures/Misc/Mixscape_Comparison_Latest.ipynb) <br> [TotalAlignmentStats_Latest.ipynb](./figures/Misc/TotalAlignmentStats_Latest.ipynb) <br> [Example Flow Cytometry Analysis.ipynb](./figures/Misc/PTEN_FBXO11_DQB1_Flow_Github.ipynb) | Extras |
| Supplemental Note | |[Supplemental Note and Genotype Calling.ipynb](./figures/SupplementaryNote/SupplementaryNote_Github.ipynb) | Fully replicatable genotyping analysis. Key files deposited in folder | 

Sample data for the PTEN, FBXO11, and DQB1 experiments is located on [Zenodo](https://zenodo.org/records/10932681). GEO records are being generated to deposit count matrices with trimmed meta data. 

File: [DNA Filtering Functions](./figures/DNA_filtering_Functions.R) has updated function on allele calling and genotyping of single cell along with extra utility functions.
