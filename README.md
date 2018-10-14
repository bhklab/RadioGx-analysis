# RadioGx Analysis

Scrips to reproduce figures. 

Manuscript: Modeling cellular response in large-scale radiogenomic databases to advance precision radiotherapy
Venkata SK. Manem, Meghan Lambie, Petr Smirnov1, Victor Kofia, Mark Freeman, Marianne Koritzinsky, Mohamed E. Abazeed, Benjamin Haibe-Kains, Scott V. Bratman

Figure2-Code.R: Fitting of dose-response data to the LQ model and concordance of radiation response across asays. (A) LQ model fit using RadioGx on the SNU-245 cholangiocarcinoma cell line (red) and SK-ES-1 Ewing sarcoma cell line (blue). The LQ model describes the fraction of cells predicted to survive (y-axis) a uniform radiation dose (x-axis) and is characterized by  and  components for each cell line. For SNU-245 and SK-ES-1, =0.14 ( Gy-1), (Gy-2)=0 and =0.45 (Gy-1),=0.02 (Gy-2), respectively. Solid curves indicate the model fit and points denote experimental data (Yard et al., 2016). (B) Histogram of AUC values calculated using the computeAUC function in RadioGx. (C) Correlation (Pearson R with standard deviation) of radioresponse results produced by the 9-day viability assay and the standard clonogenic assay according to the following indicators: SF2, SF4, SF6, SF8, and AUC. Primary data were obtained from Yard et al.(Yard et al., 2016). 

Figure3-Code.R: Concordance of SF2 and AUC. (A) Correlation between the radioresponse indicators, SF2 and AUC, across 498 cell lines. (B) Pearson correlation (with standard deviation) between SF2 and AUC across 498 cell lines based on tertiles. (C) Venn diagram illustrating the transcriptional pathways associated with radioresponse using SF2 or AUC as the response indicator. (D) False discovery rate (FDR) for each transcriptional pathway from (C) illustrating greater levels of statistical significance among pathways specific to AUC. 

Figure4-Code.R: Distinct biological underpinnings of alpha/beta derived from the LQ model. (A) Histogram of alpha/beta (Gy) values obtained from the LQ model across all cell lines. (B) Pearson correlations (with standard deviation) between AUC and the alpha and beta components of the LQ model. (C) Transcriptional pathways that are significantly associated with AUC, alpha, and/or beta.

Figure5-Code.R: Integrative analysis of radiobiological model with transcriptomic data and prediction of radioresponse under hypoxia. (A) Hypothetical illustration of cancer cell surviving fraction according to dose and oxygen partial pressure, as modeled using RadioGx. Solid curves are modeled using Equation (3) (Methods). The computed AUC values are 2.41, 2.71, 2.97 for normoxia (160 mmHg), hypoxic condition with 10 mmHg, and hypoxic condition with 5 mmHg, respectively. (B) Changes in AUC by tissue type (with minimum of 15 cell lines within RadioGx), ranked based on median AUC under normoxic conditions. Hypoxic condition: 5 mmHg. (C) Difference in ranks between the oxic and hypoxic conditions across all genes using a univariate modeling approach. 

Figure6-Code.R: Tissue specificity of molecular determinants of radioresponse. (A) The tumour types (columns) represented by a minimum of 15 cell lines were considered for analysis. A total of 281 pathways are depicted (rows) and are annotated by function. Colours designate pathways significantly associated with AUC (FDR < 5%). (B) Heterogeneity of alpha/beta ratios across various tumour types.

Figure7-Code.R: Identification of drugs and pharmacological classes with cytotoxic effects on cancer cell lines that correlate with radioresponse. Pharmacological enrichment analysis using radiation AUC as the radioresponse indicator. Pharmacological classes with statistically significant associations with radioresponse in cancer cell lines are indicated.

##############################################################################################################

> sessionInfo()
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggplot2_3.0.0        gplots_3.0.1         RColorBrewer_1.1-2   reshape_0.8.7        bindrcpp_0.2.2      
 [6] dplyr_0.7.6          VennDiagram_1.6.20   futile.logger_1.4.3  piano_1.20.0         readxl_1.1.0        
[11] RadioGx_0.0.0.0.9000 CoreGx_0.0.001       caTools_1.17.1.1     gdata_2.18.0        

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.18         lattice_0.20-35      relations_0.6-8      gtools_3.8.1         assertthat_0.2.0    
 [6] digest_0.6.16        slam_0.1-43          lsa_0.73.1           R6_2.2.2             cellranger_1.1.0    
[11] plyr_1.8.4           futile.options_1.0.1 pillar_1.3.0         rlang_0.2.2          lazyeval_0.2.1      
[16] rstudioapi_0.7       data.table_1.11.4    Matrix_1.2-14        labeling_0.3         sets_1.0-18         
[21] BiocParallel_1.14.2  stringr_1.3.1        igraph_1.2.2         munsell_0.5.0        fgsea_1.6.0         
[26] compiler_3.5.0       pkgconfig_2.0.2      BiocGenerics_0.26.0  marray_1.58.0        tidyselect_0.2.4    
[31] tibble_1.4.2         gridExtra_2.3        crayon_1.3.4         withr_2.1.2          bitops_1.0-6        
[36] SnowballC_0.5.1      gtable_0.2.0         magrittr_1.5         formatR_1.5          scales_1.0.0        
[41] KernSmooth_2.23-15   stringi_1.1.7        reshape2_1.4.3       limma_3.36.3         fastmatch_1.1-0     
[46] lambda.r_1.2.3       tools_3.5.0          Biobase_2.40.0       glue_1.3.0           purrr_0.2.5         
[51] parallel_3.5.0       yaml_2.2.0           colorspace_1.3-2     cluster_2.0.7-1      bindr_0.1.1   
