sessionInfo()
R version 4.2.1 (2022-06-23 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19042)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ROCR_1.0-11                 KernSmooth_2.23-20          fields_14.1                 viridis_0.6.2               viridisLite_0.4.1          
 [6] spam_2.9-1                  DOSE_3.22.1                 clusterProfiler_4.4.4       MeSHDbi_1.32.0              AnnotationHub_3.4.0        
[11] BiocFileCache_2.4.0         dbplyr_2.2.1                org.Hs.eg.db_3.15.0         GOSemSim_2.22.0             glmGamPoi_1.8.0            
[16] EnhancedVolcano_1.14.0      DoubletFinder_2.0.3         future_1.28.0               patchwork_1.1.2             clustree_0.5.0             
[21] ggraph_2.0.6                plotly_4.10.0               EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.20.2            AnnotationFilter_1.20.0    
[26] GenomicFeatures_1.48.4      AnnotationDbi_1.58.0        Signac_1.8.0                harmony_0.1.0               monocle3_1.2.9             
[31] SingleCellExperiment_1.18.1 SummarizedExperiment_1.26.1 GenomicRanges_1.48.0        GenomeInfoDb_1.32.4         IRanges_2.30.1             
[36] S4Vectors_0.34.0            MatrixGenerics_1.8.1        matrixStats_0.62.0          Biobase_2.56.0              BiocGenerics_0.42.0        
[41] sp_1.5-0                    SeuratObject_4.1.2          Seurat_4.2.0                dplyr_1.0.10                ggrepel_0.9.1              
[46] ggridges_0.5.4              Matrix_1.5-1                cowplot_1.1.1               ggplot2_3.3.6               Rcpp_1.0.9                 

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.3                rtracklayer_1.56.1            scattermore_0.8               tidyr_1.2.1                  
  [5] bit64_4.0.5                   irlba_2.3.5.1                 DelayedArray_0.22.0           data.table_1.14.2            
  [9] rpart_4.1.16                  KEGGREST_1.36.3               RCurl_1.98-1.9                generics_0.1.3               
 [13] callr_3.7.2                   terra_1.6-17                  usethis_2.1.6                 RSQLite_2.2.18               
 [17] shadowtext_0.1.2              RANN_2.6.1                    enrichplot_1.16.2             bit_4.0.4                    
 [21] spatstat.data_2.2-0           xml2_1.3.3                    httpuv_1.6.6                  assertthat_0.2.1             
 [25] hms_1.1.2                     promises_1.2.0.1              fansi_1.0.3                   restfulr_0.0.15              
 [29] progress_1.2.2                igraph_1.3.5                  DBI_1.1.3                     htmlwidgets_1.5.4            
 [33] spatstat.geom_2.4-0           purrr_0.3.4                   ellipsis_0.3.2                backports_1.4.1              
 [37] sparseMatrixStats_1.8.0       biomaRt_2.52.0                deldir_1.0-6                  vctrs_0.4.2                  
 [41] remotes_2.4.2                 abind_1.4-5                   cachem_1.0.6                  withr_2.5.0                  
 [45] ggforce_0.4.1                 progressr_0.11.0              checkmate_2.1.0               sctransform_0.3.5            
 [49] treeio_1.20.2                 GenomicAlignments_1.32.1      prettyunits_1.1.1             goftest_1.2-3                
 [53] cluster_2.1.3                 dotCall64_1.0-2               ape_5.6-2                     lazyeval_0.2.2               
 [57] crayon_1.5.2                  labeling_0.4.2                pkgconfig_2.0.3               tweenr_2.0.2                 
 [61] vipor_0.4.5                   nlme_3.1-157                  pkgload_1.3.0                 ProtGenerics_1.28.0          
 [65] devtools_2.4.5                rlang_1.0.6                   globals_0.16.1                lifecycle_1.0.3              
 [69] miniUI_0.1.1.1                downloader_0.4                filelock_1.0.2                ggrastr_1.0.1                
 [73] polyclip_1.10-0               lmtest_0.9-40                 aplot_0.1.8                   boot_1.3-28                  
 [77] zoo_1.8-11                    beeswarm_0.4.0                processx_3.7.0                png_0.1-7                    
 [81] rjson_0.2.21                  bitops_1.0-7                  Biostrings_2.64.1             DelayedMatrixStats_1.18.1    
 [85] blob_1.2.3                    stringr_1.4.1                 qvalue_2.28.0                 parallelly_1.32.1            
 [89] spatstat.random_2.2-0         gridGraphics_0.5-1            scales_1.2.1                  memoise_2.0.1                
 [93] magrittr_2.0.3                plyr_1.8.7                    ica_1.0-3                     zlibbioc_1.42.0              
 [97] scatterpie_0.1.8              compiler_4.2.1                BiocIO_1.6.0                  RColorBrewer_1.1-3           
[101] lme4_1.1-30                   fitdistrplus_1.1-8            Rsamtools_2.12.0              cli_3.4.1                    
[105] XVector_0.36.0                urlchecker_1.0.1              listenv_0.8.0                 pbapply_1.5-0                
[109] ps_1.7.1                      MASS_7.3-57                   mgcv_1.8-40                   tidyselect_1.2.0             
[113] stringi_1.7.8                 yaml_2.3.5                    grid_4.2.1                    fastmatch_1.1-3              
[117] tools_4.2.1                   future.apply_1.9.1            parallel_4.2.1                rstudioapi_0.14              
[121] gridExtra_2.3                 farver_2.1.1                  Rtsne_0.16                    digest_0.6.29                
[125] BiocManager_1.30.18           rgeos_0.5-9                   shiny_1.7.2                   BiocVersion_3.15.2           
[129] later_1.3.0                   RcppAnnoy_0.0.19              httr_1.4.4                    colorspace_2.0-3             
[133] XML_3.99-0.11                 fs_1.5.2                      tensor_1.5                    reticulate_1.26              
[137] splines_4.2.1                 yulab.utils_0.0.5             uwot_0.1.14                   RcppRoll_0.3.0               
[141] tidytree_0.4.1                spatstat.utils_2.3-1          graphlayouts_0.8.2            ggplotify_0.1.0              
[145] sessioninfo_1.2.2             xtable_1.8-4                  ggtree_3.4.4                  jsonlite_1.8.2               
[149] nloptr_2.0.3                  tidygraph_1.2.2               ggfun_0.0.7                   R6_2.5.1                     
[153] profvis_0.3.7                 pillar_1.8.1                  htmltools_0.5.3               mime_0.12                    
[157] glue_1.6.2                    fastmap_1.1.0                 minqa_1.2.4                   BiocParallel_1.30.3          
[161] interactiveDisplayBase_1.34.0 codetools_0.2-18              maps_3.4.0                    fgsea_1.22.0                 
[165] pkgbuild_1.3.1                utf8_1.2.2                    lattice_0.20-45               spatstat.sparse_2.1-1        
[169] tibble_3.1.8                  ggbeeswarm_0.6.0              curl_4.3.2                    leiden_0.4.3                 
[173] GO.db_3.15.0                  survival_3.3-1                munsell_0.5.0                 DO.db_2.9                    
[177] GenomeInfoDbData_1.2.8        reshape2_1.4.4                gtable_0.3.1                  spatstat.core_2.4-4 

# CODING COMPENDIUM ####
# The following set of code is a description of the analysis performed in the 
# paper entitled "enter name of paper here"
# Author Fahd Qadir FMJ Lab Tulane University, Schoool of Medicine
# Date code was written: 11/16/2022
# R version 4.2.1 (2019-12-12) 'Funny-Looking Kid'


# CODING COMPENDIUM ####
# The following set of code is a description of the analysis performed in the 
# paper entitled "enter name of paper here"
# Author Fahd Qadir FMJ Lab Tulane University, Schoool of Medicine
# Date code was written: 11/16/2022
# R version 4.2.1 (2019-12-12) 'Funny-Looking Kid'

# LOAD LIBRARIES ####
# Restart Rstudio or R

install.packages('ggplot2')
install.packages('cowplot')
install.packages('Matrix')
install.packages('ggridges')
install.packages('ggrepel')
install.packages('dplyr')
#install.packages('Seurat')
install.packages('plotly')
install.packages('clustree')
install.packages('patchwork')
install.packages('future')
install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils',
                       'HDF5Array', 'terra', 'ggrastr'))
BiocManager::install("EnhancedVolcano")
BiocManager::install("DoubletFinder")
BiocManager::install("glmGamPoi")
BiocManager::install("GOSemSim")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationHub")
BiocManager::install("GenomeInfoDb")
BiocManager::install("MeSHDbi")
BiocManager::install("clusterProfiler")
BiocManager::install("DOSE")

# install Seurat from Github (automatically updates sctransform)
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac")
BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))

devtools::install_github("satijalab/seurat", ref = "develop")
devtools::install_github("satijalab/sctransform", ref = "develop", force = TRUE)
devtools::install_github('cole-trapnell-lab/monocle3')
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
install.packages("harmony")

# Run the following code once you have Seurat installed
suppressWarnings(
  {
    library(Rcpp)
    library(ggplot2)
    library(cowplot)
    library(Matrix)
    library(ggridges)
    library(ggrepel)
    library(dplyr)
    library(Seurat)
    library(monocle3)
    library(harmony)
    library(Signac)
    library(EnsDb.Hsapiens.v86)
    library(GenomeInfoDb)
    library(plotly)
    library(clustree)
    library(patchwork)
    library(future)
    library(DoubletFinder)
    library(EnhancedVolcano)
    library(glmGamPoi)
    library(GOSemSim)
    library(org.Hs.eg.db)
    library(AnnotationHub)
    library(MeSHDbi)
    library(clusterProfiler)
    library(DOSE)
  }
)

# Set global environment parameter par-proc
#options(future.globals.maxSize = 8000 * 1024^2)
set.seed(123)

# OBJECT SETUP AND NORMALIZATION ####
# STEP 1: Load 10X data ####
{
  HP2022801.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\1_220624_Fahd_GEX1_F1_HP-20228-01\filtered_feature_bc_matrix)")
  SAMN15877725.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\2_220624_Fahd_GEX2_F2_SAMN15877725\filtered_feature_bc_matrix)")
  HP2024001.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\3_220624_Fahd_GEX3_F3_HP-20240-01\filtered_feature_bc_matrix)")
  HP2031401.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\4_220624_Fahd_GEX4_F4_HP-20314-01\filtered_feature_bc_matrix)")
  HP2105501.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\5_210406 GEX_F5_HP-21055-01\filtered_feature_bc_matrix)")
  HP2106201.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\6_210406 GEX_F6_HP-21062-01\filtered_feature_bc_matrix)")
  HP2107001.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\7_210406 GEX_F7a_HP-21070-01\filtered_feature_bc_matrix)")
  HP2107901.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\9_210714 GEX_F9a_HP-21079-01\filtered_feature_bc_matrix)")
  HP2108601.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\10_211108 GEX_F10a_HP-21086-01\filtered_feature_bc_matrix)")
  HP2108901.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\11_211108 GEX_F11a_HP-21089-01\filtered_feature_bc_matrix)")
  HP2110001.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\12_211108 GEX_F12a_HP-21100-01\filtered_feature_bc_matrix)")
  HP2121601.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\13_211108 GEX_F13a_HP-21216-01\filtered_feature_bc_matrix)")
  HP2123201.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\14_220124 GEX_F14_HP21232-01\filtered_feature_bc_matrix)")
  HP2132801.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\15_220124 GEX_F15_HP-21328-01\filtered_feature_bc_matrix)")
  HP2202101.data <- Read10X(data.dir = r"(D:\1. Sex based Study raw data\Cellranger_raw data\scRNAseq\16_220324 GEX_F16_HP-22021-01\filtered_feature_bc_matrix)")


# STEP 2: Create Seurat objects ####

  HP2022801 <- CreateSeuratObject(counts = HP2022801.data, min.features = 500)
  SAMN15877725 <- CreateSeuratObject(counts = SAMN15877725.data, min.features = 500)
  HP2024001 <- CreateSeuratObject(counts = HP2024001.data, min.features = 500)
  HP2031401 <- CreateSeuratObject(counts = HP2031401.data, min.features = 500)
  HP2105501 <- CreateSeuratObject(counts = HP2105501.data, min.features = 500)
  HP2106201 <- CreateSeuratObject(counts = HP2106201.data, min.features = 500)
  HP2107001 <- CreateSeuratObject(counts = HP2107001.data, min.features = 500)
  HP2107901 <- CreateSeuratObject(counts = HP2107901.data, min.features = 500)
  HP2108601 <- CreateSeuratObject(counts = HP2108601.data, min.features = 500)
  HP2108901 <- CreateSeuratObject(counts = HP2108901.data, min.features = 500)
  HP2110001 <- CreateSeuratObject(counts = HP2110001.data, min.features = 500)
  HP2121601 <- CreateSeuratObject(counts = HP2121601.data, min.features = 500)
  HP2123201 <- CreateSeuratObject(counts = HP2123201.data, min.features = 500)
  HP2132801 <- CreateSeuratObject(counts = HP2132801.data, min.features = 500)
  HP2202101 <- CreateSeuratObject(counts = HP2202101.data, min.features = 500)
  }

# Sample specific Metadata addition
{
  HP2022801$sample <- "HP2022801"
  SAMN15877725$sample <- "SAMN15877725"
  HP2024001$sample <- "HP2024001"
  HP2031401$sample <- "HP2031401"
  HP2105501$sample <- "HP2105501"
  HP2106201$sample <- "HP2106201"
  HP2107001$sample <- "HP2107001"
  HP2107901$sample <- "HP2107901"
  HP2108601$sample <- "HP2108601"
  HP2108901$sample <- "HP2108901"
  HP2110001$sample <- "HP2110001"
  HP2121601$sample <- "HP2121601"
  HP2123201$sample <- "HP2123201"
  HP2132801$sample <- "HP2132801"
  HP2202101$sample <- "HP2202101"


# Sex specific Metadata addition
  HP2022801$sex <- "male"
  SAMN15877725$sex <- "male"
  HP2024001$sex <- "female"
  HP2031401$sex <- "male"
  HP2105501$sex <- "female"
  HP2106201$sex <- "female"
  HP2107001$sex <- "male"
  HP2107901$sex <- "male"
  HP2108601$sex <- "female"
  HP2108901$sex <- "female"
  HP2110001$sex <- "male"
  HP2121601$sex <- "female"
  HP2123201$sex <- "male"
  HP2132801$sex <- "female"
  HP2202101$sex <- "female"

# Ancestry specific Metadata addition
  HP2022801$ancestry <- "white"
  SAMN15877725$ancestry <- "white"
  HP2024001$ancestry <- "white"
  HP2031401$ancestry <- "black"
  HP2105501$ancestry <- "white"
  HP2106201$ancestry <- "black"
  HP2107001$ancestry <- "white"
  HP2107901$ancestry <- "white"
  HP2108601$ancestry <- "white"
  HP2108901$ancestry <- "white"
  HP2110001$ancestry <- "black"
  HP2121601$ancestry <- "black"
  HP2123201$ancestry <- "black"
  HP2132801$ancestry <- "black"
  HP2202101$ancestry <- "black"

# Ancestry and sex specific UNION Metadata addition
  HP2022801$ancestry_sex <- "white_male"
  SAMN15877725$ancestry_sex <- "white_male"
  HP2024001$ancestry_sex <- "white_female"
  HP2031401$ancestry_sex <- "black_male"
  HP2105501$ancestry_sex <- "white_female"
  HP2106201$ancestry_sex <- "black_female"
  HP2107001$ancestry_sex <- "white_male"
  HP2107901$ancestry_sex <- "white_male"
  HP2108601$ancestry_sex <- "white_female"
  HP2108901$ancestry_sex <- "white_female"
  HP2110001$ancestry_sex <- "black_male"
  HP2121601$ancestry_sex <- "black_female"
  HP2123201$ancestry_sex <- "black_male"
  HP2132801$ancestry_sex <- "black_female"
  HP2202101$ancestry_sex <- "black_female"
  }

# STEP 3: Thresholding ####
# Add -MT gene percentage data as QC to cutoff for doublets 1st doublet threshold
{
  HP2022801[["percent.mt"]] <- PercentageFeatureSet(object = HP2022801, pattern = "^MT-")
  SAMN15877725[["percent.mt"]] <- PercentageFeatureSet(object = SAMN15877725, pattern = "^MT-")
  HP2024001[["percent.mt"]] <- PercentageFeatureSet(object = HP2024001, pattern = "^MT-")
  HP2031401[["percent.mt"]] <- PercentageFeatureSet(object = HP2031401, pattern = "^MT-")
  HP2105501[["percent.mt"]] <- PercentageFeatureSet(object = HP2105501, pattern = "^MT-")
  HP2106201[["percent.mt"]] <- PercentageFeatureSet(object = HP2106201, pattern = "^MT-")
  HP2107001[["percent.mt"]] <- PercentageFeatureSet(object = HP2107001, pattern = "^MT-")
  HP2107901[["percent.mt"]] <- PercentageFeatureSet(object = HP2107901, pattern = "^MT-")
  HP2108601[["percent.mt"]] <- PercentageFeatureSet(object = HP2108601, pattern = "^MT-")
  HP2108901[["percent.mt"]] <- PercentageFeatureSet(object = HP2108901, pattern = "^MT-")
  HP2110001[["percent.mt"]] <- PercentageFeatureSet(object = HP2110001, pattern = "^MT-")
  HP2121601[["percent.mt"]] <- PercentageFeatureSet(object = HP2121601, pattern = "^MT-")
  HP2123201[["percent.mt"]] <- PercentageFeatureSet(object = HP2123201, pattern = "^MT-")
  HP2132801[["percent.mt"]] <- PercentageFeatureSet(object = HP2132801, pattern = "^MT-")
  HP2202101[["percent.mt"]] <- PercentageFeatureSet(object = HP2202101, pattern = "^MT-")

# QC information before thresholding
  summary(head(HP2022801@meta.data))
  summary(head(SAMN15877725@meta.data))
  summary(head(HP2024001@meta.data))
  summary(head(HP2031401@meta.data))
  summary(head(HP2105501@meta.data))
  summary(head(HP2106201@meta.data))
  summary(head(HP2107001@meta.data))
  summary(head(HP2107901@meta.data))
  summary(head(HP2108601@meta.data))
  summary(head(HP2108901@meta.data))
  summary(head(HP2110001@meta.data))
  summary(head(HP2121601@meta.data))
  summary(head(HP2123201@meta.data))
  summary(head(HP2132801@meta.data))
  summary(head(HP2202101@meta.data))

# Visualize QC metrics as a violin plot
  VlnPlot(object = HP2022801, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = SAMN15877725, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2024001, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2031401, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2105501, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2106201, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2107001, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2107901, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2108601, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2108901, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2110001, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2121601, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2123201, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2132801, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2202101, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# RNA features + MT RNA percentage based cell thresholding 1st THRESHOLD
  HP2022801 <- subset(x = HP2022801, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  SAMN15877725 <- subset(x = SAMN15877725, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2024001 <- subset(x = HP2024001, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2031401 <- subset(x = HP2031401, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2105501 <- subset(x = HP2105501, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2106201 <- subset(x = HP2106201, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2107001 <- subset(x = HP2107001, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2107901 <- subset(x = HP2107901, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2108601 <- subset(x = HP2108601, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2108901 <- subset(x = HP2108901, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2110001 <- subset(x = HP2110001, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2121601 <- subset(x = HP2121601, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2123201 <- subset(x = HP2123201, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2132801 <- subset(x = HP2132801, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  HP2202101 <- subset(x = HP2202101, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)

# QC information after thresholding
  summary(head(HP2022801@meta.data))
  summary(head(SAMN15877725@meta.data))
  summary(head(HP2024001@meta.data))
  summary(head(HP2031401@meta.data))
  summary(head(HP2105501@meta.data))
  summary(head(HP2106201@meta.data))
  summary(head(HP2107001@meta.data))
  summary(head(HP2107901@meta.data))
  summary(head(HP2108601@meta.data))
  summary(head(HP2108901@meta.data))
  summary(head(HP2110001@meta.data))
  summary(head(HP2121601@meta.data))
  summary(head(HP2123201@meta.data))
  summary(head(HP2132801@meta.data))
  summary(head(HP2202101@meta.data))

# Visualize QC metrics post thresholding as a violin plot
  VlnPlot(object = HP2022801, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = SAMN15877725, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2024001, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2031401, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2105501, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2106201, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2107001, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2107901, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2108601, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2108901, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2110001, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2121601, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2123201, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2132801, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(object = HP2202101, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Step 4: Add cell IDs ####
# Add UNIQUE cell IDs preventing barcode overlap
  HP2022801 <- RenameCells(HP2022801, add.cell.id = "HP2022801")
  SAMN15877725 <- RenameCells(SAMN15877725, add.cell.id = "SAMN15877725")
  HP2024001 <- RenameCells(HP2024001, add.cell.id = "HP2024001")
  HP2031401 <- RenameCells(HP2031401, add.cell.id = "HP2031401")
  HP2105501 <- RenameCells(HP2105501, add.cell.id = "HP2105501")
  HP2106201 <- RenameCells(HP2106201, add.cell.id = "HP2106201")
  HP2107001 <- RenameCells(HP2107001, add.cell.id = "HP2107001")
  HP2107901 <- RenameCells(HP2107901, add.cell.id = "HP2107901")
  HP2108601 <- RenameCells(HP2108601, add.cell.id = "HP2108601")
  HP2108901 <- RenameCells(HP2108901, add.cell.id = "HP2108901")
  HP2110001 <- RenameCells(HP2110001, add.cell.id = "HP2110001")
  HP2121601 <- RenameCells(HP2121601, add.cell.id = "HP2121601")
  HP2123201 <- RenameCells(HP2123201, add.cell.id = "HP2123201")
  HP2132801 <- RenameCells(HP2132801, add.cell.id = "HP2132801")
  HP2202101 <- RenameCells(HP2202101, add.cell.id = "HP2202101")
  }

# Normalization for visualization and doublet removal
# Dooublet removal using doubletfinder 2nd THRESHOLD
{
 HP2022801 <- NormalizeData(HP2022801, verbose = TRUE)
 SAMN15877725 <- NormalizeData(SAMN15877725, verbose = TRUE)
 HP2024001 <- NormalizeData(HP2024001, verbose = TRUE)
 HP2031401 <- NormalizeData(HP2031401, verbose = TRUE)
 HP2105501 <- NormalizeData(HP2105501, verbose = TRUE)
 HP2106201 <- NormalizeData(HP2106201, verbose = TRUE)
 HP2107001 <- NormalizeData(HP2107001, verbose = TRUE)
 HP2107901 <- NormalizeData(HP2107901, verbose = TRUE)
 HP2108601 <- NormalizeData(HP2108601, verbose = TRUE)
 HP2108901 <- NormalizeData(HP2108901, verbose = TRUE)
 HP2110001 <- NormalizeData(HP2110001, verbose = TRUE)
 HP2121601 <- NormalizeData(HP2121601, verbose = TRUE)
 HP2123201 <- NormalizeData(HP2123201, verbose = TRUE)
 HP2132801 <- NormalizeData(HP2132801, verbose = TRUE)
 HP2202101 <- NormalizeData(HP2202101, verbose = TRUE)
 
# Var Feat
 HP2022801 <- FindVariableFeatures(HP2022801, selection.method = "vst", 
                                   nfeatures = 2000, verbose = TRUE)
 SAMN15877725 <- FindVariableFeatures(SAMN15877725, selection.method = "vst", 
                                      nfeatures = 2000, verbose = TRUE)
 HP2024001 <- FindVariableFeatures(HP2024001, selection.method = "vst", 
                                   nfeatures = 2000, verbose = TRUE)
 HP2031401 <- FindVariableFeatures(HP2031401, selection.method = "vst", 
                                   nfeatures = 2000, verbose = TRUE)
 HP2105501 <- FindVariableFeatures(HP2105501, selection.method = "vst", 
                                   nfeatures = 2000, verbose = TRUE)
 HP2106201 <- FindVariableFeatures(HP2106201, selection.method = "vst", 
                                   nfeatures = 2000, verbose = TRUE)
 HP2107001 <- FindVariableFeatures(HP2107001, selection.method = "vst", 
                                   nfeatures = 2000, verbose = TRUE)
 HP2107901 <- FindVariableFeatures(HP2107901, selection.method = "vst", 
                                   nfeatures = 2000, verbose = TRUE)
 HP2108601 <- FindVariableFeatures(HP2108601, selection.method = "vst", 
                                   nfeatures = 2000, verbose = TRUE)
 HP2108901 <- FindVariableFeatures(HP2108901, selection.method = "vst", 
                                   nfeatures = 2000, verbose = TRUE)
 HP2110001 <- FindVariableFeatures(HP2110001, selection.method = "vst", 
                                   nfeatures = 2000, verbose = TRUE)
 HP2121601 <- FindVariableFeatures(HP2121601, selection.method = "vst", 
                                   nfeatures = 2000, verbose = TRUE)
 HP2123201 <- FindVariableFeatures(HP2123201, selection.method = "vst", 
                                   nfeatures = 2000, verbose = TRUE)
 HP2132801 <- FindVariableFeatures(HP2132801, selection.method = "vst", 
                                   nfeatures = 2000, verbose = TRUE)
 HP2202101 <- FindVariableFeatures(HP2202101, selection.method = "vst", 
                                   nfeatures = 2000, verbose = TRUE)

# Scale data
 HP2022801 <- ScaleData(HP2022801, verbose = TRUE)
 SAMN15877725 <- ScaleData(SAMN15877725, verbose = TRUE)
 HP2024001 <- ScaleData(HP2024001, verbose = TRUE)
 HP2031401 <- ScaleData(HP2031401, verbose = TRUE)
 HP2105501 <- ScaleData(HP2105501, verbose = TRUE)
 HP2106201 <- ScaleData(HP2106201, verbose = TRUE)
 HP2107001 <- ScaleData(HP2107001, verbose = TRUE)
 HP2107901 <- ScaleData(HP2107901, verbose = TRUE)
 HP2108601 <- ScaleData(HP2108601, verbose = TRUE)
 HP2108901 <- ScaleData(HP2108901, verbose = TRUE)
 HP2110001 <- ScaleData(HP2110001, verbose = TRUE)
 HP2121601 <- ScaleData(HP2121601, verbose = TRUE)
 HP2123201 <- ScaleData(HP2123201, verbose = TRUE)
 HP2132801 <- ScaleData(HP2132801, verbose = TRUE)
 HP2202101 <- ScaleData(HP2202101, verbose = TRUE) 
 
# PCA
 HP2022801 <- RunPCA(HP2022801, verbose = TRUE, npcs = 20)
 SAMN15877725 <- RunPCA(SAMN15877725, verbose = TRUE, npcs = 20)
 HP2024001 <- RunPCA(HP2024001, verbose = TRUE, npcs = 20)
 HP2031401 <- RunPCA(HP2031401, verbose = TRUE, npcs = 20)
 HP2105501 <- RunPCA(HP2105501, verbose = TRUE, npcs = 20)
 HP2106201 <- RunPCA(HP2106201, verbose = TRUE, npcs = 20)
 HP2107001 <- RunPCA(HP2107001, verbose = TRUE, npcs = 20)
 HP2107901 <- RunPCA(HP2107901, verbose = TRUE, npcs = 20)
 HP2108601 <- RunPCA(HP2108601, verbose = TRUE, npcs = 20)
 HP2108901 <- RunPCA(HP2108901, verbose = TRUE, npcs = 20)
 HP2110001 <- RunPCA(HP2110001, verbose = TRUE, npcs = 20)
 HP2121601 <- RunPCA(HP2121601, verbose = TRUE, npcs = 20)
 HP2123201 <- RunPCA(HP2123201, verbose = TRUE, npcs = 20)
 HP2132801 <- RunPCA(HP2132801, verbose = TRUE, npcs = 20)
 HP2202101 <- RunPCA(HP2202101, verbose = TRUE, npcs = 20) 

 # UMAP
 HP2022801 <- RunUMAP(HP2022801, dims = 1:10, verbose = F)
 SAMN15877725 <- RunUMAP(SAMN15877725, dims = 1:10, verbose = F)
 HP2024001 <- RunUMAP(HP2024001, dims = 1:10, verbose = F)
 HP2031401 <- RunUMAP(HP2031401, dims = 1:10, verbose = F)
 HP2105501 <- RunUMAP(HP2105501, dims = 1:10, verbose = F)
 HP2106201 <- RunUMAP(HP2106201, dims = 1:10, verbose = F)
 HP2107001 <- RunUMAP(HP2107001, dims = 1:10, verbose = F)
 HP2107901 <- RunUMAP(HP2107901, dims = 1:10, verbose = F)
 HP2108601 <- RunUMAP(HP2108601, dims = 1:10, verbose = F)
 HP2108901 <- RunUMAP(HP2108901, dims = 1:10, verbose = F)
 HP2110001 <- RunUMAP(HP2110001, dims = 1:10, verbose = F)
 HP2121601 <- RunUMAP(HP2121601, dims = 1:10, verbose = F)
 HP2123201 <- RunUMAP(HP2123201, dims = 1:10, verbose = F)
 HP2132801 <- RunUMAP(HP2132801, dims = 1:10, verbose = F)
 HP2202101 <- RunUMAP(HP2202101, dims = 1:10, verbose = F)
}
# for (i in 1:length(pancreas.list)) {
#   pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = TRUE)
#   pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", 
#                                              nfeatures = 2000, verbose = TRUE)
# }

# Optimization
{
sweep.res <- paramSweep_v3(HP2022801) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)

barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

sweep.res <- paramSweep_v3(SAMN15877725) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

sweep.res <- paramSweep_v3(HP2024001) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

sweep.res <- paramSweep_v3(HP2031401) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

sweep.res <- paramSweep_v3(HP2105501) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

sweep.res <- paramSweep_v3(HP2106201) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

sweep.res <- paramSweep_v3(HP2107001) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

sweep.res <- paramSweep_v3(HP2107901) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

sweep.res <- paramSweep_v3(HP2108601) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

sweep.res <- paramSweep_v3(HP2108901) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

sweep.res <- paramSweep_v3(HP2110001) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

sweep.res <- paramSweep_v3(HP2121601) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

sweep.res <- paramSweep_v3(HP2123201) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

sweep.res <- paramSweep_v3(HP2132801) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

sweep.res <- paramSweep_v3(HP2202101) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

# Filtering
nExp <- round(ncol(HP2022801) * 0.04)  # expect 4% doublets
HP2022801 <- doubletFinder_v3(HP2022801, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

nExp <- round(ncol(SAMN15877725) * 0.04)  # expect 4% doublets
SAMN15877725 <- doubletFinder_v3(SAMN15877725, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

nExp <- round(ncol(HP2024001) * 0.04)  # expect 4% doublets
HP2024001 <- doubletFinder_v3(HP2024001, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

nExp <- round(ncol(HP2031401) * 0.04)  # expect 4% doublets
HP2031401 <- doubletFinder_v3(HP2031401, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

nExp <- round(ncol(HP2105501) * 0.04)  # expect 4% doublets
HP2105501 <- doubletFinder_v3(HP2105501, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

nExp <- round(ncol(HP2106201) * 0.04)  # expect 4% doublets
HP2106201 <- doubletFinder_v3(HP2106201, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

nExp <- round(ncol(HP2107001) * 0.04)  # expect 4% doublets
HP2107001 <- doubletFinder_v3(HP2107001, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

nExp <- round(ncol(HP2107901) * 0.04)  # expect 4% doublets
HP2107901 <- doubletFinder_v3(HP2107901, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

nExp <- round(ncol(HP2108601) * 0.04)  # expect 4% doublets
HP2108601 <- doubletFinder_v3(HP2108601, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

nExp <- round(ncol(HP2108901) * 0.04)  # expect 4% doublets
HP2108901 <- doubletFinder_v3(HP2108901, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

nExp <- round(ncol(HP2110001) * 0.04)  # expect 4% doublets
HP2110001 <- doubletFinder_v3(HP2110001, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

nExp <- round(ncol(HP2121601) * 0.04)  # expect 4% doublets
HP2121601 <- doubletFinder_v3(HP2121601, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

nExp <- round(ncol(HP2123201) * 0.04)  # expect 4% doublets
HP2123201 <- doubletFinder_v3(HP2123201, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

nExp <- round(ncol(HP2132801) * 0.04)  # expect 4% doublets
HP2132801 <- doubletFinder_v3(HP2132801, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

nExp <- round(ncol(HP2202101) * 0.04)  # expect 4% doublets
HP2202101 <- doubletFinder_v3(HP2202101, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
}

# Setup one metadata column
HP2022801$doublets <- HP2022801$DF.classifications_0.25_0.09_166
SAMN15877725$doublets <- SAMN15877725$DF.classifications_0.25_0.09_155
HP2107001$doublets <- HP2107001$DF.classifications_0.25_0.09_170
HP2107901$doublets <- HP2107901$DF.classifications_0.25_0.09_132
HP2024001$doublets <- HP2024001$DF.classifications_0.25_0.09_121
HP2105501$doublets <- HP2105501$DF.classifications_0.25_0.09_124
HP2108601$doublets <- HP2108601$DF.classifications_0.25_0.09_218
HP2108901$doublets <- HP2108901$DF.classifications_0.25_0.09_171
HP2031401$doublets <- HP2031401$DF.classifications_0.25_0.09_183
HP2110001$doublets <- HP2110001$DF.classifications_0.25_0.09_226
HP2123201$doublets <- HP2123201$DF.classifications_0.25_0.09_62
HP2106201$doublets <- HP2106201$DF.classifications_0.25_0.09_260
HP2121601$doublets <- HP2121601$DF.classifications_0.25_0.09_140
HP2132801$doublets <- HP2132801$DF.classifications_0.25_0.09_93
HP2202101$doublets <- HP2202101$DF.classifications_0.25_0.09_159

# Step 5: creating a list of all datasets
{
  pancreas.list <- list("HP2022801" = HP2022801, "SAMN15877725" = SAMN15877725, "HP2107001" = HP2107001, "HP2107901" = HP2107901,
                        "HP2024001" = HP2024001, "HP2105501" = HP2105501, "HP2108601" = HP2108601, "HP2108901" = HP2108901, 
                        "HP2031401" = HP2031401, "HP2110001" = HP2110001, "HP2123201" = HP2123201,
                        "HP2106201" = HP2106201, "HP2121601" = HP2121601, "HP2132801" = HP2132801, "HP2202101" = HP2202101
                         )
}

# normalize and identify variable features for each dataset independently
pancreas.list <- lapply(X = pancreas.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 2000)
pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = features)
pancreas.list <- lapply(X = pancreas.list, FUN = RunPCA, features = features)

# Perform integration (note k.anchors = 5)
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, 
                                           normalization.method = "SCT",
                                           anchor.features = features, 
                                           dims = 1:30, 
                                           reduction = "rpca", 
                                           k.anchor = 5)

#saveRDS(pancreas.anchors, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\pancreas.anchors.rds)")
#saveRDS(features, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\features.rds)")
#saveRDS(pancreas.list, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\pancreas.list.rds)")
pancreas.anchors <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\pancreas.anchors10.rds)")
features <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\features.rds)")
pancreas.combined <- readRDS(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\pancreas.combinedcorrectedSCT.rds)")
gc()
pancreas.combined <- IntegrateData(anchorset = pancreas.anchors, 
                                   normalization.method = "SCT", 
                                   dims = 1:30,
                                   verbose = TRUE,
                                   features.to.integrate = features)

# Post integration UMAP
pancreas.combined <- RunPCA(pancreas.combined, npcs = 50, verbose = TRUE)
pancreas.combined <- RunUMAP(pancreas.combined, reduction = "pca", dims = 1:50, return.model = TRUE)

# Step 9a: CLUSTER ANALYSIS ####
# Clustering, we run multiple permutations to allow clustree to analyze optimal clustering resolution.
DefaultAssay(pancreas.combined) <- "integrated"
pancreas.combined <- FindNeighbors(pancreas.combined, reduction = "pca", dims = 1:50)
pancreas.combined <- FindClusters(object = pancreas.combined, resolution = 0)
pancreas.combined <- FindClusters(object = pancreas.combined, resolution = 0.1)
pancreas.combined <- FindClusters(object = pancreas.combined, resolution = 0.2)
pancreas.combined <- FindClusters(object = pancreas.combined, resolution = 0.3)
pancreas.combined <- FindClusters(object = pancreas.combined, resolution = 0.4)
pancreas.combined <- FindClusters(object = pancreas.combined, resolution = 0.5)
pancreas.combined <- FindClusters(object = pancreas.combined, resolution = 0.6)

# Cluster-tree analysis, looking appropriate non-anomalous clustering resolution
clustree(pancreas.combined, prefix = "integrated_snn_res.")

# Based of clustree assessment choose res = 0.3
pancreas.combined <- FindClusters(pancreas.combined, resolution = 0.3)

# Alternatively build a cluster tree
DefaultAssay(object = pancreas.combined) <- "integrated"
pancreas.combined = BuildClusterTree(pancreas.combined, slot = "scale.data")
PlotClusterTree(pancreas.combined)

# Visualization
DimPlot(pancreas.combined, reduction = "umap", group.by = "ancestry_sex")
DimPlot(pancreas.combined, reduction = "umap", group.by = "doublets")

# Cluster-tree analysis, looking appropriate non-anomalous clustering resolution
clustree(pancreas.combined, prefix = "integrated_snn_res.")

# Plotting
DimPlot(pancreas.combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
DimPlot(pancreas.combined, reduction = "umap", group.by = "integrated_snn_res.0.4", label = TRUE, repel = TRUE) # this along with doublet UMAP shows cluster 7 in beta to be doublet

# Elimination of doublets
pancreas.combined.withdoublets <- pancreas.combined
# pancreas.combined <- pancreas.combined.withdoublets # OBJECT RESET RUN CAREFULLY
Idents(pancreas.combined) <- "doublets"
pancreas.combined <- subset(pancreas.combined, idents = "Singlet")

# SAVE POINT
# saveRDS(pancreas.combined, file = r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\RDS files\pancreas.combined.rds)")
# pancreas.combined <- readRDS('~/Documents/SexBasedStudy/RDSfiles/pancreas.combined.rds')

# Discovery based Plotting
DefaultAssay(pancreas.combined) <- "SCT"
FeaturePlot(object = pancreas.combined,
            features = c("GCG", "INS", "SST", "PPY", "GHRL"
                         ),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            #max.cutoff = 100,
            slot = 'counts',
            order = TRUE)

FeaturePlot(object = pancreas.combined,
            features = c("VWF", "SDS", "CD8A", "TRAC", "TPSAB1"
            ),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            #max.cutoff = 100,
            slot = 'counts',
            order = TRUE)

FeaturePlot(object = pancreas.combined,
            features = c("PDGFRA", "RGS5", "REG1A", "SPP1", "CFTR", "SOX10"
            ),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            #max.cutoff = 100,
            slot = 'counts',
            order = TRUE)

FeaturePlot(object = pancreas.combined,
            features = c("MKI67"
            ),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            #max.cutoff = 100,
            slot = 'counts',
            order = TRUE)

VlnPlot(
  object = pancreas.combined,
  features = c("TMSB4X"),
  assay = 'RNA',
  slot = 'counts',
  # cols = c("red4", "red3", "grey40", "orange", "lightgoldenrod3", "yellow4", "indianred", "orangered", "black",
  #          "royalblue2", "steelblue1", "darkcyan",
  #          "springgreen4", "green3", "darkturquoise",
  #          "purple4", "purple", "deeppink",
  #          "violetred", "violet"),
  #y.max = 3,
  pt.size = 1
)


# Observing cells
DimPlot(pancreas.combined, 
        split.by = "ancestry_sex", group.by = "celltype", 
        label = FALSE, ncol = 2,  
        cols = c("red4", "red3", "grey40", "orange", "lightgoldenrod3", "yellow4", "indianred", "orangered", "black",
                 "royalblue2", "steelblue1", "darkcyan",
                 "springgreen4", "green3", "darkturquoise",
                 "purple4", "purple", "deeppink",
                 "violetred", "violet"
))

DimPlot(pancreas.combined, 
        group.by = "celltype", 
        label = FALSE, ncol = 1,  
        cols = c("red4", "red3", "grey40", "orange", "lightgoldenrod3", "yellow4", "indianred", "orangered", "black",
                 "royalblue2", "steelblue1", "darkcyan",
                 "springgreen4", "green3", "darkturquoise",
                 "purple4", "purple", "deeppink",
                 "violetred", "violet"
        ))

markers <- FindAllMarkers(pancreas.combined, assay = "RNA",
                          logfc.threshold = 1,
                          test.use = "wilcox",
                          slot = "data",
                          min.pct = 0.5)
