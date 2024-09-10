seed_saved <- readLines("Random.seed.txt")
.Random.seed <- as.integer(seed_saved)

library(SingleCellExperiment)
library(balagan)
library(pheatmap)
library(dplyr)
library(ggpubr)
library(viridis)
library(Polychrome)

input_filename <- "./data/dataScaleSize_slide2_FOV1.csv"
fov1 <- read.csv(input_filename)

# Contruct the input file for the SingleCellExperiment. Columns were renamed to fit the function from the balagan package
Expression_data <- fov1[, 5:ncol(fov1)]
Cell_annotation <- fov1[, 1:4]
Cell_annotation$Cell_size <- Cell_annotation$cellSize
Cell_annotation$Location_Center_X <- Cell_annotation$X_cent
Cell_annotation$Location_Center_Y <- Cell_annotation$Y_cent

Gene_annotation <- colnames(Expression_data)

sce = SingleCellExperiment(
  assays = list(Raw_intensity = as.matrix(t(Expression_data))),
  colData = Cell_annotation,
  rowData = Gene_annotation,
  metadata = list(
    dimension = "2D",
    Bit_mode = 16,
    N_core = 6,
    Is_nuc_cyt = F
  )
)

sce = Count_normalization(sce, residual_normalisation = "Pearson")
sce = KNN_clustering(
  sce,
  K = 15,
  clustering_method = "Louvain",
  assay_type = "Count_normalised_intensity",
  metric = "L2"
)

message(paste0("Total number of clusters: ", length(unique(sce$label))))

# Create the heatmap
expr_data <- sce@assays@data$Count_normalised_intensity
cluster_labels <- sce$label
df <- cbind(t(expr_data), cluster_labels)
df <- as.data.frame(df)
df[, 1:(ncol(df) - 1)] <- apply(df[, 1:(ncol(df) - 1)], 2, as.numeric)
summary_data <- aggregate(. ~ cluster_labels, df, mean)
summary_data$cluster_labels <- NULL
summary_data_t <- t(summary_data)
selected_features <-
  c(
    "CD3",
    "CD4",
    "CD8",
    "CD11b",
    "CD11c",
    "CD15",
    "CD20",
    "CD21",
    "CD31",
    "CD68",
    "CD138",
    "Cytokeratin",
    "Foxp3"
  )
selected_features <-
  selected_features[selected_features %in% rownames(summary_data_t)]
colnames(summary_data_t) <- c(1:length(unique(sce$label)))

svg(paste0(input_filename, "_heatmap.svg"),
    width = 10,
    height = 8)

pheatmap(
  summary_data_t[selected_features, ],
  scale = "row",
  cluster_cols = F,
  angle_col = 45
)
dev.off()


# Create the scatter plot colored by clusters
fov1_scatter <- fov1[, c(1, 2, 3)]
fov1_scatter$cluster <- as.factor(sce$label)
color_vector <-
  setNames(as.character(palette36.colors(36)[-2]), sort(unique(fov1_scatter$cluster)))

pp1 <-
  ggscatter(
    fov1_scatter,
    x = "X_cent",
    y = "Y_cent",
    color = "cluster",
    xlab = "X_cent",
    ylab = "Y_cent",
    size = 0.5
  ) +
  scale_color_manual(values = color_vector)

ggsave(
  filename = paste0(input_filename, "_scatter.svg"),
  plot = pp1,
  width = 20,
  height = 18
)


# Run Simple sampleing analysis, and plot the number of recovered clusters versus number of sampled regions
sce$ImageNumber <- 1
Simple_sampling_analysis = Perform_sampling_analysis(
  sce,
  Selected_image = 1,
  N_times = 20,
  N_sampling_region_vector = 1:20,
  width_FOV_vector = 400,
  height_FOV_vector = 400,
  Threshold_detection_cluster = 2
)

Visualize_simple_sampling(Simple_sampling_analysis)

svg(paste0(input_filename, "_subsampling.svg"),
    width = 10,
    height = 8)
Visualize_simple_sampling(Simple_sampling_analysis)
dev.off()

# > sessionInfo()
# R version 4.4.1 (2024-06-14)
# Platform: x86_64-redhat-linux-gnu
# Running under: Red Hat Enterprise Linux 8.10 (Ootpa)
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib64/libopenblaso-r0.3.15.so;  LAPACK version 3.9.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: America/New_York
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] Polychrome_1.5.1            viridis_0.6.4               viridisLite_0.4.2           ggpubr_0.6.0                ggplot2_3.4.3               dplyr_1.1.3                 pheatmap_1.0.12             balagan_1.0.0              
# [9] spatstat_3.1-1              spatstat.linnet_3.2-1       spatstat.model_3.3-1        rpart_4.1.23                spatstat.explore_3.3-2      nlme_3.1-164                spatstat.random_3.3-1       spatstat.geom_3.3-2        
# [17] spatstat.univar_3.0-0       spatstat.data_3.1-2         SingleCellExperiment_1.22.0 SummarizedExperiment_1.30.2 Biobase_2.60.0              GenomicRanges_1.52.0        GenomeInfoDb_1.36.3         IRanges_2.34.1             
# [25] S4Vectors_0.38.1            BiocGenerics_0.46.0         MatrixGenerics_1.12.3       matrixStats_1.0.0          
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-7            deldir_1.0-9            gridExtra_2.3           rlang_1.1.1             magrittr_2.0.3          compiler_4.4.1          mgcv_1.9-1              systemfonts_1.0.4       png_0.1-8              
# [10] vctrs_0.6.3             stringr_1.5.0           pkgconfig_2.0.3         crayon_1.5.2            backports_1.4.1         XVector_0.40.0          labeling_0.4.3          magic_1.6-1             utf8_1.2.3             
# [19] tzdb_0.4.0              ragg_1.2.5              purrr_1.0.2             zlibbioc_1.46.0         goftest_1.2-3           DelayedArray_0.26.7     spatstat.utils_3.1-0    jpeg_0.1-10             tiff_0.1-12            
# [28] broom_1.0.5             parallel_4.4.1          R6_2.5.1                stringi_1.7.12          RColorBrewer_1.1-3      car_3.1-2               Rcpp_1.0.13             iterators_1.0.14        tensor_1.5             
# [37] readr_2.1.4             Matrix_1.7-0            splines_4.4.1           igraph_1.5.1            tidyselect_1.2.0        rstudioapi_0.15.0       abind_1.4-5             doParallel_1.0.17       codetools_0.2-20       
# [46] lattice_0.22-6          tibble_3.2.1            withr_2.5.0             survival_3.6-4          polyclip_1.10-4         fitdistrplus_1.1-11     pillar_1.9.0            carData_3.0-5           foreach_1.5.2          
# [55] geometry_0.5.0          generics_0.1.3          dbscan_1.2-0            RCurl_1.98-1.12         hms_1.1.3               munsell_0.5.0           scales_1.2.1            aod_1.3.3               rTensor_1.4.8          
# [64] glue_1.6.2              scatterplot3d_0.3-44    tools_4.4.1             ggsignif_0.6.4          bmp_0.3                 grid_4.4.1              tidyr_1.3.0             spatgraphs_3.4          readbitmap_0.1.5       
# [73] colorspace_2.1-0        GenomeInfoDbData_1.2.10 N2R_1.0.3               cli_3.6.1               spatstat.sparse_3.1-0   textshaping_0.3.6       fansi_1.0.4             S4Arrays_1.0.6          svglite_2.1.1          
# [82] uwot_0.1.16             zipfR_0.6-70            gtable_0.3.4            imager_1.0.2            rstatix_0.7.2           farver_2.1.1            lifecycle_1.0.3         statmod_1.5.0           MASS_7.3-60.2    
# 