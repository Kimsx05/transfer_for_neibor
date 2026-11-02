library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(textshape)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)

#####cutaneous SCC
cut1 <- Seurat::Load10X_Spatial(data.dir = "/facilities/sequencers/seq.08/230921_A01065_0249_BHLKN3DMXY/results/pipe00/Unaligned/FASTQS/spaceranger/AV5_09_D6531_15_A1_V13F06-301/outs", 
                                filename = "filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "cut1",
                                filter.matrix = TRUE, 
                                to.upper = FALSE)
cut1@meta.data$orig.ident = "cut1"
cut1 = subset(cut1, subset = nCount_Spatial > 100)
load("/data/skincan/tbudden/Cutaneous SCC 1 topic percentages.Rdata")
cut1 = AddMetaData(cut1, as.data.frame(cut1_percentages))

cut1 <- SCTransform(cut1, assay = "Spatial")
cut1 <- RunPCA(cut1, assay = "SCT")
cut1 <- FindNeighbors(cut1, dims = 1:30)
cut1 <- FindClusters(cut1, resolution = 0.6)
cut1 <- RunUMAP(cut1, reduction = "pca", dims = 1:30)

cut1@meta.data$Epithelial = cut1$Epithelial.1 + cut1$Epithelial.2 + cut1$Epithelial.3 + cut1$Epithelial.4 + cut1$Epithelial.5 + cut1$Epithelial.6 + cut1$Epithelial.7 + cut1$Epithelial.8 + cut1$Epithelial.9 + cut1$Epithelial.10
cut1@meta.data$Fibroblasts = cut1$Fibroblasts.1 + cut1$Fibroblasts.2
cut1@meta.data$Proliferative = cut1$Proliferative.1 + cut1$Proliferative.2

#Heatmap of clusters and cell signatures
signatures <- c("Proliferative", "Epithelial", "Fibroblasts", "Immune", "Plasma", "Endothelial", "Smooth.muscle")
metadata = cut1@meta.data

# Calculate average signature scores per cluster
average_scores <- metadata %>%
  group_by(seurat_clusters) %>%
  summarize(across(all_of(signatures), mean, na.rm = TRUE)) %>%
  column_to_rownames("seurat_clusters")

heatmap_data <- as.matrix(average_scores)

# Optionally, scale the data if needed
heatmap_data <- scale(heatmap_data)

pheatmap(heatmap_data,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("#2166ac", '#67a9cf','#d1e5f0', "white", '#fddbc7','#ef8a62',  "#b2182b"))(100), # Customize colors
         fontsize_row = 10,      # Adjust row font size
         fontsize_col = 10,
         scale = "row",
         border_color = NA) 

TNF_EMT = list(EMT = c("CXCL2", "SQSTM1", "PRSS2", "IER2", "IL1B", "IL1A", "FERMT2", "ZFP36", "JUN", "JUNB", "AREG", "DUSP1", "FOSB", "SERPINE1", "FOS"))
cut1 = AddModuleScore(cut1, features = TNF_EMT)
names(cut1@meta.data)[29] = "TNF_EMT"

S1P = list(S1P = c("CCN2", "CD274", "CXCL2", "EGR1", "ETS1", "FN1", "FOS", "IL6",
                   "JUN", "JUNB", "MMP9", "NCAM1", "PIEZO2", "PTGS2", "RXRA", "SERPINE1",
                   "TAGLN", "THBS1"))
cut1 = AddModuleScore(cut1, features = S1P)
names(cut1@meta.data)[30] = "S1P"

tg.chol = list(triglyceride.catab = c("ABHD5", "CAV1", "FABP1", "FABP12", "FABP2", "FABP3", "FABP4", "FABP5", "FABP6", "FABP7",
                                      "FABP9", "GPD2", "LIPE", "MGLL", "PLIN1", "PLIN3", "PNPLA4", "PNPLA5", "PPP1CA", "PPP1CB",
                                      "PPP1CC", "PRKACA", "PRKACB", "PRKACG"),
               cholesterol.synth = c("ACAT2", "ARV1", "CYP51A1", "DHCR24", "DHCR7", "EBP", "FDFT1", "FDPS", "GGPS1",
                                     "HMGCR", "HMGCS1", "HSD17B7", "IDI1", "IDI2", "LBR", "LSS", "MSMO1", "MVD",
                                     "MVK", "NSDHL", "PLPP6", "PMVK", "SC5D", "SQLE", "SREBF1", "SREBF2", "TM7SF2"))

cut1 = AddModuleScore(cut1, tg.chol)
names(cut1@meta.data)[31:32] = names(tg.chol)

#Subset fibroblasts on the tumour invasive front
cut1.fb = subset(cut1, seurat_clusters %in% c("1","7"))
cut1_fibroblast_spots <- rownames(cut1.fb@meta.data)

# Extract spatial coordinates
spatial_coords <- GetTissueCoordinates(cut1)  

# Function to calculate Euclidean distance
euclidean_distance <- function(coord1, coord2) {
  sqrt(sum((coord1 - coord2)^2))
}

# Calculate distances
distances <- sapply(spatial_coords$cell, function(spot) {
  spot_coord <- spatial_coords[spatial_coords$cell == spot, c("x", "y")]
  min_distance <- min(apply(spatial_coords[spatial_coords$cell %in% cut1_fibroblast_spots, c("x", "y")], 1, function(x) euclidean_distance(as.numeric(spot_coord), x)))
  return(min_distance)
})

# Add distances to metadata
cut1@meta.data$Distance_to_Fibroblasts <- distances

FeatureScatter(subset(cut1, Epithelial > 0.2), "Distance_to_Fibroblasts", "triglyceride.catab")



#second cutaneous
cut2 <- Seurat::Load10X_Spatial(data.dir = "/facilities/sequencers/seq.08/230921_A01065_0249_BHLKN3DMXY/results/pipe00/Unaligned/FASTQS/spaceranger/AV5_10_D02_17_A4_V13F06-301/outs", 
                                filename = "filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "cut1",
                                filter.matrix = TRUE, 
                                to.upper = FALSE)
cut2@meta.data$orig.ident = "cut2"
cut2 = subset(cut2, subset = nCount_Spatial > 100)
load("/data/skincan/tbudden/Cutaneous SCC 2 topic percentages.Rdata")
cut2 = AddMetaData(cut2, as.data.frame(cut2_percentages))

cut2 <- SCTransform(cut2, assay = "Spatial")
cut2 <- RunPCA(cut2, assay = "SCT")
cut2 <- FindNeighbors(cut2, dims = 1:30)
cut2 <- FindClusters(cut2, resolution = 0.6)
cut2 <- RunUMAP(cut2, reduction = "pca", dims = 1:30)

SpatialDimPlot(cut2, pt.size.factor = 2)

cut2@meta.data$Epithelial = cut2$Epithelial.1 + cut2$Epithelial.2 + cut2$Epithelial.3 + cut2$Epithelial.4 + cut2$Epithelial.5 + cut2$Epithelial.6 + cut2$Epithelial.7 + cut2$Epithelial.8 + cut2$Epithelial.9 + cut2$Epithelial.immune
cut2@meta.data$Stromal = cut2$Stromal.1 + cut2$Stromal.2 + cut2$Stromal.3
cut2@meta.data$Immune = cut2$Immune.1 + cut2$Immune.2
cut2@meta.data$Plasma = cut2$Plasma.1 + cut2$Plasma.2

#Heatmap of clusters and cell signatures
signatures <- c("Hair.follicle", "Epithelial", "Stromal", "Immune", "Plasma", "Melanocytes", "Muscle", "Interferon.signalling", "Endothelial")
metadata = cut2@meta.data

# Calculate average signature scores per cluster
average_scores <- metadata %>%
  group_by(seurat_clusters) %>%
  summarize(across(all_of(signatures), mean, na.rm = TRUE)) %>%
  column_to_rownames("seurat_clusters")

heatmap_data <- as.matrix(average_scores)

# Optionally, scale the data if needed
heatmap_data <- scale(heatmap_data)

pheatmap(heatmap_data,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("#2166ac", '#67a9cf','#d1e5f0', "white", '#fddbc7','#ef8a62',  "#b2182b"))(100), # Customize colors
         fontsize_row = 10,      # Adjust row font size
         fontsize_col = 10,
         scale = "row",
         border_color = NA) 


cut2 = AddModuleScore(cut2, features = TNF_EMT)
names(cut2@meta.data)[34] = "TNF_EMT"

cut2 = AddModuleScore(cut2, features = S1P)
names(cut2@meta.data)[35] = "S1P"

cut2 = AddModuleScore(cut2, tg.chol)
names(cut2@meta.data)[36:37] = names(tg.chol)


#Subset fibroblasts on the tumour invasive front
cut2.fb = subset(cut2, seurat_clusters %in% c("2","5","7"))
cut2_fibroblast_spots <- rownames(cut2.fb@meta.data)

# Extract spatial coordinates
spatial_coords <- GetTissueCoordinates(cut2)  

# Function to calculate Euclidean distance
euclidean_distance <- function(coord1, coord2) {
  sqrt(sum((coord1 - coord2)^2))
}

# Calculate distances
distances <- sapply(spatial_coords$cell, function(spot) {
  spot_coord <- spatial_coords[spatial_coords$cell == spot, c("x", "y")]
  min_distance <- min(apply(spatial_coords[spatial_coords$cell %in% cut2_fibroblast_spots, c("x", "y")], 1, function(x) euclidean_distance(as.numeric(spot_coord), x)))
  return(min_distance)
})

# Add distances to metadata
cut2@meta.data$Distance_to_Fibroblasts <- distances

FeatureScatter(subset(cut2, Epithelial > 0.2), "Distance_to_Fibroblasts", "S1P")


####CSCC single cell
GSM6735856_825_6 <- readRDS("/data/skincan/tbudden/Skin single cell/GSM6735856_825_6.rds")
GSM6735856_825_6 = UpdateSeuratObject(object = GSM6735856_825_6)
GSM6735856_825_6$subj = "UV1"
GSM6735856_825_6$age = "OLD"
GSM6735856_825_6$UV = "UV"
GSM6735856_825_6$tum.norm = "Normal"
GSM6735856_825_6[["percent.mt"]] <- PercentageFeatureSet(GSM6735856_825_6, pattern = "^MT-")

GSM6735857_865_7 <- readRDS("/data/skincan/tbudden/Skin single cell/GSM6735857_865_7.rds")
GSM6735857_865_7 = UpdateSeuratObject(object = GSM6735857_865_7)
GSM6735857_865_7$subj = "UV2"
GSM6735857_865_7$age = "OLD"
GSM6735857_865_7$UV = "UV"
GSM6735857_865_7$tum.norm = "Normal"
GSM6735857_865_7[["percent.mt"]] <- PercentageFeatureSet(GSM6735857_865_7, pattern = "^MT-")

GSM6735858_Burow_874 <- readRDS("/data/skincan/tbudden/Skin single cell/GSM6735858_Burow_874.rds")
GSM6735858_Burow_874 = UpdateSeuratObject(object = GSM6735858_Burow_874)
GSM6735858_Burow_874$subj = "UV3"
GSM6735858_Burow_874$age = "OLD"
GSM6735858_Burow_874$UV = "UV"
GSM6735858_Burow_874$tum.norm = "Normal"
GSM6735858_Burow_874[["percent.mt"]] <- PercentageFeatureSet(GSM6735858_Burow_874, pattern = "^MT-")

GSM6735862_cSCC_599 =  readRDS("/data/skincan/tbudden/Skin single cell/GSM6735862_cSCC_599.rds")
GSM6735862_cSCC_599 = UpdateSeuratObject(object = GSM6735862_cSCC_599)
GSM6735862_cSCC_599$subj = "UV4"
GSM6735862_cSCC_599$age = "OLD"
GSM6735862_cSCC_599$UV = "UV"
GSM6735862_cSCC_599$tum.norm = "Tumour"
GSM6735862_cSCC_599[["percent.mt"]] <- PercentageFeatureSet(GSM6735862_cSCC_599, pattern = "^MT-")

GSM6735863_650_8 =  readRDS("/data/skincan/tbudden/Skin single cell/GSM6735863_650_8.rds")
GSM6735863_650_8 = UpdateSeuratObject(object = GSM6735863_650_8)
GSM6735863_650_8$subj = "UV5"
GSM6735863_650_8$age = "OLD"
GSM6735863_650_8$UV = "UV"
GSM6735863_650_8$tum.norm = "Tumour"
GSM6735863_650_8[["percent.mt"]] <- PercentageFeatureSet(GSM6735863_650_8, pattern = "^MT-")

GSM6735864_719_9 =  readRDS("/data/skincan/tbudden/Skin single cell/GSM6735864_719_9.rds")
GSM6735864_719_9 = UpdateSeuratObject(object = GSM6735864_719_9)
GSM6735864_719_9$subj = "UV6"
GSM6735864_719_9$age = "OLD"
GSM6735864_719_9$UV = "UV"
GSM6735864_719_9$tum.norm = "Tumour"
GSM6735864_719_9[["percent.mt"]] <- PercentageFeatureSet(GSM6735864_719_9, pattern = "^MT-")

GSM6735865_743_6 =  readRDS("/data/skincan/tbudden/Skin single cell/GSM6735865_743_6.rds")
GSM6735865_743_6 = UpdateSeuratObject(object = GSM6735865_743_6)
GSM6735865_743_6$subj = "UV7"
GSM6735865_743_6$age = "OLD"
GSM6735865_743_6$UV = "UV"
GSM6735865_743_6$tum.norm = "Tumour"
GSM6735865_743_6[["percent.mt"]] <- PercentageFeatureSet(GSM6735865_743_6, pattern = "^MT-")

GSM6735866_800_9 =  readRDS("/data/skincan/tbudden/Skin single cell/GSM6735866_800_9.rds")
GSM6735866_800_9 = UpdateSeuratObject(object = GSM6735866_800_9)
GSM6735866_800_9$subj = "UV8"
GSM6735866_800_9$age = "OLD"
GSM6735866_800_9$UV = "UV"
GSM6735866_800_9$tum.norm = "Tumour"
GSM6735866_800_9[["percent.mt"]] <- PercentageFeatureSet(GSM6735866_800_9, pattern = "^MT-")

GSM6735859_483 =  readRDS("/data/skincan/tbudden/Skin single cell/GSM6735859_483-0.rds")
GSM6735859_483 = UpdateSeuratObject(object = GSM6735859_483)
GSM6735859_483$subj = "UV9"
GSM6735859_483$age = "OLD"
GSM6735859_483$UV = "UV"
GSM6735859_483$tum.norm = "Pre"
GSM6735859_483[["percent.mt"]] <- PercentageFeatureSet(GSM6735859_483, pattern = "^MT-")

GSM6735861_290 =  readRDS("/data/skincan/tbudden/Skin single cell/GSM6735861_290-B.rds")
GSM6735861_290 = UpdateSeuratObject(object = GSM6735861_290)
GSM6735861_290$subj = "UV10"
GSM6735861_290$age = "OLD"
GSM6735861_290$UV = "UV"
GSM6735861_290$tum.norm = "Pre"
GSM6735861_290[["percent.mt"]] <- PercentageFeatureSet(GSM6735861_290, pattern = "^MT-")
GSM6735861_290 = subset(GSM6735861_290, nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

GSM6735860_290 =  readRDS("/data/skincan/tbudden/Skin single cell/GSM6735860_290-A.rds")
GSM6735860_290 = UpdateSeuratObject(object = GSM6735860_290)
GSM6735860_290$subj = "UV11"
GSM6735860_290$age = "OLD"
GSM6735860_290$UV = "UV"
GSM6735860_290$tum.norm = "Pre"
GSM6735860_290[["percent.mt"]] <- PercentageFeatureSet(GSM6735860_290, pattern = "^MT-")

UV.skin = merge(x = GSM6735856_825_6, y = c(GSM6735857_865_7,
                                            GSM6735858_Burow_874,
                                            GSM6735862_cSCC_599,
                                            GSM6735863_650_8,
                                            GSM6735864_719_9,
                                            GSM6735865_743_6,
                                            GSM6735866_800_9,
                                            GSM6735859_483,
                                            GSM6735861_290,
                                            GSM6735860_290))

UV.skin = subset(UV.skin, nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)

UV.skin = SCTransform(UV.skin)
UV.skin = RunPCA(UV.skin, assay = "SCT")
UV.skin = RunHarmony(UV.skin, group.by.vars = "subj", assay.use = "SCT", plot_convergence = TRUE)
UV.skin = FindNeighbors(UV.skin, reduction = "harmony", dims = 1:30)
UV.skin = FindClusters(UV.skin, verbose = FALSE, resolution = 0.2)
UV.skin = RunUMAP(UV.skin, reduction = "harmony", dims = 1:30)

DimPlot(UV.skin, group.by = c("seurat_clusters", "tum.norm"), label = T)

UV.skin.markers = FindAllMarkers(UV.skin, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

UV.fb = subset(UV.skin, seurat_clusters %in% c(2,5,8))
load("~/cycle.rda")
UV.fb <- CellCycleScoring(UV.fb, 
                          g2m.features = g2m_genes, 
                          s.features = s_genes)
UV.fb = SCTransform(UV.fb, vars.to.regress = c("S.Score", "G2M.Score"))
UV.fb = RunPCA(UV.fb, assay = "SCT")
UV.fb = RunHarmony(UV.fb, group.by.vars = "subj", assay.use = "SCT", plot_convergence = TRUE)
UV.fb = FindNeighbors(UV.fb, reduction = "harmony", dims = 1:30)
UV.fb = FindClusters(UV.fb, verbose = FALSE, resolution = 0.2)
UV.fb = RunUMAP(UV.fb, reduction = "harmony", dims = 1:30)

DimPlot(UV.fb, group.by = c("seurat_clusters", "tum.norm"), label = T)
  
UV.fb.markers = FindAllMarkers(UV.fb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#remove non-fibroblast clusters
UV.fb = subset(UV.fb, seurat_clusters %in% c(7, 8, 9, 12, 14, 15), invert = T) #repeat clustering
UV.fb.markers = FindAllMarkers(UV.fb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

metadata = UV.fb@meta.data

frequency_data <- metadata %>%
  group_by(seurat_clusters, tum.norm) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(tum.norm) %>%
  mutate(frequency = count / sum(count)) %>%
  ungroup()

ggplot(frequency_data, aes(x = tum.norm, y = frequency, fill = seurat_clusters)) +
  geom_bar(stat = "identity")

cscc.fb.signatures <- UV.fb.markers %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) %>%  # Select top 10 based on log fold change
  arrange(cluster, desc(avg_log2FC)) %>%  # Ensure they are sorted
  group_split() %>%  # Split into a list
  setNames(unique(UV.fb.markers$cluster))

cscc.fb.signatures <- lapply(cscc.fb.signatures, function(df) df$gene)

cut1 = AddModuleScore(cut1, cscc.fb.signatures)
names(cut1@meta.data)[34:44] = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "FB_6", "FB_7", "FB_8", "FB_9", "FB_10")

cut2 = AddModuleScore(cut2, cscc.fb.signatures)
names(cut2@meta.data)[39:49] = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "FB_6", "FB_7", "FB_8", "FB_9", "FB_10")



fb_signatures <- c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "FB_6", "FB_7", "FB_8", "FB_9", "FB_10")
# Scale the module scores (z-score normalization)
cut1@meta.data[fb_signatures] <- scale(cut1@meta.data[fb_signatures])
expr_values <- FetchData(cut1, vars = fb_signatures)

# Calculate the global min and max
min_val <- min(expr_values, na.rm = TRUE)
max_val <- max(expr_values, na.rm = TRUE)
# Generate individual plots
plots <- SpatialFeaturePlot(subset(cut1, seurat_clusters %in% c(0,1,7,11)), c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "FB_6", "FB_7", "FB_8", "FB_9", "FB_10"), combine = F, crop = F, pt.size.factor = 2)

# Apply same color scale to each plot (e.g., white to red, limits 0 to 5)
plots <- lapply(plots, function(p) {
  p + scale_fill_viridis_c(option = "inferno", direction = 1, limits = c(min_val, max_val))
})

# Combine them again
wrap_plots(plots, ncol = 6)



# Scale the module scores (z-score normalization)
cut2@meta.data[fb_signatures] <- scale(cut2@meta.data[fb_signatures])
expr_values <- FetchData(cut2, vars = fb_signatures)

# Calculate the global min and max
min_val <- min(expr_values, na.rm = TRUE)
max_val <- max(expr_values, na.rm = TRUE)
# Generate individual plots
plots <- SpatialFeaturePlot(subset(cut2, seurat_clusters %in% c(8,10,5,6,2,7)), c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "FB_6", "FB_7", "FB_8", "FB_9", "FB_10"), combine = F, crop = F, pt.size.factor = 2)

# Apply same color scale to each plot (e.g., white to red, limits 0 to 5)
plots <- lapply(plots, function(p) {
  p + scale_fill_viridis_c(option = "inferno", direction = 1, limits = c(min_val, max_val))
})

# Combine them again
wrap_plots(plots, ncol = 6)

