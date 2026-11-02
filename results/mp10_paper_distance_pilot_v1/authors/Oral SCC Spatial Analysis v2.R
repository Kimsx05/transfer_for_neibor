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
library(patchwork)

#ORAL SCC
oral <- Seurat::Load10X_Spatial(data.dir = "/facilities/sequencers/seq.08/230921_A01065_0249_BHLKN3DMXY/results/pipe00/Unaligned/FASTQS/spaceranger/AV5_11_12882_22_V13F06-301/outs", 
                                filename = "filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "oral",
                                filter.matrix = TRUE, 
                                to.upper = FALSE)
oral@meta.data$orig.ident = "oral"
oral = subset(oral, subset = nCount_Spatial > 100)
oral = AddMetaData(oral, as.data.frame(oral1_percentages))

oral <- SCTransform(oral, assay = "Spatial")
oral <- RunPCA(oral, assay = "SCT")
oral <- FindNeighbors(oral, dims = 1:30)
oral <- FindClusters(oral, resolution = 0.6)
oral <- RunUMAP(oral, reduction = "pca", dims = 1:30)

SpatialDimPlot(object = oral, label = T)

oral@meta.data$Epithelial = oral$Epithelial.1 + oral$Epithelial.2 + oral$Epithelial.3 + oral$Epithelial.4 + oral$Epithelial.5 + oral$Epithelial.6
oral@meta.data$Fibroblasts = oral$Fibroblasts.1 + oral$Fibroblasts.2
oral@meta.data$Salivary.gland = oral$Salivary.gland.1 + oral$Salivary.gland.2

#Heatmap of clusters and cell signatures
signatures <- c("Macrophage", "Epithelial", "Fibroblasts", "Proliferating", "Plasma", "Salivary.gland")
metadata = oral@meta.data

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

#Subset fibroblasts on the tumour invasive front
oral.coords = GetTissueCoordinates(oral)
oral@meta.data$image.row = oral.coords$x
oral@meta.data$image.col = oral.coords$y

IF.fb = subset(oral, seurat_clusters %in% c("2","11","13") & image.row <18000)
IF.fb = subset(IF.fb, image.row >15000  & image.col > 19000, invert = T)
IF_fibroblast_spots <- rownames(IF.fb@meta.data)

# Extract spatial coordinates
spatial_coords <- GetTissueCoordinates(oral)  

# Function to calculate Euclidean distance
euclidean_distance <- function(coord1, coord2) {
  sqrt(sum((coord1 - coord2)^2))
}

# Calculate distances
distances <- sapply(spatial_coords$cell, function(spot) {
  spot_coord <- spatial_coords[spatial_coords$cell == spot, c("x", "y")]
  min_distance <- min(apply(spatial_coords[spatial_coords$cell %in% IF_fibroblast_spots, c("x", "y")], 1, function(x) euclidean_distance(as.numeric(spot_coord), x)))
  return(min_distance)
})

# Add distances to metadata
oral@meta.data$Distance_to_IF_Fibroblasts <- distances

IF_spots <- colnames(IF.fb)
oral@meta.data$IF_fibroblasts = colnames(oral) %in% IF_spots
oral@meta.data$IF_fibroblasts = factor(oral@meta.data$IF_fibroblasts, levels = c(TRUE, FALSE), labels = c("IF_FB", "Other"))

TNF_EMT = list(EMT = c("CXCL2", "SQSTM1", "PRSS2", "IER2", "IL1B", "IL1A", "FERMT2", "ZFP36", "JUN", "JUNB", "AREG", "DUSP1", "FOSB", "SERPINE1", "FOS"))
oral = AddModuleScore(oral, features = TNF_EMT)
names(oral@meta.data)[28] = "TNF_EMT"

S1P = list(S1P = c("CCN2", "CD274", "CXCL2", "EGR1", "ETS1", "FN1", "FOS", "IL6",
                   "JUN", "JUNB", "MMP9", "NCAM1", "PIEZO2", "PTGS2", "RXRA", "SERPINE1",
                   "TAGLN", "THBS1"))
oral = AddModuleScore(oral, features = S1P)
names(oral@meta.data)[29] = "S1P"


#plot score in epithelial clusters
SpatialFeaturePlot(subset(oral, seurat_clusters %in% c(0,1,4,6,12)), "TNF_EMT", pt.size.factor = 2.3, crop = F)
SpatialFeaturePlot(subset(oral, seurat_clusters %in% c(0,1,4,6,12)), "S1P", pt.size.factor = 2.3, crop = F)

tg.chol = list(triglyceride.catab = c("ABHD5", "CAV1", "FABP1", "FABP12", "FABP2", "FABP3", "FABP4", "FABP5", "FABP6", "FABP7",
                                      "FABP9", "GPD2", "LIPE", "MGLL", "PLIN1", "PLIN3", "PNPLA4", "PNPLA5", "PPP1CA", "PPP1CB",
                                      "PPP1CC", "PRKACA", "PRKACB", "PRKACG"),
               cholesterol.synth = c("ACAT2", "ARV1", "CYP51A1", "DHCR24", "DHCR7", "EBP", "FDFT1", "FDPS", "GGPS1",
                                     "HMGCR", "HMGCS1", "HSD17B7", "IDI1", "IDI2", "LBR", "LSS", "MSMO1", "MVD",
                                     "MVK", "NSDHL", "PLPP6", "PMVK", "SC5D", "SQLE", "SREBF1", "SREBF2", "TM7SF2"))

oral = AddModuleScore(oral, tg.chol)
names(oral@meta.data)[30:31] = names(tg.chol)


#ORAL SCC - 2
oral2 <- Seurat::Load10X_Spatial(data.dir = "/facilities/sequencers/seq.08/230921_A01065_0249_BHLKN3DMXY/results/pipe00/Unaligned/FASTQS/spaceranger/AV5_12_23CH03929_V13F06-301/outs", 
                                 filename = "filtered_feature_bc_matrix.h5",
                                 assay = "Spatial",
                                 slice = "oral2",
                                 filter.matrix = TRUE, 
                                 to.upper = FALSE)
oral2@meta.data$orig.ident = "oral2"
oral2 = subset(oral2, subset = nCount_Spatial > 100)
oral2 = AddMetaData(oral2, as.data.frame(oral2_percentages))

oral2 <- SCTransform(oral2, assay = "Spatial")
oral2 <- RunPCA(oral2, assay = "SCT")
oral2 <- FindNeighbors(oral2, dims = 1:30)
oral2 <- FindClusters(oral2, resolution = 0.6)
oral2 <- RunUMAP(oral2, reduction = "pca", dims = 1:30)

SpatialDimPlot(object = oral2, label = T, pt.size.factor = 2)

oral2@meta.data$Epithelial = oral2$Epithelial.1 + oral2$Epithelial.2 + oral2$Epithelial.3 + oral2$Epithelial.4 + oral2$Epithelial.5 + oral2$Epithelial.6 + oral2$Epithelial.7 + oral2$Proliferating.epithelial

#Heatmap of clusters and cell signatures
signatures <- c("Macrophage", "Epithelial", "Fibroblasts", "Smooth.muscle", "Plasma")
metadata = oral2@meta.data

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
oral2 = AddModuleScore(oral2, features = TNF_EMT)
names(oral2@meta.data)[21] = "TNF_EMT"

S1P = list(S1P = c("CCN2", "CD274", "CXCL2", "EGR1", "ETS1", "FN1", "FOS", "IL6",
                   "JUN", "JUNB", "MMP9", "NCAM1", "PIEZO2", "PTGS2", "RXRA", "SERPINE1",
                   "TAGLN", "THBS1"))
oral2 = AddModuleScore(oral2, features = S1P)
names(oral2@meta.data)[22] = "S1P"

oral2$tissue.group = ifelse(oral2$Fibroblasts >= 0.2, "fibroblast.enriched", ifelse(oral2$Fibroblasts <0.2 & oral2$Fibroblasts > 0, "fibroblast.low", "fibroblast.absent"))

oral2 = AddModuleScore(oral2, tg.chol)
names(oral2@meta.data)[24:25] = names(tg.chol)

##single cell sequencing for fibroblast populations
rna = read.table(gzfile("/data/skincan/tbudden/Spatial transcriptomics files/SCC spatial/GSE181919_UMI_counts.txt.gz"),sep="\t")
GSE181919_Barcode_metadata <- read.table(gzfile("/data/skincan/tbudden/Spatial transcriptomics files/SCC spatial/GSE181919_Barcode_metadata.txt.gz"),sep="\t")
GSE181919_Barcode_metadata$orig.ident = rownames(GSE181919_Barcode_metadata)
GSE181919_Barcode_metadata$orig.ident = gsub('-', '\\.', GSE181919_Barcode_metadata$orig.ident)
rownames(GSE181919_Barcode_metadata) = GSE181919_Barcode_metadata$orig.ident

meta = GSE181919_Barcode_metadata[names(rna),]
hnscc <-CreateSeuratObject(counts = rna, meta.data = meta)

hnscc[["percent.mt"]] <- PercentageFeatureSet(hnscc, pattern = "^MT-")
hnscc <- subset(hnscc, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)

hnscc <- SCTransform(hnscc, vars.to.regress = "percent.mt")
hnscc <- RunPCA(hnscc)
hnscc <- RunHarmony(hnscc, group.by.vars = "sample.id")
hnscc <- RunUMAP(hnscc, reduction = "harmony", dims = 1:30)
hnscc <- FindNeighbors(hnscc, reduction = "harmony", dims = 1:30)
hnscc <- FindClusters(hnscc, resolution = 0.2)

DimPlot(hnscc, reduction = "umap", label = T)

hnscc.markers <- FindAllMarkers(hnscc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

hnscc.fb = subset(hnscc, seurat_clusters == 0)

load("~/cycle.rda")
hnscc.fb = CellCycleScoring(hnscc.fb, g2m.features = g2m_genes, s.features = s_genes)
hnscc.fb = SCTransform(hnscc.fb, vars.to.regress = c("S.Score", "G2M.Score"))
hnscc.fb <- RunPCA(hnscc.fb)
hnscc.fb <- RunHarmony(hnscc.fb, group.by.vars = "sample.id")
hnscc.fb <- FindNeighbors(hnscc.fb, dims = 1:30, reduction = "harmony")
hnscc.fb <- FindClusters(hnscc.fb, resolution = 0.6)
hnscc.fb <- RunUMAP(hnscc.fb, dims = 1:30, reduction = "harmony")
DimPlot(hnscc.fb, label = T)

hnscc.fb.markers <- FindAllMarkers(hnscc.fb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

hnscc.fb = subset(hnscc.fb, seurat_clusters %in% c(9,10), invert = T) #repeat clustering analysis

DimPlot(hnscc.fb, group.by = "tissue.type")

metadata = hnscc.fb@meta.data

frequency_data <- metadata %>%
  group_by(seurat_clusters, tissue.type) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(tissue.type) %>%
  mutate(frequency = count / sum(count)) %>%
  ungroup()

ggplot(frequency_data, aes(x = tissue.type, y = frequency, fill = seurat_clusters)) +
  geom_bar(stat = "identity")

hnscc.fb.signatures <- hnscc.fb.markers %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC) %>%  # Select top 6 based on log fold change
  arrange(cluster, desc(avg_log2FC)) %>%  # Ensure they are sorted
  group_split() %>%  # Split into a list
  setNames(unique(hnscc.fb.markers$cluster))

hnscc.fb.signatures <- lapply(hnscc.fb.signatures, function(df) df$gene)


oral = AddModuleScore(oral, features = hnscc.fb.signatures)
names(oral@meta.data)[32:38] = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "FB_6")

oral2 = AddModuleScore(oral2, features = hnscc.fb.signatures)
names(oral2@meta.data)[26:32] = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "FB_6")

fb_signatures <- c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "FB_6")

#Oral
# Scale the module scores (z-score normalization)
oral@meta.data[fb_signatures] <- scale(oral@meta.data[fb_signatures])
expr_values <- FetchData(oral, vars = fb_signatures)

# Calculate the global min and max
min_val <- min(expr_values, na.rm = TRUE)
max_val <- max(expr_values, na.rm = TRUE)
# Generate individual plots
plots <- SpatialFeaturePlot(subset(oral, seurat_clusters %in% c("2","11","13")), c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "FB_6"), combine = F, pt.size.factor = 2)

# Apply same color scale to each plot (e.g., white to red, limits 0 to 5)
plots <- lapply(plots, function(p) {
  p + scale_fill_viridis_c(option = "inferno", direction = 1, limits = c(min_val, max_val))
})

# Combine them again
wrap_plots(plots, ncol = 4)

#Oral2
# Scale the module scores (z-score normalization)
oral2@meta.data[fb_signatures] <- scale(oral2@meta.data[fb_signatures])
expr_values <- FetchData(oral2, vars = fb_signatures)

# Calculate the global min and max
min_val <- min(expr_values, na.rm = TRUE)
max_val <- max(expr_values, na.rm = TRUE)
# Generate individual plots
plots <- SpatialFeaturePlot(subset(oral2, seurat_clusters %in% c(1,2,5,8)), c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "FB_6"), combine = F, pt.size.factor = 2)

# Apply same color scale to each plot (e.g., white to red, limits 0 to 5)
plots <- lapply(plots, function(p) {
  p + scale_fill_viridis_c(option = "inferno", direction = 1, limits = c(min_val, max_val))
})

# Combine them again
wrap_plots(plots, ncol = 4)