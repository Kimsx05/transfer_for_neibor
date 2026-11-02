library(Seurat)
library(dplyr)
library(ggplot2)
library(textshape)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
library(patchwork)
library(harmony)

#LUNG SCC
lung <- Seurat::Load10X_Spatial(data.dir = "/facilities/sequencers/seq.08/231218_A01065_0267_AHVCWWDRX3/results/pipe00/Aligned_AV7/spaceranger/AV7_03_W002415/outs", 
                                filename = "filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "lung",
                                filter.matrix = TRUE, 
                                to.upper = FALSE)
load("lung SCC topic percentages.Rdata")
lung@meta.data$orig.ident = "lung"
lung = subset(lung, subset = nCount_Spatial > 100)
lung = AddMetaData(lung, as.data.frame(lung_percentages))

lung <- SCTransform(lung, assay = "Spatial")
lung <- RunPCA(lung, assay = "SCT")
lung <- FindNeighbors(lung, dims = 1:30)
lung <- FindClusters(lung, resolution = 1.0)
lung <- RunUMAP(lung, reduction = "pca", dims = 1:30)

SpatialDimPlot(object = lung, label = T)

lung$Epithelial = lung$Epithelial.1 + lung$Epithelial.2 + lung$Epithelial.3 + lung$NRF2.activation
signatures <- c("Stromal", "Epithelial", "Alveolar", "Plasma", "Mixed", "Plueral", "Goblet")
metadata = lung@meta.data

# Calculate average signature scores per cluster
average_scores <- metadata %>%
  group_by(seurat_clusters) %>%
  summarize(across(all_of(signatures), mean, na.rm = TRUE)) %>%
  column_to_rownames("seurat_clusters")

heatmap_data <- as.matrix(average_scores)
#scale the data
heatmap_data <- scale(heatmap_data)

pheatmap(heatmap_data,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("#2166ac", '#67a9cf','#d1e5f0', "white", '#fddbc7','#ef8a62',  "#b2182b"))(100), # Customize colors
         fontsize_row = 10,      # Adjust row font size
         fontsize_col = 10,
         scale = "row",
         border_color = NA) 

#Define tumour and adjacent tissue regions
lung$Region = ifelse(lung$seurat_clusters %in% c(4,5,6,7,8,14), "Tumour", "Adjacent")
SpatialDimPlot(lung, group.by = "Region")

#Addmodule score for TG and cholesterol
tg.chol = list(triglyceride.catab = c("ABHD5", "CAV1", "FABP1", "FABP12", "FABP2", "FABP3", "FABP4", "FABP5", "FABP6", "FABP7",
                                      "FABP9", "GPD2", "LIPE", "MGLL", "PLIN1", "PLIN3", "PNPLA4", "PNPLA5", "PPP1CA", "PPP1CB",
                                      "PPP1CC", "PRKACA", "PRKACB", "PRKACG"),
               cholesterol.synth = c("ACAT2", "ARV1", "CYP51A1", "DHCR24", "DHCR7", "EBP", "FDFT1", "FDPS", "GGPS1",
                                     "HMGCR", "HMGCS1", "HSD17B7", "IDI1", "IDI2", "LBR", "LSS", "MSMO1", "MVD",
                                     "MVK", "NSDHL", "PLPP6", "PMVK", "SC5D", "SQLE", "SREBF1", "SREBF2", "TM7SF2"))

lung = AddModuleScore(lung, tg.chol)
names(lung@meta.data)[20:21] = names(tg.chol)

#Find marker genes for each cluster
lung.markers = FindAllMarkers(lung, only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)

#GSEA for each cluster
# Convert gene symbols to ENTREZ IDs
gene_conversion <- bitr(lung.markers$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Merge the conversion results with the original markers data
markers <- lung.markers %>%
  left_join(gene_conversion, by = c("gene" = "SYMBOL"))

# Filter out markers without ENTREZ IDs
markers <- markers %>% filter(!is.na(ENTREZID))

# Create a named list of ranked genes for each cluster
cluster_genes <- split(markers, markers$cluster)

# Rank genes by log fold change
ranked_genes <- lapply(cluster_genes, function(cluster) {
  genes <- cluster$ENTREZID
  scores <- cluster$avg_log2FC
  names(scores) <- genes
  return(sort(scores, decreasing = TRUE))
})

H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(H_t2g)

gsea_results_hall <- lapply(ranked_genes, function(ranked_gene_list) {
  gsea_result <- GSEA(geneList = ranked_gene_list,
                      pvalueCutoff = 0.05,
                      TERM2GENE = H_t2g)
  return(gsea_result)
})

#Plot top 10 enriched pathways in Cluster 7
dotplot(gsea_results_hall[["7"]], showCategory = 10)

TNF_EMT = list(EMT = c("CXCL2", "SQSTM1", "PRSS2", "IER2", "IL1B", "IL1A", "FERMT2", "ZFP36", "JUN", "JUNB", "AREG", "DUSP1", "FOSB", "SERPINE1", "FOS"))
S1P = list(S1P = c("CCN2", "CD274", "CXCL2", "EGR1", "ETS1", "FN1", "FOS", "IL6",
                   "JUN", "JUNB", "MMP9", "NCAM1", "PIEZO2", "PTGS2", "RXRA", "SERPINE1",
                   "TAGLN", "THBS1"))

#Oral pathways in lung sample
lung = AddModuleScore(lung, features = TNF_EMT)
names(lung@meta.data)[22] = "TNF_EMT"

lung = AddModuleScore(lung, features = S1P)
names(lung@meta.data)[23] = "S1P"

SpatialFeaturePlot(lung, c("S1P", "Pleural"), pt.size.factor = 2, crop = F)

plins = list(plins = c("PLIN2", "PLIN3", "PLIN5"))
lung = AddModuleScore(lung, plins)
names(lung@meta.data)[24] = "PLIN"

##lung scc single cell
load("/data/skincan/tbudden/Lungs/GSE148071 Lung SCC Seurat.RData")

lung.scc$orig.ident = "lung_scc"
lung.scc[["percent.mt"]] <- PercentageFeatureSet(lung.scc, pattern = "^MT-")

lung.scc = SCTransform(lung.scc, vars.to.regress = c("percent.mt"))
lung.scc <- RunPCA(lung.scc)
lung.scc <- RunHarmony(lung.scc, group.by.vars = "patient", umap.method = 'umap-learn', metric = 'correlation')
lung.scc <- FindNeighbors(lung.scc, dims = 1:30, reduction = "harmony")
lung.scc <- FindClusters(lung.scc, resolution = 0.2)
lung.scc <- RunUMAP(lung.scc, dims = 1:30, reduction = "harmony")


DimPlot(lung.scc, label = T)

markers = FindAllMarkers(lung.scc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

lung.fb = subset(lung.scc, seurat_clusters == 5)
load("~/cycle.rda")
lung.fb = CellCycleScoring(lung.fb, g2m.features = g2m_genes, s.features = s_genes)

lung.fb = SCTransform(lung.fb, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
lung.fb = RunPCA(lung.fb)
lung.fb = RunHarmony(lung.fb, "patient")

lung.fb = FindNeighbors(lung.fb, reduction = "harmony", dims = 1:30)
lung.fb = FindClusters(lung.fb, verbose = FALSE, resolution = 1)
lung.fb = RunUMAP(lung.fb, reduction = "harmony", dims = 1:30)

DimPlot(lung.fb, label = T)

lung.fb.markers = FindAllMarkers(lung.fb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

lung.fb = subset(lung.fb, seurat_clusters %in% c(5,9,10,12,14), invert = T)#repeat clustering

lung.fb = SCTransform(lung.fb, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
lung.fb = RunPCA(lung.fb)
lung.fb = RunHarmony(lung.fb, "patient")

lung.fb = FindNeighbors(lung.fb, reduction = "harmony", dims = 1:30)
lung.fb = FindClusters(lung.fb, verbose = FALSE, resolution = 0.6)
lung.fb = RunUMAP(lung.fb, reduction = "harmony", dims = 1:30)

DimPlot(lung.fb, label = T)

lung.fb.markers = FindAllMarkers(lung.fb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

lung.fb.signature <- lung.fb.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%  # Select top based on log fold change
  arrange(cluster, desc(avg_log2FC)) %>%  # Ensure they are sorted
  group_split() %>%  # Split into a list
  setNames(unique(lung.fb.markers$cluster))

lung.fb.signature <- lapply(lung.fb.signature, function(df) df$gene)

lung = AddModuleScore(lung, lung.fb.signature)
names(lung@meta.data)[25:31] = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "FB_6")

fb_signatures <- c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "FB_6")

lung@meta.data[fb_signatures] <- scale(lung@meta.data[fb_signatures])
expr_values <- FetchData(lung, vars = fb_signatures)

# Calculate the global min and max
min_val <- min(expr_values, na.rm = TRUE)
max_val <- max(expr_values, na.rm = TRUE)
# Generate individual plots
plots <- SpatialFeaturePlot(subset(lung, Stromal > 0), c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "FB_6"), combine = F, pt.size.factor = 2.2)

# Apply same color scale to each plot (e.g., white to red, limits 0 to 5)
plots <- lapply(plots, function(p) {
  p + scale_fill_viridis_c(option = "inferno", direction = 1)
})

# Combine them again
wrap_plots(plots, ncol = 4)


#LUNG SCC #2
lung2 <- Seurat::Load10X_Spatial(data.dir = "/data/skincan/tbudden/Lungs/E-MTAB-13530/P19/P19_T2", 
                                filename = "P19_T2-filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "lung2",
                                filter.matrix = TRUE, 
                                to.upper = FALSE)
load("lung2 SCC topic percentages.Rdata")
lung2@meta.data$orig.ident = "lung2"
lung2 = subset(lung2, subset = nCount_Spatial > 100)
lung2 = AddMetaData(lung2, as.data.frame(lung2_percentages))

lung2 <- SCTransform(lung2, assay = "Spatial")
lung2 <- RunPCA(lung2, assay = "SCT")
lung2 <- FindNeighbors(lung2, dims = 1:30)
lung2 <- FindClusters(lung2, resolution = 0.6)
lung2 <- RunUMAP(lung2, reduction = "pca", dims = 1:30)

SpatialDimPlot(object = lung2, label = T)

lung2@meta.data$Epithelial = lung2$Epithelial.1 + lung2$Epithelial.2 + lung2$Epithelial.3 + lung2$Epithelial.4 + lung2$Epithelial.5 + lung2$Epithelial.6 + lung2$Epithelial.7
lung2@meta.data$Mixed = lung2$Mixed.1 + lung2$Mixed.2
lung2@meta.data$Fibroblasts = lung2$Plasma.fibroblasts + lung2$Fibroblasts.1

signatures <- c("Mixed", "Epithelial", "Fibroblasts", "Endothelial", "Macrophage", "Alveolar")
metadata = lung2@meta.data

# Calculate average signature scores per cluster
average_scores <- metadata %>%
  group_by(seurat_clusters) %>%
  summarize(across(all_of(signatures), mean, na.rm = TRUE)) %>%
  column_to_rownames("seurat_clusters")

heatmap_data <- as.matrix(average_scores)
#scale the data
heatmap_data <- scale(heatmap_data)

pheatmap(heatmap_data,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("#2166ac", '#67a9cf','#d1e5f0', "white", '#fddbc7','#ef8a62',  "#b2182b"))(100), # Customize colors
         fontsize_row = 10,      # Adjust row font size
         fontsize_col = 10,
         scale = "row",
         border_color = NA) 

lung2$tumour = ifelse(lung2$seurat_clusters %in% c(0,1,3,6), "Tumour", "Stroma")

lung2 = AddModuleScore(lung2, tg.chol)
names(lung2@meta.data)[26:27] = names(tg.chol)

lung2$chol.group = ifelse(lung2$seurat_clusters %in% c(1,3), "Epithelial.high.chol", ifelse(lung2$seurat_clusters %in% c(0,6), "Epithelial.low.chol", "Stroma"))


#Find marker genes for each cluster
lung2.markers = FindAllMarkers(lung2, only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)

gene_conversion <- bitr(lung2.markers$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Merge the conversion results with the original markers data
markers <- lung2.markers %>%
  left_join(gene_conversion, by = c("gene" = "SYMBOL"))

# Filter out markers without ENTREZ IDs
markers <- markers %>% filter(!is.na(ENTREZID))

# Create a named list of ranked genes for each cluster
cluster_genes <- split(markers, markers$cluster)

# Rank genes by log fold change
ranked_genes <- lapply(cluster_genes, function(cluster) {
  genes <- cluster$ENTREZID
  scores <- cluster$avg_log2FC
  names(scores) <- genes
  return(sort(scores, decreasing = TRUE))
})

H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(H_t2g)

gsea_results_hall <- lapply(ranked_genes, function(ranked_gene_list) {
  gsea_result <- GSEA(geneList = ranked_gene_list,
                      pvalueCutoff = 0.05,
                      TERM2GENE = H_t2g)
  return(gsea_result)
})

#Plot top 10 enriched pathways in Cluster 7
dotplot(gsea_results_hall[["1"]], showCategory = 10)
dotplot(gsea_results_hall[["3"]], showCategory = 10)

lung2 = AddModuleScore(lung2, plins)
names(lung2@meta.data)[29] = "PLIN"


#Oral pathways in lung sample
lung2 = AddModuleScore(lung2, features = TNF_EMT)
names(lung2@meta.data)[30] = "TNF_EMT"

lung2 = AddModuleScore(lung2, features = S1P)
names(lung2@meta.data)[31] = "S1P"

SpatialFeaturePlot(lung2, c("TNF_EMT", "S1P"), pt.size.factor = 2)






###scRNA fb signatures
lung2 = AddModuleScore(lung2, lung.fb.signature)
names(lung2@meta.data)[32:38] = c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "FB_6")

fb_signatures <- c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "FB_6")

lung2@meta.data[fb_signatures] <- scale(lung2@meta.data[fb_signatures])
expr_values <- FetchData(lung2, vars = fb_signatures)

# Calculate the global min and max
min_val <- min(expr_values, na.rm = TRUE)
max_val <- max(expr_values, na.rm = TRUE)
# Generate individual plots
plots <- SpatialFeaturePlot(subset(lung2, Fibroblasts > 0), c("FB_0", "FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "FB_6"), combine = F, pt.size.factor = 2.2)

# Apply same color scale to each plot (e.g., white to red, limits 0 to 5)
plots <- lapply(plots, function(p) {
  p + scale_fill_viridis_c(option = "inferno", direction = 1)
})

# Combine them again
wrap_plots(plots, ncol = 4)


