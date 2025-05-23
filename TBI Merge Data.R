# Load required library 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(celldex)   # Reference datasets for annotation
library(scater)
library(SingleR)
library(pheatmap)
library(SingleCellExperiment)
library(SeuratDisk)
library(tidyr)
library(ggpubr)   # for stat_compare_means
library(svglite)  # to export SVG files
library(cowplot)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(grDevices)
library(magick)
library(ggforce)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)  # for mouse use org.Mm.eg.db ; use org.Hs.eg.db for human
library(KEGGREST)
library(STRINGdb)
library(igraph)
library(ggraph)
library(tidygraph)




# 🗂 Load individual Seurat objects
wt_sham <- Read10X(data.dir = "C:/Users/lindainesz/Desktop/WT_Sham/raw_feature_bc_matrix")
wt_d3 <- Read10X(data.dir = "C:/Users/lindainesz/Desktop/WT_D3/raw_feature_bc_matrix")
ecEXT1_sham <- Read10X(data.dir = "C:/Users/lindainesz/Desktop/Ve_Ext1_Sham/raw_feature_bc_matrix")
ecEXT1_d3 <- Read10X(data.dir = "C:/Users/lindainesz/Desktop/Ve_Ext1_D3/raw_feature_bc_matrix")

# Initialize the Seurat object with the raw (non-normalized data).
wt.sham <- CreateSeuratObject(counts = wt_sham, project = "WT Sham", min.cells = 3, min.features = 200)
wt.sham
wt.d3 <- CreateSeuratObject(counts = wt_d3, project = "WT D3", min.cells = 3, min.features = 200)
wt.d3
ecEXT1.sham <- CreateSeuratObject(counts = ecEXT1_sham, project = "ecEXT1 Sham", min.cells = 3, min.features = 200)
ecEXT1.sham
ecEXT1.d3 <- CreateSeuratObject(counts = ecEXT1_d3, project = "ecEXT1 D3", min.cells = 3, min.features = 200)
ecEXT1.d3

# 🔖 Add condition labels
wt.sham$condition <- "WT Sham"
wt.d3$condition <- "WT D3"
ecEXT1.sham$condition <- "ecEXT1 Sham"
ecEXT1.d3$condition <- "ecEXT1 D3"

# 🔄 Merge all into one object
tbi_merged_obj <- merge(wt.sham, y = list(wt.d3, ecEXT1.sham, ecEXT1.d3),
                        add.cell.ids = c("WT Sham", "WT D3", "ecEXT1 Sham", "ecEXT1 D3"),
                        project = "Merged_TBI")

# 🔬 Preprocessing: Normalize, find variable features, and scale
tbi_merged_obj <- NormalizeData(tbi_merged_obj)
tbi_merged_obj <- FindVariableFeatures(tbi_merged_obj)
tbi_merged_obj <- ScaleData(tbi_merged_obj)

# 🧬 Run PCA & UMAP
set.seed(1234) # this is so that if you run the same code, UMAP coordinates and clustering do not vary and the data is reproducible
tbi_merged_obj <- RunPCA(tbi_merged_obj, npcs = 30)
ElbowPlot(tbi_merged_obj, ndims = 50)
PC_tbiMerge_variance<- Stdev(tbi_merged_obj)
plot(cumsum(PC_tbiMerge_variance^2) / sum(PC_tbiMerge_variance^2), xlab="PCs", ylab="Cumulative Variance Explained", type="b")

# UMAP
tbi_merged_obj <- FindNeighbors(tbi_merged_obj, dims = 1:30)
tbi_merged_obj <- FindClusters(tbi_merged_obj, resolution = 0.20)
tbi_merged_obj <- RunUMAP(tbi_merged_obj, dims = 1:30)

# 🎨 Plot UMAP colored by condition
set.seed(1234)
DimPlot(tbi_merged_obj, reduction = "umap", group.by = "condition", cols = c("blue", "red", "green", "orange")) +
  ggtitle("WT & ecEXT1 Sham vs. D3") +
  theme_void()
DimPlot(tbi_merged_obj, reduction = "umap", label = TRUE, pt.size = 0.5)

set.seed(1234)


# 📌 5. Automated Cluster Annotation with SingleR

# Load reference dataset (e.g., Mouse RNA Cell Atlas)
refX <- celldex::MouseRNAseqData()

# Convert Seurat object to Single Cell Experiment

# Ensure you are using the RNA assay
DefaultAssay(tbi_merged_obj) <- "RNA"

# Convert to Single Cell Experiment with specific assay data
sceX <- as.SingleCellExperiment(tbi_merged_obj, assay = "RNA", layer = "data")
sceX <- as.SingleCellExperiment(tbi_merged_obj)

# Step 1: Extract log-normalized matrix with full cells
lognorm_matrix <- LayerData(tbi_merged_obj, assay = "RNA", layer = "data")

# Step 2: Match metadata to the matrix columns
common_cells <- intersect(colnames(lognorm_matrix), rownames(tbi_merged_obj@meta.data))
lognorm_matrix <- lognorm_matrix[, common_cells]
matched_metadata <- tbi_merged_obj@meta.data[common_cells, ]

# Step 3: Create Single Cell Experiment safely
sceX <- SingleCellExperiment(
  assays = list(logcounts = lognorm_matrix),
  colData = matched_metadata
)

# Run SingleR annotation
tbi_merge_cluster_annotation <- SingleR(
  test = sceX, 
  ref = refX, 
  labels = refX$label.main,
  clusters = sceX$seurat_clusters  # Make sure Seurat_clusters exist
)

tbi_merged_obj$SingleR_labels <- tbi_merge_cluster_annotation$labels[match(tbi_merged_obj$seurat_clusters, rownames(tbi_merge_cluster_annotation))]

DimPlot(tbi_merged_obj, reduction = "umap", 
        group.by = "SingleR_labels", label = TRUE, repel = TRUE) +
  ggtitle("rTBI: Sham vs. D3") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme_void()


# 📌 Identify Cluster Markers
# Find marker genes for each cluster
# Ensure that RNA assay layers are joined
tbi_merged_obj <- JoinLayers(tbi_merged_obj, assay = "RNA")

# Now, rerun FindAllMarkers
tbimerge_cluster_markers <- FindAllMarkers(tbi_merged_obj, 
                                           only.pos = TRUE, 
                                           min.pct = 0.25, 
                                           logfc.threshold = 0.25)

# View top marker genes per cluster
tbimerg_top_markers <- tbimerge_cluster_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
print(tbimerg_top_markers)
# Save TBI markers to file
write.csv(tbimerg_top_markers, "tbimerge_cluster_markers.csv")  
# to know where your file is save 
getwd()

# 🔍 Identify Marker Genes for Each Cluster
# This finds marker genes that define each cluster
tbi_markers <- FindAllMarkers(tbi_merged_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View top markers for each cluster
top10_tbimarkers <- tbi_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
View(top10_tbimarkers)
getwd() #To know where the data are saved


# Example marker genes for Mouse brain cells types
# marker genes for Neuron cells
marker_genes_Neuron <- c("Snap25", "Tubb3", "Rbfox3", "Syt1")
marker_genes_ExciNeuron <- c("Slc17a7", "Camk2a", "Grin1")
marker_genes_InhiNeuron <- c("Gad1", "Gad2", "Slc32a1", "Pvalb")
# marker genes for Astrocytes cells
marker_genes_Astrocytes <- c("Gfap", "Aldh1l1", "Aqp4", "Slc1a3", "Sox9")
# marker genes for Oligodendrocytes
marker_genes_Oligodendrocytes <- c("Mbp", "Plp1", "Cnp", "Sox10")
# marker genes for Oligodendrocyte Precursor Cells (OPCs)
marker_genes_OPC <- c("Pdgfra", "Cspg4", "Olig2")
# marker genes for Microglia
marker_genes_Microglia <- c("Cx3cr1", "Aif1", "Ptprc", "Trem2")
# marker genes for Endothelial Cells (EC)
marker_genes_EC <- c("Pecam1", "Cldn5", "Kdr", "Vwf")
# marker genes for Pericytes cells
marker_genes_Pericytes <- c("Pdgfrb", "Rgs5", "Des")
# marker genes for Ependymal cells
marker_genes_Ependymal <- c("Foxj1", "Tmem212", "Syt1")
# marker genes for Choroid Plexus Cells
marker_genes_Choroid <- c("Ttr", "Aqp1", "Clic6")
# marker genes for smooth muscle Cells (sm)
marker_genes_sm <- c("Acta2", "Myh11", "Tagln", "Cnn1", "Des", "Smtn","Notch3")
# marker genes for smooth muscle Cells (sm) second one
marker_genes_sm1 <- c("Acta2", "Myh11", "Tagln", "Cnn1", "Des")
# marker genes for macrophage cell
marker_genes_macrophage <- c("Cd68", "Mrc2", "Cd14", "Itgax", "Ccr2", "Cd36", "Pparg")
# marker genes for macrophage cell second one reduced 
marker_genes_macrophage1 <- c("Cd68", "Cd14", "Ccr2", "Cd36", "Pparg")
# marker genes for epithelium cells
marker_genes_epitheluim <- c("Ttr", "Foxj1", "Epcam","Krt18")


# Feature plots for selected genes
FeaturePlot(tbi_merged_obj, features = marker_genes_Neuron)
FeaturePlot(tbi_merged_obj, features = marker_genes_ExciNeuron)
FeaturePlot(tbi_merged_obj, features = marker_genes_InhiNeuron)
FeaturePlot(tbi_merged_obj, features = marker_genes_Astrocytes)
FeaturePlot(tbi_merged_obj, features = marker_genes_Microglia)
FeaturePlot(tbi_merged_obj, features = marker_genes_EC)
FeaturePlot(tbi_merged_obj, features = marker_genes_Oligodendrocytes)
FeaturePlot(tbi_merged_obj, features = marker_genes_OPC)
FeaturePlot(tbi_merged_obj, features = marker_genes_Pericytes)
FeaturePlot(tbi_merged_obj, features = marker_genes_sm)
FeaturePlot(tbi_merged_obj, features = marker_genes_Ependymal)
FeaturePlot(tbi_merged_obj, features = marker_genes_Choroid)
FeaturePlot(tbi_merged_obj, features = marker_genes_macrophage1)
FeaturePlot(tbi_merged_obj, features = marker_genes_epitheluim)

# UMAP plot with cluster labels                                                                |
DimPlot(tbi_merged_obj, reduction = "umap", label = TRUE, pt.size = 0.5)

#🏷️ Assign Cell Type Labels
# Example manual annotation
tbi_cluster_ids <- Idents(tbi_merged_obj)
tbi_cluster_annotations <- c("0" = "Astrocytes",
                             "1" = "Microglia",
                             "2" = "Endotheluim Cells",
                             "3" = "Endotheluim Cells",
                             "4" = "Microglia",
                             "5" = "Macrophages",
                             "6" = "Ependymal Cells",
                             "7" = "Microglia",
                             "8" = "Ependymal Cells",
                             "9" = "Mural Cells",
                             "10" = "Oligodendrocytes",
                             "11" = "Macrophages",
                             "13" = "Neuron",
                             "14" = "Fibroblast",
                             "15" = "Neuron")  # etc.

# Make sure the identities are set correctly
Idents(tbi_merged_obj) <- "seurat_clusters"

# Assign cell types using the full annotation list
cell_type_labels <- tbi_cluster_annotations[as.character(Idents(tbi_merged_obj))]
names(cell_type_labels) <- Cells(tbi_merged_obj)

# Add as metadata
tbi_merged_obj <- AddMetaData(tbi_merged_obj, metadata = cell_type_labels, col.name = "cell_type")
DimPlot(tbi_merged_obj, group.by = "cell_type", label = TRUE)


sum(is.na(Embeddings(tbi_merged_obj, "umap")))
tbi_merged_obj$condition <- factor(tbi_merged_obj$condition,
                                   levels = c("WT Sham", "WT D3", "ecEXT1 Sham", "ecEXT1 D3"))

DimPlot(tbi_merged_obj, reduction = "umap", group.by = "cell_type", split.by = "condition",
        label = TRUE, repel = TRUE) +
  ggtitle("TBI WT & ecEXT1 Sham vs. D3") +
  theme_void()
# Save as high-res PNG
umap_plot <- DimPlot(tbi_merged_obj, reduction = "umap", group.by = "cell_type", split.by = "condition",
                     label = TRUE, repel = TRUE) +
  ggtitle("TBI WT & ecEXT1 Sham vs. D3") +
  theme_void()
ggsave("TBI_UMAP_split.png", plot = umap_plot, width = 12, height = 6, dpi = 300)
getwd()

saveRDS(tbi_merged_obj, "tbi_merged_obj_note.rds")
tbi_merged_obj <- readRDS("tbi_merged_obj_note.rds")


# Check if NA in condition is gone you do not need to run print only keep it as a note to myself 
print(levels(tbi_merged_obj$condition))           # should NOT include "NA"
print(table(tbi_merged_obj$condition, useNA = "always"))  # should NOT include NA
# Drop NA in 'condition' if it's causing an empty facet
# Check for NA in cell_type
table(is.na(tbi_merged_obj$cell_type))  # TRUE = bad

# Subset the object to only keep complete cells
tbi_merged_obj_clean <- subset(tbi_merged_obj, 
                               subset = !is.na(cell_type) & !is.na(condition))

DimPlot(tbi_merged_obj_clean, 
        reduction = "umap", 
        group.by = "cell_type", 
        split.by = "condition", 
        label = TRUE, repel = TRUE) +
  ggtitle("TBI WT & ecEXT1 Sham vs. D3") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 20)),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )


# Save as high-res PNG
umap_plot1 <- DimPlot(tbi_merged_obj_clean, 
                      reduction = "umap", 
                      group.by = "cell_type", 
                      split.by = "condition", 
                      label = TRUE, repel = TRUE) +
  ggtitle("TBI WT & ecEXT1 Sham vs. D3") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 20)),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )
ggsave("TBI_UMAP_splitG.png", plot = umap_plot1, width = 12, height = 6, dpi = 300)
umap_split <- DimPlot(tbi_merged_obj_clean, 
                      reduction = "umap", 
                      group.by = "cell_type", 
                      split.by = "condition", 
                      label = TRUE, repel = TRUE) +
  ggtitle("TBI WT & ecEXT1 Sham vs. D3") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 20)),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

#set a 2x2 layout 
umap_layout <- DimPlot(tbi_merged_obj_clean, 
        reduction = "umap", 
        group.by = "cell_type", 
        split.by = "condition", 
        label = TRUE, 
        repel = TRUE, 
        ncol = 2) +  # 🔥 This is the key for 2 columns
  ggtitle("TBI WT & ecEXT1 Sham vs. D3") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 20)),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )
ggsave("TBI_UMAP_splitx2.png", plot = umap_layout, width = 12, height = 6, dpi = 300)


# 📊 Cluster proportion analysis per condition
tbicluster_counts <- table(tbi_merged_obj_clean$cell_type, tbi_merged_obj_clean$condition)
print(tbicluster_counts)
      
# 🧪 Chi-square test
tbichi_test <- chisq.test(tbicluster_counts)
print("Chi-Square Test:")
print(tbichi_test)

#📊 Cell proportion analysis per condition

#-------Create a table of cell counts per condition and cell type
tbi_cell_table <- table(tbi_merged_obj_clean$condition, tbi_merged_obj_clean$cell_type)

#--------Convert to data frame for ggplot
tbi_cell_df <- as.data.frame(tbi_cell_table)
colnames(tbi_cell_df) <- c("Condition", "CellType", "Count")

#-------Add proportion column
tbicell_df <- tbi_cell_df %>%
  group_by(Condition) %>%
  mutate(Percent = Count / sum(Count) * 100)

#with theme minimal used
ggplot(tbicell_df, aes(x = Condition, y = Percent, fill = CellType)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Percentage of Cells") +
  xlab("Condition") +
  ggtitle("Cell Type Composition by Condition") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

#with theme void used
ggplot(tbicell_df, aes(x = Condition, y = Percent, fill = CellType)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  ggtitle("Cell Type Composition per Condition") +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 20)),
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
    axis.text.y = element_text(),  # Show y-axis ticks
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(angle = 90, margin = margin(r = 10))  # <-- This rotates y-label
  ) +
  ylab("Percentage of Cells") +
  xlab("Condition")

# Add percentage labels on bars
ggplot(tbicell_df, aes(x = Condition, y = Percent, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +

  geom_text(aes(label = sprintf("%.1f%%", Percent)),
            position = position_dodge(width = 0.9),
            vjust = -0.25, # move text above bars
            size = 3.2) +
  
  labs(y = "Percentage of Cells", x = "Condition") +
  theme_void(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title.y = element_text(angle = 90),
    legend.position = "right"
  )

# I did Option A: Side-by-side barplot need to do now Option B: Stacked barplot (if you want to show proportions within 100%)


saveRDS(tbi_merged_obj, "tbi_merged_obj_note.rds")
saveRDS(tbi_merged_obj_clean, "tbi_merged_obj_note1.rds")
saveRDS(tbicell_df, "tbi_merged_obj_note2.rds")
tbi_merged_obj <- readRDS("tbi_merged_obj_note.rds")
tbi_merged_obj_clean <- readRDS("tbi_merged_obj_note1.rds")
tbicell_df <- readRDS("tbi_merged_obj_note2.rds")

#Option B: Stacked barplot (if you want to show proportions within 100%)
ggplot(tbicell_df, aes(x = Condition, y = Percent, fill = CellType)) +
  geom_bar(stat = "identity") +
  ylab("Percentage of Cells") +
  xlab("Condition") +
  ggtitle("Cell Type Composition per Condition") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

#Option B: Stacked barplot without gridlines added percentages in each well
ggplot(tbicell_df, aes(x = Condition, y = Percent, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Percent)),
            position = position_stack(vjust = 0.5),
            size = 2.8, color = "black") +
  labs(
    y = "Percentage of Cells",
    x = "Condition",
    title = "Cell Type Composition per Condition"
  ) +
  theme_void(base_size = 14) +
  theme(
    axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5), # <- THIS LINE
    axis.title.x = element_text(),
    axis.text.x = element_text(angle = 0, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

saveRDS(tbi_merged_obj, "tbi_merged_obj_note.rds")
saveRDS(tbi_merged_obj_clean, "tbi_merged_obj_note1.rds")
saveRDS(tbicell_df, "tbi_merged_obj_note2.rds")

tbi_merged_obj <- readRDS("tbi_merged_obj_note.rds")
tbi_merged_obj_clean <- readRDS("tbi_merged_obj_note1.rds")
tbicell_df <- readRDS("tbi_merged_obj_note2.rds")


#------------------------The code below was not used this is a model code tested DO NOT ADD IT LATER------------------------------- 
#------------------------------------------- End of the try and fail segment --------------------------------------------
#START
#violinplot generate
#✅ Step 1: List of Heparan Sulfate Genes
hs_genes <- c(
  "Ext1", "Ext2", "Extl1", "Extl2", "Extl3",   # chain elongation
  "Ndst1", "Ndst2", "Ndst3", "Ndst4",          # N-deacetylase/N-sulfotransferase
  "Hs2st1", "Hs3st1", "Hs3st2", "Hs3st3a1", "Hs3st3b1", "Hs6st1", "Hs6st2", "Hs6st3",  # sulfotransferases
  "Sulf1", "Sulf2",                            # sulfatases
  "Gpc1", "Gpc2", "Gpc3", "Gpc4", "Gpc5", "Gpc6",  # glypicans (core proteins)
  "Sdc1", "Sdc2", "Sdc3", "Sdc4"              # syndecans
)
#✅ 2. Extract Data
head(tbi_merged_obj_clean@meta.data)
# Fetch expression data + metadata
tbiexpr_data <- FetchData(tbi_merged_obj_clean, vars = c("condition", "cell_type", hs_genes))

# Reshape to long format for ggplot
tbiexpr_long <- tbiexpr_data %>%
  pivot_longer(cols = all_of(hs_genes), names_to = "Gene", values_to = "Expression") %>%
  rename(Condition = condition, CellType = cell_type)
#✅ 3 Generate Faceted Violin Plot with Stats
HSvln_plot <- ggplot(tbiexpr_long, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width") +
  facet_wrap(~Gene, ncol = 3, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("WT Sham", "WT D3"), c("ecEXT1 Sham", "ecEXT1 D3"))) +
  scale_fill_manual(values = c("WT Sham" = "#E64B35", "WT D3" = "#4DBBD5",
                               "ecEXT1 Sham" = "#00A087", "ecEXT1 D3" = "#9370DB")) +
  theme_minimal(base_size = 13) +
  labs(title = "Heparan Sulfate Gene Expression Across Conditions",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.position = "right")
print(HSvln_plot)



ggplot(tbiexpr_long, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = FALSE, width = 0.8) +
  facet_wrap(~Gene, ncol = 4) +
  stat_compare_means(
    method = "kruskal.test", 
    label = "p.signif", 
    label.y.npc = "top", 
    size = 3
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  labs(title = "Heparan Sulfate Gene Expression Across Conditions") +
  scale_fill_brewer(palette = "Set2")
#END
#------------------------The code above was not used this is a model code tested DO NOT ADD IT LATER------------------------------- 

#violinplot generate each HS genes at the time 
#✅ Step 1: List of Heparan Sulfate Genes
hs_genes_elongation <- c(
  "Ext1", "Ext2", "Extl1", "Extl2", "Extl3")  # chain elongation
hs_genes_deacetylase <- c(
  "Ndst1", "Ndst2", "Ndst3", "Ndst4")  # N-deacetylase/N-sulfotransferase
hs_genes_sulfotransferases <- c(
  "Hs2st1", "Hs3st1", "Hs3st2", "Hs3st3a1", "Hs3st3b1", "Hs6st1", "Hs6st2", "Hs6st3") # sulfotransferases
hs_genes_sulfatases <- c(
  "Sulf1", "Sulf2")   # sulfatases
hs_genes_glypicans <- c(
  "Gpc1", "Gpc2", "Gpc3", "Gpc4", "Gpc5", "Gpc6")  # glypicans (core proteins)
hs_genes_syndecans <- c(
  "Sdc1", "Sdc2", "Sdc3", "Sdc4")  # syndecans

#✅ 2. Extract Data
head(tbi_merged_obj_clean@meta.data)
# Fetch expression data + metadata
tbiexpr_data_elongation <- FetchData(tbi_merged_obj_clean, vars = c("condition", "cell_type", hs_genes_elongation))
tbiexpr_data_deacetylase <- FetchData(tbi_merged_obj_clean, vars = c("condition", "cell_type", hs_genes_deacetylase))
tbiexpr_data_sulfotransferases <- FetchData(tbi_merged_obj_clean, vars = c("condition", "cell_type", hs_genes_sulfotransferases))
tbiexpr_data_sulfatases <- FetchData(tbi_merged_obj_clean, vars = c("condition", "cell_type", hs_genes_sulfatases))
tbiexpr_data_glypicans <- FetchData(tbi_merged_obj_clean, vars = c("condition", "cell_type", hs_genes_glypicans))
tbiexpr_data_syndecans <- FetchData(tbi_merged_obj_clean, vars = c("condition", "cell_type", hs_genes_syndecans))

# Reshape to long format for ggplot
tbiexpr_hs_elongation <- tbiexpr_data_elongation %>%
  pivot_longer(cols = all_of(hs_genes_elongation), names_to = "Gene", values_to = "Expression") %>%
  rename(Condition = condition, CellType = cell_type)

tbiexpr_hs_deacetylase <- tbiexpr_data_deacetylase %>%
  pivot_longer(cols = all_of(hs_genes_deacetylase), names_to = "Gene", values_to = "Expression") %>%
  rename(Condition = condition, CellType = cell_type)

tbiexpr_hs_sulfotransferase <- tbiexpr_data_sulfotransferases %>%
  pivot_longer(cols = all_of(hs_genes_sulfotransferases), names_to = "Gene", values_to = "Expression") %>%
  rename(Condition = condition, CellType = cell_type)

tbiexpr_hs_sulfatases <- tbiexpr_data_sulfatases %>%
  pivot_longer(cols = all_of(hs_genes_sulfatases), names_to = "Gene", values_to = "Expression") %>%
  rename(Condition = condition, CellType = cell_type)

tbiexpr_hs_glypicans <- tbiexpr_data_glypicans %>%
  pivot_longer(cols = all_of(hs_genes_glypicans), names_to = "Gene", values_to = "Expression") %>%
  rename(Condition = condition, CellType = cell_type)

tbiexpr_hs_syndecans <- tbiexpr_data_syndecans %>%
  pivot_longer(cols = all_of(hs_genes_syndecans), names_to = "Gene", values_to = "Expression") %>%
  rename(Condition = condition, CellType = cell_type)


#HS Elongation
ggplot(tbiexpr_hs_elongation, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = FALSE, width = 0.8) +
  facet_wrap(~Gene, ncol = 4) +
  stat_compare_means(
    method = "kruskal.test", 
    label = "p.signif", 
    label.y.npc = "top", 
    size = 3
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  labs(title = "Heparan Sulfate Elongation Gene Expression Across Conditions") 

#HS Deacetylase
ggplot(tbiexpr_hs_deacetylase, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = FALSE, width = 0.8) +
  facet_wrap(~Gene, ncol = 4) +
  stat_compare_means(
    method = "kruskal.test", 
    label = "p.signif", 
    label.y.npc = "top", 
    size = 3
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  labs(title = "Heparan Sulfate N-deacetylase/N-sulfotransferase Gene Expression Across Conditions")

#HS Sulfotransferase
ggplot(tbiexpr_hs_sulfotransferase, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = FALSE, width = 0.8) +
  facet_wrap(~Gene, ncol = 4) +
  stat_compare_means(
    method = "kruskal.test", 
    label = "p.signif", 
    label.y.npc = "top", 
    size = 3
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  labs(title = "Heparan Sulfate Sulfotransferase Gene Expression Across Conditions")

#HS sulfatases
ggplot(tbiexpr_hs_sulfatases, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = FALSE, width = 0.8) +
  facet_wrap(~Gene, ncol = 4) +
  stat_compare_means(
    method = "kruskal.test", 
    label = "p.signif", 
    label.y.npc = "top", 
    size = 3
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  labs(title = "Heparan Sulfate Sulfatases Gene Expression Across Conditions")

#HS Glypicans
ggplot(tbiexpr_hs_glypicans, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = FALSE, width = 0.8) +
  facet_wrap(~Gene, ncol = 4) +
  stat_compare_means(
    method = "kruskal.test", 
    label = "p.signif", 
    label.y.npc = "top", 
    size = 3
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  labs(title = "Heparan Sulfate Glypicans Gene Expression Across Conditions")

#HS syndecans
ggplot(tbiexpr_hs_syndecans, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = FALSE, width = 2) +
  facet_wrap(~Gene, ncol = 4) +
  stat_compare_means(
    method = "kruskal.test", 
    label = "p.signif", 
    label.y.npc = "top", 
    size = 3
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  labs(title = "Heparan Sulfate Syndecans Gene Expression Across Conditions")
ggsave("HS_Syndecans_expression.png", dpi = 600, width = 12, height = 8)


#----------------------------- with code below I can get gene expression per condition --------------------------------
# Define your pairwise comparisons
my_comparisons_group <- list(
  c("WT Sham", "WT D3"),
  c("WT Sham", "ecEXT1 Sham"),
  c("WT Sham", "ecEXT1 D3"),
  c("WT D3", "ecEXT1 Sham"),
  c("WT D3", "ecEXT1 D3"),
  c("ecEXT1 Sham", "ecEXT1 D3")
)

hs_genes <- c(
  "Ext1", "Ext2", "Extl1", "Extl2", "Extl3",   # chain elongation
  "Ndst1", "Ndst2", "Ndst3", "Ndst4",          # N-deacetylase/N-sulfotransferase
  "Hs2st1", "Hs3st1", "Hs3st2", "Hs3st3a1", "Hs3st3b1", "Hs6st1", "Hs6st2", "Hs6st3",  # sulfotransferases
  "Sulf1", "Sulf2",                            # sulfatases
  "Gpc1", "Gpc2", "Gpc3", "Gpc4", "Gpc5", "Gpc6",  # glypicans (core proteins)
  "Sdc1", "Sdc2", "Sdc3", "Sdc4"              # syndecans
)

# Filter for HS Genes
ext1_data <- tbiexpr_long %>% filter(Gene == "Ext1")
ext2_data <- tbiexpr_long %>% filter(Gene == "Ext2")
extl1_data <- tbiexpr_long %>% filter(Gene == "Extl1")
extl2_data <- tbiexpr_long %>% filter(Gene == "Extl2")
extl3_data <- tbiexpr_long %>% filter(Gene == "Extl3")
ndst1_data <- tbiexpr_long %>% filter(Gene == "Ndst1")
ndst2_data <- tbiexpr_long %>% filter(Gene == "Ndst2")
ndst3_data <- tbiexpr_long %>% filter(Gene == "Ndst3")
ndst4_data <- tbiexpr_long %>% filter(Gene == "Ndst4")
hs2st1_data <- tbiexpr_long %>% filter(Gene == "Hs2st1")
hs3st1_data <- tbiexpr_long %>% filter(Gene == "Hs3st1")
hs3st2_data <- tbiexpr_long %>% filter(Gene == "Hs3st2")
hs3st3a1_data <- tbiexpr_long %>% filter(Gene == "Hs3st3a1")
hs3st3b1_data <- tbiexpr_long %>% filter(Gene == "Hs3st3b1")
hs6st1_data <- tbiexpr_long %>% filter(Gene == "Hs6st1")
hs6st2_data <- tbiexpr_long %>% filter(Gene == "Hs6st2")
hs6st3_data <- tbiexpr_long %>% filter(Gene == "Hs6st3")
sulf1_data <- tbiexpr_long %>% filter(Gene == "Sulf1")
sulf2_data <- tbiexpr_long %>% filter(Gene == "Sulf2")
gpc1_data <- tbiexpr_long %>% filter(Gene == "Gpc1")
gpc2_data <- tbiexpr_long %>% filter(Gene == "Gpc2")
gpc3_data <- tbiexpr_long %>% filter(Gene == "Gpc3")
gpc4_data <- tbiexpr_long %>% filter(Gene == "Gpc4")
gpc5_data <- tbiexpr_long %>% filter(Gene == "Gpc5")
gpc6_data <- tbiexpr_long %>% filter(Gene == "Gpc6")
sdc1_data <- tbiexpr_long %>% filter(Gene == "Sdc1")
sdc2_data <- tbiexpr_long %>% filter(Gene == "Sdc2")
sdc3_data <- tbiexpr_long %>% filter(Gene == "Sdc3")
sdc4_data <- tbiexpr_long %>% filter(Gene == "Sdc4")


# Plot with pairwise Wilcoxon test
# this is the original code plot I use
ggplot(ext1_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif") +
  ggtitle("Ext1 Expression Across Conditions") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )

#------------------- Everything below was me trying to find the best one and learning -----------------------------------------------------------------------------

#1st trial test
ggplot(ext1_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif") +  # use "p.format" for exact p-value labels
  ggtitle("Ext1 Expression Across Conditions") +
  ylab("Expression") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )

ggplot(ext1_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif") +  # use "p.format" for exact p-value labels
  ggtitle("Ext1 Expression Across Conditions") +
  ylab("Expression") +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, 3)) +  # or maybe up to 2. This Zoom in on the y-axis
  geom_jitter(width = 0.1, size = 0.5, alpha = 0.4) + # Add data points to show distribution
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )

# zoomed code with x and Y showing
ggplot(ext1_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.2, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  stat_compare_means(
    comparisons = my_comparisons_group,
    method = "wilcox.test",
    label = "p.signif",
    label.y = c(1.5, 1.7, 1.9, 2.1, 2.3, 2.5)  # Shift these based on new y-limits
  ) +
  scale_fill_brewer(palette = "Set2") +
  coord_cartesian(ylim = c(0.0000001, 3)) +  # <<< ZOOM HERE
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  ylab("Ext1 Expression") +
  xlab("Condition") +
  ggtitle("Ext1 Expression Across Conditions (Zoomed)")

#zoomed code with x and Y not showing
ggplot(ext1_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.2, color = NA) +
  geom_boxplot(
    width = 0.1,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  stat_compare_means(
    comparisons = my_comparisons_group,
    method = "wilcox.test",
    label = "p.signif",
    label.y = c(2.3, 2.6, 2.9, 3.2, 3.5, 3.8)
  ) +
  scale_fill_brewer(palette = "Set2") +
  coord_cartesian(ylim = c(0, 4)) +
  theme_classic(base_size = 14) +
  theme(
    axis.line.x = element_blank(),     # remove x axis line
    axis.line.y = element_blank(),     # remove y axis line
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  ylab("Expression") +
  xlab("Condition") +
  ggtitle("Ext1 Expression Across Conditions")


#log10 values
summary(ext1_data$Expression) # this give you information 


ggplot(ext1_data[ext1_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.2, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.5) +
  stat_compare_means(
    comparisons = my_comparisons_group,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE,
    label.y = c(0.5, 0.7, 0.9, 1.1, 1.3, 1.5)  # These must match the log scale
  ) +
  scale_y_log10(
    breaks = c(0.01,0.03, 0.1, 0.3, 1, 5),
    labels = scales::label_number()
  ) +
  scale_fill_brewer(palette = "Set2") +
  coord_cartesian(ylim = c(0.01, 4)) +
  theme_classic(base_size = 14) +
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  ylab("Expression (log-scaled)") +
  xlab("Condition") +
  ggtitle("Ext1 Expression Across Conditions (log10 scale)")


# You can generate a side-by-side comparison of: 
# One violin plot with raw expression (linear scale)
# One violin plot with log10-scaled expression (filtered to >0)
# This will let you visually compare how the data behaves under each scale while keeping the statistical test consistent (based on raw data).
# === Linear Scale Plot ===
ext1linear_plot <- ggplot(ext1_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif", label.y = c(3.2, 3.4, 3.6, 3.8, 4.0, 4.2)) +  # use "p.format" for exact p-value labels
  ggtitle("Ext1 Expression Across Conditions") +
  ylab("Expression") +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, 4.5)) +  # or maybe up to 2. This Zoom in on the y-axis
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )

# === Log10 Scale Plot (non-zero only) ===
ext1log_plot <- ext1_data %>% filter(Expression > 0) %>%
  ggplot(aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.2, color = NA) +
  geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.3) +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif", hide.ns = TRUE, label.y = 3) +
  scale_y_log10(breaks = c(0.1, 0.3, 1, 3), labels = scales::label_number()) +
  coord_cartesian(ylim = c(0.1, 4)) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = paste("Ext1", "- Log10 Scale"), y = "Expression (log10)") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
  )

# === Combine Side-by-Side ===
ext1combined_plot <- linear_plot + log_plot + plot_layout(ncol = 2)

# === Save as SVG or PNG ===
ggsave("Ext1_linear_vs_log10.png", dpi = 600, width = 12, height = 8)
ggsave("Extl1_linear_vs_log10.svg", plot = combined_plot, width = 12, height = 6)


#------------------- Everything above was me trying to find the best one and learning -----------------------------------------------------------------------------

#this code is really good
ggplot(ext1_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif", label.y = c(3.2, 3.4, 3.6, 3.8, 4.0, 4.2)) +  # use "p.format" for exact p-value labels
  ggtitle("Ext1 Expression Across Conditions") +
  ylab("Expression") +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, 4.5)) +  # or maybe up to 2. This Zoom in on the y-axis
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )
ggsave("HS_ext1_expression.png", dpi = 600, width = 12, height = 8)


#log10 data 
ggplot(ext1_data[ext1_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.5) +
  stat_compare_means(                           #if you do not want show the significance remove this section
    comparisons = my_comparisons_group,
    method = "wilcox.test",
    label = "p.signif",
    label.y = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)  # adjust to fit within log scale
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0.01, 70)) +  # zoom in a little
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Ext1 Expression Across Conditions",
    y = "Expression (log10)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
# Note: if this statistical values do not match between the linear and log plot you can
#📌 Option 1: Report results from the original scale
#If your data is right-skewed, log-transform for visualization only.
#Perform statistical tests on the raw (linear) data.
#📌 Option 2: Report both
#Some journals accept both:
#Stats: from raw data.
#Plot: shown in log scale for interpretability.

# EXT2 Gene

#------------------------ Code I started with------------------------------------------------------------------------------------------------------------------------
ggplot(ext2_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif") +
  ggtitle("Ext2 Expression Across Conditions") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )
ggsave("HS_ext2_expression.png", dpi = 600, width = 12, height = 8)
summary(ext2_data$Expression) # this give you information

ggplot(ext2_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif", label.y = c(1.6, 1.8, 2.0, 2.2, 2.4, 2.6)) +  # use "p.format" for exact p-value labels
  ggtitle("Ext2 Expression Across Conditions") +
  ylab("Expression") +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, 3)) +  # or maybe up to 2. This Zoom in on the y-axis
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )
#------------------------ Code I started with------------------------------------------------------------------------------------------------------------------------

#---------- Code below is what I ended up using -----------------------------
ggplot(ext2_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif", label.y = c(3.2, 3.4, 3.6, 3.8, 4.0, 4.2)) +  # use "p.format" for exact p-value labels
  ggtitle("Ext2 Expression Across Conditions") +
  ylab("Expression") +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, 4.5)) +  # or maybe up to 2. This Zoom in on the y-axis
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )

#log10 data with statistic value
ggplot(ext2_data[ext2_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.5) +
  stat_compare_means(                           #if you do not want show the significance remove this section
    comparisons = my_comparisons_group,
    method = "wilcox.test",
    label = "p.signif",
    label.y = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)  # adjust to fit within log scale
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0.01, 70)) +  # zoom in a little
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Ext2 Expression Across Conditions",
    y = "Expression (log10)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#log10 data without statistic value
ggplot(ext2_data[ext2_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.5) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0.01, 70)) +  # zoom in a little
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Ext2 Expression Across Conditions",
    y = "Expression (log10)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


#EXTL1 Genes
#Code below is the one and the one I ended up using
ggplot(extl1_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif") + # use "p.format" for exact p-value labels
  ggtitle("Extl1 Expression Across Conditions") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )

#this code is just ok not the best
ggplot(extl1_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.05, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif") +  # use "p.format" for exact p-value labels
  ggtitle("Extl1 Expression Across Conditions") +
  ylab("Expression") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )
ggsave("HS_extl1_expression.png", dpi = 600, width = 12, height = 8)


#log10 data with statistic value
ggplot(extl1_data[extl1_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.5) +
  stat_compare_means(                           #if you do not want show the significance remove this section
    comparisons = my_comparisons_group,
    method = "wilcox.test",
    label = "p.signif",
    label.y = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)  # adjust to fit within log scale
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0.01, 70)) +  # zoom in a little
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Extl1 Expression Across Conditions",
    y = "Expression (log10)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#log10 data without statistic value
ggplot(extl1_data[extl1_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.5) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0.01, 70)) +  # zoom in a little
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Extl1 Expression Across Conditions",
    y = "Expression (log10)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


#Extl2 gene
#original code
ggplot(extl2_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif") +
  ggtitle("Extl2 Expression Across Conditions") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )
ggsave("HS_extl2_expression.png", dpi = 600, width = 12, height = 8)


#the linear code I end up using  
ggplot(extl2_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif") +  # use "p.format" for exact p-value labels
  ggtitle("Extl2 Expression Across Conditions") +
  ylab("Expression") +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, 8)) +  # to 4 statistic bar do not appear, but to 8 and you get stat bar. This Zoom in on the y-axis
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )

#log10 data with statistic value
ggplot(extl2_data[extl2_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.5) +
  stat_compare_means(                           #if you do not want show the significance remove this section
    comparisons = my_comparisons_group,
    method = "wilcox.test",
    label = "p.signif",
    label.y = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)  # adjust to fit within log scale
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0.01, 70)) +  # zoom in a little
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Extl2 Expression Across Conditions",
    y = "Expression (log10)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#log10 data without statistic value
ggplot(extl2_data[extl2_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.5) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0.01, 70)) +  # zoom in a little
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Extl2 Expression Across Conditions",
    y = "Expression (log10)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


#EXTL3 Gene

ggplot(extl3_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif") +  # use "p.format" for exact p-value labels
  ggtitle("Extl3 Expression Across Conditions") +
  ylab("Expression") +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, 8)) +  # to 4 statistic bar do not appear, but to 8 and you get stat bar. This Zoom in on the y-axis
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )
ggsave("HS_extl3_expression.png", dpi = 600, width = 12, height = 8)


#log10 data with statistic value
ggplot(extl3_data[extl3_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.5) +
  stat_compare_means(                           #if you do not want show the significance remove this section
    comparisons = my_comparisons_group,
    method = "wilcox.test",
    label = "p.signif",
    label.y = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)  # adjust to fit within log scale
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0.01, 70)) +  # zoom in a little
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Extl3 Expression Across Conditions",
    y = "Expression (log10)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#log10 data without statistic value
ggplot(extl3_data[extl3_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.5) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0.01, 70)) +  # zoom in a little
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Extl3 Expression Across Conditions",
    y = "Expression (log10)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


#NDST1 Gene
ggplot(ndst1_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif") +  # use "p.format" for exact p-value labels
  ggtitle("Ndst1 Expression Across Conditions") +
  ylab("Expression") +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, 8)) +  # to 4 statistic bar do not appear, but to 8 and you get stat bar. This Zoom in on the y-axis
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )

ggsave("HS_Ndst1_expression.png", dpi = 600, width = 12, height = 8)

#log10 data with statistic value
ggplot(ndst1_data[ndst1_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.5) +
  stat_compare_means(                           #if you do not want show the significance remove this section
    comparisons = my_comparisons_group,
    method = "wilcox.test",
    label = "p.signif",
    label.y = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)  # adjust to fit within log scale
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0.01, 70)) +  # zoom in a little
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Ndst1 Expression Across Conditions",
    y = "Expression (log10)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#log10 data without statistic value
ggplot(ndst1_data[ndst1_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.5) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0.01, 7)) +  # zoom in a little
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Ndst1 Expression Across Conditions",
    y = "Expression (log10)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#NDST2
ggplot(ndst2_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif") +  # use "p.format" for exact p-value labels
  ggtitle("Ndst2 Expression Across Conditions") +
  ylab("Expression") +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, 6)) +  # 6 for statistic to appear and 2 to zoom in This Zoom in on the y-axis
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )
ggsave("HS_Ndst2_expression.png", dpi = 600, width = 12, height = 8)

#log10 data with statistic value
ggplot(ndst2_data[ndst2_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.5) +
  stat_compare_means(                           #if you do not want show the significance remove this section
    comparisons = my_comparisons_group,
    method = "wilcox.test",
    label = "p.signif",
    label.y = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)  # adjust to fit within log scale
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0.01, 50)) +  # zoom in a little
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Ndst2 Expression Across Conditions",
    y = "Expression (log10)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#log10 data without statistic value
ggplot(ndst2_data[ndst2_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.5) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0.01, 10)) +  # zoom in a little
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Ndst2 Expression Across Conditions",
    y = "Expression (log10)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


#NDST3
ggplot(ndst3_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif") +  # use "p.format" for exact p-value labels
  ggtitle("Ndst3 Expression Across Conditions") +
  ylab("Expression") +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, 0.75)) +  # 8 for statistic to appear and 0.75 to zoom in This Zoom in on the y-axis
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )
ggsave("HS_Ndst3_expression.png", dpi = 600, width = 12, height = 8)

#log10 data with statistic value
ggplot(ndst3_data[ndst3_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.5) +
  stat_compare_means(                           #if you do not want show the significance remove this section
    comparisons = my_comparisons_group,
    method = "wilcox.test",
    label = "p.signif",
    label.y = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)  # adjust to fit within log scale
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0.01, 50)) +  # zoom in a little
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Ndst3 Expression Across Conditions",
    y = "Expression (log10)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#log10 data without statistic value
ggplot(ndst3_data[ndst3_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.5) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0.1, 10)) +  # zoom in a little
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Ndst3 Expression Across Conditions",
    y = "Expression (log10)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#------------------------------------------------------------------------------------------
saveRDS(tbi_merged_obj, "tbi_merged_obj_note.rds")
saveRDS(tbi_merged_obj_clean, "tbi_merged_obj_note1.rds")
saveRDS(tbicell_df, "tbi_merged_obj_note2.rds")

tbi_merged_obj <- readRDS("tbi_merged_obj_note.rds")
tbi_merged_obj_clean <- readRDS("tbi_merged_obj_note1.rds")
tbicell_df <- readRDS("tbi_merged_obj_note2.rds")

#Ndst4
ggplot(ndst4_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif") +  # use "p.format" for exact p-value labels
  ggtitle("Ndst4 Expression Across Conditions") +
  ylab("Expression") +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, 0.5)) +  # 8 for statistic to appear and 0.5 to zoom in This Zoom in on the y-axis
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )
ggsave("HS_Ndst4_expression.png", dpi = 600, width = 12, height = 8)

#log10 data with statistic value
ggplot(ndst4_data[ndst4_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.5) +
  stat_compare_means(                           #if you do not want show the significance remove this section
    comparisons = my_comparisons_group,
    method = "wilcox.test",
    label = "p.signif",
    label.y = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)  # adjust to fit within log scale
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0.01, 70)) +  # zoom in a little
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Ndst4 Expression Across Conditions",
    y = "Expression (log10)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

#log10 data without statistic value
ggplot(ndst4_data[ndst4_data$Expression > 0, ], aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.5) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0.01, 70)) +  # zoom in a little
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Ndst4 Expression Across Conditions",
    y = "Expression (log10)",
    x = "Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
#------------------------------------------- End of the try and fail segment --------------------------------------------

# 05-08-25 after Dr Zhen's email I stop here to run dotplot, volcano plot and a DEG

#Crowded violent plot
# Plot
# ----- pull expression + metadata -----
hsexpr_df <- FetchData(
  tbi_merged_obj_clean,
  vars = c("condition", "cell_type", hs_genes)   # <- use YOUR column names here
)

# ----- convert to long tidy format -----
tbihsexpr_long <- hsexpr_df %>% 
  pivot_longer(cols = all_of(hs_genes),
               names_to   = "Gene",
               values_to  = "Expression") %>% 
  rename(Condition = condition,
         CellType  = cell_type)
# with the code we have a very crowded violent plot so I move to a dotplot
ggplot(tbihsexpr_long,
       aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", colour = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", colour = "black") +
  stat_compare_means(comparisons = my_comparisons_group,
                     method = "wilcox.test",
                     label = "p.signif",
                     size  = 3) +
  facet_grid(CellType ~ Gene, scales = "free_y", switch = "y") +
  labs(title = "Heparan-Sulfate Gene Expression Across Conditions and Cell Types",
       y = "Expression", x = "Condition") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid     = element_blank(),
    axis.text.x    = element_text(angle = 45, hjust = 1),
    strip.text     = element_text(face = "bold"),
    legend.position = "bottom"
  )

#Dotplot on tbi_merged_obj_clean and HS genes are stored in hs_genes
#code below is the best one
DotPlot(tbi_merged_obj_clean,
        features = hs_genes,
        group.by = "cell_type",
        split.by = "condition",
        cols = c("lightgrey", "pink", "salmon", "firebrick3")) +  # Use 4 colors for split conditions
  theme_void() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5),
        axis.text.y = element_text()) +
  ggtitle("Heparan Sulfate Gene Expression Across Cell Types and Conditions")

# volcano plot to look at Mmp9, Adam10, Hspe per condition
# list of HS degradation (degrad) genes in mouse
#✅ 1. Define Heparan degradation ugennes 
hs_degrad_genes <- c("Mmp9", "Adam10", "Hpse")
colnames(tbi_merged_obj_clean@meta.data)

#✅ 2. Extract Data
# Fetch expression data + metadata
exprHS_degrad_data <- FetchData(tbi_merged_obj_clean, vars = c("condition", "cell_type", hs_degrad_genes))

# Reshape to long format for ggplot
exprHS_degrad_data <- exprHS_degrad_data %>%
  pivot_longer(cols = all_of(hs_degrad_genes), names_to = "Gene", values_to = "Expression") %>%
  rename(Condition = condition, CellType = cell_type)

#✅ 3. (Optional) Filter Specific Cell Types
# Focus on a few key cell types
celltypes_to_plot <- c("Astrocytes", "Microglia", "Endothelium Cells")
exprHS_degrad_data <- exprHS_degrad_data %>% filter(CellType %in% celltypes_to_plot)

#✅ 4. Generate Faceted Violin Plot with Stats
ggplot(exprHS_degrad_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width") +
  facet_wrap(~Gene, ncol = 3, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("WT Sham", "WT D3"), c("ecEXT1 Sham", "ecEXT1 D3"))) +
  scale_fill_manual(values = c("WT Sham" = "#E64B35", "WT D3" = "#4DBBD5",
                               "ecEXT1 Sham" = "#00A087", "ecEXT1 D3" = "#9370DB")) +
  theme_minimal(base_size = 13) +
  labs(title = "Heparan Sulfate Degradation Gene Expression Across Conditions",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.position = "right")

# For No gridlines, Clean white background, Black axis lines remain use the code below
#theme_classic() or theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(color = "black"))

ggplot(exprHS_degrad_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width") +
  facet_wrap(~Gene, ncol = 3, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("WT Sham", "WT D3"), c("ecEXT1 Sham", "ecEXT1 D3"))) +
  scale_fill_manual(values = c("WT Sham" = "#E64B35", "WT D3" = "#4DBBD5",
                               "ecEXT1 Sham" = "#00A087", "ecEXT1 D3" = "#9370DB")) +
  theme_minimal(base_size = 13) +
  labs(title = "Heparan Sulfate Degradation Gene Expression Across Conditions",
       y = "Expression Level") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  )

exprHS_degrad_data <- exprHS_degrad_data %>%
  mutate(CellType = recode(as.character(CellType), !!!tbi_cluster_annotations))

exprHS_degrad_data %>%
  filter(Gene == "Hpse","Adam10","Mmp9") %>%
  ggplot(aes(x = CellType, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", alpha = 0.3) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "right"
  ) +
  labs(
    title = "Hpse Expression Across Cell Types and Conditions",
    x = "Cell Type",
    y = "Expression Level"
  )


# code below Not good
VlnPlot(
  object = tbi_merged_obj_clean,
  features = c("Hpse", "Adam10", "Mmp9"),
  group.by = "cell_type",        # or "seurat_clusters" or any cell label
  split.by = "condition",        # split by condition for each gene
  pt.size = 0                    # hide individual dots to reduce clutter
) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#----------- still working on getting the proper order ------------------------------
table(exprHS_degrad_data$cell_type)
table(exprHS_degrad_data$Gene)

# volcano plot to look at Mmp9, Adam10, Hspe per condition and per cell type
# Extract expression data for those genes
expr_mat <- GetAssayData(tbi_merged_obj_clean, slot = "data")[hs_degrad_genes, ]

# Transpose and convert to data frame
expr_df1 <- as.data.frame(t(as.matrix(expr_mat)))

# Add metadata
expr_df1$cell_type <- tbi_merged_obj_clean$cell_type
expr_df1$Condition <- tbi_merged_obj_clean$condition
expr_df1$Barcode <- colnames(tbi_merged_obj_clean)

# Pivot to long format
exprHS_degrad_data <- expr_df1 %>%
  pivot_longer(cols = all_of(hs_degrad_genes), names_to = "Gene", values_to = "Expression")


table(exprHS_degrad_data$cell_type)  # should show all cell types
table(exprHS_degrad_data$Gene)       # should show Hpse, Adam10, Mmp9


# Grouping: combine cell type and condition
exprHS_degrad_data$Group <- paste(exprHS_degrad_data$cell_type, exprHS_degrad_data$Condition, sep = " | ")

exprHS_degrad_data$Condition <- factor(
  exprHS_degrad_data$Condition,
  levels = c("WT sham", "WT D3", "ecEXT1 Sham", "ecEXT1 D3")
)

# Set Group as a factor to enforce order in plot
exprHS_degrad_data$Group <- factor(exprHS_degrad_data$Group, 
                                   levels = unique(exprHS_degrad_data %>%
                                                     arrange(cell_type, Condition) %>%
                                                     pull(Group))
)


ggplot(exprHS_degrad_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(scale = "width", trim = FALSE) +
  facet_wrap(~ Gene, ncol = 1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    legend.position = "none",
    strip.text = element_text(size = 10)
  ) +
  labs(
    title = "Expression of HS Degradation Genes by Cell Type and Condition",
    x = "Cell Type | Condition",
    y = "Expression Level"
  )


unique(exprHS_degrad_data$Condition)

table(tbi_merged_obj_clean$condition)

table(exprHS_degrad_data$Condition, useNA = "always")
#----------------------------------- still working on the code above 05/16/25------------------------------------

saveRDS(tbi_merged_obj, "tbi_merged_obj_note.rds")
saveRDS(tbi_merged_obj_clean, "tbi_merged_obj_note1.rds")
saveRDS(tbicell_df, "tbi_merged_obj_note2.rds")

tbi_merged_obj <- readRDS("tbi_merged_obj_note.rds")
tbi_merged_obj_clean <- readRDS("tbi_merged_obj_note1.rds")
tbicell_df <- readRDS("tbi_merged_obj_note2.rds")








#---------------------------------
#volcano plot
De_results_1 <- FindMarkers(
  object = tbi_merged_obj_clean,
  ident.1 = "WT D3",
  ident.2 = "WT Sham",
  group.by = "condition",
  logfc.threshold = 0,
  min.pct = 0.1
)
De_results_1$gene <- rownames(De_results_1)
De_results_1$log10_pval <- -log10(De_results_1$p_val_adj + 1e-300)  # avoid -Inf

ggplot(De_results_1, aes(x = avg_log2FC, y = log10_pval)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  labs(
    title = "Volcano Plot: WT D3 vs WT Sham",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal()

# Add columns for plotting logic
De_results_1$Significant <- De_results_1$p_val_adj < 0.01 & abs(De_results_1$avg_log2FC) > 0.25

# Basic volcano plot
ggplot(De_results_1, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significant), alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  geom_text_repel(
    data = subset(De_results_1, Significant),
    aes(label = gene),
    size = 3,
    max.overlaps = 15
  ) +
  labs(
    title = "Volcano Plot: WT D3 vs WT Sham",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Label top genes by adjusted p-value and fold change
top1_labels <- De_results_1 %>%
  filter(Significant) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)   # Label top 10 most significant genes

# Volcano plot
ggplot(De_results_1, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significant), alpha = 0.7) +
  scale_color_manual(values = c("grey", "firebrick")) +
  geom_text_repel(
    data = top1_labels,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "Volcano Plot: WT D3 vs WT Sham",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted p-value"
  ) +
  theme_minimal() +
  theme(legend.position = "right")


#To generate a volcano plot that show which one is upregulated and what is downregulated

De_results_1 <- De_results_1 %>%
  mutate(regulation = case_when(
    p_val_adj < 0.01 & avg_log2FC > 0.25  ~ "Upregulated",
    p_val_adj < 0.01 & avg_log2FC < -0.25 ~ "Downregulated",
    TRUE                                  ~ "Not Significant"
  ))

# Now create top1_labels AFTER regulation column is added
top1_labels <- De_results_1 %>%
  filter(regulation != "Not Significant") %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)


ggplot(De_results_1, aes(x = avg_log2FC, y = -log10(p_val_adj), color = regulation)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(
    values = c("Upregulated" = "firebrick3", "Downregulated" = "steelblue", "Not Significant" = "gray80")
  ) +
  geom_text_repel(
    data = top1_labels,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "Volcano Plot: WT D3 vs WT Sham",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted p-value",
    color = "Regulation"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


ggplot(De_results_1, aes(x = avg_log2FC, y = -log10(p_val_adj), color = regulation)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(
    values = c("Upregulated" = "firebrick3", "Downregulated" = "steelblue", "Not Significant" = "gray80")
  ) +
  geom_text_repel(
    data = top1_labels,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "Volcano Plot: WT D3 vs WT Sham",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted p-value",
    color = "Regulation"
  ) +
  theme_minimal()

#---------------------------$$$$$$$$$$$$$------------------------------
# volcano plot looking at heparan sulfate related genes
De_results_1 <- De_results_1 %>%
  mutate(
    regulation = case_when(
      p_val_adj < 0.01 & avg_log2FC > 0.25  ~ "Upregulated",
      p_val_adj < 0.01 & avg_log2FC < -0.25 ~ "Downregulated",
      TRUE                                  ~ "Not Significant"
    ),
    HS_gene = ifelse(gene %in% hs_genes, "HS-related", "Other")
  )

ggplot(De_results_1, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = regulation, shape = HS_gene), alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "firebrick3", "Downregulated" = "steelblue", "Not Significant" = "gray80")) +
  scale_shape_manual(values = c("HS-related" = 17, "Other" = 16)) +  # triangle vs circle
  geom_text_repel(
    data = subset(De_results_1, HS_gene == "HS-related" & p_val_adj < 0.05),
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "Volcano Plot: WT D3 vs WT Sham",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted p-value",
    color = "Regulation",
    shape = "Gene Type"
  ) +
  theme_minimal()

# this is for color depend 
De_results_1 <- De_results_1 %>%
  mutate(
    regulation = case_when(
      p_val_adj < 0.01 & avg_log2FC > 0.25  ~ "Upregulated",
      p_val_adj < 0.01 & avg_log2FC < -0.25 ~ "Downregulated",
      TRUE                                  ~ "Not Significant"
    ),
    gene_type = ifelse(gene %in% hs_genes, "HS-related", "Other"),
    color_group = case_when(
      gene_type == "HS-related" ~ "HS-related",
      regulation == "Upregulated" ~ "Upregulated",
      regulation == "Downregulated" ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

ggplot(De_results_1, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = color_group), alpha = 0.8) +
  scale_color_manual(
    values = c(
      "Upregulated" = "firebrick3",
      "Downregulated" = "steelblue",
      "Not Significant" = "gray80",
      "HS-related" = "darkorchid4"
    )
  ) +
  geom_text_repel(
    data = subset(De_results_1, color_group == "HS-related" & p_val_adj < 0.05),
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "Volcano Plot: WT D3 vs WT Sham",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted p-value",
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# ecEXT1 D3 and ecEXT1 sham
De_results_2 <- FindMarkers(
  object = tbi_merged_obj_clean,
  ident.1 = "ecEXT1 D3",
  ident.2 = "ecEXT1 Sham",
  group.by = "condition",
  logfc.threshold = 0,
  min.pct = 0.1
)
De_results_2$gene <- rownames(De_results_2)
De_results_2$log10_pval <- -log10(De_results_2$p_val_adj + 1e-300)  # avoid -Inf

# Add columns for plotting logic
De_results_2$Significant <- De_results_2$p_val_adj < 0.01 & abs(De_results_2$avg_log2FC) > 0.25

# Basic volcano plot
ggplot(De_results_2, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significant), alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  geom_text_repel(
    data = subset(De_results_2, Significant),
    aes(label = gene),
    size = 3,
    max.overlaps = 15
  ) +
  labs(
    title = "Volcano Plot: ecEXT1 D3 vs ecEXT1 Sham",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Label top genes by adjusted p-value and fold change
top2_labels <- De_results_2 %>%
  filter(Significant) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)   # Label top 10 most significant genes

# Volcano plot
ggplot(De_results_2, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significant), alpha = 0.7) +
  scale_color_manual(values = c("grey", "firebrick")) +
  geom_text_repel(
    data = top2_labels,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "Volcano Plot: ecEXT1 D3 vs ecEXT1 Sham",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted p-value"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# for up and down regulated genes
De_results_2 <- De_results_2 %>%
  mutate(regulation = case_when(
    p_val_adj < 0.01 & avg_log2FC > 0.25  ~ "Upregulated",
    p_val_adj < 0.01 & avg_log2FC < -0.25 ~ "Downregulated",
    TRUE                                  ~ "Not Significant"
  ))

# Now create top2_labels AFTER regulation column is added
top2_labels <- De_results_2 %>%
  filter(regulation != "Not Significant") %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)


ggplot(De_results_2, aes(x = avg_log2FC, y = -log10(p_val_adj), color = regulation)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(
    values = c("Upregulated" = "firebrick3", "Downregulated" = "steelblue", "Not Significant" = "gray80")
  ) +
  geom_text_repel(
    data = top2_labels,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "Volcano Plot: ecEXT1 D3 vs ecEXT1 Sham",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted p-value",
    color = "Regulation"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


#  D3 and sham
De_results_3 <- FindMarkers(
  object = tbi_merged_obj_clean,
  ident.1 = "WT D3",
  ident.2 = "WT Sham",
  ident.3 = "ecEXT1 D3",
  ident.4 = "ecEXT1 Sham",
  group.by = "condition",
  logfc.threshold = 0,
  min.pct = 0.1
)
De_results_3$gene <- rownames(De_results_3)
De_results_3$log10_pval <- -log10(De_results_3$p_val_adj + 1e-300)  # avoid -Inf

# Add columns for plotting logic
De_results_3$Significant <- De_results_3$p_val_adj < 0.01 & abs(De_results_3$avg_log2FC) > 0.25

# Basic volcano plot
ggplot(De_results_3, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significant), alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  geom_text_repel(
    data = subset(De_results_3, Significant),
    aes(label = gene),
    size = 3,
    max.overlaps = 15
  ) +
  labs(
    title = "Volcano Plot: D3 vs Sham",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Label top genes by adjusted p-value and fold change
top3_labels <- De_results_3 %>%
  filter(Significant) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)   # Label top 10 most significant genes

# Volcano plot
ggplot(De_results_3, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significant), alpha = 0.7) +
  scale_color_manual(values = c("grey", "firebrick")) +
  geom_text_repel(
    data = top3_labels,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "Volcano Plot: D3 vs Sham",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted p-value"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# for up and down regulated genes
De_results_3 <- De_results_3 %>%
  mutate(regulation = case_when(
    p_val_adj < 0.01 & avg_log2FC > 0.25  ~ "Upregulated",
    p_val_adj < 0.01 & avg_log2FC < -0.25 ~ "Downregulated",
    TRUE                                  ~ "Not Significant"
  ))

# Now create top3_labels AFTER regulation column is added
top3_labels <- De_results_3 %>%
  filter(regulation != "Not Significant") %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

ggplot(De_results_3, aes(x = avg_log2FC, y = -log10(p_val_adj), color = regulation)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(
    values = c("Upregulated" = "firebrick3", "Downregulated" = "steelblue", "Not Significant" = "gray80")
  ) +
  geom_text_repel(
    data = top3_labels,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "Volcano Plot: D3 vs Sham",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted p-value",
    color = "Regulation"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
#-------------------------------------------------------------
# working on DEG now
#WT Group
#✅ Step 1: Subset Data by Major Cell Type
celltype1_list <- SplitObject(tbi_merged_obj_clean, split.by = "cell_type")

#✅ Step 2: Run DEG Analysis per Cell Type
deg1_list <- list()
for (ct in names(celltype1_list)) {
  obj1 <- celltype1_list[[ct]]
  obj1 <- NormalizeData(obj1)  # ensure normalized
  deg1 <- FindMarkers(obj1, ident.1 = "WT D3", ident.2 = "WT Sham", group.by = "condition",
                     logfc.threshold = 0.25, min.pct = 0.1)
  deg1$gene <- rownames(deg1)
  deg1$cell_type <- ct
  deg1_list[[ct]] <- deg1
}
all_deg1 <- do.call(rbind, deg1_list)

#✅ Step 3: Filter Significant DEGs
sig_deg1 <- subset(all_deg1, p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

#✅ Step 4: Pathway Analysis (GO / KEGG)
 
# Example for one cell type:
genesM <- sig_deg1$gene[sig_deg1$cell_type == "Microglia"]
entrez_ids_1 <- mapIds(org.Mm.eg.db, keys = genesM, keytype = "SYMBOL", column = "ENTREZID")

go_res_1 <- enrichGO(gene = entrez_ids_1,
                   OrgDb = org.Mm.eg.db,
                   ont = "BP", pAdjustMethod = "BH",
                   pvalueCutoff = 0.05, readable = TRUE)

kegg_res_1 <- enrichKEGG(gene = entrez_ids_1,
                       organism = "mmu", pvalueCutoff = 0.05)

#✅ Step 5: Visualize Pathways
dotplot(go_res_1, showCategory = 13, title = "Microglia GO Biological Processes")
dotplot(kegg_res_1, showCategory = 10, title = "KEGG Pathways")

#✅ Step 6: Explore DEGs Driving Pathways
head(go_res_1@result$geneID, n = 3)  # top GO term gene lists
head(kegg_res_1@result$geneID, n = 3)  # top kegg term gene lists
#Stop here

#-------------I am still working on this one 
# workflow for all cell type

celltype2_list <- SplitObject(tbi_merged_obj_clean, split.by = "cell_type")

all_results2 <- list()

for (ct in names(celltype2_list)) {
  message("Processing cell type: ", ct)
  
  obj2 <- celltype2_list[[ct]]
  obj2 <- NormalizeData(obj2)
  
  deg2 <- tryCatch({
    FindMarkers(obj2, ident.1 = "WT D3", ident.2 = "WT Sham", group.by = "condition",
                logfc.threshold = 0.25, min.pct = 0.1)
  }, error = function(e) {
    message("Skipped: ", ct, " due to error")
    return(NULL)
  })
  
  if (!is.null(deg2) && nrow(deg2) > 0) {
    deg2$gene <- rownames(deg2)
    deg2$cell_type <- ct
    deg2 <- deg2 %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
    
    # Map to ENTREZ IDs
    entrez_ids_2 <- mapIds(org.Mm.eg.db, keys = deg2$gene,
                         keytype = "SYMBOL", column = "ENTREZID",
                         multiVals = "first")
    
    entrez_ids_2 <- na.omit(entrez_ids_2)
    
    if (length(entrez_ids_2) > 5) {
      go2 <- enrichGO(gene = entrez_ids_2,
                     OrgDb = org.Mm.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     readable = TRUE)
      
      kegg2 <- enrichKEGG(gene = entrez_ids_2,
                         organism = "mmu",
                         pvalueCutoff = 0.05)
      
      all_results2[[ct]] <- list(
        DEGs = deg2,
        GO = go2,
        KEGG = kegg2
      )
    }
  }
}

#Check result
# Example: Check Microglia GO results
dotplot(all_results2$Microglia$GO, showCategory = 10)
dotplot(all_results2$Endotheluim$GO, showCategory = 10)
dotplot(all_results2$Neuron$GO, showCategory = 10)
dotplot(all_results2$Astrocytes$GO, showCategory = 10)


go_plots1 <- list()

for (ct in names(all_results2)) {
  if (!is.null(all_results2[[ct]]$GO)) {
    go_plot1 <- dotplot(all_results2[[ct]]$GO, showCategory = 4) +
      ggtitle(paste("GO: ", ct)) +
      theme_minimal(base_size = 10)
    
    go_plots1[[ct]] <- go_plot1
  }
}
wrap_plots(go_plots1, ncol = 2)

#To investigate pathways related to heparan sulfate and Alzheimer’s disease
#🔬 1. Identify Enriched Pathways in Your Data
# Example: Find enriched terms with "heparan" or "Alzheimer"
go2@result %>% filter(grepl("heparan|Alzheimer", Description, ignore.case = TRUE))
kegg2@result %>% filter(grepl("heparan|Alzheimer", Description, ignore.case = TRUE))

#📚 2. Directly Query KEGG/GO for Related Terms

# Search KEGG for Alzheimer
kegg_alz <- keggFind("pathway", "Alzheimer")
print(kegg_alz)

# Search KEGG for Heparan sulfate
kegg_hs <- keggFind("pathway", "heparan sulfate")
print(kegg_hs)

#🧠 3. KEGG Pathway IDs to Genes
alz_genes <- keggGet("mmu05010")[[1]]$GENE

#📊 4. Visualize Overlaps
alz_gene_symbols <- c("App", "Psen1", "Mapt", "Apoe", "Bace1", "Psen2", "Csnk1e", "Mapk1", "Gsk3b", "Cdk5", "Capn1", "Nefm", "Nefh",
                      "Nefl", "Ntrk1", "Ntrk2", "Grin1", "Grin2a", "Grin2b", "Gria1", "Gria2", "Dlg4", "Syn1", "Syp",
                      "Apba1", "Apba2", "Apba3", "Apbb1", "Apbb2", "Lrp1", "Clu", "Tnf", "Nos1", "Nos2", "Bax",
                      "Bcl2", "Casp3", "Casp9", "Cycs", "Aifm1", "Ndufa9", "Ndufb8", "Uqcrc1", "Cox4i1", "Atp5a1",
                      "Atp5b", "Slc25a4", "Slc25a5", "Slc25a6", "Pink1", "Park2", "Sod1", "Sod2", "Gpx1")# expand this list

hs_gene_symbols <- c("Ext1", "Ext2", "Extl1", "Extl2", "Extl3", "Ndst1", "Ndst2", "Ndst3", "Ndst4","Hs2st1", "Hs3st1", "Hs3st2", "Hs3st3a1", "Hs3st3b1", "Hs6st1", "Hs6st2", "Hs6st3","Sulf1", "Sulf2","Gpc1", "Gpc2", "Gpc3", "Gpc4", "Gpc5", "Gpc6", "Sdc1", "Sdc2", "Sdc3", "Sdc4")

De_results_1$alz_related <- ifelse(rownames(De_results_1) %in% alz_gene_symbols, "Alzheimer gene", "Other")

topalz_labels <- De_results_1 %>%
  filter(alz_related == "Alzheimer gene", p_val_adj < 0.05) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)  # or use a different number if needed

ggplot(De_results_1, aes(x = avg_log2FC, y = -log10(p_val_adj), color = alz_related)) +
  geom_point() +
  scale_color_manual(values = c("Alzheimer gene" = "red", "Other" = "grey")) +
  geom_text_repel(
    data = topalz_labels,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

#---------------------volcano plot WT Group HS and ALZ----------------------------------------------

# Create a column indicating gene category
De_results_1$label_group <- case_when(
  De_results_1$gene %in% alz_gene_symbols ~ "Alzheimer gene",
  De_results_1$gene %in% hs_gene_symbols ~ "HS gene",
  TRUE ~ "Other"
)

# For labeling, pick significant Alzheimer or HS genes
top_labeled <- De_results_1 %>%
  filter(label_group != "Other", p_val_adj < 0.05) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 20)  # You can increase this number

# Volcano plot
ggplot(De_results_1, aes(x = avg_log2FC, y = -log10(p_val_adj), color = label_group)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(
    values = c("Alzheimer gene" = "red", "HS gene" = "blue", "Other" = "gray80")
  ) +
  geom_text_repel(
    data = top_labeled,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = " WT Volcano Plot Highlighting Alzheimer & HS Genes",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted p-value",
    color = "Gene Group"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# Add a combined group based on gene category and fold change direction
De_results_1$label_group <- case_when(
  De_results_1$gene %in% alz_gene_symbols & De_results_1$avg_log2FC > 0 ~ "Alz Up",
  De_results_1$gene %in% alz_gene_symbols & De_results_1$avg_log2FC < 0 ~ "Alz Down",
  De_results_1$gene %in% hs_gene_symbols & De_results_1$avg_log2FC > 0 ~ "HS Up",
  De_results_1$gene %in% hs_gene_symbols & De_results_1$avg_log2FC < 0 ~ "HS Down",
  TRUE ~ "Other"
)

# Label significant Alzheimer or HS genes
top_labeled <- De_results_1 %>%
  filter(label_group != "Other", p_val_adj < 0.05) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 20)

# Custom colors
custom_colors <- c(
  "Alz Up" = "firebrick2",
  "Alz Down" = "darkred",
  "HS Up" = "steelblue2",
  "HS Down" = "midnightblue",
  "Other" = "gray80"
)

# Plot
ggplot(De_results_1, aes(x = avg_log2FC, y = -log10(p_val_adj), color = label_group)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = custom_colors) +
  geom_text_repel(
    data = top_labeled,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "Volcano Plot: Alzheimer & HS Genes up and downregulated",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted p-value",
    color = "Gene Type"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


#------------still working on the gene network----------------------------------
# Use STRINGdb for protein-protein (gene-gene) interaction
# Example: Top DEGs (adjust as needed)
top_deg_genes <- De_results_1 %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>%
  pull(gene)

string_db <- STRINGdb$new(
  version = "11.5",  # or latest version
  species = 10090,   # mouse; use 9606 for human
  score_threshold = 400,
  input_directory = ""
)
mapped_genes <- string_db$map(data.frame(gene = top_deg_genes), "gene", removeUnmappedRows = TRUE)

interactions <- string_db$get_interactions(mapped_genes$STRING_id)
graph <- graph_from_data_frame(interactions[, c("from", "to")], directed = FALSE)
plot(graph, vertex.label = NA, vertex.size = 5)

graph_tbl <- as_tbl_graph(graph)

ggraph(graph_tbl, layout = "fr") +
  geom_edge_link(alpha = 0.8) +
  geom_node_point(size = 4) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_void()


saveRDS(tbi_merged_obj, "tbi_merged_obj_note.rds")
saveRDS(tbi_merged_obj_clean, "tbi_merged_obj_note1.rds")
saveRDS(tbicell_df, "tbi_merged_obj_note2.rds")

tbi_merged_obj <- readRDS("tbi_merged_obj_note.rds")
tbi_merged_obj_clean <- readRDS("tbi_merged_obj_note1.rds")
tbicell_df <- readRDS("tbi_merged_obj_note2.rds")







#





saveRDS(tbi_merged_obj, "tbi_merged_obj_note.rds")
saveRDS(tbi_merged_obj_clean, "tbi_merged_obj_note1.rds")
saveRDS(tbicell_df, "tbi_merged_obj_note2.rds")

test_obj_save <- readRDS("tbi_merged_obj_note.rds")
str(test_obj_save)


tbi_merged_obj <- readRDS("tbi_merged_obj_note.rds")
tbi_merged_obj_clean <- readRDS("tbi_merged_obj_note1.rds")
tbicell_df <- readRDS("tbi_merged_obj_note2.rds")























saveRDS(tbi_merged_obj, "tbi_merged_obj_note.rds")
saveRDS(tbi_merged_obj_clean, "tbi_merged_obj_note1.rds")
saveRDS(tbicell_df, "tbi_merged_obj_note2.rds")

tbi_merged_obj <- readRDS("tbi_merged_obj_note.rds")
tbi_merged_obj_clean <- readRDS("tbi_merged_obj_note1.rds")
tbicell_df <- readRDS("tbi_merged_obj_note2.rds")





































ggplot(gpc1_data, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width", color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_compare_means(comparisons = my_comparisons_group, method = "wilcox.test", label = "p.signif") +
  ggtitle("Gpc1 Expression Across Conditions") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )
ggsave("HS_gpc1_expression.png", dpi = 600, width = 12, height = 8)


saveRDS(tbi_merged_obj, "tbi_merged_obj_note.rds")
saveRDS(tbi_merged_obj_clean, "tbi_merged_obj_note1.rds")
saveRDS(tbicell_df, "tbi_merged_obj_note2.rds")

tbi_merged_obj <- readRDS("tbi_merged_obj_note.rds")
tbi_merged_obj_clean <- readRDS("tbi_merged_obj_note1.rds")
tbicell_df <- readRDS("tbi_merged_obj_note2.rds")























#------------------------ still working on generate generating a violin plot and heatmap------------------------------- 

























#Pathways investigation
#Heparan sulfate genes
#chain elongation
hs_genes_elongation <- c(
  "Ext1", "Ext2", "Extl1", "Extl2", "Extl3")   # chain elongation
# First make the VlnPlot separately
vln_elong <- VlnPlot(
  tbi_merged_obj_clean, 
  features = hs_genes_elongation, 
  group.by = "condition", 
  pt.size = 0.1
) + 
  ggtitle("HS-Related Elonagtion Gene Expression Across Conditions") +
  theme_minimal(base_size = 14)  # use minimal theme, not void

# Then add statistical comparison
vln_elong + 
  stat_compare_means(
    label = "p.signif", 
    method = "kruskal.test"
  )

# 📈 Violin plot of Ext1 gene across conditions
VlnPlot(tbi_merged_obj_clean, features = hs_genes_elongation, group.by = "condition", pt.size = 0.1) +
  stat_compare_means(label = "p.signif", method = "kruskal.test") +
  ggtitle("Ext1 Expression Across Conditions") +
  theme_void()








# Fetch expression data + metadata
expr_data_hs <- FetchData(tbi_merged_obj_clean, vars = c("condition", "manual_annotation", hs_genes))

# Reshape to long format for ggplot
expr_long <- expr_data %>%
  pivot_longer(cols = all_of(hs_genes), names_to = "Gene", values_to = "Expression") %>%
  rename(Condition = condition, CellType = manual_annotation)
vln_plot <- ggplot(expr_long, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(trim = FALSE, scale = "width") +
  facet_wrap(~Gene, ncol = 3, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("WT Sham", "WT D3"), c("ecEXT1 Sham", "ecEXT1 D3"))) +
  scale_fill_manual(values = c("WT Sham" = "#E64B35", "WT D3" = "#4DBBD5",
                               "ecEXT1 Sham" = "#00A087", "ecEXT1 D3" = "#9370DB")) +
  theme_minimal(base_size = 13) +
  labs(title = "Heparan Sulfate Gene Expression Across Conditions",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.position = "right")













#Heatmap generate
#✅ Step 1: List of Heparan Sulfate Genes
hs_genes <- c(
  "Ext1", "Ext2", "Extl1", "Extl2", "Extl3",   # chain elongation
  "Ndst1", "Ndst2", "Ndst3", "Ndst4",          # N-deacetylase/N-sulfotransferase
  "Hs2st1", "Hs3st1", "Hs3st2", "Hs3st3a1", "Hs3st3b1", "Hs6st1", "Hs6st2", "Hs6st3",  # sulfotransferases
  "Sulf1", "Sulf2",                            # sulfatases
  "Gpc1", "Gpc2", "Gpc3", "Gpc4", "Gpc5", "Gpc6",  # glypicans (core proteins)
  "Sdc1", "Sdc2", "Sdc3", "Sdc4"              # syndecans
)
#✅ Step 2: Calculate Average Expression by Condition
avg_exp <- AverageExpression(tbi_merged_obj_clean, features = hs_genes, assays = "RNA", slot = "data")$RNA

# Optional: scale rows for heatmap
avg_exp_scaled <- t(scale(t(avg_exp)))  # Z-score normalization by gene

pheatmap(avg_exp_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         fontsize_row = 10,
         fontsize_col = 12,
         main = "Heparan Sulfate Gene Expression by Condition")



Heatmap(
  avg_exp_scaled,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  row_names_gp = grid::gpar(fontsize = 10),
  column_names_gp = grid::gpar(fontsize = 12),
  heatmap_legend_param = list(title = "Expression"),
  name = "HS Genes",
  column_title = "Heparan Sulfate Gene Expression by Condition"   # add nice title here
)



# Suppose you have 4 conditions (columns) repeated: WT Sham, WT D3, ecEXT1 Sham, ecEXT1 D3
label_conditions <- rep(c("WT Sham", "WT D3", "ecEXT1 Sham", "ecEXT1 D3"), each = n_columns_per_condition)

# or if you know the exact conditions for each column:
conditions <- c("WT Sham", "WT Sham", ..., "WT D3", ..., "ecEXT1 Sham", ..., "ecEXT1 D3", ...)

top_anno <- HeatmapAnnotation(
  Condition = conditions,
  col = list(
    Condition = c(
      "WT Sham" = "lightblue",
      "WT D3" = "skyblue3",
      "ecEXT1 Sham" = "pink",
      "ecEXT1 D3" = "red3"
    )
  )
)

Heatmap(
  avg_exp_scaled,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  row_names_gp = grid::gpar(fontsize = 10),
  column_names_gp = grid::gpar(fontsize = 12),
  heatmap_legend_param = list(title = "Expression"),
  name = "HS Genes",
  top_annotation = top_anno,
  column_title = "Heparan Sulfate Gene Expression by Condition"
)









# 📈 Violin plot of Ext1 gene across conditions
VlnPlot(tbi_merged_obj_clean, features = hs_genes, group.by = "condition", pt.size = 0.1) +
  stat_compare_means(label = "p.signif", method = "kruskal.test") +
  ggtitle("Ext1 Expression Across Conditions") +
  theme_void()


# First make the VlnPlot separately
vln <- VlnPlot(
  tbi_merged_obj_clean, 
  features = hs_genes, 
  group.by = "condition", 
  pt.size = 0.1
) + 
  ggtitle("HS-Related Gene Expression Across Conditions") +
  theme_minimal(base_size = 14)  # use minimal theme, not void

# Then add statistical comparison
vln + 
  stat_compare_means(
    label = "p.signif", 
    method = "kruskal.test"
  )














saveRDS(tbi_merged_obj, "tbi_merged_obj_note.rds")
saveRDS(tbi_merged_obj_clean, "tbi_merged_obj_note1.rds")
saveRDS(tbicell_df, "tbi_merged_obj_note2.rds")
tbi_merged_obj <- readRDS("tbi_merged_obj_note.rds")
tbi_merged_obj_clean <- readRDS("tbi_merged_obj_note1.rds")
tbicell_df <- readRDS("tbi_merged_obj_note2.rds")



