# Scripts for DGE/DIE analysis and figures

#Setup -------
library(tidyverse)
library(wesanderson)
library(fastqcr)
library(DESeq2)
library(gProfileR)
library(ggfortify)
library(pheatmap)
library(corrplot)
library(biomaRt)
library(ggrepel)
library(readxl)
library(ComplexHeatmap)
library(colorRamp2)
library(ggVennDiagram)

#Datasets -------
# Data and count files are based on SQANTI classification files used in other scripts

barcodes <- read.csv("~/path/to/barcods")
barcodes$sex[barcodes$sex=="F"] <- "XX"
barcodes$sex[barcodes$sex=="M"] <- "XY"


agregated_flair <- flair_filter[,c(7,49:84)] # adapt to count and gene/isoform names as needed
agregated_flair <- na.omit(agregated_flair)
agregated_flair <- agregated_flair %>% group_by(associated_gene) %>% # change column names as needed
  summarise(across(BC01:BC42,~sum(.x, na.rm = TRUE))) # change column names as needed
agregated_flair <- column_to_rownames(agregated_flair, var="associated_gene") # change column names as needed

agregated_talon <- talon_filter[,c(7,49:84)] # adapt to count and gene/isoform names as needed
agregated_talon <- na.omit(agregated_talon)
agregated_talon <- agregated_talon %>% group_by(associated_gene) %>% # change column names as needed
  summarise(across(BC01:BC42,~sum(.x, na.rm = TRUE))) # change column names as needed
agregated_talon <- column_to_rownames(agregated_talon, var="associated_gene") # change column names as needed

agregated_isoseq <- isoseq_filter[,c(7,49:84)] # adapt to count and gene/isoform names as needed
agregated_isoseq <- na.omit(agregated_isoseq)
agregated_isoseq <- agregated_isoseq %>% group_by(associated_gene) %>% # change column names as needed
  summarise(across(BC01:BC42,~sum(.x, na.rm = TRUE))) # change column names as needed
agregated_isoseq <- column_to_rownames(agregated_isoseq, var="associated_gene") # change column names as needed


coldata <- barcodes
coldata$time <- factor(coldata$timepoint)
coldata$sex <- factor(coldata$sex)
coldata$condition <- factor(coldata$condition)
coldata$line <- factor(coldata$line)
rownames(coldata)<- mgl_barcodes$sample_id

# VOLCANO PLOT function --------

volcano_plot <- function(x) {
  de <- as.data.frame(x)
  

  if (!is.null(rownames(de))) {
    de$gene_symbol <- rownames(de)
  } else {
    stop("The input data frame must have row names representing gene symbols.")
  }
  
 
  de$diffexpressed <- "Neither"
  
  de$diffexpressed[de$log2FoldChange > 1.5 & de$padj < 0.05] <- "Up"
  de$diffexpressed[de$log2FoldChange < -1.5 & de$padj < 0.05] <- "Down"
  
  # Add a column for labels
  de$delabel <- NA
  de$delabel[de$diffexpressed != "Neither"] <- de$gene_symbol[de$diffexpressed != "Neither"]
  

  de <- de[order(de$padj), ]
  
  de$genelabels <- ""
  
  # Ensure we do not exceed the number of rows in de when labeling
  num_labels <- min(20, nrow(de))
  de$genelabels[1:num_labels] <- rownames(de)[1:num_labels]
  
  # Create the plot
  ggplot(data = de, 
         aes(x = log2FoldChange, 
             y = -log10(padj), 
             col = diffexpressed, 
             label = delabel)) +
    geom_point() +
    geom_text_repel(aes(label = genelabels), show.legend = FALSE) +
    scale_color_manual(values = c("Down" = "steelblue", "Neither" = "grey", "Up" = "tomato"), 
                       name = "Direction") +
    geom_vline(xintercept = c(-1.5, 1.5), col = "red") +
    geom_hline(yintercept = -log10(0.05), col = "red") +
    labs(x = "Log2(Fold Change)",
         y = "-Log10(adjusted p-value)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          panel.grid = element_blank(), 
          axis.line = element_line(color = "#5f5f5f"),
          legend.text = element_text(size = 12),
          axis.title.x = element_text(size = 18),
          axis.text.x  = element_text(size = 12),
          axis.title.y = element_text(size = 18),
          axis.text.y  = element_text(vjust = 0.5, size = 12))
}


# Dataset make function --------

dataset_make <- function(x, filename){de <- as.data.frame(x)

de$diffexpressed <- "Neither"
de$diffexpressed[de$log2FoldChange > 1.5 & de$padj < 0.05] <- "Up"
de$diffexpressed[de$log2FoldChange < -1.5 & de$padj < 0.05] <- "Down"
de$delabel <- NA
de$gene_symbol <- rownames(x)
de$delabel[de$diffexpressed != "Neither"] <- de$gene_symbol[de$diffexpressed != "Neither"]
                                      
de <- de[order(de$padj), ] 
de$Gene <- rownames(de)

de$genelabels <- ""
de$genelabels[1:20] <- rownames(de)[1:20]

filename <- de
}


# HEATMAP PLOT function--------

create_heatmap <- function(de, dds, genes, coldata, barcodes, column_split, plot_name) {
  dataset <- de  %>% filter(rownames(de) %in% genes)
  
  rld <- varianceStabilizingTransformation(dds) 
  mat <- assay(rld)[rownames(dataset), rownames(coldata)]
  based_mean <- rowMeans(mat)
  mat_scaled <- t(apply(mat, 1, scale))
  colnames(mat_scaled) <- colnames(mat)
  
  l2_val <- as.matrix(dataset[ , c("log2FoldChange")])
  colnames(l2_val) <- "logFC"
  mean <- as.matrix(dataset[ , c("baseMean")])
  colnames(mean) <- "AvgExpr"
  
  col_logFC <- colorRamp2(c(min(l2_val), 0, max(l2_val)), c("#104e8b", "#e8e8e8", "#b22222"))
  col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("#e8e8e8", "#b22222"))
  col_heat <- colorRamp2(c(min(mat_scaled), 0, max(mat_scaled)), c("#104e8b", "#e8e8e8", "#b22222"))
  
  metadata_heatmap <- barcodes %>% 
    filter(sample_id %in% colnames(mat_scaled)) %>%
    dplyr::select(sample_id, timepoint, sex, line) %>%
    arrange(sample_id) %>%
    unique
  
  ha_column <- HeatmapAnnotation(df = data.frame(Timepoint = metadata_heatmap$timepoint, 
                                                 Sex = metadata_heatmap$sex, 
                                                 Line = metadata_heatmap$line),
                                 col = list(Timepoint = c("D1" = "#000004ff", "D14" = "#721f81ff"), # adapt as needed
                                            Sex = c("XX" = "#f1605dff", "XY" = "#fcfdbfff"),
                                            Line = c("CTF007" = "#440154ff", "SCTI" = "#414487ff", 
                                                     "CTF069" = "#22a884ff", "CTM127" = "#2a788eff", 
                                                     "CTM014" = "#7ad151ff", "M3_36S" = "#fde725ff")))
  
  if (!is.null(column_split)) {
    column_split_val <- metadata_heatmap[[column_split]]
  } else {
    column_split_val <- NULL
  }
  
  plot_name <- grid.grabExpr(draw(Heatmap(mat_scaled, cluster_rows = T, cluster_columns = T,
                                          column_split = column_split_val, column_title = NULL,
                                          show_parent_dend_line = FALSE,
                                          show_row_dend = F,
                                          col = col_heat,
                                          top_annotation = ha_column,
                                          heatmap_legend_param = list(title = "Z-score")) + 
                                    Heatmap(l2_val, cluster_rows = F, name = "logFC", col = col_logFC,
                                            width = ncol(l2_val) * unit(1.5, "cm")) +
                                    Heatmap(mean, cluster_rows = F, name = "Avg Expr", col = col_AveExpr,
                                            width = ncol(mean) * unit(1.5, "cm"))))
  
  plot_name <- patchwork::wrap_elements(plot_name)
  
  return(plot_name)
}


#Differential Expression Analysis (FLAIR)-------
dds_flair <- DESeqDataSetFromMatrix(countData = agregated_flair,
                                    colData = coldata_flair,
                                    design= ~ sex + time + sex:time)

rownames(dds_flair) <- rownames(agregated_flair) #gives us which features are used in study (genes)
colnames(dds_flair) #displays samples
counts(dds_flair) #gives count table 
colData(dds_flair) #gives experimental setup

dds_flair <- dds_flair[ rowSums(DESeq2::counts(dds_flair)) > 10, ]
dds_flair <- DESeq(dds_flair)


mcols(dds_flair) <- DataFrame(mcols(dds_flair))
mcols(dds_flair)

dds_flair$time <- relevel(dds_flair$time, ref = "D1") #adapt if needed

dds_flair <- DESeq(dds_flair)
res_flair <- results(dds_flair, alpha=0.05)
res_flair %>% write.csv("flair_time_condition_results.csv")


DEresults_time_flair <- results(dds_flair, contrast = c("time", "D1" , "D14"), alpha=0.05)  #adapt if needed
DEresults_time_flair <- DEresults_time_flair[order(DEresults_time_flair$pvalue) ,]
print(DEresults_time_flair) 
DEresults_time_flair %>% write.csv("flair_deresults_time.csv") 

DEresults_sex_flair <- results(dds_flair, contrast = c("sex", "XY" , "XX"), alpha=0.05)
DEresults_sex_flair <- DEresults_sex_flair[order(DEresults_sex_flair$pvalue) ,]
print(DEresults_sex_flair) 
DEresults_sex_flair %>% write.csv("flair_deresults_sex.csv") 

##Downstream plots (FLAIR) -----

## Dispersion plot
DESeq2::plotDispEsts(dds_flair, ylim = c(1e-6, 1e1),
                     xlab = "mean of normalised counts")

DESeq2::plotMA(object = dds_flair, ylim = c(-12,12)
)

ggplot(data = as.data.frame(DEresults_time_flair), 
       aes(x = pvalue)) +
  geom_histogram(bins = 200) 

ggplot(data = as.data.frame(DEresults_sex_flair), 
       aes(x = pvalue)) +
  geom_histogram(bins = 200) 



countsNormalized_flair <- DESeq2::counts(dds_flair, normalized = TRUE)

countsNormalized %>% write.csv("flair_norm_counts.csv")

selectedGenes <- names(sort(apply(countsNormalized, 1, var),
                            decreasing = TRUE))

rld_flair <- varianceStabilizingTransformation(dds_flair) 

labels_man <- c("F_D1"="D1 XX","M_D1"="D1 XY","F_D14"="D14 XX","M_D14"="D14 XY")  #adapt if needed


pca_data_flair <- DESeq2::plotPCA(rld_flair, intgroup = c("condition","line"),returnData=TRUE) 
pca_data_flair$method <- "FLAIR"

annotated_pca_flair <- ggplot(pca_data, aes(x=PC1, y= PC2, color=condition, shape=line))+
  geom_point(size=3, guide=guide_legend(nrow=2))+
  scale_color_viridis(discrete=TRUE, labels=labels_man, guide=guide_legend(nrow=2))+
  scale_shape_manual(values =c(8,15,16,17,18,25))+
  geom_text(aes(label=substr(name, start = 1, stop = 10)),show.legend = FALSE,vjust=3,check_overlap = TRUE,size = 2)+
  labs(color="Condition",shape="Line") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"), legend.position="bottom",
        legend.text = element_text(size=12),
        title=element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.x  = element_text(size=12),
        axis.title.y = element_text(size=18),
        axis.text.y  = element_text(vjust=0.5, size=12))


#DEG/DIE (FLAIR) --------

de_flair_mgl_time <- dataset_make(DEresults_time_flair)
de_flair_mgl_time$method <- "FLAIR"
de_flair_mgl_time %>% write.csv("de_flair_mgl_time.csv")

flair_deg_time <- de_flair_mgl_time %>% filter(diffexpressed == "Up" | diffexpressed == "Down")
flair_deg_time$method <- "FLAIR"
flair_deg_time$Gene <- flair_filter$associated_gene[match(flair_deg_time$delabel, flair_filter$isoform)]
flair_deg_time$Catagory <- flair_filter$structural_category[match(flair_deg_time$delabel, flair_filter$isoform)]

flair_deg_time %>% write.csv("flair_deg_time.csv")

de_flair_sex <- dataset_make(DEresults_sex_flair)
de_flair_sex$method <- "FLAIR"

de_flair_mgl_sex %>% write.csv("deg_flair_sex.csv")

flair_deg_sex <- de_flair_mgl_sex %>% filter(diffexpressed == "Up" | diffexpressed == "Down")
flair_deg_sex$method <- "FLAIR"
flair_deg_sex$Gene <- flair_filter$associated_gene[match(flair_deg_sex$delabel, flair_filter$isoform)]
flair_deg_sex$Catagory <- flair_filter$structural_category[match(flair_deg_sex$delabel, flair_filter$isoform)]
flair_deg_sex %>% write.csv("flair_deg_sex.csv")

genes <- rownames(de_flair_time)
genes_sex <- rownames(de_flair_sex)



flair_deg_time_heatmap <- create_heatmap(flair_deg_time, dds_flair, genes, coldata_flair, barcodes, column_split = "NULL", time)

flair_deg_sex_heatmap <- create_heatmap(flair_deg_sex, dds_flair, genes_sex, coldata_flair, barcodes, column_split = "NULL", sex)


#Differential Expression Analysis (Isoseq)-------

dds_isoseq <- DESeqDataSetFromMatrix(countData = agregated_isoseq,
                                     colData = coldata,
                                     design= ~ sex + time + sex:time)

rownames(dds_isoseq) <- rownames(agregated_isoseq) #gives us which features are used in study (genes)
colnames(dds_isoseq) #displays samples
counts(dds_isoseq) #gives count table 
colData(dds_isoseq) #gives experimental setup

dds_isoseq <- dds_isoseq[ rowSums(DESeq2::counts(dds_isoseq)) > 10, ]
dds_isoseq <- DESeq(dds_isoseq)


mcols(dds_isoseq) <- DataFrame(mcols(dds_isoseq))
mcols(dds_isoseq)

dds_isoseq$time <- relevel(dds_isoseq$time, ref = "D1") #adapt if needed

dds_isoseq <- DESeq(dds_isoseq)
res_isoseq <- results(dds_isoseq, alpha=0.05)
res_isoseq %>% write.csv("isoseq_time_condition_results.csv")


DEresults_time_isoseq <- results(dds_isoseq, contrast = c("time", "D1" , "D14"), alpha=0.05) #adapt if needed
DEresults_time_isoseq <- DEresults_time_isoseq[order(DEresults_time_isoseq$pvalue) ,]
print(DEresults_time_isoseq) 
DEresults_time_isoseq %>% write.csv("isoseq_deresults_time.csv") 

DEresults_sex_isoseq <- results(dds_isoseq, contrast = c("sex", "XY" , "XX"), alpha=0.05)
DEresults_sex_isoseq <- DEresults_sex_isoseq[order(DEresults_sex_isoseq$pvalue) ,]
print(DEresults_sex_isoseq) 
DEresults_sex_isoseq %>% write.csv("isoseq_deresults_sex.csv") 


##Downstream plots (Isoseq) -----

## Dispersion plot
DESeq2::plotDispEsts(dds_isoseq, ylim = c(1e-6, 1e1),
                     xlab = "mean of normalised counts")

DESeq2::plotMA(object = dds_isoseq, ylim = c(-12,12)
)

ggplot(data = as.data.frame(DEresults_time_isoseq), 
       aes(x = pvalue)) +
  geom_histogram(bins = 200) 

ggplot(data = as.data.frame(DEresults_sex_isoseq), 
       aes(x = pvalue)) +
  geom_histogram(bins = 200) 



countsNormalized_isoseq <- DESeq2::counts(dds_isoseq, normalized = TRUE)

countsNormalized %>% write.csv("isoseq_norm_counts.csv")

selectedGenes <- names(sort(apply(countsNormalized, 1, var),
                            decreasing = TRUE))

rld_iso <- varianceStabilizingTransformation(dds_isoseq) 

pca_data_iso <- DESeq2::plotPCA(rld_iso, intgroup = c("condition","line"),returnData=TRUE) 
pca_data_iso$method <- "Isoseq"

annotated_pca_isoseq <- ggplot(pca_data, aes(x=PC1, y= PC2, color=condition, shape=line))+
  geom_point(size=3, guide=guide_legend(nrow=2))+
  scale_color_viridis(discrete=TRUE, labels=labels_man)+
  scale_shape_manual(values =c(8,15,16,17,18,25))+
  geom_text(aes(label=substr(name, start = 1, stop = 10)),show.legend = FALSE,vjust=3,check_overlap = TRUE,size = 2)+
  labs(color="Condition",shape="Line") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"), legend.position="bottom",
        legend.text = element_text(size=12),
        title=element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.x  = element_text(size=12),
        axis.title.y = element_text(size=18),
        axis.text.y  = element_text(vjust=0.5, size=12))

# DEG/DIE (Isoseq)--------
de_isoseq_time <- dataset_make(DEresults_time_isoseq)
de_isoseq_time$method <- "Isoseq"
de_isoseq_time %>% write.csv("de_isoseq_mgl_time.csv")

isoseq_deg_time <- de_isoseq_time %>% filter(diffexpressed == "Up" | diffexpressed == "Down")
isoseq_deg_time$method <- "Isoseq"
isoseq_deg_time$Gene <- isoseq_filter$associated_gene[match(isoseq_deg_time$delabel, isoseq_filter$isoform)]
isoseq_deg_time$Catagory <- isoseq_filter$structural_category[match(isoseq_deg_time$delabel, isoseq_filter$isoform)]


isoseq_deg_time %>% write.csv("isoseq_deg_time.csv")
genes <- rownames(isoseq_deg_time)


de_isoseq_sex <- dataset_make(DEresults_sex_isoseq)
de_isoseq_sex$method <- "Isoseq"

de_isoseq_sex %>% write.csv("deg_isoseq_sex.csv")

isoseq_deg_sex <- de_isoseq_sex %>% filter(diffexpressed == "Up" | diffexpressed == "Down")
isoseq_deg_sex$method <- "Isoseq"
isoseq_deg_sex$Gene <- isoseq_filter$associated_gene[match(isoseq_deg_sex$delabel, isoseq_filter$isoform)]
isoseq_deg_sex$Catagory <- isoseq_filter$structural_category[match(isoseq_deg_sex$delabel, isoseq_filter$isoform)]


isoseq_deg_sex %>% write.csv("isoseq_deg_sex.csv")
genes_sex <- rownames(isoseq_deg_sex)




isoseq_deg_time_heatmap <- create_heatmap(isoseq_deg_time, dds_isoseq, genes, coldata, barcodes, column_split = "NULL", time)

isoseq_deg_sex_heatmap <- create_heatmap(isoseq_deg_sex, dds_isoseq, genes_sex, coldata, barcodes, column_split = "NULL", sex)

#Differential Expression Analysis (TALON)-------

dds_talon <- DESeqDataSetFromMatrix(countData = talon,
                                    colData = coldata,
                                    design= ~ sex + time + sex:time)

rownames(dds_talon) <- rownames(talon) #gives us which features are used in study (genes)
colnames(dds_talon) #displays samples
counts(dds_talon) #gives count table 
colData(dds_talon) #gives experimental setup

dds_talon <- dds_talon[ rowSums(DESeq2::counts(dds_talon)) > 10, ]
dds_talon <- DESeq(dds_talon)


mcols(dds_talon) <- DataFrame(mcols(dds_talon))
mcols(dds_talon)

dds_talon$time <- relevel(dds_talon$time, ref = "D1") #adapt if needed

dds_talon <- DESeq(dds_talon)
res_talon <- results(dds_talon, alpha=0.05)
res_talon %>% write.csv("talon_time_condition_results.csv")


DEresults_time_talon <- results(dds_talon, contrast = c("time", "D1" , "D14"), alpha=0.05) #adapt if needed
DEresults_time_talon <- DEresults_time_talon[order(DEresults_time_talon$pvalue) ,]
print(DEresults_time_talon) 
DEresults_time_talon %>% write.csv("talon_deresults_time.csv") 

DEresults_sex_talon <- results(dds_talon, contrast = c("sex", "XY" , "XX"), alpha=0.05)
DEresults_sex_talon <- DEresults_sex_talon[order(DEresults_sex_talon$pvalue) ,]
print(DEresults_sex_talon) 
DEresults_sex_talon %>% write.csv("talon_deresults_sex.csv") 


##Downstream plots-----

## Dispersion plot
DESeq2::plotDispEsts(dds_talon, ylim = c(1e-6, 1e1),
                     xlab = "mean of normalised counts")

DESeq2::plotMA(object = dds_talon, ylim = c(-12,12)
)

ggplot(data = as.data.frame(DEresults_time_talon), 
       aes(x = pvalue)) +
  geom_histogram(bins = 200) 

ggplot(data = as.data.frame(DEresults_sex_talon), 
       aes(x = pvalue)) +
  geom_histogram(bins = 200) 



countsNormalized_talon <- DESeq2::counts(dds_talon, normalized = TRUE)

countsNormalized %>% write.csv("talon_norm_counts.csv")

selectedGenes <- names(sort(apply(countsNormalized, 1, var),
                            decreasing = TRUE))

rld_talon <- varianceStabilizingTransformation(dds_talon) 


pca_data_talon <- DESeq2::plotPCA(rld_talon, intgroup = c("condition","line"),returnData=TRUE) 
pca_data_talon$method <- "TALON"

pca_data <- rbind(pca_data_flair,pca_data_iso,pca_data_talon)

annotated_pca_combined <- ggplot(pca_data, aes(x=PC1, y= PC2, color=condition, shape=line, size=method))+
  geom_point(guide=guide_legend(nrow=2))+
  scale_color_viridis(discrete=TRUE, labels=labels_man)+
  scale_shape_manual(values =c(8,15,16,17,18,25))+
  geom_text(aes(label=substr(name, start = 1, stop = 10)),show.legend = FALSE,vjust=3,check_overlap = TRUE,size = 2)+
  labs(color=element_blank(),shape=element_blank(), size=element_blank()) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"), legend.position="bottom",
        legend.text = element_text(size=12),
        title=element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=12),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=12))

# DEG/DIE TALON --------
de_talon_mgl_time <- dataset_make(DEresults_time_talon)
de_talon_mgl_time$method <- "TALON"
de_talon_mgl_time %>% write.csv("de_talon_mgl_time.csv")

talon_deg_time <- de_talon_mgl_time %>% filter(diffexpressed == "Up" | diffexpressed == "Down")
talon_deg_time$method <- "TALON"
talon_deg_time$Gene <- mgl_talon_filter$associated_gene[match(talon_deg_time$delabel, talon_filter$isoform)]

talon_deg_time %>% write.csv("talon_deg_time.csv")
mgl_genes <- rownames(talon_deg_time)
talon_deg_time$Catagory <- mgl_talon_filter$structural_category[match(talon_deg_time$delabel, talon_filter$isoform)]


de_talon_sex <- dataset_make(DEresults_sex_talon)
de_talon_sex$method <- "TALON"

de_talon_mgl_sex %>% write.csv("deg_talon_sex.csv")

talon_deg_sex <- de_talon_sex %>% filter(diffexpressed == "Up" | diffexpressed == "Down")
talon_deg_sex$method <- "TALON"
talon_deg_sex$Gene <- talon_filter$associated_gene[match(talon_deg_sex$delabel, talon_filter$isoform)]
talon_deg_sex$Catagory <- talon_filter$structural_category[match(talon_deg_sex$delabel, talon_filter$isoform)]


talon_deg_sex %>% write.csv("talon_deg_sex.csv")
genes_sex <- rownames(talon_deg_sex)

talon_deg_time_heatmap <- create_heatmap(talon_deg_time, dds_talon, genes, coldata, mgl_barcodes, column_split = "NULL", time)
talon_deg_sex_heatmap <- create_heatmap(talon_deg_sex, dds_talon, genes_sex, coldata, mgl_barcodes, column_split = "NULL", sex)

# Venn diagrams --------

talon_time_down <- talon_deg_time %>% filter(diffexpressed == "Down")
isoseq_time_down <- isoseq_deg_time %>% filter(diffexpressed == "Down")
flair_time_down <- flair_deg_time %>% filter(diffexpressed == "Down")

talon_time_down$method <- "TALON"
isoseq_time_down$method <- "Isoseq"
flair_time_down$method <- "FLAIR"

down_time <- rbind(talon_time_down,isoseq_time_down,flair_time_down)

time_down_list <- list(
  FLAIR=flair_time_down$Gene,
  Isoseq=isoseq_time_down$Gene,
  TALON=talon_time_down$Gene
)

venn=Venn(time_down_list)
data=process_data(venn)

venn_set(data)
venn_region(data)

time_down_venn <- ggplot()+geom_polygon(aes(X,Y,fill=id,group=id),data=venn_regionedge(data), show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE)+
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(data), 
            show.legend = FALSE) +
  scale_color_viridis(discrete = TRUE)+
  geom_text(aes(X, Y, label = c("FLAIR","Isoseq","TALON")), 
            data = venn_setlabel(data), size=3) +
  geom_label(aes(X, Y, label = count), 
             data = venn_regionlabel(data),show.legend = FALSE) +
  coord_equal() +
  labs(title="Higher on Day 1")+
  theme(plot.title = element_text(size = 12))+
  theme_void()

talon_time_up <- talon_deg_time %>% filter(diffexpressed == "Up")
isoseq_time_up <- isoseq_deg_time %>% filter(diffexpressed == "Up")
flair_time_up <- flair_deg_time %>% filter(diffexpressed == "Up")

talon_time_up$method <- "TALON"
isoseq_time_up$method <- "Isoseq"
flair_time_up$method <- "FLAIR"

up_time <- rbind(talon_time_up,isoseq_time_up,flair_time_up)


time_up_list <- list(
  FLAIR=flair_time_up$Gene,
  Isoseq=isoseq_time_up$Gene,
  TALON=talon_time_up$Gene
)

venn=Venn(time_up_list)
data=process_data(venn)

venn_set(data)
venn_region(data)

time_up_venn <- ggplot()+geom_polygon(aes(X,Y,fill=id,group=id),data=venn_regionedge(data), show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE)+
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(data), 
            show.legend = FALSE) +
  scale_color_viridis(discrete = TRUE)+
  geom_text(aes(X, Y, label = c("FLAIR","Isoseq","TALON")), 
            data = venn_setlabel(data), size=3) +
  geom_label(aes(X, Y, label = count), 
             data = venn_regionlabel(data),show.legend = FALSE) +
  coord_equal() +
  labs(title="Higher on Day 14")+
  theme(plot.title = element_text(size = 12))+
  theme_void()


talon_sex_down <- talon_deg_sex %>% filter(diffexpressed == "Down")
isoseq_sex_down <- isoseq_deg_sex %>% filter(diffexpressed == "Down")
flair_sex_down <- flair_deg_sex %>% filter(diffexpressed == "Down")

talon_sex_down$method <- "TALON"
isoseq_sex_down$method <- "Isoseq"
flair_sex_down$method <- "FLAIR"

down_sex <- rbind(talon_sex_down,isoseq_sex_down,flair_sex_down)


sex_down_list <- list(
  FLAIR=flair_sex_down$Gene,
  Isoseq=isoseq_sex_down$Gene,
  TALON=talon_sex_down$Gene
)

venn=Venn(sex_down_list)
data=process_data(venn)

venn_set(data)
venn_region(data)

sex_down_venn <- ggplot()+geom_polygon(aes(X,Y,fill=id,group=id),data=venn_regionedge(data), show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE)+
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(data), 
            show.legend = FALSE) +
  scale_color_viridis(discrete = TRUE)+
  geom_text(aes(X, Y, label = c("FLAIR","Isoseq","TALON")), 
            data = venn_setlabel(data), size=3) +
  geom_label(aes(X, Y, label = count), 
             data = venn_regionlabel(data),show.legend = FALSE) +
  coord_equal() +
  labs(title="Higher in XY")+
  theme(plot.title = element_text(size = 12))+
  theme_void()

talon_sex_up <- talon_deg_sex %>% filter(diffexpressed == "Up")
isoseq_sex_up <- isoseq_deg_sex %>% filter(diffexpressed == "Up")
flair_sex_up <- flair_deg_sex %>% filter(diffexpressed == "Up")


sex_up_list <- list(
  FLAIR=flair_sex_up$Gene,
  Isoseq=isoseq_sex_up$Gene,
  TALON=talon_sex_up$Gene
)

talon_sex_up$method <- "TALON"
isoseq_sex_up$method <- "Isoseq"
flair_sex_up$method <- "FLAIR"

up_sex <- rbind(talon_sex_up,isoseq_sex_up,flair_sex_up)

venn=Venn(sex_up_list)
data=process_data(venn)

venn_set(data)
venn_region(data)

sex_up_venn <- ggplot()+geom_polygon(aes(X,Y,fill=id,group=id),data=venn_regionedge(data), show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE)+
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(data), 
            show.legend = FALSE) +
  scale_color_viridis(discrete = TRUE)+
  geom_text(aes(X, Y, label = c("FLAIR","Isoseq","TALON")), 
            data = venn_setlabel(data), size=3) +
  geom_label(aes(X, Y, label = count), 
             data = venn_regionlabel(data),show.legend = FALSE) +
  coord_equal() +
  labs(title="Higher in XX")+
  theme(plot.title = element_text(size = 12))+
  theme_void()

# Figures------------
deg_time <- rbind(de_flair_time,de_isoseq_time,de_talon_time)
deg_sex <- rbind(de_flair_sex,de_isoseq_sex,de_talon_sex)


colors_manual <- c("Higher in XY"="skyblue","Neither"="#5f5f5f","Higher in XX"="#f1605dff")
colors_time <- c("Higher in Day 1"="#677f51","Neither"="#5f5f5f","Higher in Day 14"="#721f81ff") # adapt as needed

deg_time$diffexpressed[deg_time$diffexpressed=="Down"] <- "Higher in Day 1"  # adapt as needed
deg_time$diffexpressed[deg_time$diffexpressed=="Up"] <- "Higher in Day 14"  # adapt as needed

deg_sex$diffexpressed[deg_sex$diffexpressed=="Down"] <- "Higher in XY"
deg_sex$diffexpressed[deg_sex$diffexpressed=="Up"] <- "Higher in XX"


deg_time <- ggplot(data = deg_time, 
       aes(x = log2FoldChange, 
           y = -log10(padj), 
           col = diffexpressed, 
           label = delabel,shape=method)) +
  geom_point() +
  geom_text_repel(aes(label = genelabels), show.legend = FALSE) +
  scale_color_manual(values = c("Higher in Day 1"="#677f51","Neither"="#5f5f5f","Higher in Day 14"="#721f81ff")) +  # adapt as needed
  geom_vline(xintercept = c(-1.5, 1.5), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  xlab(expression(Log[2](Fold~Change)))+
       ylab(expression(-Log[10](adjusted~p~value)))+
  labs(shape=element_blank(),color=element_blank()) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        panel.grid = element_blank(), 
        axis.line = element_line(color = "#5f5f5f"),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.text.x  = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        axis.text.y  = element_text(vjust = 0.5, size = 12),
        legend.position = "bottom")

deg_sex <- ggplot(data = deg_sex, 
                   aes(x = log2FoldChange, 
                       y = -log10(padj), 
                       col = diffexpressed, 
                       label = delabel,shape=method)) +
  geom_point() +
  geom_text_repel(aes(label = genelabels), show.legend = FALSE) +
  scale_color_manual(values = c("Higher in XY"="skyblue","Neither"="#5f5f5f","Higher in XX"="#f1605dff")) +  # adapt as needed
  geom_vline(xintercept = c(-1.5, 1.5), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  xlab(expression(Log[2](Fold~Change)))+
  ylab(expression(-Log[10](adjusted~p~value)))+
  labs(shape=element_blank(),color=element_blank()) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        panel.grid = element_blank(), 
        axis.line = element_line(color = "#5f5f5f"),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.text.x  = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        axis.text.y  = element_text(vjust = 0.5, size = 12),
        legend.position = "bottom")


deg_time_heatmaps <- ((flair_deg_time_heatmap) / (isoseq_deg_time_heatmap) / (talon_deg_time_heatmap)) + plot_layout(guides='keep') +
  plot_annotation(tag_levels = 'A') 

ggsave('deg_time_heatmaps.png',deg_time_heatmaps, width = 15, height = 20, bg="transparent")


deg_sex_heatmaps <- ((flair_deg_sex_heatmap) / (isoseq_deg_sex_heatmap) / (talon_deg_sex_heatmap)) + plot_layout(guides='keep') +
  plot_annotation(tag_levels = 'A') 

ggsave('deg_sex_heatmaps.png',deg_sex_heatmaps, width = 15, height = 20, bg="transparent")

deg_combo <- ((time_down_venn | time_up_venn) /
                (sex_down_venn | sex_up_venn))
venn_deg <- patchwork::wrap_elements(deg_combo)

deg_plot <- ((annotated_pca_combined | venn_deg) /(deg_time |deg_sex)) + plot_layout(guides='keep') +
  plot_annotation(tag_levels = 'A')

ggsave('deg_plot.png',deg_plot, width = 20, height = 15, bg="transparent")

