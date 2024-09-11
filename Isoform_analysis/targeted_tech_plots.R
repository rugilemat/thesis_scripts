#*** Packages---------------------
library(ggplot2)
library(ggplotify)
library(scales)
library(reshape)
library(grid)
library(gridExtra)
library(NOISeq)
library(rmarkdown)
library(htmltools)
library(DT)
library(plyr)
library(plotly)
library(dplyr)
library(ggpubr)
library(viridis)
library(patchwork)
library(Biostrings)
library(rtracklayer)
library(IRanges)
library(ComplexUpset)
library(ComplexHeatmap)
library(tidyverse)
library(here)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggcorrplot)
library(ggsankey)
library(viridis)
library(GGally)
library(ggjoy)
library(ggdist)

#General overview plots for targeted sequencing data

#*** Datasets---------------------
#Import datasets

flair <- read.delim("~/path/to/flair_sqanti_classification_file")
flair <-mgl_flair[,-1]
flair_filter <- flair %>% filter(filter_result=="Isoform")
flair_filter$method <- "flair"

flair_filter$count  <- rowSums(mgl_flair_filter[,c(49:84)]) #select relevant columns with count data
flair_filter <- mgl_flair_filter %>% filter(!is.na(count))


isoseq <- read.delim("~/path/to/isoseq_sqanti_classification_file")
isoseq <- isoseq[,-1]
isoseq_filter <- isoseq %>% filter(filter_result=="Isoform")
isoseq_filter$method <- "isoseq"

isoseq_filter$count  <- rowSums(isoseq_filter[,c(49:84)]) #select relevant columns with count data
isoseq_filter <- isoseq_filter %>% filter(!is.na(count))


talon <- read.delim("~/path/to/talon_sqanti_classification_file")
talon <- talon[,-1]
talon_filter <- talon %>% filter(filter_result=="Isoform")
talon_filter$method <- "TALON"

talon_filter$count  <- rowSums(talon_filter[,c(49:84)])
talon_filter <- talon_filter %>% filter(!is.na(count))

target_genes_lengths <- read.delim("~/path/to/targeted_panel_data", header=FALSE)
colnames(target_genes_lengths) <- c("Ensembl", "Gene", "Exons", "Length")
target_genes_lengths <- target_genes_lengths[,-5]

flair_filter <- flair_filter %>% filter(flair_filter$associated_gene %in% target_genes_lengths$Gene)
isoseq_filter <- isoseq_filter %>% filter(isoseq_filter$associated_gene %in% target_genes_lengths$Gene)
talon_filter <- talon_filter %>% filter(talon_filter$associated_gene %in% target_genes_lengths$Gene)

# Technically combined table -------------------------------------------------------

tech <- read.delim("~/path/to/gffcompare_sqanti_classification_file")
tech <- tech[,-1]
tech_filter <- tech %>% filter(filter_result=="Isoform")
tech_filter %>% filter(tech_filter$associated_gene %in% target_genes_lengths$Gene)

track <- read.delim("~/path/to/gffcompare_tracking_file", header=FALSE)
boolen <- mgl_track
boolen_sqanti <- data.frame(do.call(rbind, strsplit(as.character(boolen$V1),"\\|")))
boolen_class <- cbind(boolen_sqanti, boolen)
boolen_class <- boolen_class %>% filter(boolen_class$X1 %in% tech_filter$isoform)
boolen_class <- boolen_class %>% rename("isoform"=X1)

names(boolen_class)[8] <- "TALON"
names(boolen_class)[9] <- "FLAIR"
names(boolen_class)[10] <- "Isoseq"

merged_tech <- merge(boolen_class,merged_tech)
merged_tech$TALON <- sub("^[^|]*\\|", "", merged_tech$TALON)
merged_tech <- merged_tech_mgl %>% separate(TALON, into=c("TALON", NULL), sep="\\|", extra="drop")
merged_tech$FLAIR <- sub("^[^|]*\\|", "", merged_tech$FLAIR)
merged_tech <- merged_tech %>% separate(FLAIR, into=c("FLAIR", NULL), sep="\\|", extra="drop")
merged_tech$Isoseq <- sub("^[^|]*\\|", "", merged_tech$Isoseq)
merged_tech <- merged_tech %>% separate(Isoseq, into=c("Isoseq", NULL), sep="\\|", extra="drop")

combined_mgl <- rbind(flair_filter, isoseq_filter, talon_filter)
combined_targets <- combined %>% filter(combined$associated_gene %in% target_genes_lengths$Gene)

# Iso per gene ----------------------------------------------------------
# adapted from SQANTI3 report

if (!all(is.na(flair_filter$gene_exp))){
  isoPerGene = aggregate(flair_filter$associated_transcript,
                         by = list("associatedGene" = flair_filter$associated_gene,
                                   "novelGene" = flair_filter$novelGene,
                                   "FSM_class" = flair_filter$FSM_class,
                                   "geneExp"=flair_filter$gene_exp),
                         length)
} else {
  isoPerGene = aggregate(flair_filter$associated_transcript,
                         by = list("associatedGene" = flair_filter$associated_gene,
                                   "FSM_class" = flair_filter$FSM_class),
                         length)
}
# assign the last column with the colname "nIso" (number of isoforms)
colnames(isoPerGene)[ncol(isoPerGene)] <- "nIso"
isoPerGene$method <- "FLAIR"

if (!all(is.na(isoseq_filter$gene_exp))){
  isoPerGene_iso = aggregate(isoseq_filter$associated_transcript,
                             by = list("associatedGene" = isoseq_filter$associated_gene,
                                       "novelGene" = isoseq_filter$novelGene,
                                       "FSM_class" = isoseq_filter$FSM_class,
                                       "geneExp"=isoseq_filter$gene_exp),
                             length)
} else {
  isoPerGene_iso = aggregate(isoseq_filter$associated_transcript,
                             by = list("associatedGene" = isoseq_filter$associated_gene,
                                       "FSM_class" = isoseq_filter$FSM_class),
                             length)
}
# assign the last column with the colname "nIso" (number of isoforms)

colnames(isoPerGene_iso)[ncol(isoPerGene_iso)] <- "nIso"
isoPerGene_iso$method <- "Isoseq"

if (!all(is.na(talon_filter$gene_exp))){
  isoPerGene_talon = aggregate(talon_filter$associated_transcript,
                               by = list("associatedGene" = talon_filter$associated_gene,
                                         "novelGene" = talon_filter$novelGene,
                                         "FSM_class" = talon_filter$FSM_class,
                                         "geneExp"=talon_filter$gene_exp),
                               length)
} else {
  isoPerGene_talon = aggregate(talon_filter$associated_transcript,
                               by = list("associatedGene" = talon_filter$associated_gene,
                                         "FSM_class" = talon_filter$FSM_class),
                               length)
}
# assign the last column with the colname "nIso" (number of isoforms)

colnames(isoPerGene_talon)[ncol(isoPerGene_talon)] <- "nIso"
isoPerGene_talon$method <- "TALON"

isogene <- rbind(isoPerGene,isoPerGene_iso,isoPerGene_talon)

isogene$exons <- target_genes_lengths$Exons[match(isogene$associatedGene,target_genes_lengths$Gene)]
isogene$length <- target_genes_lengths$Length[match(isogene$associatedGene,target_genes_lengths$Gene)]
isogene_filter <- isogene %>%
  filter(isogene$associatedGene %in% target_genes_lengths$Gene)

#Correlation plots -----------------------------------
length_corr_plot <- ggscatter(isogene_filter, x = "nIso", y = "length", color="method",
                              add="reg.line",
                              conf.int = TRUE, palette = c("skyblue", "#990099","#677f51"), shape="method")+
  stat_cor(aes(color=method), label.x=1000, show_guide=FALSE, legend=element_blank())+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  scale_x_continuous(expand=c(0,0))+
  labs(y = "Transcript Length (bp)", x = "Number of Isoforms")+
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=16))

exon_corr_plot <- ggscatter(isogene_filter, x = "nIso", y = "exons", color="method",
                            add="reg.line",
                            conf.int = TRUE, palette = c("skyblue", "#990099","#677f51"), shape="method")+
  stat_cor(aes(color=method), label.x=1000, show_guide=FALSE, legend=element_blank())+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  scale_x_continuous(expand=c(0,0))+
  labs(y = "Number of exons", x = "Number of Isoforms")+
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=16))

#Calculate max iso number per gene in groups
max_iso_per_gene <- max(as.numeric(isogene_filter$nIso))
if (max_iso_per_gene > 100) {
  isogene_mgl_filter$nIsoCat <- cut(isogene_filter$nIso, breaks = c(0,1,3,5,10,50,100, max_iso_per_gene+1), labels = c("1", "2-3", "4-5", "6-10", "10-50","50-100",">100"));
} else if (max_iso_per_gene >= 100) {
  isogene_mgl_filter$nIsoCat <- cut(isogene_filter$nIso, breaks = c(0,1,3,5,10,50,100), labels = c("1", "2-3", "4-5", "6-10", "10-50","51-100"));
} else if (max_iso_per_gene >= 50) {
  isogene_mgl_filter$nIsoCat <- cut(isogene_filter$nIso, breaks = c(0,1,3,5,10,50), labels = c("1", "2-3", "4-5", "6-10", "11-50"));
} else if (max_iso_per_gene >= 10) {
  isogene_mgl_filter$nIsoCat <- cut(isogene_filter$nIso, breaks = c(0,1,3,5,10), labels = c("1", "2-3", "4-5", "6-10"));
} else if (max_iso_per_gene >= 5) {
  isogene_mgl_filter$nIsoCat <- cut(isogene_filter$nIso, breaks = c(0,1,3,5), labels = c("1", "2-3", "4-5"));
} else if (max_iso_per_gene >= 3) {
  isogene_filter$nIsoCat <- cut(isogene_filter$nIso, breaks = c(0,1,3), labels = c("1", "2-3"));
} else {
  isogene_filter$nIsoCat <- cut(isogene_filter$nIso, breaks = c(0,1), labels = c("1"));
}

colors_manual <- c("FLAIR" = "skyblue", "Isoseq" = "#990099", "TALON" = "#677f51")
colors_man <- c("flair" = "skyblue", "isoseq" = "#990099", "TALON" = "#677f51")

labels_man <- c("flair"="FLAIR","isoseq"="Isoseq","TALON"="TALON")

iso_per_gene_method <- ggplot(isogene_filter, aes(fill=method, x=nIsoCat)) + 
  geom_bar(stat="count", aes(y= (..count..)/sum(..count..)*100), position="dodge") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  scale_fill_manual(values=colors_manual)+
  labs(y = "Genes (%)", x = "Number of Isoforms", fill=element_blank())+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"), legend.position="bottom",
        legend.text = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=16))

labels_man <- c("flair"="FLAIR","isoseq"="Isoseq","TALON"="TALON")
labels_manual <- c( "full-splice_match"="Full Splice Match", "incomplete-splice_match"="Incomplete Splice Match",
                    "novel_in_catalog"="Novel in Catalogue", "novel_not_in_catalog"="Novel not in Catalogue",
                    "fusion"="Fusion","genic"="Genic")

iso_per_gene_counts <- ggplot(combined_targets, aes(fill=structural_category, x=associated_gene)) + 
  geom_bar(aes(fill=structural_category)) +
  coord_flip()+
  scale_fill_viridis(discrete=TRUE,option="D", labels=labels_manual)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  labs(x = "Gene", y = "Number of Isoforms", fill=element_blank())+
  facet_wrap(~method, labeller = as_labeller(labels_man), ncol=3,strip.position = "top")+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"),
        strip.text=element_text(size=24),
        panel.spacing.y = unit(0.75, "lines"),
        legend.text = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.text.x  = element_text(vjust = 0.5, hjust=1, size =24),
        axis.title.y = element_text(size=24),
        axis.text.y  = element_text(vjust=0.5, size=10),
        legend.position = "bottom")

category_proportion <- ggplot(combined_targets, aes(fill=structural_category, x=associated_gene)) + 
  geom_bar(position="fill",aes(y=(..count..)/sum(..count..)*100))+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  scale_fill_viridis(discrete = TRUE,option="D", labels=labels_manual)+
  labs(y = "Proportion", x = "Gene", fill=element_blank())+
  facet_wrap(~method, labeller = as_labeller(labels_man), ncol=3,strip.position = "top") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        panel.spacing.x = unit(1, "lines"),
        strip.text=element_text(size=24),
        panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"), 
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        axis.text.y  = element_text(vjust=0.5, size=10),
        legend.position = "bottom")



#Upset plot-----------
datasets = colnames(merged_tech)[8:10]

merged_tech$FLAIR <- ifelse(merged_tech$FLAIR != "-", 1, merged_tech$FLAIR)
merged_tech$Isoseq <- ifelse(merged_tech$Isoseq != "-", 1, merged_tech$Isoseq)
merged_tech$TALON <- ifelse(merged_tech$TALON != "-", 1, merged_tech$TALON)

merged_tech[datasets]=merged_tech[datasets]==1
t(head(merged_tech[datasets],3))

labels_manual <- c( "full-splice_match"="Full Splice Match", "incomplete-splice_match"="Incomplete Splice Match",
                    "novel_in_catalog"="Novel in Catalogue", "novel_not_in_catalog"="Novel not in Catalogue",
                    "fusion"="Fusion","genic"="Genic")

merged_tech <- merged_tech %>% filter(merged_tech$filter_result == "Isoform")
upset_targets <- merged_tech %>% filter(merged_tech$associated_gene %in% target_genes_lengths$Gene)


upset <- upset(upset_targets, datasets, name='Dataset', width_ratio=0.1, themes=upset_modify_themes(
  list('intersections_matrix'=theme(axis.text.y=element_text(size=16),
                                    axis.title.x=element_text(size=16)))),set_sizes = FALSE,
  base_annotations=list('Isoform number'=intersection_size(mapping=aes(fill=structural_category)) +
                          theme(panel.grid = element_blank(),
                                axis.text.y = element_text(size=16),
                                axis.title.y = element_text(size=16),
                                legend.title=element_text(size=16),
                                legend.text=element_text(size=16,hjust=0))+
                          scale_fill_viridis(discrete = TRUE,
                                             labels=labels_manual,
                                             guide=guide_legend(title = "Structural Category"))),guides="over")+
  theme(panel.grid = element_blank())

full_upset <- patchwork::wrap_elements(mgl_upset) 

#Saving tech plots-----------

combined_plot <- ((length_corr_plot / exon_corr_plot) |
                     (iso_per_gene_method)) / (full_upset) + plot_layout(guides='keep') +
  plot_annotation(tag_levels = 'A')

ggsave('overview.png',combined_plot, width = 20, height = 14, bg="transparent")


ggsave('categories_counts.png',iso_per_gene_counts, width = 20, height = 29, bg="transparent")
ggsave('categories_proportion.png',category_proportion, width = 22, height = 31, bg="transparent")


#Upsets function for conditions  -----------------------------------

generate_upset_plot <- function(dataset, plot_name) {
  
  # Calculate row sums for each condition and adapt columns as needed
  dataset$D1_F <- rowSums(dataset[, c("BC01", "BC03", "BC05", "BC13", "BC15", "BC17", "BC37", "BC39", "BC41")])
  dataset$D14_F <- rowSums(dataset[, c("BC02", "BC04", "BC06", "BC14", "BC16", "BC18", "BC38", "BC40", "BC42")])
  dataset$D1_M <- rowSums(dataset[, c("BC19", "BC21","BC23", "BC25", "BC27", "BC29", "BC31", "BC33", "BC35")])
  dataset$D14_M <- rowSums(dataset[, c("BC20", "BC22", "BC24", "BC26", "BC28", "BC30", "BC32", "BC34", "BC36")])
  

  targets <- dataset %>% filter(associated_gene %in% target_genes_lengths$Gene)
  

  dataset_cond <- c("D1_F","D14_F","D1_M","D14_M")   # Define column names for datasets and adapt as needed

  
  targets$D1_F <- ifelse(targets$D1_F != "0", 1, targets$D1_F)
  targets$D14_F <- ifelse(targets$D14_F != "0", 1, targets$D14_F)
  targets$D1_M <- ifelse(targets$D1_M != "0", 1, targets$D1_M)
  targets$D14_M <- ifelse(targets$D14_M != "0", 1, targets$D14_M)
  
 targets[dataset_cond] <- targets[dataset_cond] == 1
  

  labels_manual <- c("full-splice_match" = "Full Splice Match", "incomplete-splice_match" = "Incomplete Splice Match",
                     "novel_in_catalog" = "Novel in Catalogue", "novel_not_in_catalog" = "Novel not in Catalogue",
                     "fusion" = "Fusion", "genic" = "Genic")
  

  upset_data <- targets %>% filter(targets$associated_gene %in% target_genes_lengths$Gene)
  

  labels <- c("D1_F" = "D1 XX", "D14_F" = "D14 XX", "D1_M" = "D1 XY", "D14_M" = "D14 XY") #adapt labels as needed
  

  upset_plot <- upset(upset_data, dataset_cond, name = 'Dataset', width_ratio = 0.5, min_degree=1,
                      themes = upset_modify_themes(list('intersections_matrix' = theme(axis.text.y = element_text(size = 20),
                                                                                       axis.title.x = element_text(size = 20)))),
                      set_sizes = FALSE,
                      labeller = as_labeller(labels),
                      base_annotations = list('Isoform number' = intersection_size(mapping = aes(fill = structural_category)) +
                                                theme(panel.grid = element_blank(),
                                                      axis.text.y = element_text(size = 20),
                                                      axis.title.y = element_text(size = 20),
                                                      legend.title = element_text(size = 20),
                                                      legend.text = element_text(size = 20)) +
                                                scale_fill_viridis(discrete = TRUE, labels = labels_manual,
                                                                   guide = guide_legend(title = "Structural Category")))) +
    theme(panel.grid = element_blank(),legend.position = "bottom")
 
  plot <- patchwork::wrap_elements(upset_plot) 
  
   
  assign(plot_name, plot, envir = .GlobalEnv)
  
  
}

#Sankey function-----------------------------------

plot_sankey <- function(dataset, plot_name) {
  
  dataset <- dataset %>% filter(dataset$associated_gene %in% target_genes_lengths$Gene)
  
  
  dataset$coding_pred <- ifelse(dataset$coding != "coding", dataset$coding,
                                ifelse(dataset$predicted_NMD == "False", dataset$coding,
                                       ifelse(dataset$predicted_NMD == "True", "NMD", dataset$coding)
                                ))
  
  sank <- dataset %>% make_long(filter_result, structural_category, subcategory, coding_pred)
  
  sank$x <- as.character(sank$x)
  sank$x[sank$x == "filter_result"] <- "Transcript"
  sank$x[sank$x == "structural_category"] <- "Structural Category"
  sank$x[sank$x == "coding_pred"] <- "Coding prediction"
  sank$x[sank$x == "subcategory"] <- "Subcategory"
  sank$x <- factor(sank$x)
  
  sank$next_x <- as.character(sank$next_x)
  sank$next_x[sank$next_x == "filter_result"] <- "Transcript"
  sank$next_x[sank$next_x == "structural_category"] <- "Structural Category"
  sank$next_x[sank$next_x == "coding_pred"] <- "Coding prediction"
  sank$next_x[sank$next_x == "subcategory"] <- "Subcategory"
  sank$next_x <- factor(sank$next_x)
  
  sank$node <- recode(sank$node,
                      "genic"="Genic",
                      "novel_in_catalog" = "NIC",
                      "intron_retention" = "Intron Retention",
                      "incomplete-splice_match" = "ISM",
                      "reference_match" = "Reference Match",
                      "fusion" = "Fusion",
                      "novel_not_in_catalog" = "NNC",
                      "coding" = "Coding",
                      "non_coding" = "Non-coding",
                      "full-splice_match" = "FSM",
                      "combination_of_known_splicesites" = "Combination of Known Splice Sites or Junctions",
                      "combination_of_known_junctions" = "Combination of Known Splice Sites or Junctions",
                      "3prime_fragment" = "3' or 5' fragment",
                      "5prime_fragment" = "3' or 5' fragment",
                      "internal_fragment" = "Internal Fragment",
                      "alternative_5end" = "Alternative 5' or/and 3' end",
                      "alternative_3end" = "Alternative 5' or/and 3' end",
                      "alternative_3end5end" = "Alternative 5' or/and 3' end",
                      "mono-exon" = "Mono Exon",
                      "at_least_one_novel_splicesite" = "At Least One Novel Splice Site",
                      "multi-exon" = "Multi Exon")
  
  sank$next_node <- recode(sank$next_node,
                           "genic"="Genic",
                           "novel_in_catalog" = "NIC",
                           "intron_retention" = "Intron Retention",
                           "incomplete-splice_match" = "ISM",
                           "reference_match" = "Reference Match",
                           "fusion" = "Fusion",
                           "novel_not_in_catalog" = "NNC",
                           "coding" = "Coding",
                           "non_coding" = "Non-coding",
                           "full-splice_match" = "FSM",
                           "combination_of_known_splicesites" = "Combination of Known Splice Sites or Junctions",
                           "combination_of_known_junctions" = "Combination of Known Splice Sites or Junctions",
                           "3prime_fragment" = "3' or 5' fragment",
                           "5prime_fragment" = "3' or 5' fragment",
                           "internal_fragment" = "Internal Fragment",
                           "alternative_5end" = "Alternative 5' or/and 3' end",
                           "alternative_3end" = "Alternative 5' or/and 3' end",
                           "alternative_3end5end" = "Alternative 5' or/and 3' end",
                           "mono-exon" = "Mono Exon",
                           "at_least_one_novel_splicesite" = "At Least One Novel Splice Site",
                           "multi-exon" = "Multi Exon")
  
  sank$x <- factor(sank$x, levels = c("Transcript", "Structural Category", "Subcategory", "Coding prediction"))
  sank$next_x <- factor(sank$next_x, levels = c("Transcript", "Structural Category", "Subcategory", "Coding prediction"))
  
  dagg <- sank %>%
    dplyr::group_by(node) %>%
    tally()
  
  df_sankey2 <- merge(sank, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)
  
  sankey_plot <- ggplot(df_sankey2, aes(x = x, 
                                        next_x = next_x, 
                                        node = node, 
                                        next_node = next_node,
                                        fill = factor(node),
                                        label = paste(node, '\n', " n =", n))) +
    geom_sankey(flow.alpha = 0.5, node.color = 1, width = .3) +
    geom_sankey_label(size = 6, color = "white") +
    scale_fill_viridis(discrete = TRUE) +
    theme_sankey(base_size = 20) +
    theme(legend.position = "none",
          axis.title.x = element_blank())
  
  assign(plot_name, sankey_plot, envir = .GlobalEnv)
}


#Sankey and condition upset plots-----------------------------------

plot_sankey(flair_filter, "flair_sankey_plot")
generate_upset_plot(flair_filter,"flair_upset")


plot_sankey(isoseq_filter,"isoseq_sankey_plot")
generate_upset_plot(isoseq_filter,"isoseq_upset")

plot_sankey(talon_filter,"talon_sankey_plot")
generate_upset_plot(talon_filter,"talon_upset")


sank_plots <- (flair_sankey_plot /isoseq_sankey_plot/talon_sankey_plot) +
  plot_annotation(tag_levels = 'A')

ggsave('sank.png',sank_plots, width = 20, height = 29, bg="transparent")

upset_plots <- (flair_upset /isoseq_upset /talon_upset) +
  plot_annotation(tag_levels = 'A') + plot_layout(guides="collect") &
theme(plot.tag = element_text(size = 20), legend.position = "bottom")

ggsave('upset.png',upset_plots, width = 20, height = 29, bg="transparent")
