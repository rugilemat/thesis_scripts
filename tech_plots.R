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

#*** Global plot parameters-----------

myPalette = c("#6BAED6","#FC8D59","#78C679","#EE6A50","#969696","#66C2A4", "goldenrod1", "darksalmon", "#41B6C4","tomato3", "#FE9929")
subcat.palette = c("Alternative 3'end"='#02314d',
                   "Alternative 3'5'end"='#0e5a87',
                   "Alterantive 5'end"='#7ccdfc',
                   'Reference match'='#c4e1f2',
                   "3' fragment"='#c4531d',
                   "Internal fragment"='#e37744',  
                   "5' fragment"='#e0936e', 
                   "Comb. of annot. junctions"='#014d02',
                   "Comb. of annot. splice sites"='#379637',  
                   "Intron retention"='#81eb82', 
                   "Not comb. of annot. junctions"='#6ec091',
                   "Mono-exon by intron ret."='#4aaa72',
                   "At least 1 annot. don./accept."='#32734d',
                   "Mono-exon"='#cec2d2',
                   "Multi-exon"='#876a91')



cat.palette = c("FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                "NNC"="#EE6A50", "Genic\nGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                "Intergenic" = "darksalmon", "Genic\nIntron"="#41B6C4")

mytheme <- theme_classic(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", linewidth = 0.4),
        axis.line.y = element_line(color="black", linewidth = 0.4)) +
  theme(axis.title.x = element_text(size=20),
        axis.text.x  = element_text(size=18),
        axis.title.y = element_text(size=20),
        axis.text.y  = element_text(vjust=0.5, size=18) ) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size=20), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=22, hjust = 0.5)) +
  theme(plot.margin = unit(c(2.5,1,1,1), "cm"))


data.FSMISM <- subset(combined, structural_category %in% c('FSM', 'ISM'))
data.NICNNC <- subset(combined, structural_category %in% c("NIC", "NNC"))
data.other <- subset(combined, structural_category %in% c("Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron"))
data.FSM <- subset(combined, (structural_category=="FSM" & exons>1))
data.ISM <- subset(combined, (structural_category=="ISM" & exons>1))
data.NNC <- subset(combined, (structural_category=="NNC" & exons>1))
data.NIC <- subset(combined, (structural_category=="NIC" & exons>1))
data.GenicGenomic <- subset(combined, (structural_category=="Genic\nGenomic" & exons>1 ))
data.Antisense <- subset(combined, (structural_category=="Antisense" & exons>1))
data.Fusion <- subset(combined, (structural_category=="Fusion" & exons>1))
data.Intergenic <- subset(combined, (structural_category=="Intergenic" & exons>1))
data.GenicIntron <- subset(combined, (structural_category=="Genic\nIntron" & exons>1))



#*** Datasets---------------------
mgl_flair <- read.delim("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Experiments/Deep read seq/Targeted/MGL targeted files/SQANTI3_mgl_ctl_flair_RulesFilter_result_classification_targetgenes.txt")
mgl_flair <-mgl_flair[,-1]
mgl_flair_filter <- mgl_flair %>% filter(filter_result=="Isoform")
mgl_flair_filter$method <- "flair"
mgl_flair_filter <-  mgl_flair_filter %>%
  rename(
    "BC01" = "FL.007D1",
    "BC02" = "FL.007D14",
    "BC03" = "FL.007D1.1",
    "BC04" = "FL.007D14.1",
    "BC05" = "FL.007D1.2",
    "BC06" = "FL.007D14.2",
    "BC13" = "FL.069D1",
    "BC14" = "FL.069D14",
    "BC15" = "FL.069D1.1",
    "BC16" = "FL.069D14.1",
    "BC17" = "FL.069D1.2",
    "BC18" = "FL.069D14.2",
    "BC19" = "FL.127D1",
    "BC20" = "FL.127D14",
    "BC21" = "FL.127D1.1",
    "BC22" = "FL.127D14.1",
    "BC24" = "FL.127D14.2",
    "BC25" = "FL.014D1",
    "BC26" = "FL.014D14",
    "BC27" = "FL.014D1.1",
    "BC28" = "FL.014D14.1",
    "BC29" = "FL.014D1.2",
    "BC30" = "FL.014D14.2",
    "BC31" = "FL.36SD1",
    "BC32" = "FL.36SD14",
    "BC33" = "FL.36SD1.1",
    "BC34" = "FL.36SD14.1",
    "BC35" = "FL.36SD1.2",
    "BC36" = "FL.36SD14.2",
    "BC37" = "FL.SCTID1",
    "BC38" = "FL.SCTID14",
    "BC39" = "FL.SCTID1.1",
    "BC40" = "FL.SCTID14.1",
    "BC41" = "FL.SCTID1.2",
    "BC42" = "FL.SCTID14.2"
  )

mgl_flair_filter$count  <- rowSums(mgl_flair_filter[,c(49:83)])
mgl_flair_filter <- mgl_flair_filter %>% filter(!is.na(count))


mgl_isoseq <- read.delim("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Experiments/Deep read seq/Targeted/MGL targeted files/SQANTI3_isoseq_mgl_merged_RulesFilter_result_classification_targetgenes.txt")
mgl_isoseq <- mgl_isoseq[,-1]
mgl_isoseq_filter <-mgl_isoseq %>% filter(filter_result=="Isoform")
mgl_isoseq_filter$method <- "isoseq"
mgl_isoseq_filter <-  mgl_isoseq_filter %>%
  rename(
    "BC01" = "FL.BC01",
    "BC02" = "FL.BC02",
    "BC03" = "FL.BC03",
    "BC04" = "FL.BC04",
    "BC05" = "FL.BC05",
    "BC06" = "FL.BC06",
    "BC13" = "FL.BC13",
    "BC14" = "FL.BC14",
    "BC15" = "FL.BC15",
    "BC16" = "FL.BC16",
    "BC17" = "FL.BC17",
    "BC18" = "FL.BC18",
    "BC19" = "FL.BC19",
    "BC20" = "FL.BC20",
    "BC21" = "FL.BC21",
    "BC22" = "FL.BC22",
    "BC23" = "FL.BC23",
    "BC24" = "FL.BC24",
    "BC25" = "FL.BC25",
    "BC26" = "FL.BC26",
    "BC27" = "FL.BC27",
    "BC28" = "FL.BC28",
    "BC29" = "FL.BC29",
    "BC30" = "FL.BC30",
    "BC31" = "FL.BC31",
    "BC32" = "FL.BC32",
    "BC33" = "FL.BC33",
    "BC34" = "FL.BC34",
    "BC35" = "FL.BC35",
    "BC36" = "FL.BC36",
    "BC37" = "FL.BC37",
    "BC38" = "FL.BC38",
    "BC39" = "FL.BC39",
    "BC40" = "FL.BC40",
    "BC41" = "FL.BC41",
    "BC42" = "FL.BC42"
  )

mgl_isoseq_filter$count  <- rowSums(mgl_isoseq_filter[,c(49:84)])
mgl_isoseq_filter <- mgl_isoseq_filter %>% filter(!is.na(count))


mgl_talon <- read.delim("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Experiments/Deep read seq/Targeted/MGL targeted files/SQANTI3_mgl_ctl_talon_RulesFilter_result_classification_targetgenes.txt")
mgl_talon <- mgl_talon[,-1]
mgl_talon_filter <-mgl_talon %>% filter(filter_result=="Isoform")
mgl_talon_filter$method <- "TALON"
mgl_talon_filter <-  mgl_talon_filter %>%
  rename(
    "BC01" = "FL.BC01",
    "BC02" = "FL.BC02",
    "BC03" = "FL.BC03",
    "BC04" = "FL.BC04",
    "BC05" = "FL.BC05",
    "BC06" = "FL.BC06",
    "BC13" = "FL.BC13",
    "BC14" = "FL.BC14",
    "BC15" = "FL.BC15",
    "BC16" = "FL.BC16",
    "BC17" = "FL.BC17",
    "BC18" = "FL.BC18",
    "BC19" = "FL.BC19",
    "BC20" = "FL.BC20",
    "BC21" = "FL.BC21",
    "BC22" = "FL.BC22",
    "BC23" = "FL.BC23",
    "BC24" = "FL.BC24",
    "BC25" = "FL.BC25",
    "BC26" = "FL.BC26",
    "BC27" = "FL.BC27",
    "BC28" = "FL.BC28",
    "BC29" = "FL.BC29",
    "BC30" = "FL.BC30",
    "BC31" = "FL.BC31",
    "BC32" = "FL.BC32",
    "BC33" = "FL.BC33",
    "BC34" = "FL.BC34",
    "BC35" = "FL.BC35",
    "BC36" = "FL.BC36",
    "BC37" = "FL.BC37",
    "BC38" = "FL.BC38",
    "BC39" = "FL.BC39",
    "BC40" = "FL.BC40",
    "BC41" = "FL.BC41",
    "BC42" = "FL.BC42"
  )

mgl_talon_filter$count  <- rowSums(mgl_talon_filter[,c(49:84)])
mgl_talon_filter <- mgl_talon_filter %>% filter(!is.na(count))

target_genes_lengths <- read.delim("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Experiments/Deep read seq/target_genes_lengths.txt", header=FALSE)
colnames(target_genes_lengths) <- c("Ensembl", "Gene", "Exons", "Length")
target_genes_lengths <- target_genes_lengths[,-5]

SQANTI3_mgl <- read.delim("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Experiments/Deep read seq/Targeted/MGL targeted files/SQANTI3_mgl_tech_RulesFilter_result_classification_targetgenes.txt")
SQANTI3_mgl <- SQANTI3_mgl[,-1]
SQANTI3_mgl <-SQANTI3_mgl %>% filter(filter_result=="Isoform")
SQANTI3_mgl %>% filter(SQANTI3_mgl$associated_gene %in% target_genes_lengths$Gene)

mgl_track <- read.delim("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Experiments/Deep read seq/Targeted/MGL targeted files/mgl_tech_comp.tracking", header=FALSE)

mgl_flair_filter <- mgl_flair_filter %>% filter(mgl_flair_filter$associated_gene %in% target_genes_lengths$Gene)
mgl_isoseq_filter <- mgl_isoseq_filter %>% filter(mgl_isoseq_filter$associated_gene %in% target_genes_lengths$Gene)
mgl_talon_filter <- mgl_talon_filter %>% filter(mgl_talon_filter$associated_gene %in% target_genes_lengths$Gene)

SQANTI3_mgl %>% filter(SQANTI3_mgl$associated_gene %in% target_genes_lengths$Gene)


# Technically combined MGLs table -------------------------------------------------------

mgl_tech <- read.delim("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Experiments/Deep read seq/Targeted/MGL targeted files/SQANTI3_mgl_tech_RulesFilter_result_classification_targetgenes.txt")
mgl_tech <- mgl_tech[,-1]
mgl_tech_filter <- mgl_tech %>% filter(filter_result=="Isoform")
mgl_tech_filter %>% filter(mgl_tech_filter$associated_gene %in% target_genes_lengths$Gene)

mgl_track <- read.delim("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Experiments/Deep read seq/Targeted/MGL targeted files/mgl_tech_comp.tracking", header=FALSE)
boolen <- mgl_track
boolen_sqanti <- data.frame(do.call(rbind, strsplit(as.character(boolen$V1),"\\|")))
boolen_class <- cbind(boolen_sqanti, boolen)
boolen_class <- boolen_class %>% filter(boolen_class$X1 %in% mgl_tech_filter$isoform)
boolen_class <- boolen_class %>% rename("isoform"=X1)

names(boolen_class)[8] <- "TALON"
names(boolen_class)[9] <- "FLAIR"
names(boolen_class)[10] <- "Isoseq"

merged_tech_mgl <- merge(boolen_class,mgl_tech_filter)
merged_tech_mgl$TALON <- sub("^[^|]*\\|", "", merged_tech_mgl$TALON)
merged_tech_mgl <- merged_tech_mgl %>% separate(TALON, into=c("TALON", NULL), sep="\\|", extra="drop")
merged_tech_mgl$FLAIR <- sub("^[^|]*\\|", "", merged_tech_mgl$FLAIR)
merged_tech_mgl <- merged_tech_mgl %>% separate(FLAIR, into=c("FLAIR", NULL), sep="\\|", extra="drop")
merged_tech_mgl$Isoseq <- sub("^[^|]*\\|", "", merged_tech_mgl$Isoseq)
merged_tech_mgl <- merged_tech_mgl %>% separate(Isoseq, into=c("Isoseq", NULL), sep="\\|", extra="drop")

merged_tech_mgl$tool <- ifelse(merged_tech_mgl$TALON !="-" & merged_tech_mgl$FLAIR !="-" & merged_tech_mgl$Isoseq !="-", "Common",
                               ifelse(merged_tech_mgl$TALON !="-" & merged_tech_mgl$FLAIR =="-" & merged_tech_mgl$Isoseq =="-", "TALON",
                                      ifelse(merged_tech_mgl$TALON =="-" & merged_tech_mgl$FLAIR !="-" & merged_tech_mgl$Isoseq =="-", "FLAIR",
                                             ifelse(merged_tech_mgl$TALON =="-" & merged_tech_mgl$FLAIR =="-" & merged_tech_mgl$Isoseq !="-", "Isoseq", "Other"))))

select_mgl_flair <- mgl_flair_filter[,c(1:48,84:85)]
select_mgl_isoseq <- mgl_isoseq_filter[,c(1:48,86:87)]
select_mgl_talon <- mgl_talon_filter[,c(1:48,85:86)]

combined_mgl <- rbind(select_mgl_flair, select_mgl_isoseq, select_mgl_talon)


combined_targets_mgl <- combined_mgl %>% filter(combined_mgl$associated_gene %in% target_genes_lengths$Gene)

# Iso per gene ----------------------------------------------------------
# Make "isoPerGene" which is aggregated information by gene
#  $associatedGene - either the ref gene name or novelGene_<index>
#  $novelGene      - either "Novel Genes" or "Annotated Genes"
#  $FSM_class      - "A", "B", or "C"
#  $geneExp        - gene expression info
#  $nIso           - number of isoforms associated with this gene
#  $nIsoCat        - splicing complexity based on number of isoforms

if (!all(is.na(mgl_flair_filter$gene_exp))){
  isoPerGene_mgl = aggregate(mgl_flair_filter$associated_transcript,
                         by = list("associatedGene" = mgl_flair_filter$associated_gene,
                                   "novelGene" = mgl_flair_filter$novelGene,
                                   "FSM_class" = mgl_flair_filter$FSM_class,
                                   "geneExp"=mgl_flair_filter$gene_exp),
                         length)
} else {
  isoPerGene_mgl = aggregate(mgl_flair_filter$associated_transcript,
                         by = list("associatedGene" = mgl_flair_filter$associated_gene,
                                   "FSM_class" = mgl_flair_filter$FSM_class),
                         length)
}
# assign the last column with the colname "nIso" (number of isoforms)
colnames(isoPerGene_mgl)[ncol(isoPerGene_mgl)] <- "nIso"
isoPerGene_mgl$method <- "FLAIR"

if (!all(is.na(mgl_isoseq_filter$gene_exp))){
  isoPerGene_iso_mgl = aggregate(mgl_isoseq_filter$associated_transcript,
                             by = list("associatedGene" = mgl_isoseq_filter$associated_gene,
                                       "novelGene" = mgl_isoseq_filter$novelGene,
                                       "FSM_class" = mgl_isoseq_filter$FSM_class,
                                       "geneExp"=mgl_isoseq_filter$gene_exp),
                             length)
} else {
  isoPerGene_iso_mgl = aggregate(mgl_isoseq_filter$associated_transcript,
                             by = list("associatedGene" = mgl_isoseq_filter$associated_gene,
                                       "FSM_class" = mgl_isoseq_filter$FSM_class),
                             length)
}
# assign the last column with the colname "nIso" (number of isoforms)

colnames(isoPerGene_iso_mgl)[ncol(isoPerGene_iso_mgl)] <- "nIso"
isoPerGene_iso_mgl$method <- "Isoseq"

if (!all(is.na(mgl_talon_filter$gene_exp))){
  isoPerGene_talon_mgl = aggregate(mgl_talon_filter$associated_transcript,
                               by = list("associatedGene" = mgl_talon_filter$associated_gene,
                                         "novelGene" = mgl_talon_filter$novelGene,
                                         "FSM_class" = mgl_talon_filter$FSM_class,
                                         "geneExp"=mgl_talon_filter$gene_exp),
                               length)
} else {
  isoPerGene_talon_mgl = aggregate(mgl_talon_filter$associated_transcript,
                               by = list("associatedGene" = mgl_talon_filter$associated_gene,
                                         "FSM_class" = mgl_talon_filter$FSM_class),
                               length)
}
# assign the last column with the colname "nIso" (number of isoforms)

colnames(isoPerGene_talon_mgl)[ncol(isoPerGene_talon_mgl)] <- "nIso"
isoPerGene_talon_mgl$method <- "TALON"

isogene_mgl <- rbind(isoPerGene_mgl,isoPerGene_iso_mgl,isoPerGene_talon_mgl)

isogene_mgl$exons <- target_genes_lengths$Exons[match(isogene_mgl$associatedGene,target_genes_lengths$Gene)]
isogene_mgl$length <- target_genes_lengths$Length[match(isogene_mgl$associatedGene,target_genes_lengths$Gene)]
isogene_mgl_filter <- isogene_mgl %>%
  filter(isogene_mgl$associatedGene %in% target_genes_lengths$Gene)

#Correlation plots -----------------------------------
length_corr_plot <- ggscatter(isogene_mgl_filter, x = "nIso", y = "length", color="method",
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

exon_corr_plot <- ggscatter(isogene_mgl_filter, x = "nIso", y = "exons", color="method",
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
max_iso_per_gene_mgl <- max(as.numeric(isogene_mgl_filter$nIso))
if (max_iso_per_gene_mgl > 100) {
  isogene_mgl_filter$nIsoCat <- cut(isogene_mgl_filter$nIso, breaks = c(0,1,3,5,10,50,100, max_iso_per_gene_mgl+1), labels = c("1", "2-3", "4-5", "6-10", "10-50","50-100",">100"));
} else if (max_iso_per_gene_mgl >= 100) {
  isogene_mgl_filter$nIsoCat <- cut(isogene_mgl_filter$nIso, breaks = c(0,1,3,5,10,50,100), labels = c("1", "2-3", "4-5", "6-10", "10-50","51-100"));
} else if (max_iso_per_gene_mgl >= 50) {
  isogene_mgl_filter$nIsoCat <- cut(isogene_mgl_filter$nIso, breaks = c(0,1,3,5,10,50), labels = c("1", "2-3", "4-5", "6-10", "11-50"));
} else if (max_iso_per_gene_mgl >= 10) {
  isogene_mgl_filter$nIsoCat <- cut(isogene_mgl_filter$nIso, breaks = c(0,1,3,5,10), labels = c("1", "2-3", "4-5", "6-10"));
} else if (max_iso_per_gene_mgl >= 5) {
  isogene_mgl_filter$nIsoCat <- cut(isogene_mgl_filter$nIso, breaks = c(0,1,3,5), labels = c("1", "2-3", "4-5"));
} else if (max_iso_per_gene_mgl >= 3) {
  isogene_mgl_filter$nIsoCat <- cut(isogene_mgl_filter$nIso, breaks = c(0,1,3), labels = c("1", "2-3"));
} else {
  isogene_mgl_filter$nIsoCat <- cut(isogene_mgl_filter$nIso, breaks = c(0,1), labels = c("1"));
}

colors_manual <- c("FLAIR" = "skyblue", "Isoseq" = "#990099", "TALON" = "#677f51")
colors_man <- c("flair" = "skyblue", "isoseq" = "#990099", "TALON" = "#677f51")

labels_man <- c("flair"="FLAIR","isoseq"="Isoseq","TALON"="TALON")

iso_per_gene_method_mgl <- ggplot(isogene_mgl_filter, aes(fill=method, x=nIsoCat)) + 
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

iso_per_gene_counts <- ggplot(combined_targets_mgl, aes(fill=structural_category, x=associated_gene)) + 
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

category_proportion <- ggplot(combined_targets_mgl, aes(fill=structural_category, x=associated_gene)) + 
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
datasets = colnames(merged_tech_mgl)[8:10]

merged_tech_mgl$FLAIR <- ifelse(merged_tech_mgl$FLAIR != "-", 1, merged_tech_mgl$FLAIR)
merged_tech_mgl$Isoseq <- ifelse(merged_tech_mgl$Isoseq != "-", 1, merged_tech_mgl$Isoseq)
merged_tech_mgl$TALON <- ifelse(merged_tech_mgl$TALON != "-", 1, merged_tech_mgl$TALON)

merged_tech_mgl[datasets]=merged_tech_mgl[datasets]==1
t(head(merged_tech_mgl[datasets],3))

labels_manual <- c( "full-splice_match"="Full Splice Match", "incomplete-splice_match"="Incomplete Splice Match",
                    "novel_in_catalog"="Novel in Catalogue", "novel_not_in_catalog"="Novel not in Catalogue",
                    "fusion"="Fusion","genic"="Genic")

merged_tech_mgl <- merged_tech_mgl %>% filter(merged_tech_mgl$filter_result == "Isoform")
upset_targets_mgl <- merged_tech_mgl %>% filter(merged_tech_mgl$associated_gene %in% target_genes_lengths$Gene)


mgl_upset <- upset(upset_targets_mgl, datasets, name='Dataset', width_ratio=0.1, themes=upset_modify_themes(
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
                     (iso_per_gene_method_mgl)) / (full_upset) + plot_layout(guides='keep') +
  plot_annotation(tag_levels = 'A')

ggsave('mgl_overview.png',combined_plot, width = 20, height = 14, bg="transparent")


proportion <- (iso_per_gene_counts + plot_layout(guides="collect") &
                 theme(legend.position = "bottom")) +
  (category_proportion +theme(legend.position = "none"))+
  plot_annotation(tag_levels = 'A')


ggsave('mgl_categories_counts.png',iso_per_gene_counts, width = 20, height = 29, bg="transparent")
ggsave('mgl_categories_proportion.png',category_proportion, width = 22, height = 31, bg="transparent")


#Upsets function for conditions  -----------------------------------

generate_upset_plot <- function(dataset, plot_name) {
  
  # Calculate row sums for each condition
  dataset$D1_F <- rowSums(dataset[, c("BC01", "BC03", "BC05", "BC13", "BC15", "BC17", "BC37", "BC39", "BC41")])
  dataset$D14_F <- rowSums(dataset[, c("BC02", "BC04", "BC06", "BC14", "BC16", "BC18", "BC38", "BC40", "BC42")])
  dataset$D1_M <- rowSums(dataset[, c("BC19", "BC21", "BC25", "BC27", "BC29", "BC31", "BC33", "BC35")])
  dataset$D14_M <- rowSums(dataset[, c("BC20", "BC22", "BC24", "BC26", "BC28", "BC30", "BC32", "BC34", "BC36")])
  
  # Filter targets based on associated genes
  targets <- dataset %>% filter(associated_gene %in% target_genes_lengths$Gene)
  
  # Define column names for datasets
  dataset_cond <- c("D1_F","D14_F","D1_M","D14_M")
  
  # Convert counts to binary
  targets$D1_F <- ifelse(targets$D1_F != "0", 1, targets$D1_F)
  targets$D14_F <- ifelse(targets$D14_F != "0", 1, targets$D14_F)
  targets$D1_M <- ifelse(targets$D1_M != "0", 1, targets$D1_M)
  targets$D14_M <- ifelse(targets$D14_M != "0", 1, targets$D14_M)
  
  # Convert dataset columns to binary
  targets[dataset_cond] <- targets[dataset_cond] == 1
  
  # Define manual labels
  labels_manual <- c("full-splice_match" = "Full Splice Match", "incomplete-splice_match" = "Incomplete Splice Match",
                     "novel_in_catalog" = "Novel in Catalogue", "novel_not_in_catalog" = "Novel not in Catalogue",
                     "fusion" = "Fusion", "genic" = "Genic")
  
  # Filter targets again
  upset_data <- targets %>% filter(targets$associated_gene %in% target_genes_lengths$Gene)
  
  # Define labels for UpSet plot
  labels <- c("D1_F" = "D1 XX", "D14_F" = "D14 XX", "D1_M" = "D1 XY", "D14_M" = "D14 XY")
  
  # Create UpSet plot
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

plot_sankey(mgl_flair_filter, "flair_sankey_plot")
generate_upset_plot(mgl_flair_filter,"flair_upset")


plot_sankey(mgl_isoseq_filter,"isoseq_sankey_plot")
generate_upset_plot(mgl_isoseq_filter,"isoseq_upset")

plot_sankey(mgl_talon_filter,"talon_sankey_plot")
generate_upset_plot(mgl_talon_filter,"talon_upset")


sank_plots <- (flair_sankey_plot /isoseq_sankey_plot/talon_sankey_plot) +
  plot_annotation(tag_levels = 'A')

ggsave('mgl_sank.png',sank_plots, width = 20, height = 29, bg="transparent")

upset_plots <- (flair_upset /isoseq_upset /talon_upset) +
  plot_annotation(tag_levels = 'A') + plot_layout(guides="collect") &
theme(plot.tag = element_text(size = 20), legend.position = "bottom")

ggsave('mgl_upset.png',upset_plots, width = 20, height = 29, bg="transparent")
