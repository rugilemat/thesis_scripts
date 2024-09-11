library(tidyverse)
library(viridis)
library(gprofiler2)
library(disgenetplus2r)
library(dplyr)
library(disgenet2r)
library(enrichplot)
library(DOSE)
library(org.Hs.eg.db)
library(clusterProfiler)
library("ggupset")
library(AnnotationHub)
library("MeSHDbi")
library(meshes)


sfari_db <- read.csv("~/path/to/data")
sz_2022 <- read.csv("~/path/to/data")
eating_disorder <- read.csv("~/path/to/data")
bipolar <- read.csv("~/path/to/data")
adhd_2022 <- read.csv("~/path/to/data")
adhd_2019 <- read.csv("~/path/to/data")
substance_use <- read.csv("~/path/to/data")
depression <- read.csv("~/path/to/data")
target_genes_lengths <- read.delim("~/path/to/data", header=FALSE)
colnames(target_genes_lengths) <- c("Ensembl", "Gene", "Exons", "Length")
target_genes_lengths <- target_genes_lengths[,-5]
gProfiler <- read.csv("~/path/to/data")


sfari_targets <- sfari_db %>% filter(sfari_db$ensembl.id %in% target_genes_lengths$Ensembl)
sfari_targets$Condition <- "ASD"
sz_target <- sz_2022 %>% filter(sz_2022$Ensembl.ID %in% target_genes_lengths$Ensembl)
sz_target$Condition <- "SZ"
eating_target <- eating_disorder %>% filter(eating_disorder$GENENAME_GM %in% target_genes_lengths$Gene)
eating_target$Condition <- "ED"
bipolar_target <- bipolar %>% filter(bipolar$Ensembl.ID %in% target_genes_lengths$Ensembl)
bipolar_target$Condition <- "BP"
adhd <- adhd_2022 %>% filter(adhd_2022$GENE %in% target_genes_lengths$Ensembl)
adhd$Condition <- "ADHD"
adhd_2 <- adhd_2019 %>% filter(adhd_2019$gene_symbol %in% target_genes_lengths$Gene)
adhd_2$Condition <- "ADHD"
substance_target <- substance_use %>% filter(substance_use$GENE %in% target_genes_lengths$Ensembl)
substance_target$Condition <- "SUD"
depression_target <- depression %>% filter(depression$ENSG.ID %in% target_genes_lengths$Ensembl)
depression_target$Condition <- "Depression"

diseases <- target_genes_lengths
diseases$ASD <- sfari_targets$gene.symbol[match(diseases$Ensemb,sfari_targets$ensembl.id)]
diseases$SZ <- sz_target$Symbol.ID[match(diseases$Ensemb,sz_target$Ensembl.ID)]
diseases$ED <- eating_target$GENENAME_GM[match(diseases$Gene,eating_target$GENENAME_GM)]
diseases$ADHD <- adhd_2$gene_symbol[match(diseases$Gene,adhd_2$gene_symbol)]
diseases[216,8] <- "DCC"
diseases$SUD <- substance_target$SYMBOL[match(diseases$Ensemb,substance_target$GENE)]
diseases$Depression <- depression_target$Gene.symbol[match(diseases$Ensemb,depression_target$ENSG.ID)]

datasets = colnames(diseases)[5:10]

diseases$ASD <- ifelse(diseases$ASD != "NA", 1, 0)
diseases$SZ <- ifelse(diseases$SZ != "NA", 1, 0)
diseases$ED <- ifelse(diseases$ED != "NA", 1, 0)
diseases$ADHD <- ifelse(diseases$ADHD != "NA", 1, 0)
diseases$SUD <- ifelse(diseases$SUD != "NA", 1, 0)
diseases$Depression <- ifelse(diseases$Depression != "NA", 1, 0)
diseases <- replace_na(diseases, list(ASD=0,SZ=0,ED=0,ADHD=0,AD=0,SUD=0,Depression=0))

diseases[datasets]=diseases[datasets]==1
t(head(diseases[datasets],3))


gene_upset <- upset(diseases, datasets, name='Condition', width_ratio=0.1, themes=upset_modify_themes(
  list('intersections_matrix'=theme(axis.text.y=element_text(size=12),
                                    axis.title.x=element_text(size=12)))),set_sizes =FALSE,
  queries=list(upset_query(intersect="ADHD",color="#450C54",fill="#450C54"),
               upset_query(intersect="ASD",color="#31678D",fill="#31678D"),
               upset_query(intersect="SZ",color="#35B879",fill="#35B879"),
               upset_query(intersect="Depression",color="#FDE724",fill="#FDE724")),
 base_annotations=list('# of genes'=intersection_size() +
                          theme(panel.grid = element_blank(),
                                axis.text.y = element_text(size=12),
                                axis.title.y = element_text(size=12),
                                legend.title=element_text(size=12),
                                legend.text=element_text(size=12,hjust=0))))+
  theme(panel.grid = element_blank())

gene_upset <- patchwork::wrap_elements(gene_upset) 

labels_man <- c("GO:BP"="GO:Biological Process","GO:MF"="GO:Molecular Functions","GO:CC"="GO:Cellular Component")

GO_terms <- ggplot(gProfiler) + aes(x=negative_log10_of_adjusted_p_value, y=str_wrap(term_name), size=intersection_size, color=adjusted_p_value) + geom_point() + 
  scale_color_viridis()+
  labs(y =element_blank(), x = "-log10(p adjust)", color="p adjust",shape="Source",size="Intersection size")+
  facet_wrap(~source, labeller = as_labeller(labels_man), ncol = 3)+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 12),
        strip.text=element_text(size=12),
        panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.y  = element_text(vjust=0.5, size=12))

gene <- target_genes_lengths$Ensembl


ah <- AnnotationHub(localHub=TRUE)
hsa <- query(ah, c("MeSHDb", "Homo sapiens"))
file_hsa <- hsa[[1]]
db <- MeSHDbi::MeSHDb(file_hsa)

y <- enrichMeSH(eg$ENTREZID, MeSHDb = db, database = 'gene2pubmed', category = "C")
yx <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')
disease_network <- cnetplot(yx) + labs(size ="Size")


eg = bitr(target_genes_lengths$Ensembl, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x <- enrichDO(gene          = eg$ENTREZID,
              ont           = "DO",
              readable      = TRUE)
do_enrich <- cnetplot(x)

condition_heatmap <- heatplot(x) + 
  scale_y_discrete(expand = c(0,0))+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"),
        axis.title.x = element_text(size=12),
        axis.text.x  = element_text(size =6, angle=90),
        axis.title.y = element_text(size=12),
        axis.text.y  = element_text(vjust=0.5, size=12))


target_plot <- ((gene_upset | condition_heatmap) /( disease_network | GO_terms)) + plot_layout(guides='keep') +
  plot_annotation(tag_levels = 'A')
ggsave('targeted_overview.png',target_plot, width = 25, height = 18, bg="transparent")

