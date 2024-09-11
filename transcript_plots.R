library(ggtranscript)

gff <- "~/path/to/annotation_gff"

gff <- rtracklayer::import(gff) %>%
  as.data.frame() %>%
  tibble::as_tibble() 

gff <- gff %>% separate("gene_id", into=c("gene_id", NULL), sep="\\.", extra="drop")

gff <- gff %>%
  dplyr::filter(transcript_id %in% tech_filter$isoform)  

gff$gene_name <-tech_filter$associated_gene[match(gff$transcript_id, tech_filter$isoform)]


gene <- gff %>% filter(gene_id=="gene_name")
exons <- gene %>% dplyr::filter(type == "exon")
cds <- gene %>% dplyr::filter(type == "CDS")


labels_manual <- c( "full-splice_match"="Full Splice Match", "incomplete-splice_match"="Incomplete Splice Match",
                    "novel_in_catalog"="Novel in Catalogue", "novel_not_in_catalog"="Novel not in Catalogue",
                    "fusion"="Fusion","genic"="Genic")


colors_manual <- c( "full-splice_match"="#430154", "incomplete-splice_match"="#20908c",
                    "novel_in_catalog"="#5ec763", "novel_not_in_catalog"="#fde724",
                    "genic"="#3b528b")


transcripts <- ggplot(exons, aes(
  xstart = start,
  xend = end,
  y = transcript_id
)) +
  geom_range(aes(fill = structural_category, color=structural_category),
    height = 0.25
  ) +
  geom_range(
    data = cds, 
    aes(fill = structural_category, color=structural_category))+
  geom_intron(
    data = to_intron(exons, "transcript_id"),
    aes(strand = strand),
    arrow.min.intron.length = 2000) +
  scale_fill_manual(values=colors_manual,labels=labels_manual)+
  scale_color_manual(values=colors_manual,labels=labels_manual)+
  labs(x = "Genome Position", y=element_blank())+
  ggtitle(substitute(paste(italic("gene name"))))+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(),legend.position = "bottom",
        legend.title = element_blank(),
        panel.spacing.y = unit(0.75, "lines"),
        legend.text = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x  = element_text(vjust = 0.5, hjust=1, size =16),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=16, face="italic"),
        plot.title = element_text(size=16))



ggsave('transcripts.png',transcript, width = 11, height = 15, bg="transparent")

