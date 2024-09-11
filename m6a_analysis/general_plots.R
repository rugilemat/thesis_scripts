# Packages ----------------

library(seqLogo)
library(Biostrings)
library("plotgardener")
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("AnnotationHub")
library("ggplot2")
library("cowplot")
library("dplyr")
library(ggpubr)
library(purrr)
library(patchwork)
library("flextable")
library(ggridges)
library(ggimage)
library("magick")
library(png)
library(ggpubr)
library(cowplot)
library(ggseqlogo)
library(patchwork)
library(ggVennDiagram)
library(ComplexUpset)
library(vroom)

# Datasets ----------------
all_sites_unfiltered <- vroom("~/path/to/file")
annotated_modified <- vroom("~/path/to/file")

all_sites_unfiltered <- separate(all_sites_unfiltered, col=transcript_id, into=c("transcript_id", "gene_id", NULL), sep="\\|", extra="drop")
all_sites_unfiltered$gene_id <- sub("^[^|]*\\|", "", all_sites_unfiltered$gene_id)
all_sites_unfiltered <- separate(all_sites_unfiltered, col=transcript_id, into=c("transcript_id", NULL), sep="\\.", extra="drop")
all_sites_unfiltered <- separate(all_sites_unfiltered, col=gene_id, into=c("gene_id", NULL), sep="\\.", extra="drop")
gff_gencode <- separate(gff_gencode, col=gene_id, into=c("gene_id", NULL), sep="\\.", extra="drop")
gff_gencode <- separate(gff_gencode, col=transcript_id, into=c("transcript_id", NULL), sep="\\.", extra="drop")


all_sites_unfiltered$modification_status <- ifelse(all_sites_unfiltered$probability_modified < 0.3 & !(all_sites_unfiltered$transcript_id %in% annotated_modified$transcript_id), "unmodified",
                                                   ifelse(all_sites_unfiltered$probability_modified >=0.9,"modified","neither"))

gff_transcripts <- gff_gencode %>% filter(type=="transcript")
all_sites_unfiltered$length <- gff_transcripts$width[match(all_sites_unfiltered$transcript_id,gff_transcripts$transcript_id)]

unmodified <- unmodified %>% filter(!(transcript_id %in% m6aPerTranscript$transcript_id))

unmodified <- separate(unmodified, col=transcript_id, into=c("transcript_id", "gene_id", NULL), sep="\\.", extra="drop")
unmodified$gene_id <- sub("^[^|]*\\|", "", unmodified$gene_id)
unique_modified[unique_modified == "F_D1"] <- "D1_F"
unique_modified[unique_modified == "F_D14"] <- "D14_F"


annotated_modified[annotated_modified == "F_D1"] <- "D1_F"
annotated_modified[annotated_modified == "F_D14"] <- "D14_F"

polyA_list_filtered[polyA_list_filtered == "F_D1"] <- "D1_F"
polyA_list_filtered[polyA_list_filtered == "F_D14"] <- "D14_F"
polyA_list_filtered[polyA_list_filtered == "M_D1"] <- "D1_M"
polyA_list_filtered[polyA_list_filtered == "M_D14"] <- "D14_M"



annotated_modified$line <- ifelse(grepl("SCTI", annotated_modified$sample_id), "SCTI",
                                  ifelse(grepl("007", annotated_modified$sample_id),"CTF007",
                                         ifelse(grepl("069",annotated_modified$sample_id),"CTF069",
                                                ifelse(grepl("127", annotated_modified$sample_id),"CTM127",
                                                       ifelse(grepl("014",annotated_modified$sample_id),"CTM014","M3_36S")))))


polyA_list_filtered$line <- ifelse(grepl("SCTI", polyA_list_filtered$sample), "SCTI",
                                  ifelse(grepl("007", polyA_list_filtered$sample),"CTF007",
                                         ifelse(grepl("069",polyA_list_filtered$sample),"CTF069",
                                                ifelse(grepl("127", polyA_list_filtered$sample),"CTM127",
                                                       ifelse(grepl("014",polyA_list_filtered$sample),"CTM014","M3_36S")))))


# Sequence plot ----------------
all_sites_unfiltered$modified <- ifelse(all_sites_unfiltered$probability_modified >= 0.9, "modified", "unmodified")

all_sites_unfiltered <- separate(all_sites_unfiltered, col=transcript_id, into=c("transcript_id", "gene_id", NULL), sep="\\.", extra="drop")
all_sites_unfiltered$gene_id <- sub("^[^|]*\\|", "", all_sites_unfiltered$gene_id)

unmodified <- subset(all_sites_unfiltered, probability_modified < 0.5)
modified <- subset(all_sites_unfiltered, probability_modified >= 0.9)

all_sites_collapse <- all_sites_unfiltered %>% 
  group_by(gene_id, transcript_position, transcript_id)%>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= "|")))

transcripts_filtered <- transcripts %>% filter(grepl("\\,", sample_id))
genes_filtered <- genes %>% filter(grepl("\\,", sample_id))

combined_sites_unique <- annotated_modified %>%  group_by(gene_id, transcript_id)%>%
  dplyr::filter(length(unique(sample_id))>1) %>% 
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= "|")))

combined_sites_unique$D1_F <- ifelse(grepl("D1_F", combined_sites_unique$group_id) | grepl("F_D1", combined_sites_unique$group_id), "Yes", "-")
combined_sites_unique$D14_F <- ifelse(grepl("D14_F", combined_sites_unique$group_id) | grepl("F_D14", combined_sites_unique$group_id), "Yes", "-")
combined_sites_unique$D1_M <- ifelse(grepl("D1_M", combined_sites_unique$group_id) | grepl("M_D1", combined_sites_unique$group_id), "Yes", "-")
combined_sites_unique$D14_M <- ifelse(grepl("D14_M", combined_sites_unique$group_id) | grepl("M_D14", combined_sites_unique$group_id), "Yes", "-")

combined_sites_unique <- separate(combined_sites_unique, col=gene_id, into=c("gene_id", NULL), sep="\\.", extra="drop")



d14_f_genes <- combined_sites_unique %>% filter(D14_F=="Yes" & D1_F=="-" & D14_M=="-" & D1_M=="-")
d14_f_genes <- d14_f_genes %>% 
  group_by(gene_id)%>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= "|")))
write.csv(d14_f_genes,"d14_f_genes.csv",row.names = FALSE,quote=FALSE)

d1_f_genes <- combined_sites_unique %>% filter(D14_F=="-" & D1_F=="Yes" & D14_M=="-" & D1_M=="-")
d1_f_genes <- d1_f_genes %>% 
  group_by(gene_id)%>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= "|")))
write.csv(d1_f_genes,"d1_f_genes.csv",row.names = FALSE,quote=FALSE)

d14_m_genes <- combined_sites_unique %>% filter(D14_F=="-" & D1_F=="-" & D14_M=="Yes" & D1_M=="-")
d14_m_genes <- d14_m_genes %>% 
  group_by(gene_id)%>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= "|")))
write.csv(d14_m_genes,"d14_m_genes.csv",row.names = FALSE,quote=FALSE)

d1_m_genes <- combined_sites_unique %>% filter(D14_F=="-" & D1_F=="-" & D14_M=="-" & D1_M=="Yes")
d1_m_genes <- d1_m_genes %>% 
  group_by(gene_id)%>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= "|")))
write.csv(d1_m_genes,"d1_m_genes.csv",row.names = FALSE,quote=FALSE)

combined_sites_unanot <- combined_sites_filtered %>% filter(!(genome_pos %in% directRMDB$start | genome_pos %in% m6a_atlas$start))

combined_unmod <- all_sites_unfiltered %>% 
  group_by(transcript_position, transcript_id)%>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ",")))


all_sites_unfiltered$modified <- ifelse(all_sites_unfiltered$probability_modified >=0.9, "modified","unmodified")


all_d1_m <- subset(all_sites_unfiltered, group_id=="D1_M")
all_d14_m <- subset(all_sites_unfiltered, group_id=="D14_M")
all_d1_f <- subset(all_sites_unfiltered, group_id=="D1_F")
all_d14_f <- subset(all_sites_unfiltered, group_id=="D14_F")

modified_d1_m <- subset(annotated_modified,group_id=="D1_M")
modified_d14_m <- subset(annotated_modified, group_id=="D14_M")
modified_d1_f <- subset(annotated_modified, group_id=="D1_F")
modified_d14_f <- subset(annotated_modified, group_id=="D14_F")

transcript_d1_f <- all_d1_f %>% 
  group_by(transcript_id, modified) %>%
  summarise_all(list(~ paste0(unique(.[!is.na(.)]), collapse = ","))) %>%
  group_by(transcript_id) %>%
  summarise(modified = if_else(any(modified == "modified"), "modified", "unmodified"),
            across(everything(), ~ paste0(unique(.[!is.na(.)]), collapse = ",")))

process_dataset <- function(data, id_col, modified_col) {
  data %>%
    group_by({{ id_col }}, {{ modified_col }}) %>%
    summarise_all(list(~ paste0(unique(.[!is.na(.)]), collapse = ","))) %>%
    group_by({{ id_col }}) %>%
    summarise(modified = if_else(any({{ modified_col }} == "modified"), "modified", "unmodified"),
              across(everything(), ~ paste0(unique(.[!is.na(.)]), collapse = ",")))
}


all_d1_m_transcript <- process_dataset(all_d1_m, transcript_id, modified)
all_d1_m_gene <- process_dataset(all_d1_m, gene_id, modified)

all_d14_m_transcript <- process_dataset(all_d14_m, transcript_id, modified)
all_d14_m_gene <- process_dataset(all_d14_m, gene_id, modified)

all_d1_f_transcript <- process_dataset(all_d1_f, transcript_id, modified)
all_d1_f_gene <- process_dataset(all_d1_f, gene_id, modified)

all_d14_f_transcript <- process_dataset(all_d14_f, transcript_id, modified)
all_d14_f_gene <- process_dataset(all_d14_f, gene_id, modified)

transcripts <- rbind(all_d1_f_transcript,all_d14_f_transcript,all_d14_m_transcript,all_d1_m_transcript)
transcripts$modified <- factor(transcripts$modified, c("unmodified","modified"))
genes <- rbind(all_d1_f_gene,all_d14_f_gene,all_d14_m_gene,all_d1_m_gene)
genes$modified <- factor(genes$modified, c("unmodified","modified"))

process_modified <- function(dataset, id_col, pos_col, kmer_col) {
  dataset %>%
    group_by({{ id_col }}, {{ pos_col }}, {{ kmer_col }}) %>%
    summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ",")))
}

unmodified_edit <- process_modified(unmodified, transcript_id,transcript_position,kmer)
modified_edit <- process_modified(modified, transcript_id,transcript_position,kmer)
modified_d1_m_edit <- process_modified(modified_d1_m, transcript_id, transcript_position, kmer)
modified_d14_m_edit <- process_modified(modified_d14_m, transcript_id, transcript_position, kmer)
modified_d1_f_edit <- process_modified(modified_d1_f, transcript_id, transcript_position, kmer)
modified_d14_f_edit <- process_modified(modified_d14_f, transcript_id, transcript_position, kmer)

unmodified_kmer <- ggseqlogo(unmodified_edit$kmer, method='prob') + ggtitle("Unmodified")+
  theme(plot.title = element_text(size=40),
    axis.text.x  = element_text(vjust = 0.5, hjust=1, size = 18),axis.title.y = element_text(size=32),
    axis.text.y  = element_text(vjust=0.5, size=18))

modified_kmer <- ggseqlogo(modified_edit$kmer, method='prob') + ggtitle("All modified")+
  theme(plot.title = element_text(size=40),
        axis.text.x  = element_text(vjust = 0.5, hjust=1, size = 18),axis.title.y = element_text(size=24),
        axis.text.y  = element_text(vjust=0.5, size=18))

D14M_sites <- ggseqlogo(modified_d14_m_edit$kmer, method='prob') + ggtitle("D14 XY")+
  theme(plot.title = element_text(size=40),
        axis.text.x  = element_text(vjust = 0.5, hjust=1, size = 18),axis.title.y = element_text(size=32),
        axis.text.y  = element_text(vjust=0.5, size=18))

D1F_sites <- ggseqlogo(modified_d1_f_edit$kmer, method='prob') + ggtitle("D1 XX")+
  theme(plot.title = element_text(size=40),
        axis.text.x  = element_text(vjust = 0.5, hjust=1, size = 18),axis.title.y = element_text(size=32),
        axis.text.y  = element_text(vjust=0.5, size=18))

D1M_sites <- ggseqlogo(modified_d1_m_edit$kmer, method='prob') + ggtitle("D1 XY")+
  theme(plot.title = element_text(size=40),
        axis.text.x  = element_text(vjust = 0.5, hjust=1, size = 18),axis.title.y = element_text(size=32),
        axis.text.y  = element_text(vjust=0.5, size=18))

D14F_sites <- ggseqlogo(modified_d14_f_edit$kmer, method='prob') + ggtitle("D14 XX")+
  theme(plot.title = element_text(size=40),
        axis.text.x  = element_text(vjust = 0.5, hjust=1, size = 18),axis.title.y = element_text(size=32),
        axis.text.y  = element_text(vjust=0.5, size=18))

kmers_plots <- patchwork::wrap_elements(unmodified_kmer / modified_kmer / D14F_sites / D1F_sites / D14M_sites / D1M_sites)
ggsave('kmers.png',kmers_plots, width = 30, height = 30, bg="transparent")


# Modified proportion plot ----------------

labels_man <- c("D1_F"="D1 XX","D14_F"="D14 XX","D1_M"="D1 XY","D14_M"="D14 XY")

custom_order <- c("unmodified", "modified")
transcripts_filtered$modified <- factor(transcripts_filtered$modified, levels = custom_order)
genes_filtered$modified <- factor(genes_filtered$modified, levels = custom_order)
transcripts$modified <- factor(transcripts$modified, levels = custom_order)
genes$modified <- factor(genes$modified, levels = custom_order)


modified_transcripts <- ggplot(transcripts, aes(x=group_id, fill=group_id)) + 
  geom_bar(aes(y= (..count..)/sum(..count..)*100, alpha=modified), position="fill") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  scale_x_discrete(expand=c(0,0), labels=labels_man)+
  scale_fill_viridis(discrete=TRUE, labels=labels_man)+
  labs(y = "Proportion of modified transcripts",fill=element_blank(), alpha=element_blank())+
  guides(fill="none")+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"), legend.position="bottom",
        legend.text = element_text(size=32),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size=32),
        axis.title.y = element_text(size=32),
        axis.text.y  = element_text(vjust=0.5, size=32))

modified_genes <- ggplot(genes, aes(x=group_id, fill=group_id)) + 
  geom_bar(aes(y= (..count..)/sum(..count..)*100, alpha=modified), position="fill") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  scale_x_discrete(expand=c(0,0),labels=labels_man)+
  scale_fill_viridis(discrete=TRUE, labels=labels_man)+
  labs(y = "Proportion of modified genes", fill=element_blank(), alpha=element_blank())+
  guides(fill="none")+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"), legend.position="bottom",
        legend.text = element_text(size=32),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size=32),
        axis.title.y = element_text(size=32),
        axis.text.y  = element_text(vjust=0.5, size=32))

proportion_plots <- patchwork::wrap_elements((modified_transcripts | modified_genes) +
                                               plot_layout(guides="collect") & theme(legend.position = "bottom"))


# Upset plot ----------------

annotated_modified <- annotated_m6a_filtered %>%
    group_by(genome_pos,group_id,transcript_id) %>%
    summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ",")))
annotated_modified <- annotated_modified %>%
mutate(group_id = ifelse(group_id == "F_D1", "D1_F", group_id)) %>%
mutate(group_id = ifelse(group_id == "F_D14", "D14_F", group_id)) %>%
  group_by(genome_pos,group_id) %>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ",")))

annotated_modified$d1_m <- ifelse(annotated_modified$group_id=="D1_M", "1","0")
annotated_modified$d1_f <- ifelse(annotated_modified$group_id=="D1_F", "1","0")
annotated_modified$d14_m <- ifelse(annotated_modified$group_id=="D14_M", "1","0")
annotated_modified$d14_f <- ifelse(annotated_modified$group_id=="D14_F", "1","0")

annotated_modified <- annotated_modified %>%
  group_by(genome_pos) %>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ",")))%>%
  mutate(d1_m = ifelse(d1_m == "0,1", "1", d1_m)) %>%
  mutate(d1_f = ifelse(d1_f == "0,1", "1", d1_f)) %>%
  mutate(d1_f = ifelse(d1_f == "1,0", "1", d1_f)) %>%
  mutate(d14_f = ifelse(d14_f == "1,0", "1", d14_f)) %>%
  mutate(d14_m = ifelse(d14_f == "1,0", "1", d14_m)) %>%
  mutate(d14_m = ifelse(d14_m == "0,1", "1", d14_m)) 
  

datasets = colnames(combined_sites_filtered)[19:22]
combined_sites_filtered[datasets]=combined_sites_filtered[datasets]=="Yes"
t(head(combined_sites_filtered[datasets],3))

upset(combined_sites_filtered,datasets,name="Dataset",width_ratio = 0.1)
labels_manual <- c("D1_F"="D1 XX","D14_F"="D14 XX","D1_M"="D1 XY","D14_M"="D14 XY")


combined_sites_filtered$annotated <- ifelse(combined_sites_filtered$genome_pos %in% m6a_atlas$start | combined_sites_filtered$genome_pos %in% directRMDB$start, "Annotated","Not")

m6a_upset <- upset(combined_sites_filtered, datasets, name='Conditions', width_ratio=0.1, themes=upset_modify_themes(
  list('intersections_matrix'=theme(axis.text.y=element_text(size=18),
                                    axis.title.x=element_text(size=32)))),set_sizes = (upset_set_size(
                                      geom=geom_bar(width = 0.8, show_guide=FALSE), position="right") +
                                        theme(axis.text.x=element_text(angle=270, size=32),
                                              axis.title.x = element_text(size=40))),
  queries=list(upset_query(intersect="D1_F",color="#450C54",fill="#450C54"),
               upset_query(intersect="D1_M",color="#31678D",fill="#31678D"),
               upset_query(intersect="D14_F",color="#35B879",fill="#35B879"),
               upset_query(intersect="D14_M",color="#FDE724",fill="#FDE724")),
  labeller = as_labeller(labels_manual),
  base_annotations=list('# of m6A sites'=intersection_size() +
                          theme(panel.grid = element_blank(),
                                axis.text.y = element_text(size=32),
                                axis.title.y = element_text(size=40),
                                legend.title=element_text(size=40),
                                legend.text=element_text(size=32,hjust=0))),guides="over")+
  theme(panel.grid = element_blank())

full_upset <- patchwork::wrap_elements(m6a_upset) 

ggsave('m6a_upset.png',full_upset, width = 20, height = 15, bg="transparent")

# Venn diagram----------------
m6a_atlas <- read.delim("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Experiments/Deep read seq/Direct RNA/Full exp/m6a/Database/hg38_CL_Tech.txt")

m6a_atlas <- m6a_atlas %>%
  group_by(start,end) %>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ",")))
 
directRMDB <- read.delim("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Experiments/Deep read seq/Direct RNA/Full exp/m6a/Database/HomoSapiens/HomoSapiens_siteinfo.txt")
 
directRMDB <- directRMDB %>%
  filter(modification=="m6A") %>%
  group_by(seqnames, start,end) %>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ",")))

site_list <- list(
  microglia_m6A=combined_sites_filtered$genome_pos,
  m6A_Atlas=m6a_atlas$start,
  DirectRMDB=directRMDB$start
)

venn=Venn(site_list)
data=process_data(venn)

venn_set(data)
venn_region(data)

dataset_venn <- ggplot()+geom_polygon(aes(X,Y,fill=id,group=id),data=venn_regionedge(data), show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE)+
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(data), 
            show.legend = FALSE) +
  scale_color_viridis(discrete = TRUE)+
  geom_text(aes(X, Y, label = c("microglia\nm6A","m6A-Atlas","DirectRMDB")),size=5,
            data = venn_setlabel(data)) +
  geom_label(aes(X, Y, label = count), 
             data = venn_regionlabel(data),size=8,show.legend = FALSE) +
  coord_equal() +
  theme_void()+theme(legend.position = "none")

unique_sites_m6a_intersections <- read.csv("~/path/to/file")
d14f_highlight_intersections <- read.csv("~/path/to/file")
d1f_highlight_term_intersections <- read.csv("~/path/to/file")

female_go_terms <- rbind(d1f_highlight_term_intersections,d14f_highlight_intersections)

unique_sites <- ggplot(data=unique_sites_m6a_intersections, aes(x=negative_log10_of_adjusted_p_value, y=term_name, color=adjusted_p_value, size=intersection_size))+
  geom_point()+
  scale_color_viridis(discrete=FALSE)+
  scale_y_discrete(labels=wrap_format(45))+
  labs(y = element_blank(), color="Adjusted p value",size="Intersection size", x="-log10(p adjusted)")+
  facet_wrap(~source, nrow=3, scales="free", strip.position = "right")+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"),
        legend.position="right",
        legend.text = element_text(size=14 ),
        plot.title = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=14),
        strip.text.x = element_text(size=16))


unique_sites_fem <- ggplot(data=female_go_terms, aes(x=sample, y=term_name, color=adjusted_p_value, size=intersection_size))+
  geom_point()+
  scale_color_viridis(discrete=FALSE)+
  scale_y_discrete(labels=wrap_format(45))+
  labs(y = element_blank(), color="Adjusted p value",size="Intersection size", x="Group")+
  facet_wrap(~source, nrow=3, scales="free_y", strip.position = "right")+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"),
        legend.position="right",
        legend.text = element_text(size=14 ),
        plot.title = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=14),
        strip.text.x = element_text(size=16))

part_graph <- (dataset_venn/unique_sites)

comparison_graphs <- (part_graph|unique_sites_fem)+
  plot_annotation(tag_levels = 'A') + plot_layout(heights=c(1,1,2), widths=c(2,1))

ggsave('comparison_graphs.png',comparison_graphs, width = 20, height = 15, bg="transparent")


# Kmer plot ----------------

kmer_propportion <- ggplot(all_sites_unfiltered, aes(x=modified,y=kmer, fill=kmer)) +
  geom_bar(stat = "count",aes(y= (..count..)/sum(..count..)*100), position = "fill") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  scale_x_discrete(expand = c(0, 0))+
  labs(y = "Proportion", x="Modification status", fill=element_blank()) +
  scale_fill_viridis(discrete=TRUE) +
  theme(legend.position = "right")+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"), legend.position="bottom",
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.x  = element_text(vjust = 0.5, hjust=1, size = 12),
        axis.title.y = element_text(size=12),
        axis.text.y  = element_text(vjust=0.5, size=12))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

# Combined plot ----------------

m6a_generic <- ((
  (((((metagene / proportion_plots) | kmers_plots) +
    plot_layout(heights = c(1,1,2), widths = c(1.5,1)))
  / full_upset) + plot_layout(heights = c(4,2))))+
    plot_annotation(tag_levels = 'A')) & theme(plot.tag = element_text(size=32))

ggsave('m6a_generic.png',m6a_generic, width = 30, height = 40, bg="transparent")

# Correlation plot ----------------

m6aPerGene$nIso <- isoPerGene_talon$nIso[match(m6aPerGene$gene_name,isoPerGene_talon$associatedGene)]

length_corr_plot <- ggscatter(m6aPerTranscript, x = "site_count", y = "kbp",color="#31678d",
                              add="reg.line",
                              conf.int = TRUE)+
  stat_cor(label.x=10,label.y=50, show_guide=FALSE, size=10)+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0), limits=c(0,29))+
  labs(y = "Transcript Length (kb)", x = "# of m6A sites")+
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.text.x  = element_text(size=24),
        axis.title.y = element_text(size=24),
        axis.text.y  = element_text(vjust=0.5, size=24))

exon_corr_plot <- ggscatter(m6aPerTranscript, x = "site_count", y = "den_norm", color="#31678d",
                            add="reg.line",
                            conf.int = TRUE)+
  stat_cor(label.x=10,label.y=50, show_guide=FALSE, size=10)+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0), limits=c(0,29))+
  labs( x = "# of m6A sites", y = expression(paste("Exon density *10"^3)))+
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.text.x  = element_text(size=24),
        axis.title.y = element_text(size=24),
        axis.text.y  = element_text(vjust=0.5, size=24))

utr5_corr_plot <- ggscatter(m6aPerTranscript, x = "site_count", y = "UTR_5", color="#31678d",
                            add="reg.line",
                            conf.int = TRUE)+
  stat_cor(label.x=10,label.y=5000, show_guide=FALSE, size=10)+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0), limits=c(0,29))+
  labs( x = "# of m6A sites", y = expression(paste("5' UTR length")))+
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.text.x  = element_text(size=24),
        axis.title.y = element_text(size=24),
        axis.text.y  = element_text(vjust=0.5, size=24))

utr3_corr_plot <- ggscatter(m6aPerTranscript, x = "site_count", y = "UTR_3", color="#31678d",
                            add="reg.line",
                            conf.int = TRUE)+
  stat_cor(label.x=10,label.y=5000, show_guide=FALSE, size=10)+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0), limits=c(0,29))+
  labs( x = "# of m6A sites", y = expression(paste("3' UTR length")))+
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.text.x  = element_text(size=24),
        axis.title.y = element_text(size=24),
        axis.text.y  = element_text(vjust=0.5, size=24))

polyA_corr_plot <- ggscatter(m6aPerTranscript, x = "site_count", y = "polyA", color="#31678d",
                            add="reg.line",
                            conf.int = TRUE)+
  stat_cor(label.x=10,label.y=400, show_guide=FALSE, size=10)+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0), limits=c(0,29))+
  labs( x = "# of m6A sites", y = expression(paste("poly(A) length")))+
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.text.x  = element_text(size=24),
        axis.title.y = element_text(size=24),
        axis.text.y  = element_text(vjust=0.5, size=24))

iso_corr_plot <- ggscatter(m6aPerGene, x = "site_count", y = "nIso",color="#31678d",
                              add="reg.line",
                              conf.int = TRUE)+
  stat_cor(label.x=10,label.y=1000, show_guide=FALSE, size=10)+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0), limits=c(0,29))+
  labs(y = "# of Isoforms", x = "# of m6A sites")+
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.text.x  = element_text(size=24),
        axis.title.y = element_text(size=24),
        axis.text.y  = element_text(vjust=0.5, size=24))

corr_plots <-((iso_corr_plot | length_corr_plot) / (exon_corr_plot| polyA_corr_plot)/(utr5_corr_plot|utr3_corr_plot)) +
    plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 24))

 

sites <- setdiff(annotated_modified$genome_pos, m6a_atlas$start)
sites <- setdiff(sites,directRMDB$start)


modified <- process_modified(annotated_modified, transcript_id,genome_pos,kmer)

unique_modified <- subset(modified, !(genome_pos %in% m6a_atlas$start))
unique_modified <- subset(unique_modified, !(genome_pos %in% directRMDB$start))

unique_modified_d1_f <- combined_sites_filtered %>% filter(D1_F=="Yes"&D14_F!="Yes"&D1_M!="Yes"&D14_M!="Yes")%>%
  subset(!(genome_pos %in% m6a_atlas$start))%>%
  subset(!(genome_pos %in% directRMDB$start))

unique_modified_d14_f <- combined_sites_filtered %>% filter(D1_F!="Yes"&D14_F=="Yes"&D1_M!="Yes"&D14_M!="Yes")%>%
  subset(!(genome_pos %in% m6a_atlas$start))%>%
  subset(!(genome_pos %in% directRMDB$start))

unique_modified_d1_m <- combined_sites_filtered %>% filter(D1_F!="Yes"&D14_F!="Yes"&D1_M=="Yes"&D14_M!="Yes")%>%
  subset(!(genome_pos %in% m6a_atlas$start))%>%
  subset(!(genome_pos %in% directRMDB$start))

unique_modified_d14_m <- combined_sites_filtered %>% filter(D1_F!="Yes"&D14_F!="Yes"&D1_M!="Yes"&D14_M=="Yes")%>%
  subset(!(genome_pos %in% m6a_atlas$start))%>%
  subset(!(genome_pos %in% directRMDB$start))



unique_modified <- unique_modified %>%
  filter(transcript_type!="artifact") %>%
  group_by(genome_pos) %>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ",")))

transcripts <- unique_modified %>%
  filter(transcript_type=="protein_coding")
