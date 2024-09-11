library(data.table)
library(vroom)
library(dplyr)
library(genomation)
library(rtracklayer)
library(transPlotR)
library(IsoformSwitchAnalyzeR)
library("BSgenome.Hsapiens.UCSC.hg38")
library(factR)
library(ggbio)
library(ggtranscript)
library(viridis)
library(viridisLite)

all_comparisons <- vroom("~/path/to/xpore/data")
sex_comparisons <- vroom("~/path/to/xpore/data")
time_comparisons <- vroom("~/path/to/xpore/data")

gencode <- "~/path/to/gtf/file"
gff_gencode <- rtracklayer::import(gencode) %>%
  as.data.frame() %>%
  tibble::as_tibble() 
gff_gencode <- separate(gff_gencode, col=transcript_id, into=c("transcript_id", NULL), sep="\\.", extra="drop")


tech <- "~/path/to/gtf/file"
gff_tech <- rtracklayer::import(tech) %>%
  as.data.frame() %>%
  tibble::as_tibble()

all_comparisons$pval_F_D1_vs_F_D14_adjust <- p.adjust(all_comparisons$pval_F_D1_vs_F_D14, method = "fdr", n = length(all_comparisons$pval_F_D1_vs_F_D14))
all_comparisons$pval_F_D1_vs_M_D1_adjust <- p.adjust(all_comparisons$pval_F_D1_vs_M_D1, method = "fdr", n = length(all_comparisons$pval_F_D1_vs_M_D1))
all_comparisons$pval_F_D1_vs_M_D14_adjust <- p.adjust(all_comparisons$pval_F_D1_vs_M_D14, method = "fdr", n = length(all_comparisons$pval_F_D1_vs_M_D14))
all_comparisons$pval_F_D14_vs_M_D1_adjust <- p.adjust(all_comparisons$pval_F_D14_vs_M_D1, method = "fdr", n = length(all_comparisons$pval_F_D14_vs_M_D1))
all_comparisons$pval_F_D14_vs_M_D14_adjust <- p.adjust(all_comparisons$pval_F_D14_vs_M_D14, method = "fdr", n = length(all_comparisons$pval_F_D14_vs_M_D14))
all_comparisons$pval_M_D1_vs_M_D14_adjust <- p.adjust(all_comparisons$pval_M_D1_vs_M_D14, method = "fdr", n = length(all_comparisons$pval_M_D1_vs_M_D14))
all_comparisons$pval_F_D1_vs_all_adjust <- p.adjust(all_comparisons$pval_F_D1_vs_all, method = "fdr", n = length(all_comparisons$pval_F_D1_vs_all))
all_comparisons$pval_F_D14_vs_all_adjust <- p.adjust(all_comparisons$pval_F_D14_vs_all, method = "fdr", n = length(all_comparisons$pval_F_D14_vs_all))
all_comparisons$pval_M_D1_vs_all_adjust <- p.adjust(all_comparisons$pval_M_D1_vs_all, method = "fdr", n = length(all_comparisons$pval_M_D1_vs_all))
all_comparisons$pval_M_D14_vs_all_adjust <- p.adjust(all_comparisons$pval_M_D14_vs_all, method = "fdr", n = length(all_comparisons$pval_M_D14_vs_all))



sex_comparisons$pval_F_vs_M_adjust <- p.adjust(sex_comparisons$pval_F_vs_M, method = "fdr", n = length(sex_comparisons$pval_F_vs_M))
time_comparisons$pval_D1_vs_D14_adjust <- p.adjust(time_comparisons$pval_D1_vs_D14, method = "fdr", n = length(time_comparisons$pval_D1_vs_D14))


kmer_list <- unique(all_sites_unfiltered$kmer)


all_comparisons <- all_comparisons %>% filter(kmer %in% kmer_list)
sex_comparisons <- sex_comparisons %>% filter(kmer %in% kmer_list)
time_comparisons <- time_comparisons %>% filter(kmer %in% kmer_list)

write.csv(all_comparisons, "all_comparisons_drach_xpore_adjust_p.csv", row.names=F, quote=F)
write.csv(sex_comparisons, "sex_comparisons_drach_xpore_adjust_p.csv", row.names=F, quote=F)
write.csv(time_comparisons, "time_comparisons_drach_xpore_adjust_p.csv", row.names=F, quote=F)


all_comparisons <- separate(all_comparisons, col=id, into=c("transcript_id", "gene_id", NULL), sep="\\.", extra="drop")
all_comparisons$gene_id <- sub("^[^|]*\\|", "", all_comparisons$gene_id)

sex_comparisons <- separate(sex_comparisons, col=id, into=c("transcript_id", "gene_id", NULL), sep="\\.", extra="drop")
sex_comparisons$gene_id <- sub("^[^|]*\\|", "", sex_comparisons$gene_id)

time_comparisons <- separate(time_comparisons, col=id, into=c("transcript_id", "gene_id", NULL), sep="\\.", extra="drop")
time_comparisons$gene_id <- sub("^[^|]*\\|", "", time_comparisons$gene_id)

annotated_modified <- annotated_modified %>% separate(gene_id, into=c("gene_id", NULL), sep="\\.", extra="drop")

colnames(annotated_modified)[3]<-"position"

annotated_modified <- annotated_modified %>% 
  group_by(gene_id, genome_pos, transcript_id, position) %>% 
  dplyr::filter(length(unique(sample_id))>1) %>% 
  ungroup()

mod_columns <- c("transcript_id","gene_id","position", "kmer")
all_xpore <- semi_join(all_comparisons, annotated_modified, by=mod_columns)
all_xpore$genelabels <- NA
all_xpore$direction <- NA
all_xpore_filter <- all_xpore %>% filter(abs(diff_mod_rate_F_D1_vs_F_D14)>0.3 | abs(diff_mod_rate_F_D1_vs_M_D1)>0.3
                                         | abs(diff_mod_rate_F_D14_vs_M_D14)>0.3 | abs(diff_mod_rate_M_D1_vs_M_D14)>0.3) %>%
  filter(pval_F_D1_vs_F_D14_adjust<0.05 | pval_F_D1_vs_M_D1_adjust<0.05 | pval_F_D14_vs_M_D14_adjust<0.05 | pval_M_D1_vs_M_D14_adjust<0.05)
all_xpore_filter <- left_join(all_xpore_filter,annotated_modified, by=mod_columns)
all_xpore_filter <- all_xpore_filter[,c(-79,-80)]
all_xpore_filter <- all_xpore_filter %>% distinct(transcript_id,gene_id,genome_pos,position,kmer,.keep_all = TRUE)

female_d1_d14 <- all_xpore[,c(1:7,66,76:77)]
female_d1_d14$genelabels <- ifelse(abs(female_d1_d14$diff_mod_rate_F_D1_vs_F_D14)>0.3 & female_d1_d14$pval_F_D1_vs_F_D14_adjust<0.05, female_d1_d14$transcript_id, NA)
female_d1_d14$direction <-ifelse(female_d1_d14$diff_mod_rate_F_D1_vs_F_D14 > 0.3 & female_d1_d14$pval_F_D1_vs_F_D14_adjust<0.05, "Higher in Day 1",
                                 ifelse(female_d1_d14$diff_mod_rate_F_D1_vs_F_D14 < -0.3 & female_d1_d14$pval_F_D1_vs_F_D14_adjust<0.05,"Higher in Day 14", "Neither"))
female_d1_d14_filter <- female_d1_d14 %>% filter(abs(diff_mod_rate_F_D1_vs_F_D14)>0.3)
female_d1_d14_filter <- female_d1_d14_filter %>% filter(pval_F_D1_vs_F_D14_adjust<0.05)
female_d1_d14_filter <- left_join(female_d1_d14_filter,annotated_modified, by=mod_columns)
female_d1_d14_filter <- female_d1_d14_filter %>% distinct(transcript_id,gene_id,genome_pos,position,kmer,.keep_all = TRUE)
write.csv(female_d1_d14_filter, "female_d1_d14_filter.csv", row.names=F, quote=F)

  
male_d1_d14 <- all_xpore[,c(1:4,20:22,71)]
male_d1_d14 <- drop_na(male_d1_d14)
male_d1_d14$genelabels <- ifelse(abs(male_d1_d14$diff_mod_rate_M_D1_vs_M_D14)>0.3 & male_d1_d14$pval_M_D1_vs_M_D14_adjust<0.05, male_d1_d14$transcript_id, NA)
male_d1_d14$direction <-ifelse(male_d1_d14$diff_mod_rate_M_D1_vs_M_D14 > 0.3 & male_d1_d14$pval_M_D1_vs_M_D14_adjust<0.05, "Higher in Day 1",
                               ifelse(male_d1_d14$diff_mod_rate_M_D1_vs_M_D14 < -0.3 & male_d1_d14$pval_M_D1_vs_M_D14_adjust<0.05,"Higher in Day 14", "Neither"))
male_d1_d14_filter <- male_d1_d14 %>% filter(abs(diff_mod_rate_M_D1_vs_M_D14)>0.3)
male_d1_d14_filter <- male_d1_d14_filter %>% filter(pval_M_D1_vs_M_D14_adjust<0.05)
male_d1_d14_filter <- left_join(male_d1_d14_filter,annotated_modified, by=mod_columns)
male_d1_d14_filter <- male_d1_d14_filter %>% distinct(transcript_id,gene_id,genome_pos,position,kmer,.keep_all = TRUE)
write.csv(male_d1_d14, "male_d1_d14.csv", row.names=F, quote=F)



sex_d1 <- all_xpore[,c(1:4,8:10,67,76:77)]
sex_d1$genelabels <- ifelse(abs(sex_d1$diff_mod_rate_F_D1_vs_M_D1)>0.3 & sex_d1$pval_F_D1_vs_M_D1_adjust<0.05, sex_d1$transcript_id, NA)
sex_d1$direction <- ifelse(sex_d1$diff_mod_rate_F_D1_vs_M_D1 > 0.3 & sex_d1$pval_F_D1_vs_M_D1_adjust<0.05, "Higher in XX",
                               ifelse(sex_d1$diff_mod_rate_F_D1_vs_M_D1 < -0.3 & sex_d1$pval_F_D1_vs_M_D1_adjust<0.05,"Higher in XY", "Neither"))
sex_d1_filter <- sex_d1 %>% filter(abs(diff_mod_rate_F_D1_vs_M_D1)>0.3)
sex_d1_filter <- sex_d1_filter %>% filter(pval_F_D1_vs_M_D1_adjust<0.05)
sex_d1_filter <- left_join(sex_d1_filter,annotated_modified, by=mod_columns)
sex_d1_filter <- sex_d1_filter %>% distinct(transcript_id,gene_id,genome_pos,position,kmer,.keep_all = TRUE)
write.csv(sex_d1_filter, "sex_d1_filter.csv", row.names=F, quote=F)


sex_d14 <- all_xpore[,c(1:4,17:19,70)]
sex_d14 <- drop_na(sex_d14)
sex_d14$genelabels <- ifelse(abs(sex_d14$diff_mod_rate_F_D14_vs_M_D14)>0.3 & sex_d14$pval_F_D14_vs_M_D14_adjust<0.05, sex_d14$transcript_id, NA)
sex_d14$direction <- ifelse(sex_d14$diff_mod_rate_F_D14_vs_M_D14 > 0.3 & sex_d14$pval_F_D14_vs_M_D14_adjust<0.05, "Higher in XX",
                                                ifelse(sex_d14$diff_mod_rate_F_D14_vs_M_D14 < -0.3 & sex_d14$pval_F_D14_vs_M_D14_adjust<0.05,"Higher in XY", "Neither"))
sex_d14_filter <- sex_d14 %>% filter(abs(diff_mod_rate_F_D14_vs_M_D14)>0.3)
sex_d14_filter <- sex_d14_filter %>% filter(pval_F_D14_vs_M_D14_adjust<0.05)
sex_d14_filter <- left_join(sex_d14_filter,annotated_modified, by=mod_columns)
sex_d14_filter <- sex_d14_filter %>% distinct(transcript_id,gene_id,genome_pos,position,kmer,.keep_all = TRUE)
write.csv(sex_d14, "sex_d14.csv", row.names=F, quote=F)


sex_xpore$genelabels <- sex_xpore$label
sex_xpore$direction <- ifelse(sex_xpore$diff_mod_rate_F_vs_M > 0.3 & sex_xpore$pval_F_vs_M_adjust < 0.05, "Higher in XX",
                              ifelse(sex_xpore$diff_mod_rate_F_vs_M < -0.3 & sex_xpore$pval_F_vs_M_adjust < 0.05, "Higher in XY", "Neither"))



sex_xpore <- semi_join(sex_comparisons, annotated_modified, by=mod_columns)
sex_xpore$genelabel <- ifelse((abs(sex_xpore$diff_mod_rate_F_vs_M)>0.3) & sex_xpore$pval_F_vs_M_adjust<0.05,
                          sex_xpore$transcript_id,NA)
sex_xpore_filter <- sex_xpore %>% filter(abs(diff_mod_rate_F_vs_M)>0.3)
sex_xpore_filter$gene_name <- annotated_modified$gene_name[match(sex_xpore_filter$gene_id,annotated_modified$gene_id)]
sex_xpore_filter <- left_join(sex_xpore_filter,annotated_modified, by=mod_columns)
sex_xpore_filter <- sex_xpore_filter[,c(-43,-44)]
sex_xpore_filter <- sex_xpore_filter %>% distinct(transcript_id,gene_id,genome_pos,position,kmer,.keep_all = TRUE)
sex_xpore_filter <- sex_xpore_filter %>% filter(pval_F_vs_M_adjust<0.05)
write.csv(sex_xpore_filter, "sex_xpore_filter_significant.csv", row.names=F, quote=F)


time_xpore <- semi_join(time_comparisons, annotated_modified, by=mod_columns)
time_xpore$genelabel <- ifelse((abs(time_xpore$diff_mod_rate_D1_vs_D14)>0.3) & time_xpore$pval_D1_vs_D14_adjust<0.05,
                               time_xpore$transcript_id,NA)
time_xpore_filter <- time_xpore %>% filter(abs(diff_mod_rate_D1_vs_D14)>0.3)
time_xpore_filter$gene_name <- annotated_modified$gene_name[match(time_xpore_filter$gene_id,annotated_modified$gene_id)]
time_xpore_filter <- left_join(time_xpore_filter,annotated_modified, by=mod_columns)
time_xpore_filter <- time_xpore_filter[,c(-43,-44)]
time_xpore_filter <- time_xpore_filter %>% distinct(transcript_id,gene_id,genome_pos,position,kmer,.keep_all = TRUE)
time_xpore_filter <- time_xpore_filter %>% filter(pval_D1_vs_D14_adjust<0.05)

time_xpore$direction <- ifelse(time_xpore$diff_mod_rate_D1_vs_D14 > 0.3 & time_xpore$pval_D1_vs_D14_adjust < 0.05, "Higher in Day 1",
                              ifelse(time_xpore$diff_mod_rate_D1_vs_D14 < -0.3 & time_xpore$pval_D1_vs_D14_adjust < 0.05, "Higher in Day 14", "Neither"))

write.csv(time_xpore_filter, "time_xpore_filter_significant.csv", row.names=F, quote=F)

iso_per_gene_method <- ggplot(sex_xpore_filter, aes(fill=method, x=nIsoCat),palette = c("skyblue", "#990099","#677f51")) + 
  geom_bar(stat="count", aes(y= (..count..)), position="dodge") +
  scale_fill_manual(values=colors_manual)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  labs(y = "Genes (%)", x = "Number of Isoforms", fill=element_blank())+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"), legend.position="bottom",
        legend.text = element_text(size=24),
        axis.title.x = element_text(size=24),
        axis.text.x  = element_text(size=24),
        axis.title.y = element_text(size=24),
        axis.text.y  = element_text(vjust=0.5, size=24))

sex_xpore$genelabels <- sex_xpore$label
sex_xpore$direction <- ifelse(sex_xpore$diff_mod_rate_F_vs_M > 0.3 & sex_xpore$pval_F_vs_M_adjust < 0.05, "Higher in XX",
                              ifelse(sex_xpore$diff_mod_rate_F_vs_M < -0.3 & sex_xpore$pval_F_vs_M_adjust < 0.05, "Higher in XY", "Neither"))


# Venn--------------------

day_14_f <- female_d1_d14_filter %>% filter(direction=="Higher in Day 14")
day_14_m <- male_d1_d14_filter %>% filter(direction=="Higher in Day 14")
day_14_xpore <- time_xpore_filter %>% filter(direction=="Higher in Day 14")

day_1_f <- female_d1_d14_filter %>% filter(direction=="Higher in Day 1")
day_1_m <- male_d1_d14_filter %>% filter(direction=="Higher in Day 1")
day_1_xpore <- time_xpore_filter %>% filter(direction=="Higher in Day 1")

site_list_day_1 <- list(
  xpore_f_1=day_1_f$gene_name,
  xpore_m_1=day_1_m$gene_name,
  time_xpore_both_1=day_1_xpore$gene_name.x
)

venn_1=Venn(site_list_day_1)
data_1=process_data(venn_1)

venn_set(data_1)
venn_region(data_1)

xpore_venn <- ggplot()+geom_polygon(aes(X,Y,fill=id,group=id),data=venn_regionedge(data_1), show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE)+
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(data_1), 
            show.legend = FALSE) +
  scale_color_viridis(discrete = TRUE)+
  geom_text(aes(X, Y, label = c("XX lines","XY lines","Both")),size=5,
            data = venn_setlabel(data_1)) +
  geom_label(aes(X, Y, label = count), 
             data = venn_regionlabel(data_1),size=8,show.legend = FALSE) +
  ggtitle("Higher on Day 1")+
  coord_equal() +
  theme_void()+theme(legend.position = "none")


site_list_day_14 <- list(
  xpore_f_14=day_14_f$gene_name,
  xpore_m_14=day_14_m$gene_name,
  time_xpore_both_14=day_14_xpore$gene_name.x
)

venn=Venn(site_list_day_14)
data=process_data(venn)

venn_set(data)
venn_region(data)

xpore_venn <- ggplot()+geom_polygon(aes(X,Y,fill=id,group=id),data=venn_regionedge(data), show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE)+
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(data), 
            show.legend = FALSE) +
  scale_color_viridis(discrete = TRUE)+
  geom_text(aes(X, Y, label = c("XX lines","XY lines","Both")),size=5,
            data = venn_setlabel(data)) +
  geom_label(aes(X, Y, label = count), 
             data = venn_regionlabel(data),size=8,show.legend = FALSE) +
  ggtitle("Higher on Day 14")+
  coord_equal() +
  theme_void()+theme(legend.position = "none")

# Volcanos--------------------

colors_manual <- c("Higher in XY"="skyblue","Neither"="#5f5f5f","Higher in XX"="#f1605dff")
colors_time <- c("Higher in Day 1"="#677f51","Neither"="#5f5f5f","Higher in Day 14"="#721f81ff")

xpore_time_comp <- ggplot(data=time_xpore, aes(x=(diff_mod_rate_D1_vs_D14), y=-log10(pval_D1_vs_D14_adjust),color=direction)) +
  geom_point(size=1) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + 
  geom_vline(xintercept = c(-0.3, 0.3), linetype='dashed') + 
  geom_text_repel(aes(label = genelabel),show_guide  = FALSE)+
  scale_color_manual(values = colors_time) +
  labs(x='Differential Modification Rate', y='-Log10(adjusted p value)',
       color='Differentially Modified Site', 
       title="Differentially methylated sites\n(Day 1 vs Day 14)") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 16),
        panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=12),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=12))

xpore_sex_comp <- ggplot(data=sex_xpore, aes(x=(diff_mod_rate_F_vs_M), y=-log10(pval_F_vs_M_adjust),color=direction)) +
  geom_point(size=1) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + 
  geom_vline(xintercept = c(-0.3, 0.3), linetype='dashed') + 
  geom_text_repel(aes(label = genelabel),show_guide  = FALSE)+
  scale_color_manual(values = colors_manual) +
  labs(x='Differential Modification Rate', y='-Log10(adjusted p value)',
       color='Differentially Modified Site', 
       title="Differentially methylated sites\n(XX vs XY)") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 16),
        panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=12),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=12))


xpore_female <- ggplot(data=female_d1_d14, aes(x=-(diff_mod_rate_F_D1_vs_F_D14), y=-log10(pval_F_D1_vs_F_D14_adjust),color=direction)) +
  geom_point(size=1) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + 
  geom_vline(xintercept = c(-0.3, 0.3), linetype='dashed') + 
  geom_text_repel(aes(label = genelabels),show_guide  = FALSE)+
  scale_color_manual(values = colors_time) +
  labs(x='Differential Modification Rate', y='-Log10(adjusted p value)',
       color='Differentially Modified Site', 
       title="Differentially methylated sites\n(XX lines (Day 1 vs Day 14))") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 16),legend.position = "bottom",
        panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=12),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=12))

xpore_male <- ggplot(data=male_d1_d14, aes(x=-(diff_mod_rate_M_D1_vs_M_D14), y=-log10(pval_M_D1_vs_M_D14_adjust),color=direction)) +
  geom_point(size=1) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + 
  geom_vline(xintercept = c(-0.3, 0.3), linetype='dashed') + 
  geom_text_repel(aes(label = genelabels),show_guide  = FALSE)+
  scale_color_manual(values = colors_time) +
  labs(x='Differential Modification Rate', y='-Log10(adjusted p value)',
       color='Differentially Modified Site', 
       title="Differentially methylated sites\n(XY lines (Day 1 vs Day 14))") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 16),
        panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"),legend.position = "bottom",
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=12),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=12))


xpore_day1 <- ggplot(data=sex_d1, aes(x=(diff_mod_rate_F_D1_vs_M_D1), y=-log10(pval_F_D1_vs_M_D1_adjust),color=direction)) +
  geom_point(size=1) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + 
  geom_vline(xintercept = c(-0.3, 0.3), linetype='dashed') + 
  geom_text_repel(aes(label = genelabels),show_guide  = FALSE)+
  scale_color_manual(values = colors_manual) +
  labs(x='Differential Modification Rate', y='-Log10(adjusted p value)',
       color='Differentially Modified Site', 
       title="Differentially methylated sites\n(Day 1 (XX vs XY))") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 16),
        panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=12),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=12))

xpore_day14 <- ggplot(data=sex_d14, aes(x=(diff_mod_rate_F_D14_vs_M_D14), y=-log10(pval_F_D14_vs_M_D14_adjust),color=direction)) +
  geom_point(size=1) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + 
  geom_vline(xintercept = c(-0.3, 0.3), linetype='dashed') + 
  geom_text_repel(aes(label = genelabels),show_guide  = FALSE)+
  scale_color_manual(values = colors_manual) +
  labs(x='Differential Modification Rate', y='-Log10(adjusted p value)',
       color='Differentially Modified Site', 
       title="Differentially methylated sites\n(Day 14 (XX vs XY))") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 16),
        panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=12),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=12))

# Transcripts--------------------

FAM104A_annotation <- gff_gencode %>% 
  dplyr::filter( 
    transcript_id == "ENST00000403627" | transcript_id=="ENST00000405159"
  )

FAM104A <- sex_xpore_filter %>% filter(gene_name.y =="FAM104A")    

FAM104A_exons <- FAM104A_annotation %>% dplyr::filter(type == "exon")
FAM104A_cds <- FAM104A_annotation %>% dplyr::filter(type == "CDS")

FAM104A_transcripts <- ggplot(FAM104A_exons, aes(
  xstart = start,
  xend = end,
  y = transcript_id
)) +
  geom_range(
    fill = "#31678d",color="#31678d",
    height = 0.25
  ) +
  geom_range(
    data = FAM104A_cds,fill = "#31678d",color="#31678d")+
  geom_intron(
    data = to_intron(FAM104A_exons, "transcript_id"),
    aes(strand = strand),
    arrow.min.intron.length = 500) +
  geom_segment(aes(x=73208177, y=2.15, xend=73208177, yend=2.25), size=1.5, color="#f1605dff")+
  geom_segment(aes(x=73208177, y=0.8, xend=73208177, yend=2.15), linetype="dashed",color="#f1605dff")+
  ggtitle(substitute(paste(italic("FAM104A"))))+
  annotate("text", x=73208177, y=2.3, label= "m6A") + 
  annotate("text",x=73212000, y=2.3, label="increased methylation in XX lines")+
  labs(x = "Genome Position", y = "Transcript") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(),legend.position = "bottom",
        legend.title = element_blank(),
        strip.text=element_text(size=10),
        panel.spacing.y = unit(0.75, "lines"),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.x  = element_text(vjust = 0.5, hjust=1, size =12),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(vjust=0.5, size=12),
        plot.title = element_text(hjust = 0.5))



MARCHF7_annotation <- gff_gencode %>% 
  dplyr::filter( 
    transcript_id == "ENST00000259050"
  )

MARCHF7 <- sex_xpore_filter %>% filter(gene_name.y =="MARCHF7")    

MARCHF7_exons <- MARCHF7_annotation %>% dplyr::filter(type == "exon")
MARCHF7_cds <- MARCHF7_annotation %>% dplyr::filter(type == "CDS")

MARCHF7_transcripts <- ggplot(MARCHF7_exons, aes(
  xstart = start,
  xend = end,
  y = transcript_id
)) +
  geom_range(
    fill = "#31678d",color="#31678d",
    height = 0.25
  ) +
  geom_range(
    data = MARCHF7_cds,fill = "#31678d",color="#31678d")+
  geom_intron(
    data = to_intron(MARCHF7_exons, "transcript_id"),
    aes(strand = strand),
    arrow.min.intron.length = 500) +
  geom_segment(aes(x=159747940, y=1.27, xend=159747940, yend=1.37), size=1.5, color="skyblue")+
  geom_segment(aes(x=159747940, y=0.7, xend=159747940, yend=1.27), linetype="dashed",color="skyblue")+
  ggtitle(substitute(paste(italic("MARCHF7"))))+
  annotate("text", x=159747940, y=1.4, label= "m6A") + 
  annotate("text",x=159755940, y=1.33, label="increased methylation in XY lines")+
  labs(x = "Genome Position") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(),legend.position = "bottom",
          legend.title = element_blank(),
          strip.text=element_text(size=10),
          panel.spacing.y = unit(0.75, "lines"),
          legend.text = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.x  = element_text(vjust = 0.5, hjust=1, size =12),
          axis.title.y = element_blank(),
          axis.text.y  = element_text(vjust=0.5, size=12),
          plot.title = element_text(hjust = 0.5))


TMIGD3_annotation <- gff_gencode %>% 
  dplyr::filter( 
    transcript_id == "ENST00000369717" | transcript_id=="ENST00000369716" |
      transcript_id ==  "ENST00000442484"
  )


TMIGD3 <- sex_xpore_filter %>% filter(gene_name.y =="TMIGD3")    

TMIGD3_exons <- TMIGD3_annotation %>% dplyr::filter(type == "exon")
TMIGD3_cds <- TMIGD3_annotation %>% dplyr::filter(type == "CDS")

TMIGD3_transcripts <- ggplot(TMIGD3_exons, aes(
  xstart = start,
  xend = end,
  y = transcript_id
)) +
  geom_range(
    fill = "#31678d",color="#31678d",
    height = 0.25
  ) +
  geom_range(
    data = TMIGD3_cds, fill = "#31678d",color="#31678d")+
  geom_intron(
    data = to_intron(TMIGD3_exons, "transcript_id"),
    aes(strand = strand),
    arrow.min.intron.length = 500) +
  geom_segment(aes(x=111483576, y=3.15, xend=111483576, yend=3.25), size=1.5, color="#f1605dff")+
  geom_segment(aes(x=111483576, y=0.7, xend=111483576, yend=3.15), linetype="dashed",color="#f1605dff")+
  ggtitle(substitute(paste(italic("TMIGD3"))))+
  annotate("text", x=111483576, y=3.3, label= "m6A") + 
  annotate("text",x=111493576, y=3.2, label="increased methylation in XX lines")+
  labs(x = "Genome Position", y = "Transcript") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(),legend.position = "bottom",
        legend.title = element_blank(),
        strip.text=element_text(size=10),
        panel.spacing.y = unit(0.75, "lines"),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.x  = element_text(vjust = 0.5, hjust=1, size =12),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(vjust=0.5, size=12),
        plot.title = element_text(hjust = 0.5))


ELP6_annotation <- gff_gencode %>% 
  dplyr::filter( 
    transcript_id == "ENST00000296149" 
  )


ELP6 <- sex_xpore_filter %>% filter(gene_name.y =="ELP6")    

ELP6_exons <- ELP6_annotation %>% dplyr::filter(type == "exon")
ELP6_cds <- ELP6_annotation %>% dplyr::filter(type == "CDS")

ELP6_transcripts <- ggplot(ELP6_exons, aes(
  xstart = start,
  xend = end,
  y = transcript_id
)) +
  geom_range(
    fill = "#31678d",color="#31678d",
    height = 0.25
  ) +
  geom_range(
    data = ELP6_cds, fill = "#31678d",color="#31678d")+
  geom_intron(
    data = to_intron(ELP6_exons, "transcript_id"),
    aes(strand = strand),
    arrow.min.intron.length = 500) +
  geom_segment(aes(x=47495873, y=1.27, xend=47495873, yend=1.35), size=1.5, color="#f1605dff")+
  geom_segment(aes(x=47495873, y=0.8, xend=47495873, yend=1.27), linetype="dashed",color="#f1605dff")+
  ggtitle(substitute(paste(italic("ELP6"))))+
  annotate("text", x=47495873, y=1.37, label= "m6A") + 
  annotate("text",x=47497873, y=1.3, label="increased methylation in XX lines")+
  labs(x = "Genome Position", y = "Transcript") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(),legend.position = "bottom",
        legend.title = element_blank(),
        strip.text=element_text(size=10),
        panel.spacing.y = unit(0.75, "lines"),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.x  = element_text(vjust = 0.5, hjust=1, size =12),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(vjust=0.5, size=12),
        plot.title = element_text(hjust = 0.5))


HCLS1_annotation <- gff_gencode %>% 
  dplyr::filter( 
    transcript_id == "ENST00000314583" 
  )


HCLS1 <- sex_xpore_filter %>% filter(gene_name.y =="HCLS1")    

HCLS1_exons <- HCLS1_annotation %>% dplyr::filter(type == "exon")
HCLS1_cds <- HCLS1_annotation %>% dplyr::filter(type == "CDS")

HCLS1_transcripts <- ggplot(HCLS1_exons, aes(
  xstart = start,
  xend = end,
  y = transcript_id
)) +
  geom_range(
    fill = "#31678d",color="#31678d",
    height = 0.25
  ) +
  geom_range(
    data = HCLS1_cds, fill = "#31678d",color="#31678d")+
  geom_intron(
    data = to_intron(HCLS1_exons, "transcript_id"),
    aes(strand = strand),
    arrow.min.intron.length = 500) +
  geom_segment(aes(x=121631735, y=1.27, xend=121631735, yend=1.35), size=1.5, color="#f1605dff")+
  geom_segment(aes(x=121631735, y=0.8, xend=121631735, yend=1.27), linetype="dashed",color="#f1605dff")+
  ggtitle(substitute(paste(italic("HCLS1"))))+
  annotate("text", x=121631735, y=1.37, label= "m6A") + 
  annotate("text",x=121635735, y=1.31, label="increased methylation in XX lines")+
  labs(x = "Genome Position", y = "Transcript") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(),legend.position = "bottom",
        legend.title = element_blank(),
        strip.text=element_text(size=10),
        panel.spacing.y = unit(0.75, "lines"),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.x  = element_text(vjust = 0.5, hjust=1, size =12),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(vjust=0.5, size=12),
        plot.title = element_text(hjust = 0.5))

CNDP2_annotation <- gff_gencode %>% 
  dplyr::filter( 
    transcript_id == "ENST00000324301" 
  )


CNDP2 <- sex_xpore_filter %>% filter(gene_name.y =="CNDP2")    

CNDP2_exons <- CNDP2_annotation %>% dplyr::filter(type == "exon")
CNDP2_cds <- CNDP2_annotation %>% dplyr::filter(type == "CDS")

CNDP2_transcripts <- ggplot(CNDP2_exons, aes(
  xstart = start,
  xend = end,
  y = transcript_id
)) +
  geom_range(
    fill = "#31678d",color="#31678d",
    height = 0.25
  ) +
  geom_range(
    data = CNDP2_cds, fill = "#31678d",color="#31678d")+
  geom_intron(
    data = to_intron(CNDP2_exons, "transcript_id"),
    aes(strand = strand),
    arrow.min.intron.length = 500) +
  geom_segment(aes(x=74520721, y=1.27, xend=74520721, yend=1.35), size=1.5, color="#f1605dff")+
  geom_segment(aes(x=74520721, y=0.8, xend=74520721, yend=1.27), linetype="dashed",color="#f1605dff")+
  ggtitle(substitute(paste(italic("CNDP2"))))+
  annotate("text", x=74520721, y=1.37, label= "m6A") + 
  annotate("text",x=74517721, y=1.31, label="increased methylation in XX lines")+
  labs(x = "Genome Position", y = "Transcript") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(),legend.position = "bottom",
        legend.title = element_blank(),
        strip.text=element_text(size=10),
        panel.spacing.y = unit(0.75, "lines"),
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.x  = element_text(vjust = 0.5, hjust=1, size =12),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(vjust=0.5, size=12),
        plot.title = element_text(hjust = 0.5))

volcanos_methylation <- ((xpore_female | xpore_male)
                    /(xpore_day1|xpore_day14))+plot_annotation(tag_levels = "A")+plot_layout(guides="collect")

ggsave('volcanos_methylation.png',volcanos_methylation, width = 20, height = 15, units="in",bg="transparent")


sex_methylation <- ((xpore_sex_comp / MARCHF7_transcripts/ELP6_transcripts/HCLS1_transcripts)
                      |(CNDP2_transcripts/TMIGD3_transcripts/FAM104A_transcripts))+plot_annotation(tag_levels = "A")

ggsave('sex_methylation.png',sex_methylation, width = 25, height = 17, units="in",bg="transparent")
