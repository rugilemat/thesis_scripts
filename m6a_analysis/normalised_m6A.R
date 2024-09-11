sites_modified <- annotated_modified %>% 
  group_by(gene_id, genome_pos, transcript_id, position) %>% 
  dplyr::filter(length(unique(sample_id))>1) %>% 
  ungroup() %>%
  group_by(gene_id, gene_name,transcript_id)

unique_modified <- subset(sites_modified, !(genome_pos %in% m6a_atlas$start))
unique_modified <- subset(unique_modified, !(genome_pos %in% directRMDB$start))


unique_modified_collapse <- unique_modified %>%
  group_by(genome_pos, gene_id,transcript_position) %>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ",")))


m6aPerGene <- annotated_modified %>% 
  dplyr::group_by(gene_id, genome_pos) %>% 
  dplyr::filter(length(unique(sample_id))>1) %>% 
  dplyr::ungroup() %>%
  dplyr::group_by(gene_id, gene_name) %>%
  dplyr::summarise(site_count = n_distinct(genome_pos))

m6aPerGene_unmodified <- unmodified_collapse %>% 
  group_by(gene_id, transcript_id, transcript_position) %>% 
  dplyr::filter(length(unique(sample_id))>1) %>% 
  ungroup() %>%
  group_by(gene_id) %>%
  summarise(site_count = n_distinct(transcript_position))
m6aPerGene_unmodified <- m6aPerGene_unmodified %>% filter(!(gene_id %in% m6aPerGene$gene_id))

m6aPerGene_unmodified$gene_name <- gff_gencode$gene_name[match(m6aPerGene_unmodified$gene_id,gff_gencode$gene_id)]
m6aPerGene_unmodified$nIso <- isogene_combination$mean_iso[match(m6aPerGene_unmodified$gene_name,isogene_combination$associatedGene)]
m6aPerGene$nIso <- isogene_combination$mean_iso[match(m6aPerGene$gene_name,isogene_combination$associatedGene)]

m6aPerTranscript_sums <- annotated_modified %>% 
  group_by(gene_id, genome_pos, transcript_id, transcript_position) %>% 
  dplyr::filter(length(unique(sample_id))>1) %>% 
  ungroup() %>%
  group_by(gene_id, gene_name,transcript_id, sample_id) %>%
  summarise(sum_mod_ratio = sum(mod_ratio))
m6aPerTranscript$modified_status <- "modified"

m6aPerTranscript_sums <- annotated_modified %>% 
  dplyr::group_by(gene_id, genome_pos, transcript_id, transcript_position) %>% 
  dplyr::filter(length(unique(sample_id))>1) %>% 
  dplyr::ungroup() %>%
  dplyr::group_by(gene_name,gene_id,transcript_id,sample_id) %>% 
  dplyr::summarise(sum_mod_ratio = sum(mod_ratio))


m6aPerTranscript_sums$length <- m6aPerTranscript_filtered_normalised$length[match(m6aPerTranscript_sums$transcript_id,m6aPerTranscript_filtered_normalised$transcript_id)]
m6aPerTranscript_sums$exons <- m6aPerTranscript_filtered_normalised$exons[match(m6aPerTranscript_sums$transcript_id,m6aPerTranscript_filtered_normalised$transcript_id)]
m6aPerTranscript_sums$density <- m6aPerTranscript_sums$length/m6aPerTranscript_sums$exons
m6aPerTranscript_sums <- drop_na(m6aPerTranscript_sums)

lm_unmod=lm(sum_mod_ratio ~ length + density, data = m6aPerTranscript_sums)
m6aPerTranscript_sums$normalised <- resid(lm_unmod)

m6aPerTranscript_sums <- m6aPerTranscript_sums %>% 
  dplyr::group_by(gene_name,gene_id,transcript_id,length,exons,density) %>% 
  dplyr::summarise(mean_norm = mean(normalised))

m6aPerTranscript_hyper <- m6aPerTranscript_sums %>% filter(mean_norm >2.5)

m6aPerTranscript_umodified <- unmodified_collapse %>% 
  group_by(gene_id, transcript_id, transcript_position) %>% 
  dplyr::filter(length(unique(sample_id))>1) %>% 
  ungroup() %>%
  group_by(gene_id, transcript_id) %>%
  summarise(site_count = n_distinct(transcript_position), sum_mod_ratio = sum(mod_ratio))


m6aPerTranscript_umodified$length <- direct_tech_filter$length[match(m6aPerTranscript_umodified$transcript_id,direct_tech_filter$isoform)]
m6aPerTranscript_umodified$exons <- direct_tech_filter$exons[match(m6aPerTranscript_umodified$transcript_id,direct_tech_filter$isoform)]
m6aPerTranscript_umodified$density <- m6aPerTranscript_umodified$length/m6aPerTranscript_umodified$exons
m6aPerTranscript_umodified_filtered <- drop_na(m6aPerTranscript_umodified)
m6aPerTranscript_umodified$gene_name <- gff_gencode$gene_name[match(m6aPerTranscript_umodified$gene_id,gff_gencode$gene_id)]


lm_unmod=lm(sum_mod_ratio ~ length + density, data = m6aPerTranscript_umodified_filtered)
m6aPerTranscript_umodified_filtered$normalised <- resid(lm_unmod)

m6aPerTranscript_umodified_filtered$UTR_3 <- direct_tech_filter$UTR3_length[match(m6aPerTranscript_umodified_filtered$transcript_id,direct_tech_filter$associated_transcript)]
m6aPerTranscript_umodified_filtered$UTR_5 <- direct_tech_filter$UTR5_length[match(m6aPerTranscript_umodified_filtered$transcript_id,direct_tech_filter$associated_transcript)]
m6aPerTranscript_umodified_filtered$polyA <- polyA_list$polyA[match(m6aPerTranscript_umodified_filtered$transcript_id,polyA_list_filtered$transcript_id)]
m6aPerTranscript_umodified_filtered$modified_status <- "unmodified"
m6aPerTranscript_umodified_filtered$gene_name <- gff_gencode$gene_name[match(m6aPerTranscript_umodified_filtered$gene_id,gff_gencode$gene_id)]

direct_talon_filter <- separate(direct_talon_filter, col=isoform, into=c("isoform", NULL), sep="\\.", extra="drop")
direct_tech_filter <- separate(direct_tech_filter, col=associated_transcript, into=c("associated_transcript", NULL), sep="\\.", extra="drop")
direct_tech_filter$isoform <- ifelse(direct_tech_filter$structural_category=="full-splice_match",direct_tech_filter$associated_transcript, direct_tech_filter$isoform)
                                     
m6aPerTranscript$length <- annotated_modified$transcript_length[match(m6aPerTranscript$transcript_id,annotated_modified$transcript_id)]
m6aPerTranscript$exons <- ifelse(m6aPerTranscript$transcript_id %in% direct_talon_filter$isoform, direct_talon_filter$exons[match(m6aPerTranscript$transcript_id,direct_talon_filter$isoform)],
                               ifelse(m6aPerTranscript$transcript_id %in% direct_tech_filter$isoform, direct_tech_filter$exons[match(m6aPerTranscript$transcript_id,direct_tech_filter$isoform)],NA))
m6aPerTranscript$density <- m6aPerTranscript$length/m6aPerTranscript$exons
m6aPerTranscript_filtered <- drop_na(m6aPerTranscript)
lmTemp=lm(sum_mod_ratio ~ length + density, data = m6aPerTranscript)
m6aPerTranscript_filtered$normalised <- resid(lmTemp)

m6aPerTranscript_filtered$UTR_3 <- direct_tech_filter$UTR3_length[match(m6aPerTranscript_filtered$transcript_id,direct_tech_filter$associated_transcript)]
m6aPerTranscript_filtered$UTR_5 <- direct_tech_filter$UTR5_length[match(m6aPerTranscript_filtered$transcript_id,direct_tech_filter$associated_transcript)]
m6aPerTranscript$polyA <- polyA_list$polyA[match(m6aPerTranscript$transcript_id,polyA_list$transcript_id)]
write.csv(m6aPerTranscript_filtered,"m6aPerTranscript_filtered_normalised.csv",row.names = FALSE, quote=FALSE)
write.csv(m6aPerTranscript_umodified_filtered,"m6aPerTranscript_unmod_filtered_normalised.csv",row.names = FALSE, quote=FALSE)

annotated_modified_d14_f <- annotated_modified %>% filter(group_id=="D14_F")
annotated_modified_d1_f <- annotated_modified %>% filter(group_id=="D1_F")
annotated_modified_d14_m <- annotated_modified %>% filter(group_id=="D14_M")
annotated_modified_d1_m <- annotated_modified %>% filter(group_id=="D1_M")


all_sites_unfiltered_d14_f <- all_sites_unfiltered %>% filter(group_id=="D14_F")
all_sites_unfiltered_d1_f <- all_sites_unfiltered %>% filter(group_id=="D1_F")
all_sites_unfiltered_d14_m <- all_sites_unfiltered %>% filter(group_id=="D14_M")
all_sites_unfiltered_d1_m <- all_sites_unfiltered %>% filter(group_id=="D1_M")

m6aPerTranscript_filtered$d14_f <- ifelse(m6aPerTranscript_filtered$transcript_id %in% annotated_modified_d14_f$transcript_id, "Yes","-")
m6aPerTranscript_filtered$d1_f <- ifelse(m6aPerTranscript_filtered$transcript_id %in% annotated_modified_d1_f$transcript_id, "Yes","-")
m6aPerTranscript_filtered$d14_m <- ifelse(m6aPerTranscript_filtered$transcript_id %in% annotated_modified_d14_m$transcript_id, "Yes","-")
m6aPerTranscript_filtered$d1_m <- ifelse(m6aPerTranscript_filtered$transcript_id %in% annotated_modified_d1_m$transcript_id, "Yes","-")

m6aPerTranscript_umodified_filtered$d14_f <- ifelse(m6aPerTranscript_umodified_filtered$transcript_id %in% all_sites_unfiltered_d14_f$transcript_id, "Yes","-")
m6aPerTranscript_umodified_filtered$d1_f <- ifelse(m6aPerTranscript_umodified_filtered$transcript_id %in% all_sites_unfiltered_d1_f$transcript_id, "Yes","-")
m6aPerTranscript_umodified_filtered$d14_m <- ifelse(m6aPerTranscript_umodified_filtered$transcript_id %in% all_sites_unfiltered_d14_m$transcript_id, "Yes","-")
m6aPerTranscript_umodified_filtered$d1_m <- ifelse(m6aPerTranscript_umodified_filtered$transcript_id %in% all_sites_unfiltered_d1_m$transcript_id, "Yes","-")


m6aPerTranscript$kbp <- (m6aPerTranscript$length)/1000
m6aPerTranscript$den_norm <- (m6aPerTranscript$density)/1000


modified_isoforms <- rbind(m6aPerTranscript_filtered,m6aPerTranscript_umodified_filtered)
modified_isoforms$kbp <- (modified_isoforms$length)/1000
modified_isoforms$den_norm <- (modified_isoforms$density)/1000

modified_gene <- rbind(m6aPerGene,m6aPerGene_unmodified)
