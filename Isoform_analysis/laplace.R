#Laplace entropy analysis

library(SplicingFactory)
library(treemap)

barcodes$condition <- ifelse(barcodes$condition =="F_D14","XX D14", ## adapt as needed
                                 ifelse(barcodes$condition =="F_D1","XX D1", ## adapt as needed
                                        ifelse(barcodes$condition =="M_D14","XY D14", ## adapt as needed
                                               ifelse(barcodes$condition =="M_D1","XY D1",NA)))) ## adapt as needed


diversity <- talon_filter[,c(7,49:57,59:61)] ## adapt as needed
cols <- colnames(diversity)
diversity <- na.omit(diversity)

genes <- diversity[, 1]
readcounts <- diversity[, -1]

tokeep <- rowSums(readcounts > 5) > 5
readcounts <- readcounts[tokeep, ]
genes      <- genes[tokeep]

laplace_entropy_flair <- calculate_diversity(readcounts, genes,  method = "laplace",
                                       norm = TRUE, verbose = TRUE)


laplace_entropy_talon <- calculate_diversity(readcounts, genes,  method = "laplace",
                                             norm = TRUE, verbose = TRUE)


laplace_entropy_isoseq <- calculate_diversity(readcounts, genes,  method = "laplace",
                                             norm = TRUE, verbose = TRUE)

laplace_data_flair <- cbind(assay(laplace_entropy_flair),
                      Gene = rowData(laplace_entropy_flair)$genes)

# Reshape data.frame
laplace_data_flair <- pivot_longer(laplace_data_flair, -Gene, names_to = "sample",
                             values_to = "entropy")

# Add sample type information
laplace_data_flair$sample_type <- barcodes$condition[match(laplace_data_flair$sample,barcodes$sample_id)]

# Filter genes with NA entropy values
laplace_data_flair <- drop_na(laplace_data_flair)

laplace_data_flair$Gene <- paste0(laplace_data_flair$Gene, "_", laplace_data_flair$sample_type)
laplace_data_flair$diversity <-  "Normalized Laplace entropy"
laplace_data_flair$method <- "FLAIR"
laplace_data_flair$sample <- paste0(laplace_data_flair$method, "_", laplace_data_flair$sample)


laplace_entropy <- calculate_diversity(readcounts, genes,  method = "laplace",
                                       norm = TRUE, verbose = TRUE)

laplace_data_talon <- cbind(assay(laplace_entropy_talon),
                          Gene = rowData(laplace_entropy_talon)$genes)

# Reshape data.frame
laplace_data_talon <- pivot_longer(laplace_data_talon, -Gene, names_to = "sample",
                                 values_to = "entropy")

# Add sample type information
laplace_data_talon$sample_type <- barcodes$condition[match(laplace_data_talon$sample,barcodes$sample_id)]

# Filter genes with NA entropy values
laplace_data_talon <- drop_na(laplace_data_talon)

laplace_data_talon$Gene <- paste0(laplace_data_talon$Gene, "_", laplace_data_talon$sample_type)
laplace_data_talon$diversity <-  "Normalized Laplace entropy"
laplace_data_talon$method <- "TALON"
laplace_data_talon$sample <- paste0(laplace_data_talon$method, "_", laplace_data_talon$sample)

laplace_entropy <- calculate_diversity(readcounts, genes,  method = "laplace",
                                       norm = TRUE, verbose = TRUE)

laplace_data_iso <- cbind(assay(laplace_entropy_isoseq),
                      Gene = rowData(laplace_entropy_isoseq)$genes)

# Reshape data.frame
laplace_data_iso <- pivot_longer(laplace_data_iso, -Gene, names_to = "sample",
                             values_to = "entropy")

# Add sample type information
laplace_data_iso$sample_type <- barcodes$condition[match(laplace_data_iso$sample,barcodes$sample_id)]

# Filter genes with NA entropy values
laplace_data_iso <- drop_na(laplace_data_iso)

laplace_data_iso$Gene <- paste0(laplace_data_iso$Gene, "_", laplace_data_iso$sample_type)
laplace_data_iso$diversity <-  "Normalized Laplace entropy"
laplace_data_iso$method <- "Isoseq"
laplace_data_iso$sample <- paste0(laplace_data_iso$method, "_", laplace_data_iso$sample)

laplace_data <- rbind(laplace_data_flair,laplace_data_iso,laplace_data_talon)


colors_manual <- c("FLAIR" = "skyblue", "Isoseq" = "#990099", "TALON" = "#677f51")

# Plot diversity data
gene_div <- ggplot() +
  geom_density(data = laplace_data, alpha = 0.3,
               aes(x = entropy, group = sample, color = method)) +
  scale_color_manual(values=c("FLAIR" = "skyblue", "Isoseq" = "#990099", "TALON" = "#677f51")) +
  labs(x = "Diversity values",
       y = "Density", color=element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  scale_x_continuous(expand=c(0,0))+
  facet_wrap(~sample_type,strip.position = "top")+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"),
        legend.position = "bottom",
          strip.text=element_text(size=16),
          panel.spacing.y = unit(0.75, "lines"),
          legend.text = element_text(size=16),
          axis.title.x = element_text(size=16),
          axis.text.x  = element_text(vjust = 0.5, hjust=1, size =16),
          axis.title.y = element_text(size=16),
          axis.text.y  = element_text(vjust=0.5, size=16))


# Mean entropy calculation across samples for each gene/sample type combination
laplace_entropy_mean <- aggregate(laplace_data$entropy,
                                  by = list(laplace_data$Gene,laplace_data$method), mean)
colnames(laplace_entropy_mean)[2] <- "method"
colnames(laplace_entropy_mean)[3] <- "mean_entropy"

laplace_entropy_mean <- as_tibble(laplace_entropy_mean)

# Add sample type information
laplace_entropy_mean <- laplace_entropy_mean %>%
  mutate(sample_type = case_when(
    grepl("_XX D14", .[[1]]) ~ "XX D14", ## adapt as needed
    grepl("_XX D1", .[[1]]) ~ "XX D1", ## adapt as needed
    grepl("_XY D14", .[[1]]) ~ "XY D14", ## adapt as needed
    grepl("_XY D1", .[[1]]) ~ "XY D1", ## adapt as needed
    TRUE ~ NA_character_
  ))
# Add diversity type column
laplace_entropy_mean$diversity <-  "Normalized Laplace entropy"

mean_div <- ggplot() +
  geom_violin(data = laplace_entropy_mean, aes(x = sample_type, y = mean_entropy,
                                               fill = method),
              alpha = 0.6) +
  scale_fill_manual(values=colors_manual) +
  coord_flip() +
  labs(x = "Samples",
       y = "Diversity",
       fill=element_blank())+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"),
        legend.position = "bottom",
        strip.text=element_text(size=16),
        panel.spacing.y = unit(0.75, "lines"),
        legend.text = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x  = element_text(vjust = 0.5, hjust=1, size =16),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=16))

## adapt as needed
colData(laplace_entropy_flair) <- cbind(colData(laplace_entropy_flair),
                                        sample_type = ifelse(laplace_entropy_flair$samples %in% c("007D1", "069D1", "SCTID1"),# "BC13","BC15","BC17","BC37","BC39","BC41"),
                                                             "XX",
                                                             ifelse(laplace_entropy_flair$samples %in% c("007D14", "069D14", "SCTID14"),#"BC14","BC16","BC18","BC38","BC40","BC42"),
                                                                    "XX",ifelse(laplace_entropy_flair$samples %in% c("36SD1", "127D1", "014D1"),#"BC25","BC27","BC29","BC31","BC33","BC35"),
                                                                                 "XY","XY"))),
                                  method="FLAIR")

## adapt as needed
colData(laplace_entropy_talon) <- cbind(colData(laplace_entropy_talon),
                                        sample_type = ifelse(laplace_entropy_talon$samples %in% c("007D1", "069D1", "SCTID1"),# "BC13","BC15","BC17","BC37","BC39","BC41"),
                                                             "XX",
                                                             ifelse(laplace_entropy_talon$samples %in% c("007D14", "069D14", "SCTID14"),#"BC14","BC16","BC18","BC38","BC40","BC42"),
                                                                    "XX",ifelse(laplace_entropy_talon$samples %in% c("36SD1", "127D1", "014D1"),#"BC25","BC27","BC29","BC31","BC33","BC35"),
                                                                                 "XY","XY"))),
                                        method="TALON")

## adapt as needed
colData(laplace_entropy_isoseq) <- cbind(colData(laplace_entropy_isoseq),
                                        sample_type = ifelse(laplace_entropy_isoseq$samples %in% c("007D1", "069D1", "SCTID1"),# "BC13","BC15","BC17","BC37","BC39","BC41"),
                                                             "XX",
                                                             ifelse(laplace_entropy_isoseq$samples %in% c("007D14", "069D14", "SCTID14"),#"BC14","BC16","BC18","BC38","BC40","BC42"),
                                                                    "XX",ifelse(laplace_entropy_isoseq$samples %in% c("36SD1", "127D1", "014D1"),#"BC25","BC27","BC29","BC31","BC33","BC35"),
                                                                                "XY","XY"))),
                                        method="Isoseq")


# Calculate significant entropy changes
entropy_significance_time <- calculate_difference(x = laplace_entropy_isoseq, samples = "sample_type",
                                             control = "D1", ## adapt as needed
                                             method = "mean", test = "wilcoxon",
                                             verbose = TRUE)


entropy_significance_time$label <- apply(entropy_significance_time[, c(4, 7)], 1,
                                    function(x) ifelse(abs(x[1]) >= 0.1 & x[2] < 0.05,
                                                       "significant", "non-significant"))

entropy_significance_time$mean <- apply(entropy_significance_time[, c(2, 3)], 1,
                                   function(x) (x[1] + x[2]) / 2)

entropy_significance_time$delabel <- entropy_significance_time$genes
entropy_significance_time <- entropy_significance_time[order(entropy_significance_time$adjusted_p_values), ] 
entropy_significance_time$genelabels <- ""
entropy_significance_time$genelabels[1:10] <- entropy_significance_time$genes[1:10]

entropy_significance_time %>% write.csv("isoseq_laplace_sex.csv")

entropy_significance_time_flair <- entropy_significance_time
entropy_significance_time_flair$method <- "FLAIR"

entropy_significance_time_talon <- entropy_significance_time
entropy_significance_time_talon$method <- "TALON"

entropy_significance_time_isoseq <- entropy_significance_time
entropy_significance_time_isoseq$method <- "Isoseq"

entropy_significance_time <- rbind(entropy_significance_time_flair,entropy_significance_time_isoseq,entropy_significance_time_talon)
entropy_significance_sex <- rbind(entropy_significance_time_flair,entropy_significance_time_isoseq,entropy_significance_time_talon)


sex <-ggplot(entropy_significance_sex, aes(x = mean_difference,
                                                   y = -log10(adjusted_p_values),
                                                   color = label,
                                                   label=genelabels,
                                                  shape=method)) +
  geom_point() +
  geom_text_repel(aes(label = genelabels),show_guide  = FALSE) +
  scale_color_manual(values = c("grey", "red"), guide = "none") +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_vline(xintercept = c(0.1, -0.1)) +
  theme_minimal() +
  labs(title = "Normalized Laplace entropy (XY vs XX)",
       x = "Mean difference of entropy values",
       y = "-Log10(adjusted p-value)",shape=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 16),
        panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"),
        legend.text = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=16),
        legend.position = "bottom"
  )



time <-ggplot(entropy_significance_time, aes(x = mean_difference,
                                 y = -log10(adjusted_p_values),
                                 color = label,
                                 label=genelabels,shape=method)) +
  geom_point() +
  geom_text_repel(aes(label = genelabels),show_guide  = FALSE) +
  scale_color_manual(values = c("grey", "red"), guide = "none") +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_vline(xintercept = c(0.1, -0.1)) +
  theme_minimal() +
  labs(title = "Normalized Laplace entropy (Day 1 vs Day 14)",
       x = "Mean difference of entropy values",
       y = "-Log10(adjusted p-value)",shape=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 16),
        panel.grid = element_blank(), axis.line = element_line(color ="#5f5f5f"),
        legend.text = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(vjust=0.5, size=16),
        legend.position = "bottom"
  )

diversity_plot <- ((gene_div | mean_div) /
                    (time |  sex)) + plot_annotation(tag_levels = 'A')

ggsave('combined_diversity.png',diversity_plot, width = 20, height = 14, bg="transparent")


