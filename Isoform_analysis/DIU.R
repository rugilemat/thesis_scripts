library(genomation)
library(rtracklayer)
library(IsoformSwitchAnalyzeR)
library("BSgenome.Hsapiens.UCSC.hg38")

### Adapt as needed, use count and gtf files from SQANTI classification

flair_counts <- flair_counts[,c(1,49:84)]

flair_counts <- flair_counts %>% 
  distinct(.keep_all = TRUE)
flair_counts <- na.omit(flair_counts)
flair_counts <- column_to_rownames(flair_counts, var="isoform")


talon_counts <- talon_counts[,c(1,49:84)]

talon_counts <- talon_counts %>% 
  distinct(.keep_all = TRUE)
talon_counts <- na.omit(talon_counts)
talon_counts <- column_to_rownames(talon_counts, var="isoform")

isoseq_counts <- isoseq_counts[,c(1,49:84)]

isoseq_counts <- isoseq_counts %>% 
  distinct(.keep_all = TRUE)
isoseq_counts <- na.omit(isoseq_counts)
isoseq_counts <- column_to_rownames(isoseq_counts, var="isoform")

### Adapt as needed
myDesign <- data.frame(
    sampleID = colnames(counts),
    condition =   factor(c("D1","D14","D1","D14","D1","D14",
                           "D1","D14","D1","D14","D1","D14",
                           "D1","D14","D1","D14","D1","D14",
                           "D1","D14","D1","D14","D1","D14",
                           "D1","D14","D1","D14","D1","D14",
                           "D1","D14","D1","D14","D1","D14")))
  
  myDesign$sex <- factor(c("F","F","F","F","F","F",
                           "F","F","F","F","F","F",
                           "M","M","M","M","M","M",
                           "M","M","M","M","M","M",
                           "M","M","M","M","M","M",
                           "F","F","F","F","F","F"))
  
  myDesign$rep <- factor(c("B1","B1","B2","B2","B3","B3",
                           "B1","B1","B2","B2","B3","B3",
                           "B1","B1","B2","B2","B3","B3",
                           "B1","B1","B2","B2","B3","B3",
                           "B1","B1","B2","B2","B3","B3",
                           "B1","B1","B2","B2","B3","B3"))

### Create switchAnalyzeRlist FLAIR --------------
flairSwitchList <- importRdata(
  isoformCountMatrix   = flair,
  designMatrix         = myDesign,
  isoformExonAnnoation = "~/path/to/gtf/file", removeNonConvensionalChr=TRUE, ignoreAfterSpace = TRUE,
  isoformNtFasta       = "~/path/to/fasta/file")

flair_mgl_list <- cbind(flairSwitchList$isoformFeatures$isoform_id,  flairSwitchList$isoformFeatures$gene_id)
write.table(flair_list, 'flair_list.txt', sep = '\t', quote = FALSE)

flairSwitchList <- preFilter(
  switchAnalyzeRlist = flairSwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE)

flairSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = flairSwitchList,
  reduceToSwitchingGenes=FALSE
)

bsg <- BSgenome.Hsapiens.UCSC.hg38

flairSwitchListAnalyzed <- addORFfromGTF(flairSwitchListAnalyzed,"~/path/to/gtf/file", overwriteExistingORF=TRUE)
flairSwitchListAnalyzed <- analyzeNovelIsoformORF(flairSwitchListAnalyzed, genomeObject = bsg, analysisAllIsoformsWithoutORF = TRUE)
flairSwitchListAnalyzed <- analyzeORF(flairSwitchListAnalyzed, genomeObject = bsg)

flairSwitchListAnalyzed <- extractSequence(
  flairSwitchListAnalyzed,
  removeLongAAseq=FALSE,
  alsoSplitFastaFile=FALSE,
  # onlySwitchingGenes=FALSE,
  pathToOutput = "~/path/to/output",
  writeToFile=FALSE # to avoid output when running this example data
)


flairSwitchListAnalyzed <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = flairSwitchListAnalyzed, 
  #  dIFcutoff                 = 0.1,   # Cutoff for defining switch size 
  n                         = 10,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = FALSE,  # Because ORF was predicted de novo
  pathToCPC2resultFile = "~/path/to/flair/CPC2",
  pathToPFAMresultFile = "~/path/to/flair/PFAM",
  pathToSignalPresultFile = "~/path/to/flair/Signal",
  pathToIUPred2AresultFile = "~/path/to/flair/IUPred2",
  pathToDeepTMHMMresultFile = "~/path/to/flair/DeepTMHMM",
  pathToDeepLoc2resultFile = "~/path/to/flair/DeepLoc2",
  outputPlots               = FALSE,
  quiet = FALSE)


gene_flair_details_all <- extractTopSwitches(flairSwitchListAnalyzed, filterForConsequences = FALSE, n=711)
write.table(gene_flair_details_all, 'gene_time_flair_details.txt', quote = FALSE)

extractSwitchSummary(flairSwitchListAnalyzed, filterForConsequences = TRUE)


flair_time_summary <- extractConsequenceSummary(
  flairSwitchListAnalyzed,
  consequencesToAnalyze='all',
  plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = FALSE,
  returnResult = TRUE    # enables analysis of fraction of significant features
)
flair_time_summary$method <- "FLAIR"



flair_conc <- extractConsequenceEnrichment(
  flairSwitchListAnalyzed,
  consequencesToAnalyze='all',
  analysisOppositeConsequence = TRUE,
  returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
)

write.table(flair_conc, 'flair_conc_time.txt', sep = '\t', quote = FALSE)

flair_time_summary_splice <- extractSplicingSummary(
  flairSwitchListAnalyzed,
  asFractionTotal = FALSE,
  plotGenes=FALSE,
  returnResult = TRUE
)
flair_time_summary_splice$method <- "FLAIR"

splicingflairEnrichment <- extractSplicingEnrichment(
  flairSwitchListAnalyzed,
  splicingToAnalyze='all',
  returnResult=TRUE,
  returnSummary=TRUE
)

write.table(splicingflairEnrichment, 'DTU_SplicingEnrichment_flair_time.txt', sep = '\t', quote = FALSE)

flair_consequences <- extractConsequenceGenomeWide(flairSwitchListAnalyzed,  featureToExtract = 'all') 
write.table(flair_mgl_consequences, 'consequences_flair_time.txt', sep = '\t', quote = FALSE)

flair_splicing <- extractSplicingGenomeWide(flairSwitchListAnalyzed)
write.table(flair_splicing, 'splicing_flair_time.txt', sep = '\t', quote = FALSE)

table(flairSwitchListAnalyzed$AlternativeSplicingAnalysis$ES)
table(flairSwitchListAnalyzed$AlternativeSplicingAnalysis$MES)
table(flairSwitchListAnalyzed$AlternativeSplicingAnalysis$MEE)
table(flairSwitchListAnalyzed$AlternativeSplicingAnalysis$IR)
table(flairSwitchListAnalyzed$AlternativeSplicingAnalysis$A5)
table(flairSwitchListAnalyzed$AlternativeSplicingAnalysis$A3)
table(flairSwitchListAnalyzed$AlternativeSplicingAnalysis$ATSS)
table(flairSwitchListAnalyzed$AlternativeSplicingAnalysis$ATTS)

bioMechanismeAnalysis <- analyzeSwitchConsequences(
  flairSwitchListAnalyzed, 
  consequencesToAnalyze = c('tss','tts','intron_structure'),
  showProgress = FALSE
)$switchConsequence # only the consequences are interesting here
write.table(bioMechanismeAnalysis, 'biomechanism_flair_time.txt', sep = '\t', quote = FALSE)


### Create switchAnalyzeRlist TALON --------------
talonSwitchList <- importRdata(
  isoformCountMatrix   = talon,
  designMatrix         = myDesign,
  isoformExonAnnoation = "~/path/to/gtf/file", removeNonConvensionalChr=TRUE, ignoreAfterSpace = TRUE,
  isoformNtFasta       = "~/path/to/fasta/file")

talon_list <- cbind(talonSwitchList$isoformFeatures$isoform_id,  talonSwitchList$isoformFeatures$gene_id)
write.table(talon_list, 'talon_list.txt', sep = '\t', quote = FALSE)

talonSwitchList <- preFilter(
  switchAnalyzeRlist = talonSwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE)

talonSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = talonSwitchList,
  reduceToSwitchingGenes=FALSE
)

bsg <- BSgenome.Hsapiens.UCSC.hg38

talonSwitchListAnalyzed <- addORFfromGTF(talonSwitchListAnalyzed,"~/path/to/gff", overwriteExistingORF=TRUE)
talonSwitchListAnalyzed <- analyzeNovelIsoformORF(talonSwitchListAnalyzed, genomeObject = bsg, analysisAllIsoformsWithoutORF = TRUE)
talonSwitchListAnalyzed <- analyzeORF(talonSwitchListAnalyzed, genomeObject = bsg)

talonSwitchListAnalyzed <- extractSequence(
  talonSwitchListAnalyzed,
  removeLongAAseq=FALSE,
  alsoSplitFastaFile=FALSE,
  pathToOutput = "~/path/to/output",
  writeToFile=TRUE # to avoid output when running this example data
)

talonSwitchListAnalyzed <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = talonSwitchListAnalyzed, 
  #  dIFcutoff                 = 0.1,   # Cutoff for defining switch size 
  n                         = 10,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = FALSE,  # Because ORF was predicted de novo
  pathToCPC2resultFile = "~/path/to/talon/CPC2",
  pathToPFAMresultFile = "~/path/to/talon/PFAM",
  pathToSignalPresultFile = "~/path/to/talon/Signal",
  pathToIUPred2AresultFile = "~/path/to/talon/IUPred2",
  pathToDeepTMHMMresultFile = "~/path/to/talon/DeepTMHMM",
  pathToDeepLoc2resultFile = "~/path/to/talon/DeepLoc2",
  outputPlots               = FALSE,
  quiet = FALSE)


gene_talon_details_all <- extractTopSwitches(talonSwitchListAnalyzed, filterForConsequences = FALSE, n=711)
write.table(gene_talon_details_all, 'gene_time_talon_details.txt', quote = FALSE)

extractSwitchSummary(talonSwitchListAnalyzed, filterForConsequences = TRUE)

talon_time_summary <- extractConsequenceSummary(
  talonSwitchListAnalyzed,
  consequencesToAnalyze='all',
  plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = FALSE,
  returnResult = TRUE    # enables analysis of fraction of significant features
)
talon_time_summary$method <- "TALON"

talon_conc <- extractConsequenceEnrichment(
  talonSwitchListAnalyzed,
  consequencesToAnalyze='all',
  analysisOppositeConsequence = TRUE,
  returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
)

write.table(talon_conc, 'talon_conc_time.txt', sep = '\t', quote = FALSE)

talon_time_summary_splice <- extractSplicingSummary(
  talonSwitchListAnalyzed,
  asFractionTotal = FALSE,
  plotGenes=FALSE,
  returnResult = TRUE
)
talon_time_summary_splice$method <- "TALON"

splicingtalonEnrichment <- extractSplicingEnrichment(
  talonSwitchListAnalyzed,
  splicingToAnalyze='all',
  returnResult=TRUE,
  returnSummary=TRUE
)

write.table(splicingtalonEnrichment, 'DTU_SplicingEnrichment_talon_time.txt', sep = '\t', quote = FALSE)

talon_consequences <- extractConsequenceGenomeWide(talonSwitchListAnalyzed,  featureToExtract = 'all') 
write.table(talon_consequences, 'consequences_talon_time.txt', sep = '\t', quote = FALSE)

talon_splicing <- extractSplicingGenomeWide(talonSwitchListAnalyzed)
write.table(talon_splicing, 'splicing_talon_time.txt', sep = '\t', quote = FALSE)

table(talonSwitchListAnalyzed$AlternativeSplicingAnalysis$ES)
table(talonSwitchListAnalyzed$AlternativeSplicingAnalysis$MES)
table(talonSwitchListAnalyzed$AlternativeSplicingAnalysis$MEE)
table(talonSwitchListAnalyzed$AlternativeSplicingAnalysis$IR)
table(talonSwitchListAnalyzed$AlternativeSplicingAnalysis$A5)
table(talonSwitchListAnalyzed$AlternativeSplicingAnalysis$A3)
table(talonSwitchListAnalyzed$AlternativeSplicingAnalysis$ATSS)
table(talonSwitchListAnalyzed$AlternativeSplicingAnalysis$ATTS)

bioMechanismeAnalysis <- analyzeSwitchConsequences(
  talonSwitchListAnalyzed, 
  consequencesToAnalyze = c('tss','tts','intron_structure'),
  showProgress = FALSE
)$switchConsequence # only the consequences are interesting here
write.table(bioMechanismeAnalysis, 'biomechanism_talon_time.txt', sep = '\t', quote = FALSE)


### subset to those with differences
bioMechanismeAnalysis <- bioMechanismeAnalysis[which(bioMechanismeAnalysis$isoformsDifferent),]



### Create switchAnalyzeRlist Isoseq --------------
isoseqSwitchList <- importRdata(
  isoformCountMatrix   = isoseq,
  designMatrix         = myDesign,
  isoformExonAnnoation = "~/path/to/gff/file", removeNonConvensionalChr=TRUE, ignoreAfterSpace = TRUE,
  isoformNtFasta       = "~/path/to/fasta/file")

isoseq_list <- cbind(isoseqSwitchList$isoformFeatures$isoform_id,  isoseqSwitchList$isoformFeatures$gene_id)
write.table(isoseq_mgl_list, 'isoseq_list.txt', sep = '\t', quote = FALSE)

isoseqSwitchList <- preFilter(
  switchAnalyzeRlist = isoseqSwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE)

isoseqSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = isoseqSwitchList,
  reduceToSwitchingGenes=FALSE
)

bsg <- BSgenome.Hsapiens.UCSC.hg38

isoseqSwitchListAnalyzed <- analyzeORF(isoseqSwitchListAnalyzed, genomeObject = bsg)

isoseqSwitchListAnalyzed <- extractSequence(
  isoseqSwitchListAnalyzed,
  removeLongAAseq=FALSE,
  alsoSplitFastaFile=FALSE,
  pathToOutput = "~/path/to/output",
  writeToFile=TRUE # to avoid output when running this example data
)

isoseqSwitchListAnalyzed <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = isoseqSwitchListAnalyzed, 
  #  dIFcutoff                 = 0.1,   # Cutoff for defining switch size 
  n                         = 10,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = FALSE,  # Because ORF was predicted de novo
  pathToCPC2resultFile = "~/path/to/isoseq/CPC2",
  pathToPFAMresultFile = "~/path/to/isoseq/PFAM",
  pathToSignalPresultFile = "~/path/to/isoseq/Signal",
  pathToIUPred2AresultFile = "~/path/to/isoseq/IUPred2",
  pathToDeepTMHMMresultFile = "~/path/to/isoseq/DeepTMHMM",
  pathToDeepLoc2resultFile = "~/path/to/isoseq/DeepLoc2",
  outputPlots               = FALSE,
  quiet = FALSE)


gene_isoseq_details_all <- extractTopSwitches(isoseqSwitchListAnalyzed, filterForConsequences = FALSE, n=711)
write.table(gene_isoseq_details_all, 'gene_time_isoseq_details.txt', quote = FALSE)

extractSwitchSummary(isoseqSwitchListAnalyzed, filterForConsequences = TRUE)

isoseq_time_summary <- extractConsequenceSummary(
  isoseqmglSwitchListAnalyzed,
  consequencesToAnalyze='all',
  plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = FALSE,
  returnResult = TRUE    # enables analysis of fraction of significant features
)
isoseq_time_summary$method <- "Isoseq"

isoseq_conc <- extractConsequenceEnrichment(
  isoseqSwitchListAnalyzed,
  consequencesToAnalyze='all',
  analysisOppositeConsequence = TRUE,
  returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
)

write.table(isoseq_conc, 'isoseq_conc_time.txt', sep = '\t', quote = FALSE)


isoseq_time_summary_splice <- extractSplicingSummary(
  isoseqSwitchListAnalyzed,
  asFractionTotal = FALSE,
  plotGenes=FALSE,
  returnResult = TRUE
)
isoseq_time_summary_splice$method <- "Isoseq"


splicingisoseqEnrichment <- extractSplicingEnrichment(
  isoseqSwitchListAnalyzed,
  splicingToAnalyze='all',
  returnResult=TRUE,
  returnSummary=TRUE
)

write.table(splicingisoseqEnrichment, 'DTU_SplicingEnrichment_isoseq_time.txt', sep = '\t', quote = FALSE)

isoseq_consequences <- extractConsequenceGenomeWide(isoseqSwitchListAnalyzed,  featureToExtract = 'all') 
write.table(isoseq_mgl_consequences, 'consequences_isoseq_time_mgl.txt', sep = '\t', quote = FALSE)

isoseq_splicing <- extractSplicingGenomeWide(isoseqSwitchListAnalyzed)
write.table(isoseq_splicing, 'splicing_isoseq_time.txt', sep = '\t', quote = FALSE)

table(isoseqSwitchListAnalyzed$AlternativeSplicingAnalysis$ES)
table(isoseqSwitchListAnalyzed$AlternativeSplicingAnalysis$MES)
table(isoseqSwitchListAnalyzed$AlternativeSplicingAnalysis$MEE)
table(isoseqSwitchListAnalyzed$AlternativeSplicingAnalysis$IR)
table(isoseqSwitchListAnalyzed$AlternativeSplicingAnalysis$A5)
table(isoseqSwitchListAnalyzed$AlternativeSplicingAnalysis$A3)
table(isoseqSwitchListAnalyzed$AlternativeSplicingAnalysis$ATSS)
table(isoseqSwitchListAnalyzed$AlternativeSplicingAnalysis$ATTS)

bioMechanismeAnalysis <- analyzeSwitchConsequences(
  isoseqSwitchListAnalyzed, 
  consequencesToAnalyze = c('tss','tts','intron_structure'),
  showProgress = FALSE
)$switchConsequence # only the consequences are interesting here
write.table(bioMechanismeAnalysis, 'biomechanism_isoseq_time_mgl.txt', sep = '\t', quote = FALSE)




switch_cons <- rbind(isoseq_time_summary,talon_time_summary,flair_time_summary)
splicing_cons <- rbind(isoseq_time_summary_splice,talon_time_summary_splice,flair_time_summary_splice)

