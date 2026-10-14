##### DADA2 pipeline #####

# Author: Eva Tanneau on the 05.03.2024
# inspired by Janine Haueisen DADA2 script
# from raw data to creation of the phyloseq object

# load packages
library(dada2)           # version 1.38.0
library(Rcpp)            # version 1.1.1
library("BiocParallel")  # version 1.44.0
library(ShortRead)       # version 1.68.0
library(Biostrings)      # version 2.78.0
library(phyloseq)        # version 1.54.0
library(stringr)         # version 1.6.0
library(pander)          # version 0.6.6


# raw reads available in: https://www.ebi.ac.uk/ena/browser/view/PRJEB127283
# Download into a new folder "Data/raw_reads"

# read data
path.1 <- "Data/raw_reads"

fnFs.1 <- sort(list.files(path.1, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs.1 <- sort(list.files(path.1, pattern = "_R2_001.fastq.gz", full.names = TRUE))


# create directory for filtering
fnFs.1.filtR <- file.path("Data/raw_reads/filtR", 
                          basename(fnFs.1)) 
fnRs.1.filtR <- file.path("Data/raw_reads/filtR", 
                          basename(fnRs.1))

# filter
filterAndTrim(fnFs.1, fnFs.1.filtR, fnRs.1, fnRs.1.filtR, matchIDs=TRUE, maxN=0, multithread = FALSE)

# plot quality profile
plotQualityProfile(fnFs.1.filtR[1:10])
plotQualityProfile(fnRs.1.filtR[1:10])

# plotQualityProfile(fnFs.1.filtR, aggregate = T)
# plotQualityProfile(fnRs.1.filtR, aggregate = T)


# ITS primer removal
FWD <- "TCACGCGCAGGCTTGGTCATTTAGAGGAAGTAA"  ## change accordingly
REV <- "CGTACTGTGGAGAGCTGCGTTCTTCATCGATGC"  ## change accordingly

# To ensure we have the right primers, and the correct orientation of the primers on the reads, we will verify the presence and orientation of these primers in the data.
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# Count how often primers appear in the forward and reverse reads - considering all possible primer orientations - for a subset of files only.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.1.filtR[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.1.filtR[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.1.filtR[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.1.filtR[[1]]))

# Download and provide the path to cutadapt: https://github.com/marcelm/cutadapt/releases
# For UNIX
# cutadapt <- "cutadapt/bin/cutadapt"
# For Windows
cutadapt <- "cutadapt.exe"

system2(cutadapt, args = "--version") # Run shell commands from R

path.1.cut <- file.path(path.1, "cutadapt_2")
if(!dir.exists(path.1.cut)) dir.create(path.1.cut)
fnFs.1.cut <- file.path(path.1.cut, basename(fnFs.1))
fnRs.1.cut <- file.path(path.1.cut, basename(fnRs.1))

# reverse complement
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

R1.flags <- paste("-g", FWD, "-a", REV.RC) 
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# run cutadapt
for(i in seq_along(fnFs.1)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-j", 4, # multithreading
                             "-o", fnFs.1.cut[i], "-p", fnRs.1.cut[i], # output files
                             fnFs.1.filtR[i], fnRs.1.filtR[i])) # input files
}



# check the primer sequences presence/absence after
# before
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.1.filtR), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.1.filtR), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.1.filtR), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.1.filtR))
# after
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.1.cut), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.1.cut), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.1.cut), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.1.cut))

# they have all been removed

# save files as cut format
cutFs.1 <- sort(list.files(path.1.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs.1 <- sort(list.files(path.1.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))


# Extract sample names, assuming filenames have format:
# 20_61_CON_01_L_R2_L001_R1_001.fastq.gz
# we need them later to extract the files that remain after quality filtering
get.sample.name <- function(fname) strsplit(basename(fname), "_L001")[[1]][1]
sample.names.1 <- unname(sapply(cutFs.1, get.sample.name))
head(sample.names.1)


# set path ofr dataset filtering
filtFs.1 <- file.path(path.1.cut, "filtered", basename(cutFs.1))
filtRs.1 <- file.path(path.1.cut, "filtered", basename(cutRs.1))


# quality triming and csv
# trunLen = the position from where we set (Default 0 (no truncation). Truncate reads after truncLen bases. Reads shorter than this are discarded.)
# truncQ = Default 2. Truncate reads at the first instance of a quality score less than or equal to truncQ.
out.1 <- filterAndTrim(cutFs.1, filtFs.1, cutRs.1, filtRs.1, maxN = 0, maxEE = c(2, 2), 
                       truncQ = 2, minLen = 200, maxLen = 400, rm.phix = TRUE, compress = TRUE, multithread = FALSE)


# ==== REMOVE ?????? ===== #
# # less strict filter (test from mols 2025)
# out.1 <- filterAndTrim(cutFs.1, filtFs.1, cutRs.1, filtRs.1, maxN = 0, maxEE = c(2, 2), 
#                        truncQ = 2, minLen = 30, maxLen = 400, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
# ========================

head(out.1)
write.csv(out.1, file="Data/raw_reads/Reads_in-out_quality_filtering.csv", sep = ",")

# check if the file exists after filtering
exists.1 <- file.exists(filtFs.1)


# error learning
errF.1 <- learnErrors(filtFs.1[exists.1], nbases= 1e10, randomize = TRUE, multithread = TRUE)
errR.1 <- learnErrors(filtRs.1[exists.1], nbases= 1e10, randomize = TRUE, multithread = TRUE)

# save data error
saveRDS(errF.1, "Data/raw_reads/error_rate_F_test.rds")
saveRDS(errR.1, "Data/raw_reads/error_rate_R_test.rds")

# plot errors
plotErrors(errF.1, nominalQ=TRUE)
plotErrors(errR.1, nominalQ=TRUE)

# Dereplicate identical reads
derepFs.1 <- derepFastq(filtFs.1[exists.1], verbose=TRUE)
derepRs.1 <- derepFastq(filtRs.1[exists.1], verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs.1) <- sample.names.1[exists.1]
names(derepRs.1) <- sample.names.1[exists.1]

# Sample inference
dadaFs.1 <- dada(derepFs.1, err = errF.1, multithread = TRUE)
saveRDS(dadaFs.1, "Data/raw_reads/dada_F.rds")

dadaRs.1 <- dada(derepRs.1, err = errR.1, multithread = TRUE)
saveRDS(dadaRs.1, "Data/raw_reads/dada_R.rds")


# merge pairs
mergers.1 <- mergePairs(dadaFs.1, derepFs.1, dadaRs.1, derepRs.1, verbose=TRUE)

# Inspect the merger data.frame from the first sample of lib1
head(mergers.1[[1]])

# save pairs
saveRDS(mergers.1, "Data/raw_reads/merged-reads.rds")


# generate sequence table
seqtab.1 <- makeSequenceTable(mergers.1)

# save sequence table
saveRDS(seqtab.1, "Data/raw_reads/seqtab.1.rds")


# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab.1, method="consensus", multithread=TRUE, verbose=TRUE)

# Save final sequence table as .rds file
saveRDS(seqtab.nochim, "Data/seqtab.nonchim.rds")
write.csv(seqtab.nochim, "Data/seqtab.nonchim.csv")


# inspect sequence length distribution
length.dis <- data.frame(table(nchar(getSequences(seqtab.nochim))))

ggplot(length.dis, aes(Var1, Freq)) + geom_col() + ggtitle("Sequence length distribution") + xlab("length in bp") + ylab("frequency") + theme(axis.text.x = element_text(size=5, angle=90))


# assign taxonomy
# Download UNITE database: https://unite.ut.ee
unite.ref <- "Data/sh_general_release_dynamic_29.11.2022.fasta" 

taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

write.csv(taxa.print, "Data/taxonomy_ITS.csv")


# get sample names of all samples that made it through the filters and turn into data.frame
samples.out <- rownames(seqtab.nochim)
samples.out

sam <- rownames(seqtab.nochim)
sam <- gsub(  "_S.*", "", sam)
sam <- gsub(  ".*22Sep", "", sam)
sam <- gsub(  "1117", "22Sep1117", sam)
sam

rownames(seqtab.nochim) <- sam

samples.out.frame <- data.frame(samples.out)
colnames(samples.out.frame) <- c("Sample")

names <- as.data.frame(rownames(seqtab.nochim))
write.csv(names, file="Data/raw_reads/names_samples.csv", row.names = F)
sample.det <- read.csv("Data/metadata.csv")

# merge samples.out and sample detail table - to only keep the details for the files that passed all filters
samdf <- merge(samples.out.frame,sample.det, by = "Sample")
samdf <- sample.det
rownames(samdf) <- sample.det$Sample

# write.csv(samdf, "/home/tanneau/Documents/Hohenschulen/Hohenschulen_revised_datasets/FWHoh_2023_ITS/condition_table_all_libs.csv")
# samdf <- read.csv("/home/tanneau/Documents/Hohenschulen/Hohenschulen_revised_datasets/FWHoh_2023_ITS/condition_table_all_libs.csv")
# turn 61 and 62 in Block from numeric into characters (better for handling with ggplot)
samdf$Plot <- as.character(samdf$Plot)
samdf$Year <- as.character(samdf$Year)

# rownames(seqtab.nochim) <- samples.out

# construct phyloseq object from DADA2 output
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

# save phyloseq object as .rds file
saveRDS(ps, "Data/ps_not_pruned.rds")
ps <- readRDS("Data/ps_not_pruned.rds")
ps


# change the sequence name in ASV number
n_seqs.ITS <- seq(ntaxa(ps))
len_n_seqs.ITS <- nchar(max(n_seqs.ITS))
taxa_names(ps) <- paste0("ASV", formatC(n_seqs.ITS))
taxa_names(ps)



# extract the otu table
mat <- otu_table(ps, taxa_are_rows = FALSE)

write.csv(mat, "Data/ASV_table.csv", row.names = TRUE, quote = FALSE)


# extract taxonomy table
tax <- tax_table(ps)
tax <- as.data.frame(tax)
tax <- rownames_to_column(tax, var = "ASV")

write.table(tax, "Data/taxonomy.txt", quote = FALSE, sep = "\t", row.names = FALSE)


# extract ASV sequence
asv_seq = data.frame(ASV = taxa_names(ps), Sequence = getSequences(seqtab.nochim))

write.table(asv_seq, "Data/ASV_sequences.txt", quote = FALSE, sep = "\t", row.names = FALSE)








##### ——— Data analysis ——— #####
# Authors: Sima Tayyebi Fallah and Marco Alexandre Guerreiro


##### Loading libraries #####
library(data.table)  # version 1.18.0
library(dplyr)       # version 1.1.4
library(ggplot2)     # version 4.0.3 
library(ggpmisc)     # version 0.6.3
library(vegan)       # version 2.7-2
library(ggpubr)      # version 0.6.2
library(reshape2)    # version 1.4.5
library(tidyr)       # version 1.3.2
library(tibble)      # version 3.3.1
library(biomformat)  # version 1.38.3
library(decontam)    # version 1.30.0 
library(ggvenn)      # version 0.1.19 
library(cowplot)     # version 1.2.0
library(betapart)    # version 1.6.1
library(ggtree)      # version 4.0.4 
library(ggtreeExtra) # version 1.20.0
library(tidytree)    # version 0.4.7
library(ape)         # version 5.8-1 
library(wesanderson) # version 0.3.7


##### ——— NGS data analysis ——— #####

##### Load NGS data #####
asv_table = fread("Data/ASV_table.csv")
asv_tax   = fread("Data/taxonomy.txt")
asv_seq   = fread("Data/ASV_sequences.txt", header = T)
metadata  = fread("Data/metadata.csv")



##### Data formatting #####
asv_table = asv_table %>% dplyr::rename(sample = V1)

asv_table$sample = gsub("220900001117-DS[0-9]+_", "", asv_table$sample)
asv_table$sample = gsub("_S[0-9]+", "", asv_table$sample)

##### Decontam #####
neg = subset(metadata, Treatment == "Neg-ctrl")$Sample
neg = ifelse( metadata$Sample %in% neg, TRUE, FALSE)

contam = isContaminant(as.matrix(asv_table[,-1]), neg = neg)

contam_asv = as.data.frame(contam) %>% filter(contaminant == TRUE) %>% rownames()



##### Data filtering ##### 

# Remove negative and positive controls
controls = subset(metadata, Treatment == "Neg-ctrl" | Treatment == "Pos-ctrl")$Sample
asv_table = subset(asv_table, !sample %in% controls)

# Remove contaminant ASVs
asv_table = select(all_of(asv_table), -contam_asv)

# Remove ASVs that are 0 in every sample
asv_0 = as.data.frame(colSums(asv_table[,-1]))
asv_0$ASV = rownames(asv_0)
asv_0 = asv_0 %>% dplyr::rename(reads = 1)
asv_0 = subset(asv_0, reads > 0)$ASV

asv_table = asv_table %>% select(sample, asv_0)

rm(asv_0, controls, contam_asv, neg)



##### Data summary #####  

# ASVs
total_asvs_n = ncol(asv_table[,-"sample"])
cat("ASVs: ", total_asvs_n, "\n")

# reads
total_reads = sum(asv_table[,-"sample"])
cat("reads: ", total_reads, "\n")

# samples
total_samples_n = nrow(asv_table)
cat("Samples: ", total_samples_n)

# Shared and exclusive ASVs per treatment
samples_treated = metadata$Sample[metadata$Treatment == "treated"]
samples_untreated = metadata$Sample[metadata$Treatment == "untreated"]

asvs_treated = asv_table %>% 
  subset(sample %in% samples_treated) %>% 
  select(-sample) %>% 
  select(where(~ sum(.x) > 0)) %>% 
  names()

asvs_untreated = asv_table %>% 
  subset(sample %in% samples_untreated) %>% 
  select(-sample) %>% 
  select(where(~ sum(.x) > 0)) %>% 
  names()

asvs_shared = intersect(asvs_treated, asvs_untreated)
asvs_exclusive_treated = setdiff(asvs_treated, asvs_untreated)
asvs_exclusive_untreated = setdiff(asvs_untreated, asvs_treated)

asvs_shared_reads = asv_table %>% 
  select(asvs_shared) %>% 
  sum() / total_reads * 100

asvs_treated_reads = asv_table %>% 
  select(asvs_exclusive_treated) %>% 
  sum() / total_reads * 100

asvs_untreated_reads = asv_table %>% 
  select(asvs_exclusive_untreated) %>% 
  sum() / total_reads * 100

cat("ASVs shared: ", length(asvs_shared), "(", asvs_shared_reads, "%)\n")
cat("ASVs exclusive treated: ", length(asvs_exclusive_treated), "(", asvs_treated_reads, "%)\n")
cat("ASVs exclusive untreated: ", length(asvs_exclusive_untreated), "(", asvs_untreated_reads, "%)\n")


reads = asv_table %>% 
  mutate(reads = rowSums(across(-sample))) %>%
  select(sample, reads)

asvs = asv_table %>% 
  mutate(asvs = rowSums(across(-sample) > 0)) %>%
  select(sample, asvs)

metadata = merge(metadata, reads, by.x = "Sample", by.y = "sample", all = T)
metadata = merge(metadata, asvs, by.x = "Sample", by.y = "sample", all = T)



##### Rarefication #####
min( rowSums(select(asv_table, -sample)) )

set.seed(123)
asv_rarefied = rrarefy(as.data.frame(select(asv_table, -sample)), sample = 60000)
asv_rarefied = as.data.frame(asv_rarefied)
asv_rarefied$sample = asv_table$sample
asv_rarefied = asv_rarefied %>% dplyr::relocate(sample)

write.csv(asv_rarefied, "Data/ASV_table_rarefied.csv", row.names = FALSE)
asv_rarefied = fread("Data/ASV_table_rarefied.csv")


##### standardize by total ##### 
asv_rel_ab = decostand(select(asv_rarefied, -sample), method = "total")
asv_rel_ab$sample = asv_rarefied$sample
asv_rel_ab = merge(asv_rel_ab, select(metadata, c(Sample, Treatment, Year)), by.x = "sample", by.y = "Sample")

asv_rel_ab_2020 = subset(asv_rel_ab, Year == "2020")
asv_rel_ab_2021 = subset(asv_rel_ab, Year == "2021")

asv_rel_ab_dist_2020 = vegdist(select(asv_rel_ab_2020, -c(sample, Treatment, Year)), method = "bray", diag = T)
asv_rel_ab_dist_2021 = vegdist(select(asv_rel_ab_2021, -c(sample, Treatment, Year)), method = "bray", diag = T)


##### NMDS ##### 
asv_nmds_2020 = metaMDS(asv_rel_ab_dist_2020, distance = "bray")
asv_nmds_2021 = metaMDS(asv_rel_ab_dist_2021, distance = "bray")

asv_nmds_data.scores_2020 = as.data.frame(scores(asv_nmds_2020))
asv_nmds_data.scores_2020$Sample = asv_rel_ab_2020$sample
asv_nmds_data.scores_2020$Treatment = asv_rel_ab_2020$Treatment
asv_nmds_data.scores_2020$Year = asv_rel_ab_2020$Year

asv_nmds_data.scores_2021 = as.data.frame(scores(asv_nmds_2021))
asv_nmds_data.scores_2021$Sample = asv_rel_ab_2021$sample
asv_nmds_data.scores_2021$Treatment = asv_rel_ab_2021$Treatment
asv_nmds_data.scores_2021$Year = asv_rel_ab_2021$Year

asv_nmds_df = rbind(asv_nmds_data.scores_2020, asv_nmds_data.scores_2021)




##### ——— Isolates and experimental data analysis ——— #####
##### Load isolates and experimental data #####
isolates = read.delim("Data/Isolates.txt")
inhibit = read.delim("Data/Inhibition_assays.txt")


##### Data formatting #####
isolates = subset(isolates, Identification != "Unknown")

isolates$Treatment = gsub("Fluconazole", "treated", isolates$Treatment)
isolates$Treatment = gsub("DMSO", "untreated", isolates$Treatment)

isolates$Identification = gsub("Penicillium bialowiezense", "Penicillium aff. bialowiezense", isolates$Identification)
isolates$Identification = gsub("Penicillium brevicompactum", "Penicillium aff. brevicompactum", isolates$Identification)
isolates$Identification = gsub("Penicillium chrysogenum", "Penicillium aff. chrysogenum", isolates$Identification)

isolates$Replicate = paste0("Rep_", isolates$Replicate)



inhibit$Group = gsub("[0-9]+", "", inhibit$Group)
inhibit$Group = gsub("\\.", "", inhibit$Group)

inhibit$Species = gsub( "P. ", "Penicillium ", inhibit$Species)
inhibit$Species = gsub( "Penicillium bialowiezense", "Penicillium aff. bialowiezense", inhibit$Species)
inhibit$Species = gsub( "Penicillium brevicompactum", "Penicillium aff. brevicompactum", inhibit$Species)
inhibit$Species = gsub( "Penicillium chrysogenum", "Penicillium aff. chrysogenum", inhibit$Species)

inhibit_df = inhibit %>% select( -Group ) %>% 
  pivot_wider( 
    names_from = Concentration,
    values_from = Value,
    values_fill = 0 )

inhibit_df = inhibit_df %>% group_by(Isolate) %>% dplyr::mutate(control = mean(`0`))

inhibit_df$`3` = 1 - ( inhibit_df$`3` / inhibit_df$control )
inhibit_df$`16` = 1 - ( inhibit_df$`16` / inhibit_df$control )
inhibit_df$`32` = 1 - ( inhibit_df$`32` / inhibit_df$control )
inhibit_df$`96` = 1 - ( inhibit_df$`96` / inhibit_df$control )
inhibit_df$`256` = 1 - ( inhibit_df$`256` / inhibit_df$control )
inhibit_df$`625` = 1 - ( inhibit_df$`625` / inhibit_df$control )

inhibit_df = select(inhibit_df, -c(`0`, Replicate, control) )

inhibit_df = melt(inhibit_df, variable.name = "Concentration")

treated = unique( subset(isolates_venn_df, Treatment == "treated")$Identification )
untreated = unique( subset(isolates_venn_df, Treatment == "untreated")$Identification )

inhibit_df$Source_Group[inhibit_df$Species %in% untreated & !(inhibit_df$Species %in% treated) ] <- "untreated"
inhibit_df$Source_Group[inhibit_df$Species %in% treated] <- "treated + untreated"




##### ——— MAIN FIGURES AND TABLES ——— #####


#### Figure 1 ####

asv_diversity = rbind(
  data.frame(test = "Shannon index", value = vegan::diversity(select(asv_rel_ab, -c(sample,  Treatment, Year)), groups = asv_rel_ab$sample, index = "shannon"), sample = asv_rel_ab$sample) 
  ,
  data.frame(test = "Simpson index", value = vegan::diversity(select(asv_rel_ab, -c(sample,  Treatment, Year)), groups = asv_rel_ab$sample, index = "simpson"), sample = asv_rel_ab$sample) 
  ,
  data.frame(test = "Observed ASVs (N)", value = vegan::specnumber(select(asv_rel_ab, -c(sample,  Treatment, Year)), groups = asv_rel_ab$sample), sample = asv_rel_ab$sample)  
)

asv_diversity = merge(asv_diversity, select(metadata, c(Sample, Treatment, Year)), by.x = "sample", by.y = "Sample")


Fig1 = ggplot(asv_diversity, aes(x = Treatment, y = value)) + 
  geom_boxplot( aes(colour = Treatment )) + 
  geom_jitter(size = 0.8, alpha = 0.7, aes(colour = Treatment), position=position_jitterdodge(jitter.width = 0.7)) +
  facet_grid(test ~ Year, scales = "free_y", switch = "y") +
  ylab(NULL) +
  geom_pwc(label = "p.adj",
           label.size = 2.5,
           bracket.nudge.y = 0.15,
           vjust = -0.1) +
  scale_y_continuous(expand = expansion(mult = 0.17)) +
  scale_colour_manual(values = c("#C52E19", "#AC9765")) +
  theme(axis.text.y = element_text(colour = "black", size = 8),
        axis.text.x = element_text(colour = "black", size = 8),
        axis.title.x = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 10), 
        text = element_text(size = 10),
        legend.text = element_text(size = 10, colour ="black"), 
        legend.position = "bottom", 
        legend.title = element_text(size = 10, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.placement = "outside",
        strip.position = "left") 

Fig1

pdf("Figures/Fig1.pdf", height = 5, width = 3, useDingbats = FALSE)
Fig1
dev.off()

tiff("Figures/Fig1.tiff", height = 5, width = 3, units = "in", compression = "lzw+p", res = 360)
Fig1
dev.off()



#### Figure 2 ####
Fig2 = ggplot(asv_nmds_df, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 2.5, aes( colour = Treatment), alpha = 0.8) + 
  theme(axis.text.y = element_text(colour = "black", size = 8),
        axis.text.x = element_text(colour = "black", size = 8),
        axis.title.x = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 10), 
        text = element_text(size = 10),
        legend.text = element_text(size = 10, colour ="black"), 
        legend.position = "right", 
        legend.title = element_text(size = 10, colour = "black"), 
        legend.key=element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) + 
  labs(x = "NMDS1", colour = "Treatment", y = "NMDS2")  + 
  scale_colour_manual(values = c("#C52E19", "#AC9765")) +
  facet_wrap(Year ~ ., scales = "free") +
  scale_y_continuous(expand = expansion(mult = 0.1))+
  scale_x_continuous(expand = expansion(mult = 0.1)) 

Fig2

pdf("Figures/Fig2.pdf", height = 2.5, width = 6.2, useDingbats = FALSE)
Fig2
dev.off()

tiff("Figures/Fig2.tiff", height = 2.5, width = 6.2, units = "in", compression = "lzw+p", res = 360)
Fig2
dev.off()


## PERMANOVA
PERM2020 = adonis2(select(asv_rel_ab_2020, -c(sample, Treatment, Year)) ~ Treatment, 
                   data = subset(metadata, Year == "2020" & ( Treatment == "treated" | Treatment == "untreated" ) ), 
                   permutations = 9999, 
                   method = "bray")
PERM2020

PERM2021 = adonis2(select(asv_rel_ab_2021, -c(sample, Treatment, Year)) ~ Treatment, 
                   data = subset(metadata, Year == "2021" & ( Treatment == "treated" | Treatment == "untreated" ) ), 
                   permutations = 9999, 
                   method = "bray")
PERM2021


## Distance to centroids

# Calculate distance to treatment centroids
disp_2020 <- betadisper(asv_rel_ab_dist_2020, subset(metadata, Year == "2020" & ( Treatment == "treated" | Treatment == "untreated" ))$Treatment )
disp_2021 <- betadisper(asv_rel_ab_dist_2021, subset(metadata, Year == "2021" & ( Treatment == "treated" | Treatment == "untreated" ))$Treatment )

# Distances of individual samples to their treatment centroid
dist_centroid_2020 <- disp_2020$distances
dist_centroid_2021 <- disp_2021$distances

# Add treatment information
disp_2020_df <- data.frame(
  Sample = asv_rel_ab_2020$sample,
  Treatment = subset(metadata, Year == "2020" & ( Treatment == "treated" | Treatment == "untreated" ))$Treatment,
  Distance = dist_centroid_2020)

disp_2021_df <- data.frame(
  Sample = asv_rel_ab_2021$sample,
  Treatment = subset(metadata, Year == "2021" & ( Treatment == "treated" | Treatment == "untreated" ))$Treatment,
  Distance = dist_centroid_2021)

anova(disp_2020)
anova(disp_2021)




#### Table 1 #### 

# 2020
asv_tax_filtered = subset(asv_tax, ASV %in% names(asv_rel_ab))

asv_mean_2020_treated = asv_rel_ab_2020 %>% subset(Treatment == "treated")
asv_mean_2020_untreated = asv_rel_ab_2020 %>% subset(Treatment == "untreated")

asv_mean_2021_treated = asv_rel_ab_2021 %>% subset(Treatment == "treated")
asv_mean_2021_untreated = asv_rel_ab_2021 %>% subset(Treatment == "untreated")

asv_mean_2020_treated = asv_mean_2020_treated %>% select(-sample, -Year, -Treatment) %>% colMeans() %>% as.data.frame() %>% dplyr::rename(abundance = ".")
asv_mean_2020_untreated = asv_mean_2020_untreated %>% select(-sample, -Year, -Treatment) %>% colMeans() %>% as.data.frame() %>% dplyr::rename(abundance = ".")
asv_mean_2020_treated$Treatment = "treated"
asv_mean_2020_untreated$Treatment = "untreated"
asv_mean_2020_treated$ASV = rownames(asv_mean_2020_treated)
asv_mean_2020_untreated$ASV = rownames(asv_mean_2020_untreated)

asv_mean_2020_treated = merge(asv_mean_2020_treated, select(asv_tax_filtered, c(ASV, Genus)), by = "ASV", all = T, sort = FALSE)
asv_mean_2020_untreated = merge(asv_mean_2020_untreated, select(asv_tax_filtered, c(ASV, Genus)), by = "ASV", all = T, sort = FALSE)

asv_mean_2020_treated = asv_mean_2020_treated %>% mutate(
  Genus = case_when( grepl("gen_Incertae_sedis", Genus) ~ "Other", TRUE ~ Genus ) )
asv_mean_2020_untreated = asv_mean_2020_untreated %>% mutate(
  Genus = case_when( grepl("gen_Incertae_sedis", Genus) ~ "Other", TRUE ~ Genus ) )

asv_mean_2020_treated = asv_mean_2020_treated %>% group_by(Treatment, Genus) %>% select(-ASV) %>% dplyr::summarize(abundance = sum(abundance))
asv_mean_2020_untreated = asv_mean_2020_untreated %>% group_by(Treatment, Genus) %>% select(-ASV) %>% dplyr::summarize(abundance = sum(abundance))



Top10_2020_treated = asv_mean_2020_treated %>% subset(Genus != "Other") %>% na.omit() %>% arrange(desc(abundance)) %>% dplyr::slice_head(n = 10)
Top10_2020_untreated = asv_mean_2020_untreated %>% subset(Genus != "Other") %>% na.omit() %>% arrange(desc(abundance)) %>% dplyr::slice_head(n = 10)

asv_mean_2020_treated$Genus[ !(asv_mean_2020_treated$Genus %in% Top10_2020_treated$Genus) ] <- "Other"
asv_mean_2020_untreated$Genus[ !(asv_mean_2020_untreated$Genus %in% Top10_2020_untreated$Genus) ] <- "Other"

asv_mean_2020 = rbind(asv_mean_2020_treated, asv_mean_2020_untreated)
rm(asv_mean_2020_treated, asv_mean_2020_untreated)

asv_mean_2020$Year = "2020"

# 2021
asv_mean_2021_treated = asv_rel_ab_2021 %>% subset(Treatment == "treated")
asv_mean_2021_untreated = asv_rel_ab_2021 %>% subset(Treatment == "untreated")

asv_mean_2021_treated = asv_rel_ab_2021 %>% subset(Treatment == "treated")
asv_mean_2021_untreated = asv_rel_ab_2021 %>% subset(Treatment == "untreated")

asv_mean_2021_treated = asv_mean_2021_treated %>% select(-sample, -Year, -Treatment) %>% colMeans() %>% as.data.frame() %>% dplyr::rename(abundance = ".")
asv_mean_2021_untreated = asv_mean_2021_untreated %>% select(-sample, -Year, -Treatment) %>% colMeans() %>% as.data.frame() %>% dplyr::rename(abundance = ".")
asv_mean_2021_treated$Treatment = "treated"
asv_mean_2021_untreated$Treatment = "untreated"
asv_mean_2021_treated$ASV = rownames(asv_mean_2021_treated)
asv_mean_2021_untreated$ASV = rownames(asv_mean_2021_untreated)

asv_mean_2021_treated = merge(asv_mean_2021_treated, select(asv_tax_filtered, c(ASV, Genus)), by = "ASV", all = T, sort = FALSE)
asv_mean_2021_untreated = merge(asv_mean_2021_untreated, select(asv_tax_filtered, c(ASV, Genus)), by = "ASV", all = T, sort = FALSE)

asv_mean_2021_treated = asv_mean_2021_treated %>% mutate(
  Genus = case_when( grepl("gen_Incertae_sedis", Genus) ~ "Other", TRUE ~ Genus ) )
asv_mean_2021_untreated = asv_mean_2021_untreated %>% mutate(
  Genus = case_when( grepl("gen_Incertae_sedis", Genus) ~ "Other", TRUE ~ Genus ) )

asv_mean_2021_treated = asv_mean_2021_treated %>% group_by(Treatment, Genus) %>% select(-ASV) %>% dplyr::summarize(abundance = sum(abundance))
asv_mean_2021_untreated = asv_mean_2021_untreated %>% group_by(Treatment, Genus) %>% select(-ASV) %>% dplyr::summarize(abundance = sum(abundance))


Top10_2021_treated = asv_mean_2021_treated %>% subset(Genus != "Other") %>% na.omit() %>% arrange(desc(abundance)) %>% dplyr::slice_head(n = 10)
Top10_2021_untreated = asv_mean_2021_untreated %>% subset(Genus != "Other") %>% na.omit() %>% arrange(desc(abundance)) %>% dplyr::slice_head(n = 10)

asv_mean_2021_treated$Genus[ !(asv_mean_2021_treated$Genus %in% Top10_2021_treated$Genus) ] <- "Other"
asv_mean_2021_untreated$Genus[ !(asv_mean_2021_untreated$Genus %in% Top10_2021_untreated$Genus) ] <- "Other"

asv_mean_2021 = rbind(asv_mean_2021_treated, asv_mean_2021_untreated)
rm(asv_mean_2021_treated, asv_mean_2021_untreated)

asv_mean_2021$Year = "2021"


asv_mean_df = rbind(asv_mean_2020, asv_mean_2021)
asv_mean_df$Genus = gsub("g__", "", asv_mean_df$Genus)

asv_mean_df_Top10 = asv_mean_df %>% 
  group_by(Treatment, Genus, Year) %>% 
  dplyr::summarize(abundance = sum(abundance)) %>% 
  arrange(desc(abundance))


table1 = asv_mean_df %>% 
  group_by(Treatment, Genus, Year) %>% 
  dplyr::summarize(abundance = (sum(abundance) ) * 100) %>% 
  ungroup() %>% 
  dplyr::mutate(Treatment_Year = paste0(Treatment, "_", Year)) %>% 
  select(-c(Treatment, Year)) %>% 
  pivot_wider(
    names_from = Treatment_Year,
    values_from = abundance,
    values_fill = 0
  )


#### Figure 3 ####

colourCount = length(unique(asv_mean_df_Top10$Genus))
getPalette = colorRampPalette(RColorBrewer::brewer.pal(8,"Pastel1"))

Fig3 = ggplot(asv_mean_df_Top10, aes ( x = Treatment, y = abundance, fill = factor(Genus, levels = c(sort(setdiff(Genus, "Other")), "Other")) )) +
  geom_col(color = "black", linewidth = 0.3) + 
  scale_fill_manual(name = "Genus", values = getPalette(colourCount)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = 0.01) +
  labs(x = "Treatment", y = "Averaged relative abundance") + 
  guides(fill = guide_legend(ncol = 1)) +
  facet_wrap(Year ~ .) +
  theme(axis.text.x = element_text(angle = 0, size = 8, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour ="black"), 
        legend.title = element_text(size = 8, colour = "black"),
        legend.key.size = unit(4, "mm"), 
        panel.background = element_blank(), 
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

Fig3

pdf("Figures/Fig3.pdf", height = 4, width = 5, useDingbats = FALSE)
Fig3
dev.off()

tiff("Figures/Fig3.tiff", height = 4, width = 5, units = "in", compression = "lzw+p", res = 360)
Fig3
dev.off()






#### Figure 4 ####

isolates_venn_df = isolates %>% dplyr::count(Treatment, Identification) 
isolates_venn_df = isolates_venn_df %>% group_by(Treatment) %>% mutate( perc = n/sum(n) )

isolates_venn_list = list( untreated = unique( subset(isolates_venn_df, Treatment == "untreated")$Identification ),
                           treated = unique( subset(isolates_venn_df, Treatment == "treated")$Identification ) )

venn_plot = ggvenn(isolates_venn_list, auto_scale = TRUE, show_percentage = FALSE, stroke_color = c("#AC9765", "#C52E19"), fill_alpha = 1, fill_color = c("#AC9765", "#C52E19"), stroke_size = 0, padding = 0.2, text_color = "black", set_name_size = 3)


colourCount = length(unique(isolates_venn_df$Identification))
getPalette = colorRampPalette(RColorBrewer::brewer.pal(8,"Pastel1"))

isolates_N = ggplot( isolates_venn_df, aes ( x = Treatment, y = n, fill = Identification) ) +
  geom_col(colour = "black", linewidth = 0.25)  +  
  scale_fill_manual(name = "Identification", values = getPalette(colourCount)) +
  scale_y_continuous(expand = 0.01) +
  labs(x = "Treatment", y = "Isolate count (N)") + 
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(angle = 0, size = 8, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour ="black", face = "italic"), 
        legend.title = element_text(size = 8, colour = "black"),
        legend.key.size = unit(4, "mm"), 
        panel.background = element_blank(), 
        strip.background = element_blank() )

isolates_perc = ggplot( isolates_venn_df, aes ( x = Treatment, y = perc, fill = Identification) ) +
  geom_col(colour = "black", linewidth = 0.25) +
  scale_fill_manual(name = "Identification", values = getPalette(colourCount)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = 0.01) +
  labs(x = "Treatment", y = "Relative abundance",) + 
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(angle = 0, size = 8, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour ="black", face = "italic"), 
        legend.title = element_text(size = 8, colour = "black"),
        legend.key.size = unit(4, "mm"), 
        panel.background = element_blank(), 
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))



Fig4 = plot_grid(venn_plot, isolates_perc, align = "h", rel_widths = c(0.4, 1), labels = "AUTO", label_size = 12 )


pdf("Figures/Fig4.pdf", height = 2.7, width = 5.5, useDingbats = FALSE)
Fig4
dev.off()


tiff("Figures/Fig4.tiff", height = 2.7, width = 5.5, units = "in", compression = "lzw+p", res = 360)
Fig4
dev.off()




#### Figure 5 ####


# Diversity - isolates
isolates_div_df = isolates %>% group_by(Treatment) %>% dplyr::count(Treatment, Identification, Replicate)  %>% 
  pivot_wider( names_from = Replicate, values_from = n  ) 

isolate_counts <- isolates %>%
  dplyr::count(Treatment, Replicate, Identification) %>%
  pivot_wider(
    names_from = Identification,
    values_from = n,
    values_fill = 0
  )

isolate_counts$isolates_n = select(isolate_counts, unique(isolates$Identification)) %>% rowSums()

isolate_counts$sample = paste0(isolate_counts$Treatment, "_", isolate_counts$Replicate)

isolate_counts_PA = vegan::decostand(select(isolate_counts, -c(Replicate, Treatment, sample)), method = "pa")
isolate_counts_PA$sample = isolate_counts$sample
isolate_counts_PA$Replicate = isolate_counts$Replicate
isolate_counts_PA$Treatment = isolate_counts$Treatment


isolate_diversity = rbind(
  data.frame(test = "Shannon index", value = vegan::diversity(select(isolate_counts, -c(Replicate, Treatment, sample, isolates_n)), groups = isolate_counts$sample, index = "shannon"), sample = isolate_counts$sample) 
  ,
  data.frame(test = "Simpson index", value = vegan::diversity(select(isolate_counts, -c(Replicate, Treatment, sample, isolates_n)), groups = isolate_counts$sample, index = "simpson"), sample = isolate_counts$sample) 
  ,
  data.frame(test = "Isolates (N)", value = isolate_counts$isolates_n, sample = isolate_counts$sample)
  ,
  data.frame(test = "Species (N)", value = vegan::specnumber(select(isolate_counts_PA, -c(Replicate, Treatment, sample)), groups = isolate_counts$sample), sample = isolate_counts$sample)
)

isolate_diversity = merge(isolate_diversity, select(isolate_counts, c(sample, Treatment)), by.x = "sample", by.y = "sample")


isolate_diversity$test <- factor(
  isolate_diversity$test,
  levels = c("Isolates (N)", "Species (N)", "Shannon index", "Simpson index")
)


isolate_diversity %>% 
  select(test, Treatment, value) %>% 
  group_by(test, Treatment) %>% 
  summarise(mean_value = mean(value), .groups = "drop")



Fig5 = ggplot(isolate_diversity, aes(x = Treatment, y = value)) + 
  geom_boxplot( aes(colour = Treatment )) + 
  geom_jitter(size = 0.8, alpha = 0.7, aes(colour = Treatment), position=position_jitterdodge(jitter.width = 0.7)) +
  facet_wrap(. ~ test, scales = "free", strip.position = "left", nrow = 1 ) +
  ylab(NULL) +
  geom_pwc(label = "p.adj",
           label.size = 2.5,
           bracket.nudge.y = 0.15,
           vjust = -0.1) +
  scale_y_continuous(expand = expansion(mult = 0.08)) +
  scale_colour_manual(values = c("#C52E19", "#AC9765")) +
  theme(axis.text.y = element_text(colour = "black", size = 8),
        axis.text.x = element_text(colour = "black", size = 8),
        axis.title.x = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 10), 
        text = element_text(size = 10),
        legend.text = element_text(size = 10, colour ="black"), 
        legend.position = "bottom", 
        legend.title = element_text(size = 10, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        panel.spacing = unit(0.3, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", margin = margin(r = 0)),
        strip.placement = "outside",
        strip.position = "left")
Fig5

pdf("Figures/Fig5.pdf", height = 3, width = 6.2, useDingbats = FALSE)
Fig5
dev.off()

tiff("Figures/Fig5.tiff", height = 3, width = 6.2, units = "in", compression = "lzw+p", res = 360)
Fig5
dev.off()


#### Figure 6 ####
Fig6 = inhibit_df %>% 
  subset(Concentration == "625") %>% 
  ggplot(aes( x = Source_Group, y = value, colour = Source_Group )) +
  geom_boxplot() +
  geom_jitter(size = 1, alpha = 0.7, aes(colour = Source_Group), position=position_jitterdodge(jitter.width = 0.7)) +
  ylab("Relative inhibition") +
  xlab("Occurrence") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = 0.03, limits = c(0,1), breaks = seq(0, 1, 0.2)) +
  geom_pwc(label = "p.adj",
           label.size = 2.5,
           bracket.nudge.y = 0.07,
           vjust = -0.1) +
  scale_colour_manual(name = "Occurrence", values = c("#C52E19", "#AC9765")) +
  theme(axis.text.y = element_text(colour = "black", size = 8),
        axis.text.x = element_text(colour = "black", size = 8),
        axis.title.x = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 10), 
        text = element_text(size = 10),
        legend.text = element_text(size = 10, colour ="black"), 
        legend.position = "right", 
        legend.title = element_text(size = 10, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        panel.spacing = unit(0.3, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", margin = margin(r = 0)),
        strip.placement = "outside",
        strip.position = "left")

Fig6

pdf("Figures/Fig6.pdf", height = 2.5, width = 4.5, useDingbats = FALSE)
Fig6
dev.off()

tiff("Figures/Fig6.tiff", height = 2.5, width = 4.5, units = "in", compression = "lzw+p", res = 360)
Fig6
dev.off()


#### Figure 7 ####

colourCount = length(unique(inhibit_df$Concentration))
getPalette = colorRampPalette(RColorBrewer::brewer.pal(10,"Pastel1"))

Fig7 = inhibit_df %>% 
  subset(Concentration != "0") %>% 
  subset(Species == "Aspergillus fumigatus") %>% 
  mutate( facet_label = paste0( "bolditalic('", Species, "')~bold('", Isolate, "')" ) ) %>% 
  ggplot(aes( x = Concentration, y = value, fill = Concentration )) +
  geom_boxplot() +
  geom_jitter(size = 1, alpha = 0.7, colour = "black", position=position_jitterdodge(jitter.width = 0.7)) +
  ylab("Relative inhibition") +
  xlab("Concentration (µg/mL)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = 0.01, limits = c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_fill_manual(name = "Concentration (µg/mL)", values = getPalette(colourCount)) +
  theme(axis.text.y = element_text(colour = "black", size = 8),
        axis.text.x = element_text(colour = "black", size = 8),
        axis.title.x = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 10), 
        text = element_text(size = 10),
        legend.text = element_text(size = 10, colour ="black"), 
        legend.position = "none", 
        legend.title = element_text(size = 10, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        panel.spacing = unit(0.3, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic")) +
  facet_wrap(facet_label ~ ., labeller = label_parsed)

Fig7

pdf("Figures/Fig7.pdf", height = 3, width = 6.2, useDingbats = FALSE)
Fig7
dev.off()

tiff("Figures/Fig7.tiff", height = 3, width = 5.5, units = "in", compression = "lzw+p", res = 360)
Fig7
dev.off()




##### ——— SUPPLEMENTARY FIGURES ——— #####

#### Supplementary Figure 1 ####

asvs_venn_list = list( untreated = asvs_untreated,
                       treated = asvs_treated )

SupFig1 = ggvenn(asvs_venn_list, auto_scale = TRUE, show_percentage = FALSE, stroke_color = c("#AC9765", "#C52E19"), fill_alpha = 0.7, fill_color = c("#AC9765", "#C52E19"), stroke_size = 2, padding = 0.2, text_color = "black", set_name_size = 3) +
  geom_text(aes( x = -1, y = -0.15, label = paste0( "(", format(asvs_untreated_reads, digits = 2), "%)" ) ) ) +
  geom_text(aes( x = 0, y = -0.15, label = paste0(  "(", format(asvs_shared_reads, digits = 4), "%)" ) ) ) +
  geom_text(aes( x = 1, y = -0.15, label = paste0(  "(", format(asvs_treated_reads, digits = 2), "%)" ) ) ) 

gt <- ggplotGrob(SupFig1)

panel_id <- which(gt$layout$name == "panel")
panel <- gt$grobs[[panel_id]]

venn <- panel$children[["geom_venn"]]

text_ids <- names(venn$children)[
  sapply(venn$children, inherits, what = "text")
]

for (i in text_ids) {
  venn$children[[i]]$gp$fontface <- "bold"
}

panel$children[["geom_venn"]] <- venn
gt$grobs[[panel_id]] <- panel


pdf("Figures/SupFig1.pdf", height = 4, width = 4, useDingbats = FALSE)
grid::grid.draw(gt)
dev.off()

tiff("Figures/SupFig1.tiff", height = 4, width = 4, units = "in", compression = "lzw+p", res = 360)
grid::grid.draw(gt)
dev.off()



#### Supplementary Figure 2 ####

# 2020
asv_tax_filtered = subset(asv_tax, ASV %in% names(asv_rel_ab))

asv_mean_2020_treated = asv_rel_ab_2020 %>% subset(Treatment == "treated")
asv_mean_2020_untreated = asv_rel_ab_2020 %>% subset(Treatment == "untreated")

asv_mean_2021_treated = asv_rel_ab_2021 %>% subset(Treatment == "treated")
asv_mean_2021_untreated = asv_rel_ab_2021 %>% subset(Treatment == "untreated")

asv_mean_2020_treated = asv_mean_2020_treated %>% select(-sample, -Year, -Treatment) %>% colMeans() %>% as.data.frame() %>% dplyr::rename(abundance = ".")
asv_mean_2020_untreated = asv_mean_2020_untreated %>% select(-sample, -Year, -Treatment) %>% colMeans() %>% as.data.frame() %>% dplyr::rename(abundance = ".")
asv_mean_2020_treated$Treatment = "treated"
asv_mean_2020_untreated$Treatment = "untreated"
asv_mean_2020_treated$ASV = rownames(asv_mean_2020_treated)
asv_mean_2020_untreated$ASV = rownames(asv_mean_2020_untreated)

asv_mean_2020_treated = merge(asv_mean_2020_treated, select(asv_tax_filtered, c(ASV, Genus)), by = "ASV", all = T, sort = FALSE)
asv_mean_2020_untreated = merge(asv_mean_2020_untreated, select(asv_tax_filtered, c(ASV, Genus)), by = "ASV", all = T, sort = FALSE)

asv_mean_2020_treated = asv_mean_2020_treated %>% mutate(
  Genus = case_when( grepl("gen_Incertae_sedis", Genus) ~ "Other", TRUE ~ Genus ) )
asv_mean_2020_untreated = asv_mean_2020_untreated %>% mutate(
  Genus = case_when( grepl("gen_Incertae_sedis", Genus) ~ "Other", TRUE ~ Genus ) )

asv_mean_2020_treated = asv_mean_2020_treated %>% group_by(Treatment, Genus) %>% select(-ASV) %>% dplyr::summarize(abundance = sum(abundance))
asv_mean_2020_untreated = asv_mean_2020_untreated %>% group_by(Treatment, Genus) %>% select(-ASV) %>% dplyr::summarize(abundance = sum(abundance))



Top30_2020_treated = asv_mean_2020_treated %>% subset(Genus != "Other") %>% na.omit() %>% arrange(desc(abundance)) %>% dplyr::slice_head(n = 30)
Top30_2020_untreated = asv_mean_2020_untreated %>% subset(Genus != "Other") %>% na.omit() %>% arrange(desc(abundance)) %>% dplyr::slice_head(n = 30)

asv_mean_2020_treated$Genus[ !(asv_mean_2020_treated$Genus %in% Top30_2020_treated$Genus) ] <- "Other"
asv_mean_2020_untreated$Genus[ !(asv_mean_2020_untreated$Genus %in% Top30_2020_untreated$Genus) ] <- "Other"

asv_mean_2020 = rbind(asv_mean_2020_treated, asv_mean_2020_untreated)
rm(asv_mean_2020_treated, asv_mean_2020_untreated)

asv_mean_2020$Year = "2020"




# 2021

asv_mean_2021_treated = asv_rel_ab_2021 %>% subset(Treatment == "treated")
asv_mean_2021_untreated = asv_rel_ab_2021 %>% subset(Treatment == "untreated")

asv_mean_2021_treated = asv_rel_ab_2021 %>% subset(Treatment == "treated")
asv_mean_2021_untreated = asv_rel_ab_2021 %>% subset(Treatment == "untreated")

asv_mean_2021_treated = asv_mean_2021_treated %>% select(-sample, -Year, -Treatment) %>% colMeans() %>% as.data.frame() %>% dplyr::rename(abundance = ".")
asv_mean_2021_untreated = asv_mean_2021_untreated %>% select(-sample, -Year, -Treatment) %>% colMeans() %>% as.data.frame() %>% dplyr::rename(abundance = ".")
asv_mean_2021_treated$Treatment = "treated"
asv_mean_2021_untreated$Treatment = "untreated"
asv_mean_2021_treated$ASV = rownames(asv_mean_2021_treated)
asv_mean_2021_untreated$ASV = rownames(asv_mean_2021_untreated)

asv_mean_2021_treated = merge(asv_mean_2021_treated, select(asv_tax_filtered, c(ASV, Genus)), by = "ASV", all = T, sort = FALSE)
asv_mean_2021_untreated = merge(asv_mean_2021_untreated, select(asv_tax_filtered, c(ASV, Genus)), by = "ASV", all = T, sort = FALSE)

asv_mean_2021_treated = asv_mean_2021_treated %>% mutate(
  Genus = case_when( grepl("gen_Incertae_sedis", Genus) ~ "Other", TRUE ~ Genus ) )
asv_mean_2021_untreated = asv_mean_2021_untreated %>% mutate(
  Genus = case_when( grepl("gen_Incertae_sedis", Genus) ~ "Other", TRUE ~ Genus ) )

asv_mean_2021_treated = asv_mean_2021_treated %>% group_by(Treatment, Genus) %>% select(-ASV) %>% dplyr::summarize(abundance = sum(abundance))
asv_mean_2021_untreated = asv_mean_2021_untreated %>% group_by(Treatment, Genus) %>% select(-ASV) %>% dplyr::summarize(abundance = sum(abundance))


Top30_2021_treated = asv_mean_2021_treated %>% subset(Genus != "Other") %>% na.omit() %>% arrange(desc(abundance)) %>% dplyr::slice_head(n = 30)
Top30_2021_untreated = asv_mean_2021_untreated %>% subset(Genus != "Other") %>% na.omit() %>% arrange(desc(abundance)) %>% dplyr::slice_head(n = 30)

asv_mean_2021_treated$Genus[ !(asv_mean_2021_treated$Genus %in% Top30_2021_treated$Genus) ] <- "Other"
asv_mean_2021_untreated$Genus[ !(asv_mean_2021_untreated$Genus %in% Top30_2021_untreated$Genus) ] <- "Other"

asv_mean_2021 = rbind(asv_mean_2021_treated, asv_mean_2021_untreated)
rm(asv_mean_2021_treated, asv_mean_2021_untreated)

asv_mean_2021$Year = "2021"


asv_mean_df = rbind(asv_mean_2020, asv_mean_2021)
asv_mean_df$Genus = gsub("g__", "", asv_mean_df$Genus)

asv_mean_df_Top30 = asv_mean_df %>% 
  group_by(Treatment, Genus, Year) %>% 
  dplyr::summarize(abundance = sum(abundance)) %>% 
  arrange(desc(abundance))

colourCount = length(unique(asv_mean_df_Top30$Genus))
getPalette = colorRampPalette(RColorBrewer::brewer.pal(8,"Pastel1"))

SupFig2 = ggplot(asv_mean_df_Top30, aes ( x = Treatment, y = abundance, fill = factor(Genus, levels = c(sort(setdiff(Genus, "Other")), "Other")) )) +
  geom_col(color = "black", linewidth = 0.25) + 
  scale_fill_manual(name = "Genus", values = getPalette(colourCount)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = 0.01) +
  labs(x = "Treatment", y = "Averaged relative abundance") + 
  guides(fill = guide_legend(ncol = 2)) +
  facet_wrap(Year ~ .) +
  theme(axis.text.x = element_text(angle = 0, size = 8, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour ="black"), 
        legend.title = element_text(size = 8, colour = "black"),
        legend.key.size = unit(4, "mm"), 
        panel.background = element_blank(), 
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

SupFig2

pdf("Figures/SupFig2.pdf", height = 4, width = 6.2, useDingbats = FALSE)
SupFig2
dev.off()

tiff("Figures/SupFig2.tiff", height = 4, width = 6.2, units = "in", compression = "lzw+p", res = 360)
SupFig2
dev.off()


#### Supplementary Figure 3 ####

tree1 = ape::read.tree("Data/Non-Penicillium_tree.nwk")

tree_mid <- phangorn::midpoint(tree1)

tree_mid$node.label[tree_mid$node.label < 0.7 | is.na(tree_mid$node.label)] <- as.numeric("")

tree_mid$node.label = format(as.numeric(tree_mid$node.label)*100, digits = 1)
tree_mid$node.label[tree_mid$node.label == " NA"] <- ""

tree_mid$tip.label = gsub( "_ITS_region_from_TYPE_material" , "_T", tree_mid$tip.label)

SupFig3 <- ggtree(tree_mid) +
  geom_treescale( x = 0, y = -1, fontsize = 3, linesize = 1 ) +
  geom_tiplab( data = function(x) subset(x, !grepl("^SF", label)), fontface = "plain", size = 3, hjust = -0.01 ) +
  geom_tiplab( data = function(x) subset(x, grepl("^SF", label)), fontface = "bold", size = 3, hjust = -0.05 ) +
  geom_nodelab( size = 2.5, hjust = 1.5, vjust = -0.5 ) +
  xlim(0, 1.6)

SupFig3


pdf("Figures/SupFig3.pdf", height = 7, width = 6.2, useDingbats = FALSE)
SupFig3
dev.off()

tiff("Figures/SupFig3.tiff", height = 7, width = 6.2, units = "in", compression = "lzw+p", res = 360)
SupFig3
dev.off()



#### Supplementary Figure 4 ####

tree2 = ape::read.tree("Data/Penicillium_ITS_tree.nwk")

tree2_rooted = root(tree2, outgroup = "NR_135337.1_Aspergillus_glaucus_NRRL_116_T", resolve.root = TRUE)

tree2_rooted$node.label[tree2_rooted$node.label < 0.7 | is.na(tree2_rooted$node.label)] <- as.numeric("")

tree2_rooted$node.label = format(as.numeric(tree2_rooted$node.label)*100, digits = 1)
tree2_rooted$node.label[tree2_rooted$node.label == " NA"] <- ""

SupFig4 <- ggtree(tree2_rooted) +
  geom_treescale( x = 0, y = -1, fontsize = 3, linesize = 1 ) +
  geom_tiplab( data = function(x) subset(x, !grepl("^SF", label)), fontface = "plain", size = 3, hjust = -0.01 ) +
  geom_tiplab( data = function(x) subset(x, grepl("^SF", label)), fontface = "bold", size = 3, hjust = -0.05 ) +
  geom_nodelab( size = 2.5, hjust = 1.5, vjust = -0.5 ) +
  xlim(0, 0.3)

SupFig4


pdf("Figures/SupFig4.pdf", height = 7.5, width = 6.2, useDingbats = FALSE)
SupFig4
dev.off()

tiff("Figures/SupFig4.tiff", height = 7.5, width = 6.2, units = "in", compression = "lzw+p", res = 360)
SupFig4
dev.off()



#### Supplementary Figure 5 ####
tree3 = ape::read.tree("Data/Penicillium_beta_tubulin_tree.nwk")

tree3_rooted = root(tree3, outgroup = "LC589346.1_Aspergillus_glaucus_NRRL_116_T", resolve.root = TRUE)
tree3_rooted$tip.label = gsub( "'" , "", tree3_rooted$tip.label)

tree3_rooted$node.label[tree3_rooted$node.label < 0.7 | is.na(tree3_rooted$node.label)] <- as.numeric("")

tree3_rooted$node.label = format(as.numeric(tree3_rooted$node.label)*100, digits = 1)
tree3_rooted$node.label[tree3_rooted$node.label == " NA"] <- ""

SupFig5 <- ggtree(tree3_rooted) +
  geom_treescale( x = 0, y = -1, fontsize = 3, linesize = 1 ) +
  geom_tiplab( data = function(x) subset(x, !grepl("^SF", label)), fontface = "plain", size = 3, hjust = -0.01 ) +
  geom_tiplab( data = function(x) subset(x, grepl("^SF", label)), fontface = "bold", size = 3, hjust = -0.05 ) +
  geom_nodelab( size = 2.5, hjust = 1.5, vjust = -0.5 ) +
  xlim(0, 1.5)

SupFig5


pdf("Figures/SupFig5.pdf", height = 6, width = 6.2, useDingbats = FALSE)
SupFig5
dev.off()

tiff("Figures/SupFig5.tiff", height = 6, width = 6.2, units = "in", compression = "lzw+p", res = 360)
SupFig5
dev.off()


#### Supplementary Figure 6 ####

colourCount = length(unique(inhibit_df$Species))
getPalette = colorRampPalette(RColorBrewer::brewer.pal(50,"Pastel1"))


vlines <- inhibit_df %>%
  filter(Concentration != "0") %>%
  distinct(Source_Group, Species) %>%
  group_by(Source_Group) %>%
  arrange(Species) %>%
  mutate(xintercept = row_number() + 0.5) %>%
  ungroup()

SupFig6 = inhibit_df %>% 
  subset(Concentration != "0" ) %>% 
  ggplot(aes( x = Species, y = value, fill = Concentration )) +
  geom_boxplot() +
  geom_jitter(size = 1, alpha = 0.7, aes(fill = Concentration), position=position_jitterdodge(jitter.width = 0.4)) +
  ylab("Relative inhibition") +
  xlab("Species") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = 0.01, limits = c(0,1), breaks = seq(0, 1, 0.2)) +
  scale_fill_manual(name = "Concentration (µg/mL)", values = getPalette(colourCount)) +
  theme(axis.text.y = element_text(colour = "black", size = 8),
        axis.text.x = element_text(colour = "black", size = 8),
        axis.title.x = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 10), 
        text = element_text(size = 10),
        legend.text = element_text(size = 10, colour ="black"), 
        legend.position = "bottom", 
        legend.title = element_text(size = 10, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        panel.spacing = unit(0.3, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.grid.minor.x = element_line(colour = "red")) +
  facet_wrap(Source_Group ~ ., scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10), expand = 0) +
  geom_vline( data = vlines, aes(xintercept = xintercept),  colour = "grey70", linewidth = 0.4, inherit.aes = FALSE ) +
  guides(  fill = guide_legend(nrow = 1) )

SupFig6


pdf("Figures/SupFig6.pdf", height = 5, width = 12, useDingbats = FALSE)
SupFig6
dev.off()

tiff("Figures/SupFig6.tiff", height = 5, width = 12, units = "in", compression = "lzw+p", res = 360)
SupFig6
dev.off()


