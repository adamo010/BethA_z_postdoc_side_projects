#Hello! Today I am doing 16s analysis on samples from the Mississippi microbiome project. 

#this is based on 08.27.21_mayo_microbiome_dada2_attempt1_clean.R
#and 10.27.22_mississippi_mayo_microbiome_data2.R

######step 1: load packages
#need to download first, but only need to do this once, fortunately. 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.13")
library(dada2); packageVersion("dada2") #this prints the package version that you're using. Here, I get 1.20.0.
library(stringr) #for editing strings; a tidyverse creation
library("ggplot2")
library("readxl")
library("dplyr")
library("tidyverse")
library("ggpubr")
library("ape")
library("lme4")
library("gplots")
library("plotly")
library("tidyr")
library("vegan")
library("data.table")
library("stringr")

rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/Blekhman_lab_external_consulting_projects/Mississippi_CRC_microbiome_project_2021/") 

######step 2: set the path to the fastq files
path <- "/Users/adamo010/Documents/Blekhman_lab_external_consulting_projects/Mississippi_CRC_microbiome_project_2021/"
head(list.files(path)) #print a few file names in the file path

#Now we read in the names of the fastq files, and perform some string manipulation to get matched lists of the forward and reverse fastq files.

######step 3: Get matched lists of the forward and reverse fastq.gz files: look at the file names to identify the pattern
forward_reads <- sort(list.files(path, pattern="_R1_001_trimmed.fastq.gz", full.names = TRUE))
reverse_reads <- sort(list.files(path, pattern="_R2_001_trimmed.fastq.gz", full.names = TRUE))

forward_reads[[1]]; reverse_reads[[1]] #this prints off the first element in each list (b/c R isn't 0- indexed like my old friend python)

######step 4:  Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(forward_reads), "_R1_001_trimmed.fastq.gz"), `[`, 1) #this creates a list called sample_names based on fnFs, which has _R1_001.fastq.gz trimmed

######step 5: inspect visual quality of the reads
#first, do forward reads
plotQualityProfile(forward_reads[95:102])
#somehow, from this, we're supposed to decide how many reads to trim. For illumina forward reads of 250bp, usually trim the last 10.
#since these data are 300bp reads, trim the top 15, to keep 285.
#now, look at reverse reads
plotQualityProfile(reverse_reads[1:8])
#well... we'll try removing 100 reads, to keep 200. This is recommended in

#######step 6: Place filtered files in filtered/ subdirectory
#first, create output files for our filtering
filtered_forward_reads <- file.path(path, "filtered", paste0(sample.names, "_F_trimmed_filtered.fastq.gz"))
filtered_reverse_reads <- file.path(path, "filtered", paste0(sample.names, "_R_trimmed_filtered.fastq.gz"))
filtered_forward_reads2 <- paste0(sample.names, "_F_trimmed_filtered.fastq.gz")
filtered_reverse_reads2 <- paste0(sample.names, "_R_trimmed_filtered.fastq.gz")

#create output directory
path2 <- "/Users/adamo010/Documents/Blekhman_lab_external_consulting_projects/Mississippi_CRC_microbiome_project_2021/filtered/"

#######step 7: filter and trim, using default parameters
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
                              rm.phix=TRUE, minLen=175, truncLen=c(250,200), 
                              compress=TRUE, matchIDs = TRUE, truncQ=2)

plotQualityProfile(filtered_forward_reads[1:6])
plotQualityProfile(filtered_reverse_reads[1:6])
#I think I want to trim more than 250 from the forward reads. Maybe 230. Lots are still messy.
#ditto for the reverse reads. Trim an extra 20bp from those too
filtered_out2 <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
                              rm.phix=TRUE, minLen=175, truncLen=c(230,180), 
                              compress=TRUE, matchIDs = TRUE, truncQ=2)

#that looks better.

#######step 8: errors!
#this takes a while. 
path <- "/Users/adamo010/Documents/Blekhman_lab_external_consulting_projects/Mississippi_CRC_microbiome_project_2021/filtered/"
errF <- learnErrors(filtered_forward_reads, multithread=TRUE)
errR <- learnErrors(filtered_reverse_reads, multithread=TRUE)
#error! Not all provided files exist. So, double check.
#table(file.exists(filtFs2)) #only four files exist
#whatever the fuck that is. 
#plot:
plotErrors(errF, nominalQ=TRUE) 
#looks nice.

#######step 9: INFERRING ASVs: aka sample inference
#note that we're rolling denoising into this step.
dadaFs <- dada(filtered_forward_reads, err=errF, pool= "pseudo", multithread=TRUE) 
#specify files (filtered_forward_reads), errors (errF), turn ON pooling to detect rare
#variants (pool=pseudo), and turn on multithreading to get things to run faster. 
dadaRs <- dada(filtered_reverse_reads, err=errR, pool= "pseudo", multithread=TRUE)
#now, we look at the results:
dadaFs[[1]]
dadaFs[[2]]
dadaFs[[3]]
dadaFs[[4]]
dadaRs[[1]]
dadaRs[[2]]
dadaRs[[3]]
dadaRs[[4]]
#all right, we have 424 sequence variants inferred from 18899 unique sequences in sample 1.

#######step 10: merge paired reads
mergers <- mergePairs(dadaFs, filtered_forward_reads, dadaRs, filtered_reverse_reads, verbose=TRUE)
#by default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical 
#to each other in the overlap region (but these conditions can be changed via function arguments).
#if you did the trimming right, most of the pairs should merge. In my case, I think I only lost 1-2% of reads here.
head(mergers[[1]])
#from tutorial: The mergers object is a list of data.frames from each sample. Each data.frame contains the merged $sequence, 
#its $abundance, and the indices of the $forward and $reverse sequence variants that were merged. Paired reads that did not exactly overlap 
#were removed by mergePairs, further reducing spurious output.
#if you did the trimming right, most of the 

#######step 11: make an ASV TABLE!!!
#so exciting...
seqtab <- makeSequenceTable(mergers)
# Inspect distribution of sequence lengths
dim(seqtab) #have 7942 ASVs across 102 samples 
table(nchar(getSequences(seqtab)))
#from tutorial: The sequence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants. 

#######step 12: remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#identified 324 bimeras out of 735 input sequences.
#according to the tutorial, this is not uncommon; a lot of ASVs might be removed, but not very many reads. 
#so, check how many reads are non-chimeric
sum(seqtab.nochim)/sum(seqtab)
#great, 96% of reads remain. This means our upstream processing was successful.
#if a lot of reads are removed, it likely means primer sequences were not removed.

#######step 13: sanity check: let's do a summary of our reads. 
#I have no idea what's going on here. But it generates a summary table, so let's just run it.
getN <- function(x) sum(getUniques(x)) #write a function
track <- cbind(filtered_out2, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)) #stick a bunch of data together
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names2
head(track)
#the goal here is to make sure that we're keeping the majority of raw reads, and that no single step has a large drop in # reads
#save this breakdown as a csv
write.csv(as.data.frame(track), file="10.29.21_DADA2_analysis_Mississippi_data_read_count_tracking.csv")


#######step 14: assigning taxonomy (newish version)
#downloading DECIPHER-formatted SILVA v138 reference
#download.file(url="http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", destfile="SILVA_SSU_r138_2019.RData")
#loading reference taxonomy object (whatever the fuck that is)
#loading DECIPHER
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER", force= TRUE)
library(DECIPHER)
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("/Users/adamo010/SILVA_SSU_r138_2019.RData")
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors. Still takes FOREVER
#######step 15: evaluate accuracy: can't find where I should have done this

#######step 16: creating outputs
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim) #this is the unmanageble names
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character") #this is to make the ASV names useful and not some horrific seq 
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
  } #this creates asv_headers, but doesn't actually change the column names. 
#the tutorial seems to have forgotten to do that.
#first, we need to create a 2nd ASV header thing
asv_headers2 <- vector(dim(seqtab.nochim)[2], mode="character") #this is to make the ASV names useful and not some horrific seq 
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers2[i] <- paste("ASV", i, sep="_") #use ASV instead of >ASV
} #this creates asv_headers, but doesn't actually change the column names. 
colnames(seqtab.nochim) <- c(asv_headers2)

#making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "03.11.22_DADA2_analysis_Mississippi_data_microbiome_ASVs.fa")

#count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "03.11.22_DADA2_analysis_Mississippi_data_microbiome_ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

#tax table:
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim) #smash everything together?
row.names(taxid) <- sub(">", "", asv_headers) #this renames ASVs to something useful
write.table(taxid, "03.11.22_DADA2_analysis_Mississippi_data_microbiome_ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

#######step 17: removing contaminants
library(decontam)
colnames(asv_tab) # print column names to used to find the blanks
#have 102 samples; #1-96 are real, 97-101 are blanks, 102 is a mock community
vector_for_decontam <- c(rep(FALSE, 96), rep(TRUE, 5), rep(FALSE, 1)) #define filtering vector to id contams
contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam) #create new dataframe of contaminants
table(contam_df$contaminant) # identified 18 as contaminants (print as TRUE in this table)
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ]) #get vector of contaminant ids
contam_asv_table <- taxid[row.names(taxid) %in% contam_asvs, ] #convert to a taxonomy table
write.table(contam_asv_table, "03.11.22_DADA2_analysis_Mississippi_data_contaminant_ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

#NEW 03.14.22: actually remove the contaminant sequences
#making new fasta file
contam_indices <- which(asv_fasta %in% paste0(">", contam_asv_table))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_no_contam <- asv_fasta[- dont_want]
#making new count table
asv_tab_no_contam <- asv_tab[!row.names(asv_tab) %in% contam_asvs, ]
#making new taxonomy table
asv_tax_no_contam <- taxid[!row.names(taxid) %in% contam_asvs, ]
#and now writing them out to files
write(asv_fasta_no_contam, "03.14.22_Mississippi_data_microbiome_ASVs_no_contams.fa")
write.table(asv_tab_no_contam, "03.14.22_Mississippi_data_microbiome_ASVs_counts_no_contam.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam, "03.14.22_Mississippi_data_microbiome_ASVs_taxonomy_no_contam.tsv",
            sep="\t", quote=F, col.names=NA)

#######step 18 NEW as of 03.11.22: handing off to phyloseq
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

#here(finally) is our metadata
metadata <- read_csv("/Users/adamo010/Documents/Blekhman_lab_external_consulting_projects/Mississippi_CRC_microbiome_project_2021/Metadata_and_extra_stuff/Mississippi_seqsonly_metadata.csv")
samples.out <- rownames(seqtab.nochim) 
#the next few lines are from the tutorial- extracting metadata from the sample names.
#subject <- sapply(strsplit(samples.out, "D"), `[`, 1) #edit me here
#gender <- substr(subject,1,1) #edit me here
#subject <- substr(subject,2,999) #edit me here
#day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
#samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
#samdf$When <- "Early"
#samdf$When[samdf$Day>100] <- "Late"
#rownames(metadata) <- samples.out

#one thing I need to do is edit the sample names. been getting an error that says "component sample names do not match"
#remove the "_F_trimmed_filtered.fastq.gz" from all rownames
#sample.names is the vector we want to use to replace these
row.names(seqtab.nochim) <- sample.names #nice. this replaces the row names. 

#double check that the sample names in the metadata match what is in the table
#we want Seq_sample in metadata to be row name
metadata <- column_to_rownames(metadata, var = "Seq_sample")
rownames(metadata) #okay good, now sample names in table match metadata

#ding ding ding
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxid))
ps <- prune_samples(sample_names(ps) != "H2O_", ps) # Remove H20 blanks
ps <- prune_samples(sample_names(ps) != "UMGC_Mock_S102", ps) # Remove mock community. Weird that it had to be so specific.
#there's still an extra sample here... what is it? one of the blanks...
sample_names(ps) #prints sample names in phyloseq object ps

#store the DNA sequences of our ASVs in the refseq slot of the phyloseq object, and then rename 
#our taxa to a short string. That way, the short new taxa names will appear in tables and plots, 
#and we can still recover the DNA sequences corresponding to each ASV as needed with refseq(ps).
#FUCK IT I DON"T CARE
#dna <- Biostrings::DNAStringSet(taxa_names(ps))
#names(dna) <- taxa_names(ps)
#ps <- merge_phyloseq(ps, dna)
#taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps #prints info about phyloseq object ps

#######step 18a NEW as of 03.14.22: uploading
rm(list=ls())
setwd("/Users/adamo010/Documents/Blekhman_lab_external_consulting_projects/Mississippi_CRC_microbiome_project_2021/") 

count_tab <- read.table("03.14.22_Mississippi_data_microbiome_ASVs_counts_no_contam.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")
#drop blank samples; last 6 samples in this dataset (five water blanks plus mock community)
count_tab <- count_tab[1:(length(count_tab)-6)]
#import taxonomy
tax_tab <- as.matrix(read.table("03.14.22_Mississippi_data_microbiome_ASVs_taxonomy_no_contam.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))
#import metadata
sample_info_tab <- read_csv("/Users/adamo010/Documents/Blekhman_lab_external_consulting_projects/Mississippi_CRC_microbiome_project_2021/Metadata_and_extra_stuff/Mississippi_seqsonly_metadata.csv")
#make a deseq2 object
library(DESeq2)
#deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~type) #okay, troublesho0ot here.

#or, just make a phyloseq object.
sample_info_tab <- column_to_rownames(sample_info_tab, var = "Seq_sample")
rownames(metadata) #okay good, now sample names in table match metadata
rownames(sample_info_tab)
colnames(count_tab); #need to edit these... #####START HERE
#convert column names of count_tab to vector. Then, apply below code. 
sample.names <- sapply(strsplit(colnames(count_tab), "_F_trimmed_filtered.fastq.gz"), `[`, 1) #this creates a list called sample_names based on fnFs, which has _R1_001.fastq.gz trimmed
#then, replace column names with edited code below:
colnames(count_tab) <- sample.names #nice. this replaces the column names. 

ps_reloaded <- phyloseq(otu_table(count_tab, taxa_are_rows=TRUE), 
               sample_data(sample_info_tab), 
               tax_table(tax_tab)) #troubleshoot here
#neat.
#one thing I want to do is subset the samples too. Because we have such chaos in the number of samples we have, I want to be able to only look at
#paired samples, for example. Or cancer samples.

sample_variables(ps_reloaded) #print possible filters
ps_pairedonly = subset_samples(ps_reloaded, paired=='y') #DOUBLE == REQUIRED HERE!!!
ps_canceronly = subset_samples(ps_reloaded, Cancer_or_other=='cancer')
#this filtering did not work.

#######step 19 running phyloseq: looking at alpha diversity
plot_richness(ps_reloaded, x="Sample_description", measures=c("Shannon", "Simpson"), color="Sample_description",
              title="Tumor v.normal alpha diversity, all samples")
plot_richness(ps_pairedonly, x="Sample_description", measures=c("Shannon", "Simpson"), color="Sample_description",
              title="Tumor v.normal alpha diversity, paired samples only")
plot_richness(ps_canceronly, x="Sample_description", measures=c("Shannon", "Simpson"), color="Sample_description",
              title="Tumor v.normal alpha diversity, cancer samples only")

plot_richness(ps_reloaded, x="Ethnicity", measures=c("Shannon", "Simpson"), color="Ethnicity",
              title="Alpha diversity by ethnicity, all samples")
plot_richness(ps_pairedonly, x="Ethnicity", measures=c("Shannon", "Simpson"), color="Ethnicity",
              title="Alpha diversity by ethnicity, paired samples only")
plot_richness(ps_canceronly, x="Ethnicity", measures=c("Shannon", "Simpson"), color="Ethnicity",
              title="Alpha diversity by ethnicity, cancer samples only")

plot_richness(ps_reloaded, x="Cancer_or_other", measures=c("Shannon", "Simpson"), color="Cancer_or_other",
              title="Alpha diversity by diagnosis, all samples")
plot_richness(ps_reloaded, x="Age_at_dx", measures=c("Shannon", "Simpson"), color="Ethnicity",
              title="Alpha diversity by diagnosis, all samples")
plot_richness(ps_canceronly, x="Age_at_dx", measures=c("Shannon", "Simpson"), color="Ethnicity",
              title="Alpha diversity by age at diagnosis, cancer samples only")
plot_richness(ps_canceronly, x="Stage", measures=c("Shannon", "Simpson"), color="Ethnicity",
              title="Alpha diversity by stage, cancer samples only")

#######step 20: ordinate
# Transform data to proportions as appropriate for Bray-Curtis distances
ps_reloaded.prop <- transform_sample_counts(ps_reloaded, function(otu) otu/sum(otu))
ps_pairedonly.prop <- transform_sample_counts(ps_pairedonly, function(otu) otu/sum(otu))
ps_canceronly.prop <- transform_sample_counts(ps_canceronly, function(otu) otu/sum(otu))

ord.nmds.bray.reloaded <- ordinate(ps_reloaded.prop, method="NMDS", distance="bray") #oops.no convergance
ord.nmds.bray.pairedonly <- ordinate(ps_pairedonly.prop, method="NMDS", distance="bray") #oops.no convergance
ord.nmds.bray.canceronly <- ordinate(ps_canceronly.prop, method="NMDS", distance="bray") #oops.no convergance

ord.PCoA.bray.reloaded <- ordinate(ps_reloaded.prop, method="PCoA", distance="bray") #oops.no convergance
ord.PCoA.bray.pairedonly <- ordinate(ps_pairedonly.prop, method="PCoA", distance="bray") #oops.no convergance
ord.PCoA.bray.canceronly <- ordinate(ps_canceronly.prop, method="PCoA", distance="bray") #oops.no convergance

plot_ordination(ps_reloaded.prop, ord.PCoA.bray.reloaded, color="Sample_description", title="Bray PCOA by Description, all samples")
plot_ordination(ps_pairedonly.prop, ord.PCoA.bray.pairedonly, color="Sample_description", title="Bray PCOA by Description, paired samples only")
plot_ordination(ps_canceronly.prop, ord.PCoA.bray.canceronly, color="Sample_description", title="Bray PCOA by Description, cancer samples only")

plot_ordination(ps_reloaded.prop, ord.PCoA.bray.reloaded, color="Gender", title="Bray PCOA by Gender, all samples")
plot_ordination(ps_pairedonly.prop, ord.PCoA.bray.pairedonly, color="Gender", title="Bray PCOA by Gender, paired samples only")
plot_ordination(ps_canceronly.prop, ord.PCoA.bray.canceronly, color="Gender", title="Bray PCOA by Gender, cancer samples only")

plot_ordination(ps_reloaded.prop, ord.PCoA.bray.reloaded, color="Ethnicity", title="Bray PCOA by Ethnicity, all samples")
plot_ordination(ps_pairedonly.prop, ord.PCoA.bray.pairedonly, color="Ethnicity", title="Bray PCOA by Ethnicity, paired samples only")
plot_ordination(ps_canceronly.prop, ord.PCoA.bray.canceronly, color="Ethnicity", title="Bray PCOA by Ethnicity, cancer samples only")

plot_ordination(ps_reloaded.prop, ord.PCoA.bray.reloaded, color="Location", title="Bray PCOA by Location, all samples")
plot_ordination(ps_pairedonly.prop, ord.PCoA.bray.pairedonly, color="Location", title="Bray PCOA by Location, paired samples only")
plot_ordination(ps_canceronly.prop, ord.PCoA.bray.canceronly, color="Location", title="Bray PCOA by Location, cancer samples only")

plot_ordination(ps_reloaded.prop, ord.PCoA.bray.reloaded, color="Cancer_or_other", title="Bray PCOA by Location, all samples")

#meh

##########step 21: plot taxa
#filter by top 20 taxa
top20_reloaded <- names(sort(taxa_sums(ps_reloaded), decreasing=TRUE))[1:20]
top20_pairedonly <- names(sort(taxa_sums(ps_pairedonly), decreasing=TRUE))[1:20]
top20_canceronly <- names(sort(taxa_sums(ps_canceronly), decreasing=TRUE))[1:20]

#transform to relative abundance
standf_reloaded <- transform_sample_counts(ps_reloaded, function(OTU) OTU/sum(OTU))
standf_pairedonly <- transform_sample_counts(ps_pairedonly, function(OTU) OTU/sum(OTU))
standf_canceronly <- transform_sample_counts(ps_canceronly, function(OTU) OTU/sum(OTU))

#prune
standf_reloaded <- prune_taxa(top20_reloaded, standf_reloaded)
standf_pairedonly <- prune_taxa(top20_pairedonly, standf_pairedonly)
standf_canceronly <- prune_taxa(top20_canceronly, standf_canceronly)

#now, graph

plot_bar(standf_reloaded, x="Sample_description", fill="phylum") + geom_bar(stat="identity") + facet_wrap(~Location, scales="free_x")
plot_bar(standf_reloaded, x="Ethnicity", fill="phylum") + geom_bar(stat="identity") + facet_wrap(~Cancer_or_other, scales="free_x")
plot_bar(standf_canceronly, x="Ethnicity", fill="phylum") + 
  geom_bar(stat="identity") +
  ggtitle("Phylum abundance by ethnicity, cancer samples only")
plot_bar(standf_reloaded, x="sequenced_samples", fill="phylum") + 
  #geom_bar(stat="identity") +
  ggtitle("Phylum abundance by sample, all samples")

#try this agian, without filtering.
#transform to relative abundance
standf_reloaded2 <- transform_sample_counts(ps_reloaded, function(OTU) OTU/sum(OTU))
standf_pairedonly2 <- transform_sample_counts(ps_pairedonly, function(OTU) OTU/sum(OTU))
standf_canceronly2 <- transform_sample_counts(ps_canceronly, function(OTU) OTU/sum(OTU))

plot_bar(standf_reloaded2, x="Sample_description", fill="phylum") + geom_bar(stat="identity") + facet_wrap(~Location, scales="free_x")
plot_bar(standf_reloaded2, x="Ethnicity", fill="phylum") + geom_bar(stat="identity") + facet_wrap(~Cancer_or_other, scales="free_x")
plot_bar(standf_canceronly2, x="Ethnicity", fill="phylum") + 
  geom_bar(stat="identity") +
  ggtitle("Phylum abundance by ethnicity, cancer samples only")
plot_bar(standf_reloaded2, x="sequenced_samples", fill="phylum") + 
  #geom_bar(stat="identity") +
  ggtitle("Phylum abundance by sample, all samples")

#non-filtering sucks

#filter by top 100 taxa
top100_reloaded <- names(sort(taxa_sums(ps_reloaded), decreasing=TRUE))[1:100]
top100_pairedonly <- names(sort(taxa_sums(ps_pairedonly), decreasing=TRUE))[1:100]
top100_canceronly <- names(sort(taxa_sums(ps_canceronly), decreasing=TRUE))[1:100]

#transform to relative abundance
standf_reloaded3 <- transform_sample_counts(ps_reloaded, function(OTU) OTU/sum(OTU))
standf_pairedonly3 <- transform_sample_counts(ps_pairedonly, function(OTU) OTU/sum(OTU))
standf_canceronly3 <- transform_sample_counts(ps_canceronly, function(OTU) OTU/sum(OTU))

#prune
standf_reloaded3 <- prune_taxa(top100_reloaded, standf_reloaded3)
standf_pairedonly3 <- prune_taxa(top100_pairedonly, standf_pairedonly3)
standf_canceronly3 <- prune_taxa(top100_canceronly, standf_canceronly3)

#graph
plot_bar(standf_reloaded3, x="Sample_description", fill="phylum") + geom_bar(stat="identity") + facet_wrap(~Location, scales="free_x")
plot_bar(standf_reloaded3, x="Ethnicity", fill="phylum") + geom_bar(stat="identity") + facet_wrap(~Cancer_or_other, scales="free_x")
plot_bar(standf_canceronly3, x="Ethnicity", fill="phylum") + 
  geom_bar(stat="identity") +
  ggtitle("Phylum abundance by ethnicity, cancer samples only")
plot_bar(standf_reloaded3, x="sequenced_samples", fill="phylum") + 
  geom_bar(stat="identity") +
  aes(x = reorder(sequenced_samples, -Abundance), y = Abundance) +
  ggtitle("Phylum abundance by sample, all samples") 
plot_bar(standf_canceronly3, x="sequenced_samples", fill="phylum") + 
  geom_bar(stat="identity") +
  aes(x = reorder(sequenced_samples, -Abundance), y = Abundance) +
  ggtitle("Phylum abundance by sample, cancer samples only")

plot_bar(standf_pairedonly3, x="sequenced_samples", fill="phylum") + 
  geom_bar(stat="identity") +
  aes(x = reorder(sequenced_samples, -Abundance), y = Abundance) +
  ggtitle("Phylum abundance by sample, cancer samples only") +
  facet_wrap(~Sample_description, scales="free_x")


#why is the abundance weird? shouldn't it be proportional? seems like abundances are being added across samples
#I hate these graphs. 


#########step 22: OKAY, save everything and... hope to find something new and excting later.

#######################################junk
#this standardizes abundances to median sequencing depth
total_reloaded = median(sample_sums(ps_reloaded))
total_pairedonly = median(sample_sums(ps_pairedonly))
total_canceronly = median(sample_sums(ps_canceronly))

standf_reloaded = function(x, t=total_reloaded) round(t * (x / sum(x)))
standf_pairedonly = function(x, t=total_pairedonly) round(t * (x / sum(x)))
standf_canceronly = function(x, t=total_canceronly) round(t * (x / sum(x)))



gps_reloaded = transform_sample_counts(ps_reloaded, standf_reloaded)
gps_pairedonly = transform_sample_counts(ps_pairedonly, standf_pairedonly)
gps_canceronly = transform_sample_counts(ps_canceronly, standf_canceronly)

#now, just use top 20 taxa
top20_reloaded <- names(sort(taxa_sums(ps_reloaded), decreasing=TRUE))[1:20]
top20_pairedonly <- names(sort(taxa_sums(ps_pairedonly), decreasing=TRUE))[1:20]
top20_canceronly <- names(sort(taxa_sums(ps_canceronly), decreasing=TRUE))[1:20]

ps_reloaded.top20 <- prune_taxa(top20_reloaded, gps_reloaded)
ps_pairedonly.top20 <- prune_taxa(top20_pairedonly, gps_pairedonly)
ps_canceronly.top20 <- prune_taxa(top20_canceronly, gps_canceronly)



