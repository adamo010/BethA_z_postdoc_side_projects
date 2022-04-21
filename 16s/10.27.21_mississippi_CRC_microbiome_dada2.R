#Hello! Today I am doing 16s analysis on samples from the Mississippi microbiome project. 

#this is based on 08.27.21_mayo_microbiome_dada2_attempt1_clean.R
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/Mississippi_microbiome_project_2021/") 


######step 1: load packages
#need to download first, but only need to do this once, fortunately. 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.13")
library(dada2); packageVersion("dada2") #this prints the package version that you're using. Here, I get 1.20.0.
library(stringr) #for editing strings; a tidyverse creation

######step 2: set the path to the fastq files
path <- "/Users/adamo010/Documents/Mississippi_microbiome_project_2021/"
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

#create output directory
path2 <- "/Users/adamo010/Documents/Mississippi_microbiome_project_2021/filtered/"

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
path <- "/Users/adamo010/Documents/Mississippi_microbiome_project_2021/filtered/"
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
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors. Still takes FOREVER
#######step 15: evaluate accuracy: can't find where I should have done this

#######step 16: creating outputs
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character") #this is to make the ASV names useful and not some horrific seq 
for (i in 1:dim(seqtab.nochim)[2]) {asv_headers[i] <- paste(">ASV", i, sep="_")}

#making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "10.29.21_DADA2_analysis_Mississippi_data_microbiome_ASVs.fa")

#count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "10.29.21_DADA2_analysis_Mississippi_data_microbiome_ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

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
write.table(taxid, "10.29.21_DADA2_analysis_Mississippi_data_microbiome_ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

#######step 17: removing contaminants
library(decontam)
colnames(asv_tab) # print column names to used to find the blanks
#have 102 samples; #1-96 are real, 97-101 are blanks, 102 is a mock community
vector_for_decontam <- c(rep(FALSE, 96), rep(TRUE, 5), rep(FALSE, 1)) #define filtering vector to id contams
contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam) #create new dataframe of contaminants
table(contam_df$contaminant) # identified 18 as contaminants (print as TRUE in this table)
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ]) #get vector of contaminant ids
contam_asv_table <- taxid[row.names(taxid) %in% contam_asvs, ] #convert to a taxonomy table
write.table(contam_asv_table, "10.29.21_DADA2_analysis_Mississippi_data_contaminant_ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)
