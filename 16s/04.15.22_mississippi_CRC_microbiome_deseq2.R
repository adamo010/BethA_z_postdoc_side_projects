library("DESeq2")
library("dada2")
library("stringr")
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
library("viridis")

#Today I am going to run deseq2 on the Mississippi CRC microbiome data

#######step 1: prep work
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/Blekhman_lab_external_consulting_projects/Mississippi_CRC_microbiome_project_2021/") 

##############step 1: import counts data and metadata
countData <- read.table("03.14.22_Mississippi_data_microbiome_ASVs_counts_no_contam.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")
#edit countData column headers to be useful
sample.names <- sapply(strsplit(colnames(countData), "_F_trimmed_filtered.fastq.gz"), `[`, 1) #this creates a list called sample_names based on fnFs, which has _R1_001.fastq.gz trimmed
colnames(countData) <- sample.names #nice. this replaces the column names. 
#NEW: add and merge in taxonomic data. Doesn't really help, but is an option. 
tax_info <- read.table(file="03.11.22_DADA2_analysis_Mississippi_data_microbiome_ASVs_taxonomy.tsv", sep="\t", header=TRUE) 
countData2 <- tibble::rownames_to_column(countData, "X") #set rownames to a column called X
countData2 <- dplyr::left_join(countData2, tax_info, by= "X") #merge based on column called X
countData2 <- column_to_rownames(countData2, "X")

#import metadata
metaData <- read_csv("/Users/adamo010/Documents/Blekhman_lab_external_consulting_projects/Mississippi_CRC_microbiome_project_2021/Metadata_and_extra_stuff/Mississippi_seqsonly_metadata.csv")
#may need to set patient_blind_id as factor
metaData$Patient_ID <- as.factor(metaData$Patient_ID)
#drop everything in metadata that's not specifically necessary
metaData2 = subset(metaData, select= c(Patient_ID, Seq_sample, Sample_description, Ethnicity, Cancer_or_other, paired))
metaData2 <- column_to_rownames(metaData2, 'Seq_sample') #set Seq_sample to rownames

##############step 2: filter data to only include paired, cancer, ethnicity=Black/white data
#start by filtering metadata, then use that to filter countdata
metaData_filt_paired <- filter(metaData2, Cancer_or_other == "cancer" & paired=="y" & Ethnicity != "Hispanic")
metaData_filt_paired_ids <- rownames(metaData_filt_paired) #extract row names as vector
countData_filt <- countData[, which((names(countData) %in% metaData_filt_paired_ids)==TRUE)] #then, filter countdata by the list of metadata rownames

##############step 3: construct deseq2 object:
dds_paired <- DESeqDataSetFromMatrix(countData=countData_filt, 
                              colData=metaData_filt_paired, 
                              design = ~factor(Patient_ID) + Sample_description) #this part indicates that we have paired data.
#got a warning message that: In DESeqDataSet(se, design = design, ignoreRank)::some variables in design formula are characters, converting to factors"
#which should be fine
#ALSO: factor levels were dropped which had no samples
#hmmm, there was an error estimating size factors due to every gene containing at least one zero. I guess that means we need pseudocounts?
dds_from_deseq <- DESeq(dds_paired, test = "Wald", sfType = "poscounts") #poscounts is a solution to this. it changes the estimatesizefactors parameter

#look at the results
res <- results(dds_from_deseq)
#head(results(dds_from_deseq, tidy=TRUE)) #let's look at the results table
#summary(res) #summary of results for differential ASV abundance
res <- res[order(res$padj),] #re-order by p-value
head(res)

#pull out significant results
sigtab = res[which(res$padj < 0.05), ]
#import taxonomic info
tax_info <- read.table(file="03.11.22_DADA2_analysis_Mississippi_data_microbiome_ASVs_taxonomy.tsv", sep="\t", header=TRUE) 
sigtab_df <- as(sigtab, "data.frame") #convert to dataframe
sigtab_df <- tibble::rownames_to_column(sigtab_df, "X") #set rownames to a column called X
sigtab_df <- dplyr::left_join(sigtab_df, tax_info, by= "X")

#now we can plot:
require("ggrepel")
theme_set(theme_bw())
ggplot(sigtab_df, aes(x=genus, y=log2FoldChange, color=order)) + geom_point(size=5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.4, vjust=0.5)) +
  scale_y_continuous(name = "log2FoldChange in taxon abundance\n normal(base) vs tumor", limits=c(0,30)) +
  scale_x_discrete(name = "Genus") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = expression("Differentially abundant taxa, normal vs tumor"),
       subtitle = expression("Design: Patient_ID + Sample_description")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  geom_label_repel(aes(label = signif(padj, digits=3)), color="black", box.padding= 0.5, size=3, xlim = c(NA, Inf)) 
  
#plot individual taxa
#graph the top six taxa. Most of them are fuso! So that's reassuring. Just for looking
par(mfrow=c(2,3))
plotCounts(dds_from_deseq, gene="ASV_254", intgroup="Sample_description")
plotCounts(dds_from_deseq, gene="ASV_641", intgroup="Sample_description")
plotCounts(dds_from_deseq, gene="ASV_1087", intgroup="Sample_description")
plotCounts(dds_from_deseq, gene="ASV_830", intgroup="Sample_description")
plotCounts(dds_from_deseq, gene="ASV_391", intgroup="Sample_description")
plotCounts(dds_from_deseq, gene="ASV_212", intgroup="Sample_description")

##############step 4: make a PCA plot:
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation
#getting an error here too: " it is recommended to use varianceStabilizingTransformation directly"
#see solution here: https://www.biostars.org/p/456209/: decrease nsub

vsdata <- vst(dds_from_deseq, blind=FALSE, nsub=500) #nsub default is 1000; 500 is the biggest working value
pca1 <-plotPCA(vsdata, intgroup="Sample_description") #using the DESEQ2 plotPCA fxn 
pca2 <- plotPCA(vsdata, intgroup=c("Sample_description", "Ethnicity"))
pca3 <- plotPCA(vsdata, intgroup=c("Ethnicity"))
pca4 <- plotPCA(vsdata, intgroup=c("Patient_ID"))

pca1 + labs(title = expression("PCA of taxon abundance"),
              subtitle = expression("Design: Patient_ID + Sample_description")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_fill_brewer(palette = "Set1") +  scale_color_brewer(palette = "Set1") 

pca2 + labs(title = expression("PCA of taxon abundance"),
            subtitle = expression("Design: Patient_ID + Sample_description")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_fill_brewer(palette = "Set1") +  scale_color_brewer(palette = "Set1")

pca3 + labs(title = expression("PCA of taxon abundance"),
            subtitle = expression("Design: Patient_ID + Sample_description")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_fill_viridis(discrete=TRUE) +  scale_color_viridis(discrete=TRUE)

pca4 + labs(title = expression("PCA of taxon abundance, by patient"),
            subtitle = expression("Design: Patient_ID + Sample_description")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), legend.position="below") +
  scale_fill_viridis(discrete=TRUE) +  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=Patient_ID))

#############################################################################################################################
#a new design

##############step 2: filter data to only include paired, cancer, ethnicity=Black/white data
#start by filtering metadata, then use that to filter countdata
metaData_filt_cancer <- filter(metaData2, Cancer_or_other == "cancer" & Ethnicity != "Hispanic")
metaData_filt_cancer_ids <- rownames(metaData_filt_cancer) #extract row names as vector
countData_filt_cancer <- countData[, which((names(countData) %in% metaData_filt_cancer_ids)==TRUE)] #then, filter countdata by the list of metadata rownames

##############step 3: construct deseq2 object:
dds_cancer <- DESeqDataSetFromMatrix(countData=countData_filt_cancer, 
                                     colData=metaData_filt_cancer, 
                                     design = ~Sample_description + Ethnicity + Sample_description:Ethnicity)
#got a warning message that: In DESeqDataSet(se, design = design, ignoreRank)::some variables in design formula are characters, converting to factors"
#which should be fine
#ALSO: factor levels were dropped which had no samples
#hmmm, there was an error estimating size factors due to every gene containing at least one zero. I guess that means we need pseudocounts?
dds_from_deseq_cancer <- DESeq(dds_cancer, test = "Wald", sfType = "poscounts") #poscounts is a solution to this. it changes the estimatesizefactors parameter

#look at the results
res_cancer <- results(dds_from_deseq_cancer)
#head(results(dds_from_deseq, tidy=TRUE)) #let's look at the results table
#summary(res) #summary of results for differential ASV abundance
res_cancer <- res_cancer[order(res_cancer$padj),] #re-order by p-value
head(res_cancer)

#pull out significant results
sigtab_cancer = res_cancer[which(res_cancer$padj < 0.05), ]
#import taxonomic info
tax_info <- read.table(file="03.11.22_DADA2_analysis_Mississippi_data_microbiome_ASVs_taxonomy.tsv", sep="\t", header=TRUE) 
sigtab_df_cancer <- as(sigtab_cancer, "data.frame") #convert to dataframe
sigtab_df_cancer <- tibble::rownames_to_column(sigtab_df_cancer, "X") #set rownames to a column called X
sigtab_df_cancer <- dplyr::left_join(sigtab_df_cancer, tax_info, by= "X")

#now we can plot:
require("ggrepel")
theme_set(theme_bw())
ggplot(sigtab_df_cancer, aes(x=genus, y=log2FoldChange, color=order)) + geom_point(size=4) + 
  theme(axis.text.x = element_text(vjust=0.5, hjust = 0.5, angle=90), panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(name = "log2FoldChange in taxon abundance\n normal(base) vs tumor") +
  scale_x_discrete(name = "Genus", expand=c(-1.55,1.55)) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = expression("Differentially abundant taxa, tumor vs normal"),
       subtitle = expression("Sample_description + Ethnicity + Sample_description*Ethnicity")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  geom_label_repel(aes(label = signif(padj, digits=3)), color="black", box.padding= 0.6, size=2.5, xlim = c(NA, Inf),
                   point.padding=0.3) 

#plot individual taxa
#graph the top six taxa. Most of them are fuso! So that's reassuring. Just for looking
par(mfrow=c(2,3))
plotCounts(dds_from_deseq_cancer, gene="ASV_85", intgroup="Sample_description")
plotCounts(dds_from_deseq_cancer, gene="ASV_179", intgroup="Sample_description")
plotCounts(dds_from_deseq_cancer, gene="ASV_48", intgroup="Sample_description")
plotCounts(dds_from_deseq_cancer, gene="ASV_543", intgroup="Sample_description")
plotCounts(dds_from_deseq_cancer, gene="ASV_469", intgroup="Sample_description")
plotCounts(dds_from_deseq_cancer, gene="ASV_143", intgroup="Sample_description")

##############step 4: make a PCA plot:
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation
#getting an error here too: " it is recommended to use varianceStabilizingTransformation directly"
#see solution here: https://www.biostars.org/p/456209/: decrease nsub

vsdata <- vst(dds_from_deseq_cancer, blind=FALSE, nsub=500) #nsub default is 1000; 500 is the biggest working value
pca1 <-plotPCA(vsdata, intgroup="Sample_description") #using the DESEQ2 plotPCA fxn 
pca2 <- plotPCA(vsdata, intgroup=c("Sample_description", "Ethnicity"))
pca3 <- plotPCA(vsdata, intgroup=c("Ethnicity"))
pca4 <- plotPCA(vsdata, intgroup=c("Patient_ID"))

pca1 + labs(title = expression("PCA of taxon abundance"),
            subtitle = expression("Design: Sample_description + Ethnicity + Sample_description:Ethnicity")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_fill_brewer(palette = "Set1") +  scale_color_brewer(palette = "Set1") 

pca2 + labs(title = expression("PCA of taxon abundance"),
            subtitle = expression("Design: Sample_description + Ethnicity + Sample_description:Ethnicity")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_fill_brewer(palette = "Set1") +  scale_color_brewer(palette = "Set1")

pca3 + labs(title = expression("PCA of taxon abundance"),
            subtitle = expression("Design: Sample_description + Ethnicity + Sample_description:Ethnicity")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_fill_viridis(discrete=TRUE) +  scale_color_viridis(discrete=TRUE)

pca4 + labs(title = expression("PCA of taxon abundance, by patient"),
            subtitle = expression("Design: Patient_ID + Sample_description")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), legend.position="below") +
  scale_fill_viridis(discrete=TRUE) +  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=Patient_ID))

#############################################################################################################################
#a new design by ethnicity

##############step 2: filter data 
#start by filtering metadata, then use that to filter countdata
metaData_filt_Black <- filter(metaData2, Cancer_or_other == "cancer" & Ethnicity == "Black")
metaData_filt_Black_ids <- rownames(metaData_filt_Black) #extract row names as vector
countData_filt_Black <- countData[, which((names(countData) %in% metaData_filt_Black_ids)==TRUE)] #then, filter countdata by the list of metadata rownames

##############step 3: construct deseq2 object:
dds_Black <- DESeqDataSetFromMatrix(countData=countData_filt_Black, 
                                     colData=metaData_filt_Black, 
                                     design = ~factor(Patient_ID) + Sample_description)
#got a warning message that: In DESeqDataSet(se, design = design, ignoreRank)::some variables in design formula are characters, converting to factors"
#which should be fine
#ALSO: factor levels were dropped which had no samples
#hmmm, there was an error estimating size factors due to every gene containing at least one zero. I guess that means we need pseudocounts?
dds_from_deseq_Black <- DESeq(dds_Black, test = "Wald", sfType = "poscounts") #poscounts is a solution to this. it changes the estimatesizefactors parameter

#look at the results
res_Black <- results(dds_from_deseq_Black)
#head(results(dds_from_deseq, tidy=TRUE)) #let's look at the results table
#summary(res) #summary of results for differential ASV abundance
res_Black <- res_Black[order(res_Black$padj),] #re-order by p-value
head(res_Black)

#pull out significant results
sigtab_Black = res_Black[which(res_Black$padj < 0.05), ]
#import taxonomic info
tax_info <- read.table(file="03.11.22_DADA2_analysis_Mississippi_data_microbiome_ASVs_taxonomy.tsv", sep="\t", header=TRUE) 
sigtab_df_Black <- as(sigtab_Black, "data.frame") #convert to dataframe
sigtab_df_Black <- tibble::rownames_to_column(sigtab_df_Black, "X") #set rownames to a column called X
sigtab_df_Black <- dplyr::left_join(sigtab_df_Black, tax_info, by= "X")

#now we can plot:
require("ggrepel")
theme_set(theme_bw())
ggplot(sigtab_df_Black, aes(x=genus, y=log2FoldChange, color=phylum)) + geom_point(size=4) + 
  theme(axis.text.x = element_text(vjust=0.5, hjust = 0.5), panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(name = "log2FoldChange in taxon abundance\n normal(base) vs tumor") +
  scale_x_discrete(name = "Genus") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = expression("Differentially abundant taxa, tumor vs normal,\nBlack patients only"),
       subtitle = expression("Patient_ID + Sample_description")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  geom_label_repel(aes(label = signif(padj, digits=3)), color="black", box.padding= 0.6, size=2.5, xlim = c(NA, Inf),
                   point.padding=0.3) 

#plot individual taxa
#graph the top six taxa. Most of them are fuso! So that's reassuring. Just for looking
par(mfrow=c(1,2))
plotCounts(dds_from_deseq_Black, gene="ASV_254", intgroup="Sample_description")
plotCounts(dds_from_deseq_Black, gene="ASV_317", intgroup="Sample_description")

##############step 4: make a PCA plot:
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation
#getting an error here too: " it is recommended to use varianceStabilizingTransformation directly"
#see solution here: https://www.biostars.org/p/456209/: decrease nsub
vsdata <- vst(dds_from_deseq_Black, blind=FALSE, nsub=500) #nsub default is 1000; 500 is the biggest working value
pca1 <-plotPCA(vsdata, intgroup="Sample_description") #using the DESEQ2 plotPCA fxn 
pca4 <- plotPCA(vsdata, intgroup=c("Patient_ID"))

pca1 + labs(title = expression("PCA of taxon abundance, Black patients only"),
            subtitle = expression("Design: SPatient_ID + Sample_description")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_fill_brewer(palette = "Set1") +  scale_color_brewer(palette = "Set1") 

pca4 + labs(title = expression("PCA of taxon abundance, by patient, Black patients only"),
            subtitle = expression("Design: Patient_ID + Sample_description")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), legend.position="below") +
  scale_fill_viridis(discrete=TRUE) +  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=Patient_ID))

#############################################################################################################################
#a new design by ethnicity

##############step 2: filter data 
#start by filtering metadata, then use that to filter countdata
metaData_filt_White <- filter(metaData2, Cancer_or_other == "cancer" & Ethnicity == "White")
metaData_filt_White_ids <- rownames(metaData_filt_White) #extract row names as vector
countData_filt_White <- countData[, which((names(countData) %in% metaData_filt_White_ids)==TRUE)] #then, filter countdata by the list of metadata rownames

##############step 3: construct deseq2 object:
dds_White <- DESeqDataSetFromMatrix(countData=countData_filt_White, 
                                    colData=metaData_filt_White, 
                                    design = ~factor(Patient_ID) + Sample_description)
#got a warning message that: In DESeqDataSet(se, design = design, ignoreRank)::some variables in design formula are characters, converting to factors"
#which should be fine
#ALSO: factor levels were dropped which had no samples
#hmmm, there was an error estimating size factors due to every gene containing at least one zero. I guess that means we need pseudocounts?
dds_from_deseq_White <- DESeq(dds_White, test = "Wald", sfType = "poscounts") #poscounts is a solution to this. it changes the estimatesizefactors parameter

#look at the results
res_White <- results(dds_from_deseq_White)
#head(results(dds_from_deseq, tidy=TRUE)) #let's look at the results table
#summary(res) #summary of results for differential ASV abundance
res_White <- res_White[order(res_White$padj),] #re-order by p-value
head(res_White)

#pull out significant results
sigtab_White = res_White[which(res_White$padj < 0.05), ]
#import taxonomic info
tax_info <- read.table(file="03.11.22_DADA2_analysis_Mississippi_data_microbiome_ASVs_taxonomy.tsv", sep="\t", header=TRUE) 
sigtab_df_White <- as(sigtab_White, "data.frame") #convert to dataframe
sigtab_df_White <- tibble::rownames_to_column(sigtab_df_White, "X") #set rownames to a column called X
sigtab_df_White <- dplyr::left_join(sigtab_df_White, tax_info, by= "X")

#now we can plot:
require("ggrepel")
theme_set(theme_bw())
ggplot(sigtab_df_White, aes(x=genus, y=log2FoldChange, color=phylum)) + geom_point(size=4) + 
  theme(axis.text.x = element_text(vjust=0.5, hjust = 0.5), panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(name = "log2FoldChange in taxon abundance\n normal(base) vs tumor") +
  scale_x_discrete(name = "Genus") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = expression("Differentially abundant taxa, tumor vs normal,\nwhite patients only"),
       subtitle = expression("Patient_ID + Sample_description")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  geom_label_repel(aes(label = signif(padj, digits=3)), color="black", box.padding= 0.6, size=2.5, xlim = c(NA, Inf),
                   point.padding=0.3) 

#plot individual taxa
#graph the top six taxa. Most of them are fuso! So that's reassuring. Just for looking
par(mfrow=c(1,1))
plotCounts(dds_from_deseq_White, gene="ASV_212", intgroup="Sample_description")
#plotCounts(dds_from_deseq_White, gene="ASV_317", intgroup="Sample_description")

##############step 4: make a PCA plot:
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation
#getting an error here too: " it is recommended to use varianceStabilizingTransformation directly"
#see solution here: https://www.biostars.org/p/456209/: decrease nsub
vsdata <- vst(dds_from_deseq_White, blind=FALSE, nsub=300) #nsub default is 1000; 500 is the biggest working value
pca1 <-plotPCA(vsdata, intgroup="Sample_description") #using the DESEQ2 plotPCA fxn 
pca4 <- plotPCA(vsdata, intgroup=c("Patient_ID"))

pca1 + labs(title = expression("PCA of taxon abundance, white patients only"),
            subtitle = expression("Design: SPatient_ID + Sample_description")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_fill_brewer(palette = "Set1") +  scale_color_brewer(palette = "Set1") 

pca4 + labs(title = expression("PCA of taxon abundance, by patient, white patients only"),
            subtitle = expression("Design: Patient_ID + Sample_description")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), legend.position="below") +
  scale_fill_viridis(discrete=TRUE) +  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=Patient_ID))

#############################################################################################################################
#a new design to test controls

##############step 2: filter data 
#start by filtering metadata, then use that to filter countdata
metaData_filt_all <- filter(metaData2, paired=="y")
metaData_filt_all_ids <- rownames(metaData_filt_all) #extract row names as vector
countData_filt_all <- countData[, which((names(countData) %in% metaData_filt_all_ids)==TRUE)] #then, filter countdata by the list of metadata rownames

##############step 3: construct deseq2 object:
dds_all <- DESeqDataSetFromMatrix(countData=countData_filt_all, 
                                    colData=metaData_filt_all, 
                                    design = ~Sample_description + Ethnicity + Cancer_or_other + Sample_description:Cancer_or_other)
#got a warning message that: In DESeqDataSet(se, design = design, ignoreRank)::some variables in design formula are characters, converting to factors"
#which should be fine
#ALSO: factor levels were dropped which had no samples
#hmmm, there was an error estimating size factors due to every gene containing at least one zero. I guess that means we need pseudocounts?
dds_from_deseq_all <- DESeq(dds_all, test = "Wald", sfType = "poscounts") #poscounts is a solution to this. it changes the estimatesizefactors parameter

#look at the results
res_all <- results(dds_from_deseq_all)
#head(results(dds_from_deseq, tidy=TRUE)) #let's look at the results table
#summary(res) #summary of results for differential ASV abundance
res_all <- res_all[order(res_all$padj),] #re-order by p-value
head(res_all)

#pull out significant results
sigtab_all = res_all[which(res_all$padj < 0.05), ]
#import taxonomic info
tax_info <- read.table(file="03.11.22_DADA2_analysis_Mississippi_data_microbiome_ASVs_taxonomy.tsv", sep="\t", header=TRUE) 
sigtab_df_all <- as(sigtab_all, "data.frame") #convert to dataframe
sigtab_df_all <- tibble::rownames_to_column(sigtab_df_all, "X") #set rownames to a column called X
sigtab_df_all <- dplyr::left_join(sigtab_df_all, tax_info, by= "X")

#now we can plot:
require("ggrepel")
theme_set(theme_bw())
ggplot(sigtab_df_all, aes(x=genus, y=log2FoldChange, color=phylum)) + geom_point(size=4) + 
  theme(axis.text.x = element_text(vjust=0.5, hjust = 0.5), panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(name = "log2FoldChange in taxon abundance\n normal(base) vs tumor") +
  scale_x_discrete(name = "Genus") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = expression("Differentially abundant taxa, tumor vs normal,\nwhite patients only"),
       subtitle = expression("Patient_ID + Sample_description")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  geom_label_repel(aes(label = signif(padj, digits=3)), color="black", box.padding= 0.6, size=2.5, xlim = c(NA, Inf),
                   point.padding=0.3) 

#plot individual taxa
#graph the top six taxa. Most of them are fuso! So that's reassuring. Just for looking
par(mfrow=c(1,1))
plotCounts(dds_from_deseq_White, gene="ASV_212", intgroup="Sample_description")
#plotCounts(dds_from_deseq_White, gene="ASV_317", intgroup="Sample_description")

##############step 4: make a PCA plot:
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation
#getting an error here too: " it is recommended to use varianceStabilizingTransformation directly"
#see solution here: https://www.biostars.org/p/456209/: decrease nsub
vsdata <- vst(dds_from_deseq_all, blind=FALSE, nsub=500) #nsub default is 1000; 500 is the biggest working value
pca1 <-plotPCA(vsdata, intgroup="Sample_description") #using the DESEQ2 plotPCA fxn 
pca2 <- plotPCA(vsdata, intgroup=c("Sample_description", "Ethnicity"))
pca5 <- plotPCA(vsdata, intgroup="Cancer_or_other") 
pca6 <- plotPCA(vsdata, intgroup=c("Cancer_or_other", "Sample_description"))
pca7 <- plotPCA(vsdata, intgroup=c("Ethnicity", "Sample_description"))

pca1 + labs(title = expression("PCA of taxon abundance"),
            subtitle = expression("Design: Sample_description + Ethnicity + Sample_description:Ethnicity")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_fill_brewer(palette = "Set1") +  scale_color_brewer(palette = "Set1") 

pca2 + labs(title = expression("PCA of taxon abundance"),
            subtitle = expression("Design: Sample_description + Ethnicity + Sample_description:Ethnicity")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_fill_brewer(palette = "Set1") +  scale_color_brewer(palette = "Set1")

pca5 + labs(title = expression("PCA of taxon abundance, sample type"),
            subtitle = expression("Design: Sample_description + Ethnicity + Cancer_or_other + Sample_description:Cancer_or_other")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_fill_viridis(discrete=TRUE) +  scale_color_viridis(discrete=TRUE) 

pca6 + labs(title = expression("PCA of taxon abundance, sample type"),
            subtitle = expression("Design: Sample_description + Ethnicity + Cancer_or_other + Sample_description:Cancer_or_other")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_fill_viridis(discrete=TRUE) +  scale_color_viridis(discrete=TRUE) 

pca7 + labs(title = expression("PCA of taxon abundance, sample type"),
            subtitle = expression("Design: Sample_description + Ethnicity + Cancer_or_other + Sample_description:Cancer_or_other")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_fill_viridis(discrete=TRUE) +  scale_color_viridis(discrete=TRUE) 


