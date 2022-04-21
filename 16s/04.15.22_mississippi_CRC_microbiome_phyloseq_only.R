#Hello! Today I am doing 16s analysis on samples from the Mississippi microbiome project. 
######step 1: load packages
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
library("phyloseq")
library("Biostrings")

#######step 1: prep work
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/Blekhman_lab_external_consulting_projects/Mississippi_CRC_microbiome_project_2021/") 

######step 2: set the path to the fastq files
path <- "/Users/adamo010/Documents/Blekhman_lab_external_consulting_projects/Mississippi_CRC_microbiome_project_2021/"
head(list.files(path)) #print a few file names in the file path

count_tab <- read.table("03.14.22_Mississippi_data_microbiome_ASVs_counts_no_contam.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")
#drop blank samples; last 6 samples in this dataset (five water blanks plus mock community)
count_tab <- count_tab[1:(length(count_tab)-6)]
#import taxonomy
tax_tab <- as.matrix(read.table("03.14.22_Mississippi_data_microbiome_ASVs_taxonomy_no_contam.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))
#import metadata
sample_info_tab <- read_csv("/Users/adamo010/Documents/Blekhman_lab_external_consulting_projects/Mississippi_CRC_microbiome_project_2021/Metadata_and_extra_stuff/Mississippi_seqsonly_metadata.csv")

# make a phyloseq object.
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

#######Step 3: running phyloseq: looking at alpha diversity
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

########NEW analysis: look at alpha diversity within ethnicity
#remove Hispanic samples from dataset; not helpful and screws with statistics
ps_canceronly_bw <- subset_samples(ps_canceronly, Ethnicity!="Hispanic")

#run stats
canceronly_rich = estimate_richness(ps_canceronly_bw) #export dataframe of alpha diversity metrics
pairwise.wilcox.test(canceronly_rich$Shannon, sample_data(ps_canceronly_bw)$Sample_description)
#should just do those within the graphs.

#graph! one graph per metric. 
alpha1 <- plot_richness(ps_canceronly_bw, x="Sample_description", measures=c("Shannon"), color="Sample_description") +
  geom_boxplot(color="black", aes(fill= sample_data(ps_canceronly_bw)$Sample_description)) +
  geom_point(color="black") +
  theme_bw() +
  scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") +
  scale_y_continuous(name = "Shannon") +
  theme(legend.position='none', axis.title.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("tumor", "normal")), label = "p.signif", symnum.args = symnum.args) +
  facet_grid(~Ethnicity, scales="free")  +
  theme(strip.background = element_rect(color="black", size=3, linetype="blank"), strip.placement = "outside",
        strip.text.x = element_text(size = 10, color = "black", face = "bold"))
alpha2 <- plot_richness(ps_canceronly_bw, x="Sample_description", measures=c("Chao1"), color="Sample_description") +
  geom_boxplot(color="black", aes(fill= sample_data(ps_canceronly_bw)$Sample_description)) +
  geom_point(color="black") +
  theme_bw() +
  scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") +
  scale_y_continuous(name = "Chao1") +
  theme(legend.position='none', axis.title.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("tumor", "normal")), label = "p.signif", symnum.args = symnum.args) +
  facet_grid(~Ethnicity, scales="free")  +
  theme(strip.background = element_rect(color="black", size=3, linetype="blank"), strip.placement = "outside",
        strip.text.x = element_text(size = 10, color = "black", face = "bold"))
alpha3 <- plot_richness(ps_canceronly_bw, x="Sample_description", measures=c("Simpson"), color="Sample_description") +
  geom_boxplot(color="black", aes(fill= sample_data(ps_canceronly_bw)$Sample_description)) +
  geom_point(color="black") +
  theme_bw() +
  scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") +
  scale_y_continuous(name = "Simpson") +
  theme(legend.position='none', axis.title.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("tumor", "normal")), label = "p.signif", symnum.args = symnum.args) +
  facet_grid(~Ethnicity, scales="free_x") +
  theme(strip.background = element_rect(color="black", size=3, linetype="blank"), strip.placement = "outside",
        strip.text.x = element_text(size = 10, color = "black", face = "bold"))
#combine graphs
ggarrange(alpha1, alpha2, alpha3, 
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)

########NEW analysis: within tumor samples, compare Black and white patients
#only look at tumor samples
ps_tumoronly_bw <- subset_samples(ps_canceronly, Sample_description=="tumor")

#should just do those within the graphs.

#graph! one graph per metric. 
alpha4 <- plot_richness(ps_tumoronly_bw, x="Ethnicity", measures=c("Shannon"), color="Ethnicity") +
  geom_boxplot(color="black", aes(fill= sample_data(ps_tumoronly_bw)$Ethnicity)) +
  geom_point(color="black") +
  theme_bw() +
  scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") +
  #labs(title = expression("Shannon diversity by ethnicity, tumor samples only")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  #scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_y_continuous(name = "Shannon") +
  theme(legend.position='none', axis.title.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Black", "White")), label = "p.signif", symnum.args = symnum.args) +
  theme(strip.background = element_rect(color="black", size=3, linetype="blank"), strip.placement = "outside",
        strip.text.x = element_text(size = 10, color = "black", face = "bold"))
alpha5 <- plot_richness(ps_tumoronly_bw, x="Ethnicity", measures=c("Chao1"), color="Ethnicity") +
  geom_boxplot(color="black", aes(fill= sample_data(ps_tumoronly_bw)$Ethnicity)) +
  geom_point(color="black") +
  theme_bw() +
  scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") +
  #labs(title = expression("Chao1 diversity by ethnicity, tumor samples only")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  #scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_y_continuous(name = "Chao1") +
  theme(legend.position='none', axis.title.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Black", "White")), label = "p.signif", symnum.args = symnum.args) +
  theme(strip.background = element_rect(color="black", size=3, linetype="blank"), strip.placement = "outside",
        strip.text.x = element_text(size = 10, color = "black", face = "bold"))
alpha6 <- plot_richness(ps_tumoronly_bw, x="Ethnicity", measures=c("Simpson"), color="Ethnicity") +
  geom_boxplot(color="black", aes(fill= sample_data(ps_tumoronly_bw)$Ethnicity)) +
  geom_point(color="black") +
  theme_bw() +
  scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") +
  #labs(title = expression("Simpson diversity by ethnicity, tumor samples only")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  #scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_y_continuous(name = "Simpson") +
  theme(legend.position='none', axis.title.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Black", "White")), label = "p.signif", symnum.args = symnum.args) +
  theme(strip.background = element_rect(color="black", size=3, linetype="blank"), strip.placement = "outside",
        strip.text.x = element_text(size = 10, color = "black", face = "bold"))
#combine graphs
ggarrange(alpha4, alpha5, alpha6, 
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)

#######step 4: ordinate
# Transform data to proportions as appropriate for Bray-Curtis distances
ps_reloaded.prop <- transform_sample_counts(ps_reloaded, function(otu) otu/sum(otu))
ps_pairedonly.prop <- transform_sample_counts(ps_pairedonly, function(otu) otu/sum(otu))
ps_canceronly.prop <- transform_sample_counts(ps_canceronly, function(otu) otu/sum(otu))
ps_canceronly_bw.prop <- transform_sample_counts(ps_canceronly_bw, function(otu) otu/sum(otu))

ord.nmds.bray.reloaded <- ordinate(ps_reloaded.prop, method="NMDS", distance="bray") #oops.no convergance
ord.nmds.bray.pairedonly <- ordinate(ps_pairedonly.prop, method="NMDS", distance="bray") #oops.no convergance
ord.nmds.bray.canceronly <- ordinate(ps_canceronly.prop, method="NMDS", distance="bray") #oops.no convergance
ord.nmds.bray.canceronlybw <- ordinate(ps_canceronly_bw.prop, method="NMDS", distance="bray") #oops.no convergance

ord.PCoA.bray.reloaded <- ordinate(ps_reloaded.prop, method="PCoA", distance="bray") 
ord.PCoA.bray.pairedonly <- ordinate(ps_pairedonly.prop, method="PCoA", distance="bray") 
ord.PCoA.bray.canceronly <- ordinate(ps_canceronly.prop, method="PCoA", distance="bray")
ord.PCoA.bray.canceronlybw <- ordinate(ps_canceronly_bw.prop, method="PCoA", distance="bray")

plot_ordination(ps_reloaded.prop, ord.PCoA.bray.reloaded, color="Sample_description", title="Bray PCOA by Description, all samples")
plot_ordination(ps_pairedonly.prop, ord.PCoA.bray.pairedonly, color="Sample_description", title="Bray PCOA by Description, paired samples only")
plot_ordination(ps_canceronly.prop, ord.PCoA.bray.canceronly, color="Sample_description", title="Bray PCOA by Description, cancer samples only")
plot_ordination(ps_canceronly_bw.prop, ord.PCoA.bray.canceronlybw, color="Sample_description", title="Bray PCOA by Description, cancer samples only")


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

##########step 5: plot taxa
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


#########step 6: OKAY, save everything and... hope to find something new and exciting later.

#########step 7: now it is later! Let's go to DESeq2
library("DESeq2")
#let's make a couple of deseq2 objects with different formulae
diagdds = phyloseq_to_deseq2(kostic, ~ DIAGNOSIS)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")



