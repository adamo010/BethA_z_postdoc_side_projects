#The goal of this analysis is to analyze differential gene expression in RNAseq data from mice.
#see 08.14.20 meeting with Ran and Gargi re: RNAseq data in Evernote for an explanation
#see 08.19.20 Getting into RNA-seq analysis: MSI overview in Evernote for the kallisto notes
#see 09.15.20 RNA-seq analysis: post-kallisto in Evernote for notes on this process.

#####################step 0: packages and dependencies##############
#what does DESeq2 need?

library("DESeq2")
#here are some other packages that I might need:
library("ggplot2")
library("readxl")
library("dplyr")
library("tidyverse")
library("ggpubr")
library("ape")
library("gplots")
library("plotly")
library("tidyr")
library("vegan")
library("VennDiagram")
library("rhdf5")
library("tximportData")

#also need rhdf5 package: install that too.
#out of time order, but see the 09.16.20_tximport_vignette.R

#"First, we locate the directory containing the files.

dir <- system.file("extdata", package = "tximportData")
list.files(dir)

