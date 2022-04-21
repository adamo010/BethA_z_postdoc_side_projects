#lets's work through a sample first, if I can
#see https://bioconductor.riken.jp/packages/3.4/bioc/vignettes/tximport/inst/doc/tximport.html

#step 1: locate the directory containing the files. For the tutorial we use system.file, but for
#most uses we'd use a file path e.g "/my/files/here"
dir <- system.file("extdata", package = "tximportData")

#what files are present in the directory?
list.files(dir)

#step 2: create a vector pointing to the kallisto files (whatever the fuck that means)
#read in a table that contains sample_ids (presumably this is samples.txt) and combining this with dir
#and abundance.tsv (now where the fuck is this?)

samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
samples
#soooo. it turns out samples has info like pop/center/assay/sample/experiment/run. Presumably this info
#I have to create myself?

#all right, this is the next step in the tutorial but I don't know what's actually going on here.
files <- file.path(dir, "kallisto", samples$run, "abundance.tsv")
#according to StackOverflow, the $ is for extracting specific information. 
#E.g if x <- list(a=1, b=2, c=3); printing x$b returns a value 2.
#the files object contains the following from file.path: directory, "kallisto" object, the column "run"
#from the samples dataframe, and the tsv file "abundance.

names(files) <- paste0("sample", 1:6)
all(file.exists(files))
#some other dataframe shit.
#well, that last command was supposed to return TRUE (or at least it did in the tutorial), but it returned FALSE
#for me. We'll see if it works.
#NOPE. Falls apart when we try to use tximport

#Side note: I think I should just run through the tutorial and build some of this shit by hand/in python.

#first make a data.frame called tx2gene with two columns: 1) transcript ID and 2) gene ID. 
#The column names do not matter but this column order must be used. The transcript ID must be 
#the same one used in the abundance files.
#Creating this tx2gene data.frame can be accomplished from a TxDb object and the 
#select function from the AnnotationDbi package (which I have???)
#nope. Had to install that too.

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]  # tx ID, then gene ID
#df and tx2gene are the same data, just with the column order reversed

#okay, so it seems like what this is doing is pusing the TxDb.Hsapiens.UCSC.hg19.knownGene library/database
#to generate a lookup table of gene ids and their txnames (transcript names?)- maybe there's a way to export and 
#save this table, since it seems like it would be useful for any transcriptome analysis.

#then the tutorial has some other info about using an Ensembl transcriptome- not sure if/what that is.
#skip for now.

library(tximport)
library(readr)

txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
#Did not work; need all(file.exists(files)) to return TRUE. What is going on?

############################################################
#okay, let's run through this stupid tutorial AGAIN 

#"First, we locate the directory containing the files. (Here we use system.file to locate 
#the package directory, but for a typical use, we would just provide a path, e.g. "/path/to/dir".)"

library(tximportData)
dir <- system.file("extdata", package = "tximportData")
list.files(dir)

#"Next, we create a named vector pointing to the quantification files. We will create a vector 
#of filenames first by reading in a table that contains the sample IDs, and then combining this
#with dir and "quant.sf.gz". (We gzipped the quantification files to make the data package smaller, 
#this is not a problem for R functions that we use to import the files.)"

samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
samples

files <- file.path(dir, "kallisto", samples$run, "quant.sf.gz")
#okay this is so odd; the file.path directory is a string combining all of those components? dir + kallisto + samples$run + quant.sf.gz
#I don't get it. Maybe this is why things are showing up as false? There's no files called quant.sf.gz in those directories?
names(files) <- paste0("sample", 1:6)
all(file.exists(files))
#this part is still not working. This all(file.exists(files)) command is supposed to return TRUE
#but keeps returning FALSE

#try a modification of this command so it matches what's really in the directories:
files2 <- file.path(dir, "kallisto", samples$run, "abundance.h5")
files3 <- file.path(dir, "kallisto", samples$run, "abundance.tsv.gz") #try this because gene names aren't matching up; did not help

names(files) <- paste0("sample", 1:6)
names(files2) <- paste0("sample", 1:6)
all(file.exists(files2))
#AHA! files2 works (i.e. the all command returns True). The CODE IN THE VIGNETTE IS WRONG.

#moving on...
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]  # tx ID, then gene ID
#so far so good.

tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))
head(tx2gene)

library(tximport)
library(readr)

txi <- tximport(files2, type = "kallisto", tx2gene = tx2gene)

#okay, now this is a helpful error.
#Error in .local(object, ...) : 
#"None of the transcripts in the quantification files are present
#in the first column of tx2gene. Check to see that you are using
#the same annotation for both."

txi <- tximport(files2, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar= TRUE)
#okay, neither of the possible fixes worked. I wonder if trying the ensembl approach will help...

#oh yep, I googled one of the example ids in the file and it is deffo an ensembl id.
#who wrote this terrible tutorial?

#all right, switching gears to the ensembl approach (vs the library(TxDb.Hsapiens.UCSC.hg19.knownGene) approach)

#as far as I can tell, the steps are as follows:
#1: call the library
#2: assign the library to a shortcut name (e.g. txdb2)
#3: create a referenceable dataframe??? (I don't know enough R for this)
#4: swap the columns

#step 1: call the library (after installing it)
#library(ensembldb)
library(EnsDb.Hsapiens.v86) #use this one, have to do less filtering

#step 2: assign the library to a shortcut name
edb <- EnsDb.Hsapiens.v86

#the next step in the tximport vignette is "the transcripts function can be used 
#with return.type="DataFrame", in order to obtain something like the df object 
#constructed in the code chunk above. See the ensembldb package vignette for 
#more details."- because, you know, that makes sense.

#this is from the Ensembl vignette (I DO NOT TRUST vignettes anymore, but...)
Tx <- transcripts(edb, return.type="DataFrame")
#I'm not sure what to filter by here. Let's try nothing.
Tx #print the DF so we know which column to pull out.
#so, it seems like the return.type function combines steps 2 and 3 (maybe??)

#step 4: get tx id, then gene ID into a single column
ebd_tx2gene <- Tx[c("tx_id", "gene_id")]
ebd_tx2gene

#I miss python.

#Okay, ebd_tx2gene should now function like tx2gene

library(tximport)
library(readr)

txi_ebd <- tximport(files2, type = "kallisto", tx2gene = ebd_tx2gene, ignoreTxVersion = TRUE)
#Okay,this is SO unhelpful; the actual variable name is tx2gene, and in the tutorial you create a database called tx2gene
#that you then assign to the function tx2gene. What the fuck.

#if I use the ignoreTxVersion = True option here, I get "transcripts missing from tx2gene: 2768" but it still runs.

#I think I have the wrong ensembl version. Which is... fixable, but a pain in the ass.
#how do I figure out which ensembl version to use? Or is that really the problem?

#characterize the identity of the missing transcripts: see https://support.bioconductor.org/p/123134/
#basically, I want all the transcript IDs that are in ebd_tx2gene but not txi_ebd
head(ebd_tx2gene)
head(txi_ebd)

#something is fishy here. ebd_tx2gene should be a long list of TxName and geneids.  Currently the gene_ids are... incorrect. 
#this is what happens when we use R. no one explains what's going on so I can't see when steps aren't working.

###################################RESTART THE TUTORIAL#####################################
library(tximportData)
dir <- system.file("extdata", package = "tximportData")
list.files(dir)

samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
samples

files <- file.path(dir, "kallisto", samples$run, "abundance.tsv")
names(files) <- paste0("sample", 1:6)
all(file.exists(files))

files2 <- file.path(dir, "kallisto", samples$run, "abundance.h5")
names(files2) <- paste0("sample", 1:6)
all(file.exists(files2))

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]  # tx ID, then gene ID

tx2gene #has two column headings: Txname and GeneID

#now, I want to make something similar with ensembldb

#step 1: call the library (after installing it)
#library(ensembldb)
library(EnsDb.Hsapiens.v86) #use this one, have to do less filtering

#step 2: assign the library to a shortcut name
edb <- EnsDb.Hsapiens.v86

#this is from the Ensembl vignette (I DO NOT TRUST vignettes anymore, but...)
Tx <- transcripts(edb, return.type="DataFrame")
#I'm not sure what to filter by here. Let's try nothing.
Tx #print the DF so we know which column to pull out.

#something about this is weird. There's only six rows in Tx, and there should be thousands (I assume).

#found another tutorial online about this (https://support.bioconductor.org/p/109092/). Install AnnotationHub and try to sort this out.
library(AnnotationHub)
ah <- AnnotationHub()
query(ah, "EnsDb.Hsapiens") #okay, this tells me what all the databases for H. Sapiens are. Let's try 87.
edb <- ah[["AH53211"]] #load Ensembl 87 EnsDb for Homo Sapiens

## Now get the transcripts table as a data frame
txs <- transcripts(edb, return.type = "DataFrame")
txs

#okay, NOW we're cooking.
#get tx id, then gene ID into a single column
ebd_tx2gene <- txs[c("tx_id", "gene_id")]
ebd_tx2gene

#YES, now we're onto 

#Okay, ebd_tx2gene should now function like tx2gene

library(tximport)
library(readr)

txi_ebd <- tximport(files2, type = "kallisto", tx2gene = ebd_tx2gene, ignoreTxVersion = TRUE)
#still get the Error in .local(object, ...) without the ignoreTxVersion=TRUE option
#weird. Still get 2768 transcripts missing.

#what are those transcripts? Found in files2 but not in ebd_tx2gene. How do I figure that out?
#well, the hard thing is that files2 is a vector of file names for transcript-level abundances. It's not like comparing to 2 data frames.

#I miss python.

#the tutorial at (https://support.bioconductor.org/p/123134/) has instructions for reading in a tsv, but I'm using h5 files so I have to do extra work. 
#library(rhdf5) I started with this and it didn't work, so I'm trying a different package
install.packages("hdf5r")
#can't install this without hdf5
install.packages("hdf5")
#Oh for fuck's sake: Warning in install.packages : package 'hdf5' is not available (for R version 4.0.2)

