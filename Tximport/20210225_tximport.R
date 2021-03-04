#' ---
#' title: "Import transcript abundances from kallisto using tximport"
#' author: "Chelsea Herdman"
#' date: "February 25th, 2021"
#' output: github_document
#' ---
#'
#' In order to perform differential expression analysis using DESeq2 on the 
#' estimated transcript abundances produced by kallisto, we used the tximport
#' package in R.
#'
#' We followed the vignette provided by Michael I. Love, Charlotte Soneson and 
#' Mark D. Robinson found [here](http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#salmon_with_inferential_replicates)
#' 
#' #### Prepare transcript to gene dataframe for gene-level summarization
#' 
#' Load required libraries.
#' 
#+ libraries, message=FALSE, error=FALSE, warning=FALSE
library(tximport)
library(data.table)
library(ensembldb)
library(RMariaDB)
library(here)
library(rhdf5)

#' Fetch gene annotations from Ensembl and create one to one 
#' transcript id/gene id data frame.

txdb = makeTxDbFromEnsembl(organism="Danio rerio",
                           release=89,
                           server="ensembldb.ensembl.org",
                           username="anonymous", password=NULL, port=0L,
                           tx_attrib=NULL)

k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)

length(unique(tx2gene$TXNAME))

length(unique(tx2gene$GENEID))

#' #### Prepare named vector of files
#' Use kallisto abundance.tsv files. These are listed in alphabetical order in the 
#' kallisto_out directory so ensure naming associates accurately and that names 
#' associate with sample info table that will be used for DESeq2.

file_path_vec = list.files(path=here("kallisto"),
                   pattern="abundance.tsv",
                   recursive = TRUE,
                   full.names = TRUE)
all(file.exists(file_path_vec))

file_tab = data.table(file_path=file_path_vec,
                      sample_id=basename(dirname(file_path_vec)))

files = file_tab$file_path
names(files) = file_tab$sample_id

#'
#' #### Run tximport for gene level estimation
#'
#' ##### Bias corrected counts without an offset
txi.kallisto.biascorr = tximport(files, type= "kallisto", tx2gene = tx2gene, 
                                 ignoreTxVersion = TRUE, 
                                 countsFromAbundance="lengthScaledTPM")
txi.kallisto.biascorr$counts[1:6, 1:6]

counts_tab = as.data.table(txi.kallisto.biascorr$counts)
counts_tab$ensembl_gene_id = rownames(txi.kallisto.biascorr$counts)
setcolorder(counts_tab, c("embryonic_heart_0", "embryonic_heart_1", "embryonic_heart_2",
                          "adult_heart_0", "adult_heart_1", "adult_heart_2",
                          "adult_muscle_0", "adult_muscle_1", "adult_muscle_2"))

#' Save the counts table in order to use as input for DESeq2 (*see DESeq2 folder*)
fwrite(counts_tab, file=here("Tximport", "20210225_Wangetal_counts_fromtximport_biascorrected.txt"), sep="\t")

