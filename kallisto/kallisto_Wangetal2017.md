
# Estimate Transcript Abundances using Kallisto

### Prepare kallisto index

Download fasta files for coding and non-coding annotated sequences from [ensembl](http://www.ensembl.org/info/data/ftp/index.html) and create the kallisto index

```
kallisto index \
    -i GRCz10.r89.cdna.all.ncrna.kallisto_index \
    Danio_rerio.GRCz10.cdna.all.fa.gz
    Danio_rerio.GRCz10.ncrna.fa.gz
```
### Download fastq files from this dataset

Downloaded fastq files from [Gene Expression Omnibus (Accession #GSE85416)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85416), (Wang et al. 2017, Shih et al. 2015).

### Run kallisto

Run kallisto version 0.46.1 using GRCz10 release 89 and the following arguments.

```
/usr/local/kallisto/kallisto quant \
-i /data3/cherdman/genome/GRCz10.r89.cdna.all.ncrna.kallisto_index \
-o embryonic_heart_0 \
-b 100 \
-t 14 \
../fastq/SRR4017367.1_1.fastq.gz ../fastq/SRR4017367.1_2.fastq.gz >> \ ./embryonic_heart_0/kallisto_quant_embryonic_heart_0_log.txt 2>&1
```
