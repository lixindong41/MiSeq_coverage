
### Install third-party tools for the pipeline:
sudo apt-get install bowtie2  
sudo apt-get install samtools  
sudo apt-get install bedtools

### Install R packages for the pipeline:
install.packages("gridBase")  

### Reference genome and Bowtie indexing libraries:
The reference genome is downloaded from: http://www.ncbi.nlm.nih.gov/nuccore/K03455.1  
Create bowtie2 index:  
bowtie2-build ./ref/HXB2.fasta ./ref/HXB2  

### Data set:
MiSeq data set is downloaded from: http://www.ncbi.nlm.nih.gov/sra/?term=SRR961514  
Split paired end SRA file into fastq files:  
fastq-dump --split-files SRR961514.sra
