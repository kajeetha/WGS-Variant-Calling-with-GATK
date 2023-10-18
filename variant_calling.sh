#!/bin/bash

#script to call germline variants in a human WGS paired end reads following the GATK best practices workflow. 

if false
then
#download the data

wget -P ~/VC/reads/ ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
wget -P ~/VC/reads/ ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz

echo "Run Prep Files..."

#download the reference file
wget -P ~/VC/ref/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip ~/VC/ref/hg38.fa.gz

index the fasta file
samtools faidx ~/VC/reads/hg38.fa

#create a dictionary important for re-alignment and re-calibration 
gatk CreateSequenceDictionary R= ~/VC/reads/hg38.fa O= /Users/kajeethasarvananthan/Desktop/VC/reads/hg38.dict

#download the known variant sites for the base quality score recalibration step and the index file of it
wget -P ~/VC/reads/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P ~/VC/reads/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

fi
#directories 
ref="~/VC/reads/hg38.fa"
known_sites="~/VC/reads/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="~/VC/aligned_reads"
reads="~/VC/reads"
results="~/VC/results"
data="~/VC/data"


#Step 1: Run fastqc

echo "STEP 1: QC- Run fastqc"

fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/

#no trimming of adaptors or for quality needed. 

echo "STEP 2: Map to reference using BWA-MEM"
#BWA index for reference
bwa index ${ref}


#BWA alignment and provide a read group since this will be missing from the sam file(this field is tab delimited)
#in Illumina data the read group IDs are composed of the flowcell name, lane number(unique identifiers)
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam

#if there is not enough RAM you can pipe the results directly into a BAM file but we want to look at the sam file 

#bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz | samtools sort -@ 4 -o ${aligned_reads}/SRR062634.sorted.bam

#flagstat to look at bitwise data of pairs of reads, how many mapped, how many duplicates and how many supplments

echo "Mark duplicate reads and flag them(ignore them)"

gatk MarkDuplicateSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam

#correct base quality scores through recalibration (based on ML model with known variants (stored in file called known sites)). A table will be output in the data folder.

gatk Baserecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam - R ${ref} --known-sites ${known_sites} -O {data}/recal_data.table

echo "adjust the base quality score"

gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} -bqsr-recal-file {data}/recal_data.table - O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam

echo "Collect Alignment & Insert Size Metrics"
gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/alignment_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf

#can go to the aligned reads values and combine alignment and insert metrics in a file using multiqc(tells us the number of reads that passes the threshold and a way to validate library construction)

#variants with analysis ready reads. 

echo "Collect Alignment & Insert Size Metrics" 
#an index vcf file is also created
gatk Haplotypecaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf

#extract SNVs and indels

gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf
