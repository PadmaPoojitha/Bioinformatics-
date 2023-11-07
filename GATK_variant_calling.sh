
#! /bin/sh.

#Variant Analysis: 

#Data: 
#1000 genomes project - Phase 3 HG00096
#Paired end whole genome sequencing reads

#Download Data: done 
wget -P ~/Desktop/Bioinformatics/WGS/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz

wget -P ~/Desktop/Bioinformatics/WGS/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz

echo "Run Prep files..."

#Download Reference File: done
wget -P ~/Desktop/Bioinformatics/WGS https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

gunzip ~Desktop/Bioinformatics/WGS/hg38.fa.gz

# index ref - .fai file: done
samtools faidx ~/Desktop/Bioinformatics/WGS/hg38.fa

# ref dict - .dict file:
java -jar picard.jar CreateSequenceDictionary R=/Users/poojithaalla/Desktop/Bioinformatics/WGS/hg38.fa O=/Users/poojithaalla/Desktop/Bioinformatics/WGS/hg38.dict

# download known sites files for BQSR from GATK resource bundle: done 

wget -P ~/Desktop/Bioinformatics/WGS/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf

wget -P ~/Desktop/Bioinformatics/WGS/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

#Setting up directories: done
ref="/Users/poojithaalla/Desktop/Bioinformatics/WGS/hg38.fa"
known_sites="/Users/poojithaalla/Desktop/Bioinformatics/WGS/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/Users/poojithaalla/Desktop/Bioinformatics/WGS/aligned_reads"
reads="/Users/poojithaalla/Desktop/Bioinformatics/WGS/reads"
results="/Users/poojithaalla/Desktop/Bioinformatics/WGS/results"
data="/Users/poojithaalla/Desktop/Bioinformatics/WGS/data"

#FASTQC: done 
fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}

#BWA-MEM: done
bwa index ${ref}

bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam

samtools flagstat ${aligned_reads}/SRR062634.paired.sam 

#99.66% - no duplicates 

#Mark & flag duplicate reads: done  

gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dd_reads.bam

samtools flagstat ${aligned_reads}/SRR062634_sorted_dd_reads.bam

#Base quality recalibration: Done 

gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dd_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table

gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dd_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_bqsr_reads.bam

#Collect alignement & insert metric: running 
gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt

gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf

#QC for insert and alignment summary metrics: (go into the directory)
multiqc . 

#Call Variants: 
gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_bqsr_reads.bam -O ${results}/raw_variants.vcf


#SNPs and indels: 
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/snps.vcf

gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/indels.vcf

#Filter SNPs: 
gatk VariantFiltration -R ${ref} -V ${results}/snps.vcf -O ${results}/filtered_snps.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 60.0" -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 4.0" -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" -genotype-filter-expression "DP < 10" -genotype-filter-name "DP_filter" -genotype-filter-expression "GQ < 10" -genotype-filter-name "GQ_filter"

#Filter Indels:
gatk VariantFiltration -R ${ref} -V ${results}/indels.vcf -O ${results}/filtered_indels.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 200.0" -filter-name "SOR_filter" -filter "SOR > 10.0" -genotype-filter-expression "DP < 10" -genotype-filter-name "DP_filter" -genotype-filter-expression "GQ < 10" -genotype-filter-name "GQ_filter" 

#Retaining Passed Variants: 
gatk SelectVariants --exclude-filtered -V ${results}/filtered_snps.vcf -O ${results}/final_snps.vcf

gatk SelectVariants --exclude-filtered -V ${results}/filtered_indels.vcf -O ${results}/final_indels.vcf

#Removing Genotype filter failed variants: (make sure to be in the directory)
cat final_snps.vcf|grep -v -E "DP_filter|GQ_filter" > final_snps_filteredGT.vcf

cat final_indels.vcf| grep -v -E "DP_filter|GQ_filter" > final_indels_filteredGT.vcf

#Annotation: 

gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download

gatk Funcotator --variant ${results}/final_snps_filteredGT.vcf --reference ${ref} --ref-version hg38 --data-sources-path /Users/poojithaalla/funcotator_dataSources.v1.7.20200521g --output ${results}/funcanto_snps.vcf --output-file-format VCF

gatk Funcotator --variant ${results}/final_indels_filteredGT.vcf --reference ${ref} --ref-version hg38 --data-sources-path /Users/poojithaalla/funcotator_dataSources.v1.7.20200521g --output ${results}/funcanto_indels.vcf --output-file-format VCF


#Converting to tab-delimited file: 
gatk VariantsToTable -V ${results}/funcanto_snps.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION -O ${results}/output_snps.table

gatk VariantsToTable -V ${results}/funcanto_indels.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION -O ${results}/output_indels.table

cat funcanto_snps.vcf | grep "Funcotation fields are: " | sed 's/|/\t/g' > snps_variants.txt

cat output_snps.table | cut -f 5 | grep "NBPF1" | sed 's/|/\t/g' >> snps_variants.txt

cat funcanto_indels.vcf | grep "Funcotation fields are: " | sed 's/|/\t/g' > indels_variants.txt

cat output_indels.table | cut -f 5 | grep "NBPF1" | sed 's/|/\t/g' >> indels_variants.txt


