# UV_mutagenesis_Bgh-K1-

#trimming the reads

java -jar ~/trimmomatic-0.38.jar PE -threads 8 -trimlog ~/log-file ~/F.fastq.gz ~/R.fastq.gz ~/F_paired.fq ~/F_unpaired.fq ~/R_paired.fq ~/R_unpaired.fq ILLUMINACLIP:Trimmomatic-0.38/adapters/TruSeq3-PE.fa:5:30:10 SLIDINGWINDOW:3:18 LEADING:6 TRAILING:6 MINLEN:60

#map all K1, UV2 and UV8 reads to the DH14 genome (BWA)

~/bwa-0.7.17/bwa mem -t 12 -M -R '@RG\tID:sample_1\tLB:sample_1\tPL:ILLUMINA\tPM:HISEQ\tSM:sample_1' ~/Bgh_genome/bgh_dh14_v4.fa ~/F_paired.fq ~/B_paired.fq > ~/BWA-K1/UV2/UV8-DH14.sam

#sam to bam

~/samtools-1.9/samtools view -b -S ~/BWA-K1/UV2/UV8-DH14.sam > ~/BWA-K1/UV2/UV8-DH14.bam

#sort bam

java -jar ~/picardcloud.jar AddOrReplaceReadGroups I=~/BWA-K1/UV2/UV8-DH14.bam O=~/BWA-K1/UV2/UV8-DH14.bam_PicardSort.bam SO=coordinate RGID=sample_dh14 RGLB=sample_dh14-Bgh RGPL=ILLUMINA RGSM=dh14 RGPU=Mi

#mark duplicate
java -jar /home/mb297167/tools/Picard/picard/build/libs/picardcloud.jar MarkDuplicates INPUT=~/BWA-K1/UV2/UV8-DH14.bam_PicardSort.bam OUTPUT=~/BWA-K1/UV2/UV8-DH14.bam_PicardSort_dedup.bam METRICS_File=metric.txt

#index

java -jar ~/picardcloud.jar BuildBamIndex I=~/BWA-K1/UV2/UV8-DH14_PicardSort_dedup.bam

#call raw variants (GATK)

java -jar ~/GenomeAnalysisTK.jar -T HaplotypeCaller -ploidy 1 -R ~/Bgh_genome/bgh_dh14_v4.fa -I ~/BWA-K1/UV2/UV8-DH14_PicardSort_dedup.bam -o ~/BWA-K1/UV2/UV8-DH14_PicardSort_dedup_raw_variants_GATK.vcf

#call raw SNPs (GATK)

java -jar ~/GenomeAnalysisTK.jar -T SelectVariants -R ~/Bgh_genome/bgh_dh14_v4.fa -V ~/BWA-K1/UV2/UV8 DH14_PicardSort_dedup_raw_variants_GATK.vcf -selectType SNP -o ~/BWA-K1/UV2/UV8 DH14_PicardSort_dedup_GATK_raw_snps.vcf

#call raw INDELs (GATK)

java -jar ~/GenomeAnalysisTK.jar -T SelectVariants -R ~/Bgh_genome/bgh_dh14_v4.fa -V ~/BWA-K1/UV2/UV8 DH14_PicardSort_dedup_raw_variants_GATK.vcf -selectType INDEL -o ~/BWA-K1/UV2/UV8 DH14_PicardSort_dedup_GATK_raw_INDELs.vcf

#filter SNPs (GATK)

java -jar GenomeAnalysisTK.jar -T VariantFiltration -R ~/Bgh_genome/bgh_dh14_v4.fa -V ~/BWA-K1/UV2/UV8 DH14_PicardSort_dedup_GATK_raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o ~/BWA-K1/UV2/UV8 DH14_PicardSort_dedup_GATK_filter_snps.vcf

#filter INDELs (GATK)



#call raw variants (Freebayes)

#filter variants (Freebayes)

#call raw variants (Mpileup)

#filter variants (Mpileup)

#Unique and common variants





