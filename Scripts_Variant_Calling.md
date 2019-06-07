UV_mutagenesis_Bgh-K1

#trimming the reads (Trimmomatic)

```
java -jar ~/trimmomatic-0.38.jar PE -threads 8 -trimlog ~/log-file ~/F.fastq.gz ~/R.fastq.gz ~/F_paired.fq ~/F_unpaired.fq ~/R_paired.fq ~/R_unpaired.fq ILLUMINACLIP:Trimmomatic-0.38/adapters/TruSeq3-PE.fa:5:30:10 SLIDINGWINDOW:3:18 LEADING:6 TRAILING:6 MINLEN:60
```
#map all K1, UV2 and UV8 reads to the DH14 genome (BWA)
```
~/bwa-0.7.17/bwa mem -t 12 -M -R '@RG\tID:sample_1\tLB:sample_1\tPL:ILLUMINA\tPM:HISEQ\tSM:sample_1' ~/Bgh_genome/bgh_dh14_v4.fa ~/F_paired.fq ~/B_paired.fq > ~/BWA-K1/UV2/UV8-DH14.sam
```
#sam to bam
```
~/samtools-1.9/samtools view -b -S ~/BWA-K1/UV2/UV8-DH14.sam > ~/BWA-K1/UV2/UV8-DH14.bam
```
#sort bam
```
java -jar /picardcloud.jar AddOrReplaceReadGroups I=/BWA-K1/UV2/UV8-DH14.bam O=~/BWA-K1/UV2/UV8-DH14.bam_PicardSort.bam SO=coordinate RGID=sample_dh14 RGLB=sample_dh14-Bgh RGPL=ILLUMINA RGSM=dh14 RGPU=Mi
```
#mark duplicate
```
java -jar /home/mb297167/tools/Picard/picard/build/libs/picardcloud.jar MarkDuplicates INPUT=/BWA-K1/UV2/UV8-DH14.bam_PicardSort.bam OUTPUT=/BWA-K1/UV2/UV8-DH14.bam_PicardSort_dedup.bam METRICS_File=metric.txt
```
#index
```
java -jar /picardcloud.jar BuildBamIndex I=/BWA-K1/UV2/UV8-DH14_PicardSort_dedup.bam
```
#Create Realignment Targets (GATK) #This is the first step in a two-step process of realigning around indels

```
java -jar ~/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ~/Bgh_genome/bgh_dh14_v4.fa -I ~/K1DH14mapped-n-sort.markdup-RG.bam -o ~/K1DH14mapped-n-sort-markdup-RG-realignment_targets.list
```
#Realign Indels (GATK) #This step performs the realignment around the indels which were identified in the previous step (the ‘realignment targets’)
```
java -jar ~/GenomeAnalysisTK.jar -T IndelRealigner -R ref -I ~/dedup_reads.bam -targetIntervals ~/realignment_targets.list -o ~/realigned_reads.bam
```
#Base Quality Score Recalibration BQSR (GATK)
```
java -jar ~/GenomeAnalysisTK.jar -T BaseRecalibrator -R ~/Bgh_genome/bgh_dh14_v4.fa -I ~/realigned_reads.bam -knownSites ~/filtered_snps.vcf -knownSites ~/filtered_indels.vcf -o ~/recal_data.table
```
#Apply BQSR (GATK)
```
java -jar ~/GenomeAnalysisTK.jar -T PrintReads -R ~/Bgh_genome/bgh_dh14_v4.fa -I ~/realigned_reads.bam -BQSR recal_data.table -o recal_reads.bam
```
#call raw variants from dedup.bam, realigned.bam or recal_data.bam (GATK)
```
java -jar ~/GenomeAnalysisTK.jar -T HaplotypeCaller -ploidy 1 -R ~/Bgh_genome/bgh_dh14_v4.fa -I ~/BWA-K1/UV2/UV8-DH14_PicardSort_dedup.bam -o ~/BWA-K1/UV2/UV8-DH14_PicardSort_dedup_raw_variants_GATK.vcf
```
#call raw SNPs (GATK)
```
java -jar ~/GenomeAnalysisTK.jar -T SelectVariants -R ~/Bgh_genome/bgh_dh14_v4.fa -V ~/BWA-K1/UV2/UV8 DH14_PicardSort_dedup_raw_variants_GATK.vcf -selectType SNP -o ~/BWA-K1/UV2/UV8 DH14_PicardSort_dedup_GATK_raw_snps.vcf
```
#call raw INDELs (GATK)
```
java -jar ~/GenomeAnalysisTK.jar -T SelectVariants -R ~/Bgh_genome/bgh_dh14_v4.fa -V ~/BWA-K1/UV2/UV8 DH14_PicardSort_dedup_raw_variants_GATK.vcf -selectType INDEL -o ~/BWA-K1/UV2/UV8 DH14_PicardSort_dedup_GATK_raw_INDELs.vcf
```
#filter SNPs (GATK)
```
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R ~/Bgh_genome/bgh_dh14_v4.fa -V ~/BWA-K1/UV2/UV8 DH14_PicardSort_dedup_GATK_raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o ~/BWA-K1/UV2/UV8 DH14_PicardSort_dedup_GATK_filter_snps.vcf
```
#filter INDELs (GATK)
```
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R ~/Bgh_genome/bgh_dh14_v4.fa -V ~/BWA-UV2-DH14_PicardSort.dedup_GATKNORMAL_raw_INDEL.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o ~/BWA-UV2-DH14_PicardSort.dedup_GATKNORMAL_filtered_INDEL.vcf
```
#call raw variants (Freebayes)
```
~/freebayes -f ~/Bgh_genome/bgh_dh14_v4.fa -p 1 ~/BWA-UV2-DH14_PicardSort.dedup.bam -u > ~/BWA UV2_DH14_PicardSort.dedup._raw_variant_freebayes.vcf
```
#filter variants (Freebayes)
```
~/vcffilter -f "QUAL > 20 & QUAL / AO > 10 & SAF > 1 & SAR > 1 & RPR > 1 & RPL > 1" ~/BWA UV2_DH14_PicardSort.dedup._raw_variant_freebayes.vcf > ~/BWA UV2_DH14_PicardSort.dedup_vcffiltered.freebayes.vcf
```
#call raw variants (Mpileup)
```
~/samtools mpileup -g -f ~/Bgh_genome/bgh_dh14_v4.fa ~/BWA-UV2-DH14_PicardSort.dedup.bam > ~/mpileup_UV2_raw.bcf ~/bcftools call -c -v --output-type b --ploidy 1 ~/mpileup_UV2_raw.bcf > ~/mpileup_UV2_var.bcf ~/bcftools view ~/mpileup_UV2_var.bcf > ~/mpileup_final_UV2.vcf
```
#filter variants (Mpileup)
```
java -jar ~/SnpSift.jar filter "( QUAL >= 20 && DP > 3 && MQ > 50 )" /mpileup_final_UV2.vcf > //mpileup_final_UV2.vcf_SNPsift.vcf
```
#remove INDELs (Freebayes and mpileup)
```
~/vcftools --vcf ~/BWA UV2_DH14_PicardSort.dedup_vcffiltered.freebayes.vcf --remove-indels --recode --recode-INFO-all --out ~/BWA UV2_DH14_PicardSort.dedup_vcffiltered.freebayes.snpsonly.vcf
```
#Unique and common variants (bcftools) #I could also reduce the false ''unique'' positive by switching filtering and finding the unique variants steps, especially in freebayes and mpileup where the filtered variants are deleted from the vcf file; if i filter the variants before finding the unique variants, i will have more false positive (filtered in K1 and not in UV).
```
~/bcftools isec -p ~/isec_common_Unique_K1UV2 ~/K1.vcf.gz ~/UV.vcf.gz
```
#Predict the effect
```
java -jar ~/snpEff.jar eff DH14 /work/mb297167/UV_BWA/BWA-UV2-DH14_PicardSort.dedup.vcf > /work/mb297167/UV_BWA/BWA-UV2-DH14_PicardSort.dedup.annotated.vcf
```
