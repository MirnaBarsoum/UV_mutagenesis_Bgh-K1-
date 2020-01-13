#subset for single copy genes (the list of sigle copy genes (Frantzeskakis et al., 2018)in DH14 was used

````
grep -Fwf /home/mb297167/One_copy_genes /home/mb297167/proteincoding_single_copy /work/mb297167/DH14/raw_variant_DH14_freebayes.annotated.singlecopy.vcf > /home/mb297167/proteincoding_single_copy /work/mb297167/DH14/raw_variant_DH14_freebayes.annotated.singlecopy_proteincoding.vcf
grep -Fwf /home/mb297167/proteincoding_single_copy /work/mb297167/DH14/raw_variant_DH14_freebayes.annotated.singlecopy.vcf > /home/mb297167/proteincoding_single_copy /work/mb297167/DH14/raw_variant_DH14_freebayes.annotated.singlecopy_proteincoding.vcf
````
the position (Pos), the depth at each position (DP) and the read supporting the alternative allele/mutation (AO) were pulled out from the vcf using GATK tool
java -jar /home/mb297167/tools/gatk-4.0.12.0/GenomeAnalysisTK.jar -R /home/mb297167/temp_mirna/Bgh_genome/bgh_dh14_v4.fa -T VariantsToTable -V /work/mb297167/DH14/raw_variant_DH14_freebayes.annotated.singlecopy.proteincoding.vcf -F POS -F AO -F DP /work/mb297167/DH14/raw_variant_DH14_freebayes.annotated.singlecopy.proteincoding.vcf.table

the data was then converted to excel file, the percentage of reads supporting each of the mutation compared to the total number of reads at each position were calculated, and a PivotTable was generated, from which the frequency histogram was made