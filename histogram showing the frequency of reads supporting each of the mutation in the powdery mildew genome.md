#subset for single copy genes (the list of sigle copy genes (Frantzeskakis et al., 2018)in DH14 was used

````
grep -Fwf ~/One_copy_genes ~/raw_variant_DH14_freebayes.annotated.singlecopy.vcf > ~/raw_variant_DH14_freebayes.annotated.one_copy_genes.vcf
grep -Fwf ~/proteincoding_single_copy ~/raw_variant_DH14_freebayes.annotated.one_copy_genes.vcf > ~/raw_variant_DH14_freebayes.annotated.one_copy_genes_proteincoding.vcf
````
#the position (Pos), the depth at each position (DP) and the read supporting the alternative allele/mutation (AO) were pulled out from the vcf using GATK tool
````
java -jar ~/GenomeAnalysisTK.jar -R ~/Bgh_genome/bgh_dh14_v4.fa -T VariantsToTable -V ~/raw_variant_DH14_freebayes.annotated.one_copy_genes_proteincoding.vcf -F POS -F AO -F DP ~/raw_variant_DH14_freebayes.annotated.one_copy_genes_proteincoding.vcf.table
````

#the data was then converted to excel file, the percentage of reads supporting each of the mutation compared to the total number of reads at each position were calculated, and a PivotTable was generated, from which the frequency histogram was made
