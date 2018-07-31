plink=/home/ubuntu/gwas/old_gwas/tools/plink/plink

for i in {1..22}; do 
echo $i;
$plink --vcf phase3/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --make-bed --out by_chr/ALL.chr${i};
done
