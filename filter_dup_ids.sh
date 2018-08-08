#script to filter bfiles obtained from downloaded VCFs
plink=/home/ubuntu/gwas/old_gwas/tools/plink/plink
path=/home/ubuntu/gwas/old_gwas/sega/ld_real/by_chr
for i in {1..22}; do 
echo $i;
cut -f 2 ${path}/ALL.chr${i}.bim | sort | uniq -d > ${i}.dups;
grep -v 'rs' ${path}/ALL.chr${i}.bim | cut -f 2 >> ${i}.dups;
$plink --bfile ${path}/ALL.chr${i}  --exclude ${i}.dups --make-bed --out ${path}/ALL.filt.chr${i};
done
