pop=$1

for chr in {1..22}
do
	vcftools --gzvcf ../1_download_vcf/1000G.chr"$chr".vcf.gz --keep "$pop".txt --snps ../2_calculate_In/1000G_2000AIMs.txt --IMPUTE --out "$pop"_chr"$chr"_2000AIMs
	rm "$pop"_chr"$chr"_2000AIMs.impute.hap.indv "$pop"_chr"$chr"_2000AIMs.log
done
