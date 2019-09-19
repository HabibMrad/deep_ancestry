for ((chr=22; chr>=1; chr--))
do
	# remove non-phased, non-snps, non-biallelic, monomorphic, no PASS, and missing (no flag for this, but no missing values in 1000G assumed)
	x-terminal-emulator -e "bcftools view -O z -o 1000G.chr"$chr".vcf.gz -v snps -c 1:minor -p -m 2 -M 2 -f PASS --threads 4 ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ; rm ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" &
	# bcftools view -O z -o 1000G.chr"$chr".vcf.gz -v snps -c 1:minor -p -m 2 -M 2 -f PASS --threads 4 ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
done

# slower, does the same (except that it removes the info column and outputs a '.' instead)
# vcftools --gzvcf ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --remove-indels --mac 1 --min-alleles 2 --max-alleles 2 --max-missing 1 --phased --remove-filtered-all --recode --stdout | gzip -c > 1000G."$chr".vcf.gz
