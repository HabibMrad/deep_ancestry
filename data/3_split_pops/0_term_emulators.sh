# # create lists of inds belonging to each population
# for popname in `awk '{print $1}' ../0_1000G_info/pop_names_unrelated |sort | uniq`; do grep $popname ../0_1000G_info/pop_names_unrelated | awk '{print $2}' > $popname.txt ; done

# ls -v *txt > 0_pops

# extract from the vcf file one population a time
# output in the IMPUTE format
# only non-monomorphic (among all 1000G populations), phased (--impute), biallelic (--impute) snps (--remove-indels) with no-missingness (--impute)
# only the AIMs (SNPs whose MAF are higher within a population and lower in the overall dataset)
for pop in `sed 's/\.txt//g' 0_pops`
do
	x-terminal-emulator -e "bash 1_extract_AIMs.sh "$pop"" &
done
