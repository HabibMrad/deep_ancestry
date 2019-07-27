for pop in `sed 's/\.txt//g' 0_pops`
do
	cat `ls -v "$pop"*hap` > "$pop"_2000AIMs.impute.hap
done
