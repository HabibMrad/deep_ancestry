for chr in {1..2}
do
	x-terminal-emulator -e "python3 calc_In_VCF_1000G.py -chr "$chr" -thres 0.2 -exclude ASW ACB MXL PUR CLM" &
done
