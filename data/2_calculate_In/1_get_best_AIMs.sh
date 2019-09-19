# best AIMs
cat *in | sort -k2,2r | head -n 100

# AIMs per chromosome
wc -l *in

# print file for next step
cat *in | sort -k2,2r | head -n 2000 | awk '{print$1}' > 1000G_2000AIMs.txt
cat *in | sort -k2,2r | head -n 4000 | awk '{print$1}' > 1000G_4000AIMs.txt
