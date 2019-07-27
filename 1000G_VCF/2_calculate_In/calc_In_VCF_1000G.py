from utils import informativeness_assignment
from gzip import open as gzopen
import argparse


parser = argparse.ArgumentParser(description='Find AIMs in 1000G VCF files')
parser.add_argument("-chr", type=int, default=21, help='Chromosome number')
parser.add_argument("-thres", type=float, default=0.5, help='Informativeness (In) threshold value for accepting a SNP')
parser.add_argument("-inds", type=list, default=False, help='Python list of the population to which each individual belongs, in the same order as in the VCF files. E.g. ["GBR", "GBR", "GBR", ...]')
parser.add_argument("-exclude", type=str, nargs='*', default=False, help='Python list of the populations to exclude e.g. YRI CEU GBR')
parser.add_argument("-skip", type=int, default=0, help='Skip the first m SNPs')
parser.add_argument("-path", type=str, default="../1_download_vcf/", help='Path to VCF files')
args = parser.parse_args()

chr = args.chr
threshold = args.thres
vcf_name = f"{args.path}1000G.chr{chr}.vcf.gz"
output = f"1000G.chr{chr}.in"

if args.inds is False:
	inds = ["GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR",
		"GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR",
		"GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR",
		"GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN",
		"FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR",
		"GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR",
		"GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "GBR", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN",
		"FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN",
		"FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN",
		"FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN",
		"FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN",
		"FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS",
		"CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS",
		"CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS",
		"CHS", "CHS", "CHS", "PUR", "PUR", "PUR", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS",
		"CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS",
		"CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "PUR", "PUR", "PUR", "PUR", "CHS", "CHS", "CHS",
		"CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS",
		"CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "CHS", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR",
		"PUR", "PUR", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX",
		"CDX", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR",
		"PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR",
		"PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM",
		"CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "PUR", "PUR", "PUR", "PUR",
		"PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR",
		"PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM",
		"CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR",
		"GBR", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM",
		"CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "PUR", "PUR", "PUR", "PUR",
		"PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "PUR", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM",
		"CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM", "CLM",
		"CLM", "CLM", "CLM", "CLM", "CLM", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS",
		"IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "CLM", "CLM", "CLM",
		"PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PJL", "PJL", "PJL", "PJL", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV",
		"IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS",
		"IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS",
		"IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS",
		"IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS",
		"IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "GBR", "GBR", "GBR", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX",
		"CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "KHV",
		"KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV",
		"KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV",
		"KHV", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "PEL", "PEL", "ACB", "ACB", "ACB", "ACB", "ACB",
		"PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL",
		"PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "ACB", "ACB", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL",
		"PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "ACB", "ACB", "ACB", "ACB", "ACB", "PEL", "PEL", "PEL", "PEL",
		"PEL", "PEL", "PEL", "ACB", "ACB", "ACB", "ACB", "ACB", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV",
		"KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "ACB", "ACB", "ACB", "ACB", "KHV", "KHV", "KHV",
		"KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV",
		"KHV", "KHV", "PEL", "PEL", "ACB", "PEL", "PEL", "PEL", "ACB", "ACB", "ACB", "KHV", "KHV", "KHV", "KHV", "KHV",
		"KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "KHV", "ACB", "ACB", "PEL", "PEL",
		"PEL", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX",
		"CDX", "CDX", "CDX", "CDX", "CDX", "GBR", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS", "IBS",
		"IBS", "IBS", "IBS", "CDX", "PEL", "PEL", "ACB", "ACB", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL",
		"PEL", "PEL", "PEL", "ACB", "ACB", "ACB", "ACB", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "PEL", "ACB",
		"ACB", "ACB", "PEL", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB",
		"PEL", "PEL", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX",
		"CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX", "CDX",
		"CDX", "CDX", "CDX", "CDX", "CDX", "ACB", "ACB", "PEL", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB",
		"ACB", "GWD", "GWD", "GWD", "GWD", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "PJL", "PJL",
		"PJL", "PJL", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "KHV", "KHV", "KHV", "KHV", "ACB", "ACB", "ACB",
		"ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "ACB", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "ACB", "ACB",
		"GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "PJL", "PJL", "PJL", "PJL", "PJL", "GWD", "GWD", "GWD",
		"GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "PJL", "PJL", "PJL",
		"PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "PJL", "PJL", "PJL",
		"PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD",
		"PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD",
		"GWD", "GWD", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "GWD", "GWD",
		"GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD",
		"GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD",
		"GWD", "GWD", "GWD", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN",
		"ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "GWD", "GWD", "BEB", "BEB", "BEB", "BEB", "PJL", "PJL", "PJL", "PJL",
		"PJL", "PJL", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "GWD", "MSL", "MSL", "MSL", "MSL",
		"MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL",
		"MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN",
		"ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN",
		"ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN",
		"ESN", "ESN", "MSL", "MSL", "MSL", "MSL", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "GWD", "GWD", "GWD", "GWD",
		"GWD", "GWD", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN",
		"ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "MSL",
		"MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL",
		"MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL",
		"MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "PJL", "PJL", "PJL", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN", "ESN",
		"ESN", "GWD", "GWD", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL", "MSL",
		"MSL", "MSL", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "PJL",
		"PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "STU", "STU", "STU", "STU", "STU", "PJL", "PJL", "PJL", "PJL",
		"PJL", "PJL", "PJL", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU",
		"STU", "STU", "STU", "STU", "STU", "STU", "PJL", "PJL", "PJL", "PJL", "PJL", "PJL", "STU", "ITU", "ITU", "ITU",
		"ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "STU", "STU", "STU", "STU", "STU", "ITU", "STU", "STU",
		"STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "PJL", "PJL", "PJL", "ITU", "ITU", "ITU",
		"ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU",
		"ITU", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB",
		"BEB", "BEB", "BEB", "BEB", "BEB", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU",
		"STU", "STU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU",
		"ITU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "BEB", "BEB",
		"BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB",
		"BEB", "BEB", "BEB", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "ITU", "ITU", "ITU", "ITU",
		"ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "STU", "STU", "STU", "STU", "STU", "STU", "STU", "STU",
		"ITU", "ITU", "STU", "STU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "STU", "STU",
		"STU", "STU", "STU", "STU", "STU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "ITU", "STU", "ITU", "ITU",
		"ITU", "ITU", "ITU", "ITU", "ITU", "STU", "STU", "STU", "STU", "ITU", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB",
		"BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "BEB",
		"BEB", "BEB", "BEB", "BEB", "BEB", "BEB", "ITU", "ITU", "ITU", "ITU", "ITU", "STU", "ITU", "ITU", "ITU", "ITU",
		"ITU", "ITU", "ITU", "STU", "STU", "ITU", "ITU", "ITU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU",
		"CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU",
		"CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU",
		"CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU",
		"CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU",
		"CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU",
		"CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "CEU", "YRI", "YRI", "YRI", "YRI", "YRI",
		"YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "CHB", "CHB",
		"CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB",
		"CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB",
		"CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB",
		"CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB",
		"CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB",
		"CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB", "CHB",
		"CHB", "CHB", "CHB", "CHB", "CHB", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI",
		"YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI",
		"YRI", "YRI", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT",
		"JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT",
		"JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT",
		"JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT",
		"JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK",
		"LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT",
		"JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT",
		"JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "JPT", "YRI", "YRI", "YRI", "YRI",
		"YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI",
		"YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI",
		"YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI",
		"YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK",
		"LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK",
		"LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK",
		"LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK",
		"LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK",
		"LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "LWK", "ASW", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL",
		"MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "ASW", "ASW", "ASW",
		"ASW", "ASW", "ASW", "ASW", "ASW", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL",
		"MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL",
		"MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL", "MXL",
		"MXL", "MXL", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW",
		"ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW",
		"ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "ASW",
		"ASW", "ASW", "ASW", "ASW", "ASW", "ASW", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI",
		"TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI",
		"TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI",
		"TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI",
		"TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI",
		"TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI",
		"TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI", "TSI",
		"TSI", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH",
		"GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH",
		"GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH",
		"GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH",
		"GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH",
		"GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH",
		"GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH", "GIH"]
else:
	inds = args.inds

if args.exclude == False:
	# include only CHS, YRI, and CEU (representing only the 2 pure ancestral components from East Asia, Africa and Europe, respectively, in 1000G 2015 paper)
	pops_exclude = ['BEB', 'ITU', 'MXL', 'PUR', 'ACB', 'TSI', 'MSL', 'PEL', 'GWD', 'GIH', 'ASW', 'CHB', 'CLM', 'FIN', 'KHV', 'JPT', 'IBS', 'GBR', 'LWK', 'CDX', 'STU', 'ESN', 'PJL']
else:
	pops_exclude = []
	for pop in args.exclude:
		pops_exclude.append(pop)

pops = dict()
for pop in sorted(set(inds)):
	if pop not in pops_exclude:
		pops[pop] = []
num_pops = len(pops)

# open vcf file and output file
with gzopen(vcf_name, "rt") as vcf:
	with open(output, "w") as out:

		print(f"\n\n###########################################################################\n")
		print(f"Find AIMs in VCF files by calculating each SNP's informativeness assignment\n")
		print(f"###########################################################################\n\n")

		# skip headers and, when done, rewind to first snp line
		line = vcf.readline()
		while line[0] == "#":
			pos = vcf.tell()
			line = vcf.readline()
		vcf.seek(pos)

		print(f"Processing chromosome {chr}...\n")

		# loop through lines (snps) to obtain genotypes

		m = args.skip
		# skip first m snps, if specified
		if m > 0:
			for i, snp in enumerate(vcf, 0):
				if i == m-1:
					break

		for i, snp in enumerate(vcf, m+1):
			line = snp.split()
			if i % 10000 == 0:
				print(f"\tReached {line[2]}, {i}th SNP of chromosome {chr}. Continue processing...\n")
			genotypes = line[9:]
			# calculate In
			In = informativeness_assignment(genotypes, inds, pops, num_pops)
			if In >= threshold:
				# write to output the good ones
				rs = f"{line[2]}"
				out.write(f"{rs} {In}\n")
				out.flush()
				print(f"\t{i}th SNP, {rs}, has a In = {In}")

		print(f"Finished processing chromosome {chr}! Results written to {output}\n")
