import utils	# functions file, in this folder
from numpy import savetxt, loadtxt, transpose
from argparse import ArgumentParser
from random import choice, sample, seed
from contextlib import ExitStack
from os.path import isfile, getsize
from sys import exit


# parse arguments
parser = ArgumentParser(description='Mandatory commands: file path and at least a pop_probs name')

parser.add_argument('-path', type=str, default="./", help='Path of the infile name, whose name format is assumed to be [path][pop][suffix].(t)hap')
parser.add_argument('-pops', type=str, nargs='*', help='List of population names to sample from')
parser.add_argument('-suffix', type=str, default='_2000AIMs.impute', help='File suffix after population name, e.g. .impute')
parser.add_argument('-train', type=int, default=1, help='Number of simulated individuals for ANN training, a single one by default')
parser.add_argument('-test', type=int, default=0, help='Number of simulated individuals for ANN testing, none by default')
parser.add_argument('-set', type=float, default=0.5, help='Fraction of the input file that will be used for simulating the individuals for the training, the remaining fraction will be used for the simulations for the testing. 0.5, i.e. the 1st and 2nd halves of each population file by default')
parser.add_argument('-rc', type=int, nargs='*', default=77, help='Number of recombinations per individual, 77 by default. If -admix is selected, a 2nd integer is allowed (see -admix)')
parser.add_argument('-mut', type=int, default=0, help='Number of mutations per individual, none by default')
parser.add_argument('-mode', type=str, default='geno', help='geno (default) indicates the genotype (e.g. 0, 1, 1, 2 constitutes 4 individuals), 01 indicates to output the phased genotypes in 2-digit format (e.g. the previous individuals are 0 0 0 1 1 0 1 1), while the option 0123 outputs them in single-digit format (e.g. the previous individuals are coded as 0 1 2 3)')
parser.add_argument('-admix', type=float, nargs='*', default=False, help='List of floats indicating the proportion of chunks coming from the individual of each pop, in the order of -pops. If two floats per population are provided, they are interpreted as lower and upper bounds whence randomly pick a proportion (all proportions will add up to 1). There can be an additional second integer in -rc, in which case the two ints are also interpreted as bounds for the number of recombination events')
parser.add_argument('-seed', type=int, default=False, help='Seed number')
parser.add_argument('-out', type=str, default="simulated_genotypes/simulated_genotypes", help='Path and file name for the simulated genotypes, by default simulated_genotypes/simulated_genotypes(.train/test.gz)')
# parser.add_argument('-thread', type=int, default=1, help='Number of threads to use')
parser.add_argument('-v', action='store_true', help='Output includes the recombinant haplotypes, original parents haplotypes, their column number in the original .hap infile, and the cut points of the crossovers. '
														'WARNING: '
														'This option will default to a single simulated individual. '
														'Select this if you want to check the crossover procedure. '
														'This output is not valid for the ANN training nor testing')
args = parser.parse_args()

if args.seed is not False:
	seed(args.seed)

# infile name
# Example: "path/to/file/YRI.impute[.(t)hap]"
filenames = [(f"{args.path}{population}{args.suffix}", population) for population in args.pops]

# encode population one-hot probabilities for cost function (dict),
# and list pops for admixture mode
dict_prob_pops = dict()
list_pops = []
for i, pop in enumerate(args.pops):
	list_pops.append(pop)

	probs = ['0' if j != i else '1' for j, pop in enumerate(args.pops)]
	dict_prob_pops[pop] = probs

num_pops = len(list_pops)

# check that there already are all transposed impute (.thap) files, and if so, simulate offspring
transposed_file_exist = False
while transposed_file_exist is False:
	try:
		# open haplotypes files
		print(f"\nOpening {args.path}{args.pops}{args.suffix}.thap files...\n")

		haplotypes_files = [(ExitStack().enter_context(open(f"{filename[0]}.thap")), filename[1])
							if getsize(f"{filename[0]}.thap") > 0 else exit(f"Error: At least one .thap file ({filename[0]}.thap) exists but is empty. Please delete all empty .thap files before running the program.")
							for filename in filenames]

		# don't repeat this loop, since all the transposed impute files exist
		transposed_file_exist = True

		# delimite haplotypes to use for each subset (training or testing) for each population
		limit_hap = utils.find_limit_hap(haplotypes_files, args.set)

		# create and store matrices for each set
		haplotypes_matrices_training = utils.create_haplotypes_matrices(haplotypes_files, "training", limit_hap)
		haplotypes_matrices_testing = utils.create_haplotypes_matrices(haplotypes_files, "testing", limit_hap)

		# check for empty files
		if args.train > 0:
			for pop in haplotypes_matrices_training:
				if len(haplotypes_matrices_training[pop]) < 4:
					exit("Error: There are less than 4 training haplotypes in at least one .thap file.")
		if args.test > 0:
			for pop in haplotypes_matrices_testing:
				if len(haplotypes_matrices_testing[pop]) < 4:
					exit("Error: There are less than 4 testing haplotypes in at least one .thap file.")
		# check that the number of loci be the same in all populations
		try:
			if len(set([len(haplotypes_matrices_training[pop][0]) for pop in args.pops])) > 1:
				exit("Error: The number of loci is not the same in all populations.")
		except IndexError:
			if len(set([len(haplotypes_matrices_testing[pop][0]) for pop in args.pops])) > 1:
				exit("Error: The number of loci is not the same in all populations.")

		# check whether we have one or two values in -rc
		try:
			length_rc = len(args.rc)
			# value/s is/are provided
			if length_rc == 1:
				# single value -rc, specified
				fixed_rc = True
				num_rcs = args.rc[0]
			elif length_rc == 2:
				if args.admix is not False:
					# -rc contains an interval, only for admixture mode
					fixed_rc = False
					num_rcs = max(args.rc)
				else:
					# -rc has an interval but the -admix mode is not selected
					exit("Error: Make sure to not provide more than two integers in -rc when -admix is in use, and only one if the -admix mode is not in use.")
		except TypeError:
			# single value -rc, the default value
			fixed_rc = True
			num_rcs = args.rc

		# if nº mutations or recombinations excede m, raise error
		try:
			num_snps= len(haplotypes_matrices_training[args.pops[0]][0])	# num. loci
		except IndexError:
			num_snps= len(haplotypes_matrices_testing[args.pops[0]][0])
		if args.mut > num_snps or num_rcs > num_snps - 2:
			exit("Error: The specified nº of mutations excedes the nº or loci, and/or the nº of recombination events excedes the nº of loci - 2.")

		if args.admix is not False:
			# ADMIXTURE MODE
			if len(args.admix) == num_pops:
				# fixed proportion of chunks from each pop
				fixed_props = True
			elif len(args.admix) == 2 * num_pops:
				# variable prop. of chunks from each pop
				fixed_props = False
			else:
				exit("Error: the values provided in -admix are not equal or double the number of -pops")

		# open output files (check if they exist already to append to them)
		if args.train == 0 and args.test == 0:
			exit("Error: There must be specified > 0 simulations for, at least, either the training or testing set.")

		trainingFileExists, testingFileExists = isfile(f"{args.out}.train.gz"), isfile(f"{args.out}.test.gz")

		if trainingFileExists:
			# append to existing output file
			try:
				out_training, out_pops_training, out_info_training, old_num_sims = utils.optionalRead(args.out, args.train, "train")
			except TypeError:
				pass
		else:
			try:
				out_training, out_pops_training, out_info_training = utils.optionalWrite(args.out, args.train, "train")
			except TypeError:
				pass

		if testingFileExists:
			# append to existing output file
			try:
				out_testing, out_pops_testing = utils.optionalRead(args.out, args.test, "test")
			except TypeError:
				pass
		else:
			try:
				out_testing, out_pops_testing = utils.optionalWrite(args.out, args.test, "test")
			except TypeError:
				pass

		# # number of simulated individuals to create (= number of loop iterations)
		# if args.v:
		# 	# print header
		# 	out_training.write("# Number of loci\n".encode('utf-8'))
		# 	if args.admix is False:
		# 		out_training.write("# Non-admixture mode\n".encode('utf-8'))
		# 	else:
		# 		out_training.write("# Admixture mode\n".encode('utf-8'))
		# 	# only one simulation
		# 	total_sim = 1
		# else:
		# 	# no header, all specified simulations
		total_sim = args.train + args.test

		num_pops = len(args.pops)

		try:
			# print number of loci in training
			out_info_training.write(f"{num_snps}\n".encode('utf-8'))
			# print number of pops (= number of output neurons)
			out_info_training.write(f"{num_pops}\n".encode('utf-8'))
			# print num.of simulated offspring
			if trainingFileExists:
				# sum to number of simulations in original file
				sum_sims = args.train + int(old_num_sims)
				out_info_training.write(f"{sum_sims}\n".encode('utf-8'))
			else:
				out_info_training.write(f"{args.train}\n".encode('utf-8'))

			out_info_training.close()
		except NameError:
			# no training file requested
			pass

		# loop and keep the count of the number of simulated individuals
		num_sim_inds = 0

		# start loop of simulations
		print(f"Creating {args.train} training and {args.test} testing simulations...\n")

		while num_sim_inds < total_sim:
			if num_sim_inds != 0 and num_sim_inds % 10000 == 0:
				print(f"\tReached {num_sim_inds} simulations, continue...\n")

			if args.admix is False:
				# no admixture
				# in each iteration, choose a single random population for creating non-admixed offspring
				population = choice(args.pops)
			else:
				# ADMIXTURE MODE
				if fixed_rc is not True:
					# determine num of recombinations in this iteration
					num_rcs = sample(range(args.rc[0], args.rc[1] + 1), 1)[0]

			# store random cuts
			random_cuts = utils.make_cuts(num_rcs, num_snps)

			# check whether it is still training or already testing set
			if num_sim_inds < args.train:
				# still training
				if args.admix is False:
					# no admixture mode
					# store random ints to pick 4 parent haplotypes
					random_indexes = utils.get_random_indexes(len(haplotypes_matrices_training[population]))
					# store random haplotypes
					parent_haplotypes = utils.get_random_parent_haps(haplotypes_matrices_training[population], random_indexes)
					# store 2 random recombinant haplotypes into the offspring individual:
					haplotypes_offspring = utils.crossovers(parent_haplotypes, random_indexes, random_cuts)
				else:
					# ADMIXTURE mode
					# store random haplotypes from ancestral populations
					ancestral_haplotypes = utils.ancestral_haps(haplotypes_matrices_training)
					# create admixed haplotypes
					haplotypes_offspring, pop_probs = utils.crossovers_admix(list_pops, args.admix, ancestral_haplotypes, random_cuts, num_rcs, num_snps)
			else:
				# already testing
				if args.admix is False:
					# no admixture mode
					# store randonumints to pick 4 parent haplotypes
					random_indexes = utils.get_random_indexes(len(haplotypes_matrices_testing[population]))
					# store random haplotypes
					parent_haplotypes = utils.get_random_parent_haps(haplotypes_matrices_testing[population], random_indexes)
					# store 2 random recombinant haplotypes into the offspring individual:
					haplotypes_offspring = utils.crossovers(parent_haplotypes, random_indexes, random_cuts)
				else:
					# ADMIXTURE mode
					# store random haplotypes from ancestral populations
					ancestral_haplotypes = utils.ancestral_haps(haplotypes_matrices_testing)
					# create admixed haplotypes
					haplotypes_offspring, pop_probs = utils.crossovers_admix(list_pops, args.admix, ancestral_haplotypes, random_cuts, num_rcs, num_snps)

			# store phased genotypes of offspring
			if args.mode == 'geno':
				simulated_genotypes_offspring = utils.simulated_genotypes_geno(haplotypes_offspring, args.mut)
			elif args.mode == '01':
				simulated_genotypes_offspring = utils.simulated_genotypes_01(haplotypes_offspring, args.mut)
			else:
				simulated_genotypes_offspring = utils.simulated_genotypes_0123(haplotypes_offspring, args.mut)

			# output the offspring genotype
			if args.v is False:
				# assign probabilities depending on population, for cost function
				if args.admix is False:
					pop_probs = dict_prob_pops[population]

				if num_sim_inds < args.train:
					# training phased genotypes for input and pop name for loss function
					out_training.write(f"{' '.join([i for i in map(str, simulated_genotypes_offspring)])}\n".encode('utf-8'))
					out_pops_training.write(f"{' '.join(pop_probs)}\n".encode('utf-8'))
				else:
					# testing phased genotypes for input and pop name for loss function
					out_testing.write(f"{' '.join([i for i in map(str, simulated_genotypes_offspring)])}\n".encode('utf-8'))
					out_pops_testing.write(f"{' '.join(pop_probs)}\n".encode('utf-8'))
				num_sim_inds += 1

			elif args.v is not False and args.train == 0:
				exit("Error: Verbose mode with 0 training simulations. Exiting...")
			# if verbose, output the phased genotype and pop info in separated lines, plus the recomb haplotype format, original parents haplotypes, their column number in the infile and the cut points
			else:
				# population
				out_training.write("\n# mode\n".encode('utf-8'))
				out_training.write(f"{args.mode}\n\n".encode('utf-8'))

				# population
				out_training.write("# population\n".encode('utf-8'))
				out_training.write(f"{pop_probs}\n\n".encode('utf-8'))

				# phased genotypes
				out_training.write("# phased/unphased genotypes for feeding the ANN, obtained from the recombinant haplotypes\n".encode('utf-8'))
				out_training.write(f"{' '.join([i for i in map(str, simulated_genotypes_offspring)])}\n\n".encode('utf-8'))

				# recomb_hap1
				out_training.write("# recombinant haplotype coming from parent 1\n".encode('utf-8'))
				out_training.write(f"{' '.join([i for i in map(str, haplotypes_offspring['parent1_recomb'][0])])}\n".encode('utf-8'))
				# recomb_hap2
				out_training.write("# recombinant haplotype coming from parent 2\n".encode('utf-8'))
				out_training.write(f"{' '.join([i for i in map(str, haplotypes_offspring['parent2_recomb'][0])])}\n\n".encode('utf-8'))

				# in the impute infile, the index (column number) in which is located the parent1 hap1
				out_training.write("# in the impute infile, the index (column number, 0-based) in which is located the parent1 hap1\n".encode('utf-8'))
				index_parent1_hap1 = haplotypes_offspring['parent1_recomb'][1][0]
				out_training.write(f"{index_parent1_hap1}\n".encode('utf-8'))
				# parent1 hap1
				out_training.write("# original haplotype 1 of parent1\n".encode('utf-8'))
				out_training.write(f"{' '.join([i for i in map(str, haplotypes_matrices_training[pop_probs][index_parent1_hap1])])}\n".encode('utf-8'))
				# in the impute infile, the index (column number) in which is located the parent1 hap2
				out_training.write("# index parent1 hap2\n".encode('utf-8'))
				index_parent1_hap2 = haplotypes_offspring['parent1_recomb'][1][1]
				out_training.write(f"{index_parent1_hap2}\n".encode('utf-8'))
				# parent1 hap2
				out_training.write("# original haplotype 2 of parent1\n".encode('utf-8'))
				if num_rcs > 0:
					out_training.write(f"{' '.join([i for i in map(str, haplotypes_matrices_training[pop_probs][index_parent1_hap2])])}\n".encode('utf-8'))
				else:
					out_training.write(f"No recombination\n".encode('utf-8'))
				# cut positions parent1
				out_training.write("# cut positions parent1\n".encode('utf-8'))
				out_training.write(f"{haplotypes_offspring['parent1_recomb'][2]}\n".encode('utf-8'))

				# in the impute infile, the index (column number) in which is located the parent2 hap1
				out_training.write("# in the impute infile, the index (column number, 0-based) in which is located the parent2 hap1\n".encode('utf-8'))
				index_parent2_hap1 = haplotypes_offspring['parent2_recomb'][1][0]
				out_training.write(f"{index_parent2_hap1}\n".encode('utf-8'))
				# parent2 hap1
				out_training.write("# original haplotype 1 of parent2\n".encode('utf-8'))
				out_training.write(f"{' '.join([i for i in map(str, haplotypes_matrices_training[pop_probs][index_parent2_hap1])])}\n".encode('utf-8'))
				# in the impute infile, the index (column number) in which is located the parent2 hap2
				out_training.write("# index parent2 hap2\n".encode('utf-8'))
				index_parent2_hap2 = haplotypes_offspring['parent2_recomb'][1][1]
				out_training.write(f"{index_parent2_hap2}\n".encode('utf-8'))
				# parent2 hap2
				out_training.write("# original haplotype 2 of parent2\n".encode('utf-8'))
				if num_rcs > 0:
					out_training.write(f"{' '.join([i for i in map(str, haplotypes_matrices_training[pop_probs][index_parent2_hap2])])}\n".encode('utf-8'))
				else:
					out_training.write(f"No recombination\n".encode('utf-8'))
				# cut positions parent2
				out_training.write("# cut positions parent2\n".encode('utf-8'))
				out_training.write(f"{haplotypes_offspring['parent2_recomb'][2]}\n".encode('utf-8'))

				num_sim_inds += 1

		# close files
		try:
			out_training.close()
		except NameError:
			pass
		try:
			out_pops_training.close()
		except NameError:
			pass
		try:
			out_testing.close()
		except NameError:
			pass
		try:
			out_pops_testing.close()
		except NameError:
			pass

		print(f"\tFinished simulating individuals! Results written to {args.out}\n")


	# if there is any no transposed impute (.thap) file yet, create it and write it out for this iteration and the following
	except FileNotFoundError:
		# make list of pops for which there is no haplotypes (.thap) file yet, just impute file
		for file in filenames:
			try:
				open(f"{file[0]}.thap", "r")
			except FileNotFoundError:
				# read files in .hap format and create haplotype files (.thap)
				# store file into numpy array and transpose
				haplotypes_matrix = loadtxt(f"{file[0]}.hap", delimiter=" ").transpose()
				# save to text to reuse it
				savetxt(f"{file[0]}.thap", haplotypes_matrix, fmt='%0.0f', delimiter=' ', newline='\n')
		# now return to the top of the while loop and create the simulated offspring
