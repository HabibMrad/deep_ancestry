from random import randint, sample, uniform, choice, seed
from operator import add
import gzip
from linecache import getline


def optionalWrite(outfilename, num_sims, mode):
	# not opening the training or testing files if they are not requested
	if num_sims > 0:
		simulated_genotypes_file = gzip.open(f"{outfilename}.{mode}.gz", "wb")
		pops_file = gzip.open(f"{outfilename}.pops.{mode}.gz", "wb")
		if mode == "train":
			info_file = gzip.open(f"{outfilename}.info.train.gz", "wb")
			return simulated_genotypes_file, pops_file, info_file
		else:
			return simulated_genotypes_file, pops_file


def optionalRead(outfilename, num_sims, mode):
	# not opening the existing training or testing files if they are not requested
	if num_sims > 0:
		simulated_genotypes_file = gzip.open(f"{outfilename}.{mode}.gz", "ab")
		pops_file = gzip.open(f"{outfilename}.pops.{mode}.gz", "ab")
		if mode == "train":
			# get previous number of simulations
			old_info = gzip.open(f"{outfilename}.info.train.gz", "rt")
			old_num_sims = old_info.readlines()[2]
			# write it new
			info_file = gzip.open(f"{outfilename}.info.train.gz", "wb")
			return simulated_genotypes_file, pops_file, info_file, old_num_sims
		else:
			return simulated_genotypes_file, pops_file


def find_limit_hap(list_tuples_of_file_and_pop, fraction):
	# find subsets haplotypes for each subset at each population
	limit_haps = {}
	for haplotypes_file in list_tuples_of_file_and_pop:
		pop = haplotypes_file[1]
		num_haps = len(haplotypes_file[0].readlines())

		limit_haps[pop] = int(num_haps * fraction)
	return limit_haps


def create_haplotypes_matrices(list_tuples_of_file_and_pop, subset, limit_hap):
	# create matrices
	haplotypes_matrices = {}

	for haplotypes_file in list_tuples_of_file_and_pop:
		haplotypes_file[0].seek(0)
		pop = haplotypes_file[1]
		haplotypes_matrices[pop] = []
		if subset is "training":
			for i, haplotype in enumerate(haplotypes_file[0]):
				if i < limit_hap[pop]:
					# this takes 9 secs, cause it has to split the haplotypes into loci
					haplotypes_matrices[haplotypes_file[1]].append([locus for locus in map(int, haplotype.split())])
				else:
					# training subset reached
					break
		else:
			for i, haplotype in enumerate(haplotypes_file[0]):
				if i >= limit_hap[pop]:
					# this takes 9 secs, cause it has to split the haplotypes into loci
					haplotypes_matrices[haplotypes_file[1]].append([locus for locus in map(int, haplotype.split())])
				else:
					# testing subset not reached yet
					pass
	return haplotypes_matrices


def get_random_indexes(length_hap_matrix):
	# gets 4 ints from 0 to the number of haplotypes - 1 for then picking 4 random parent haplotypes
	hap0, hap1, hap2, hap3 = 0, 0, 0, 0
	# do not allow repeating same index
	while len({hap0, hap1, hap2, hap3}) < 4:
		hap0 = randint(0, length_hap_matrix - 1)
		hap1 = randint(0, length_hap_matrix - 1)
		hap2 = randint(0, length_hap_matrix - 1)
		hap3 = randint(0, length_hap_matrix - 1)
	return {'hap0': hap0, 'hap1': hap1, 'hap2': hap2, 'hap3': hap3}


def get_random_parent_haps(haplotype_matrix, random_integers):
	# pick 4 random parent haplotypes
	# notice that the same haplotype can be picked more than once
	parent_haplotypes = dict()
	parent_haplotypes['hap0'] = haplotype_matrix[random_integers['hap0']]
	parent_haplotypes['hap1'] = haplotype_matrix[random_integers['hap1']]
	parent_haplotypes['hap2'] = haplotype_matrix[random_integers['hap2']]
	parent_haplotypes['hap3'] = haplotype_matrix[random_integers['hap3']]

	return parent_haplotypes


def make_cuts(num_cuts, m):
	# get x random cut points for parent1 and 2
	cut_points_dict = dict()

	for parent in ['parent1', 'parent2']:
		# get unique random points (not the 1st and last snps), and append to the parent key
		cut_points = sample(range(1, m), num_cuts)

		# make sure the bp where the first cut takes place is located BEFORE the second, etc.
		cut_points_dict[parent] = sorted(cut_points)

	return cut_points_dict


def crossovers(haplotypes, index_haps, random_cuts):
	# simulate crossovers for each pair of parent haplotypes A and B, and then pick the a+b+a... recombinant haplotype from each parent
	# assign parent1 the hap0 and hap1 haplotypes in the list etc.
	dict_parents = {'parent1': {'haplotypes': (haplotypes['hap0'], haplotypes['hap1']), 'index_haps': (index_haps['hap0'], index_haps['hap1'])}, 'parent2': {'haplotypes': (haplotypes['hap2'], haplotypes['hap3']), 'index_haps': (index_haps['hap2'], index_haps['hap3'])}}

	# get x crossovers in each pair of parent genomes
	offspring = dict()
	# if more than 1 recombination events specified
	if len(random_cuts['parent1']) not in [0, 1]:
		for parent in ['parent1', 'parent2']:
			# crossovers
			recombinant = []
			for i, cut in enumerate(random_cuts[parent]):
				# not first nor last cuts:
				if i not in [0, len(random_cuts[parent])-1]:
					if i%2 is 1:
						recombinant = recombinant + dict_parents[parent]['haplotypes'][1][previous_cut_pos:cut]
						previous_cut_pos = cut
					elif i%2 is 0:
						recombinant = recombinant + dict_parents[parent]['haplotypes'][0][previous_cut_pos:cut]
						previous_cut_pos = cut
				# first cut
				elif i is 0:
					recombinant = recombinant + dict_parents[parent]['haplotypes'][0][:cut]
					previous_cut_pos = cut
				# last cut
				else:
					if i%2 is 1:
						recombinant = recombinant + dict_parents[parent]['haplotypes'][1][previous_cut_pos:cut] + dict_parents[parent]['haplotypes'][0][cut:]
					elif i%2 is 0:
						recombinant = recombinant + dict_parents[parent]['haplotypes'][0][previous_cut_pos:cut] + dict_parents[parent]['haplotypes'][1][cut:]
			# store the recombinant chromosome a+b+a...
			offspring[f"{parent}_recomb"] = (recombinant, (dict_parents[parent]['index_haps'][0], dict_parents[parent]['index_haps'][1]), random_cuts[parent])
	# -rc 1
	elif len(random_cuts['parent1']) is 1:
		for parent in ['parent1', 'parent2']:
			cut = random_cuts[parent][0]
			recombinant = dict_parents[parent]['haplotypes'][0][:cut] + dict_parents[parent]['haplotypes'][1][cut:]
			offspring[f"{parent}_recomb"] = (recombinant, (dict_parents[parent]['index_haps'][0], dict_parents[parent]['index_haps'][1]), random_cuts[parent])
	# -rc 0
	else:
		for parent in ['parent1', 'parent2']:
			offspring[f"{parent}_recomb"] = (dict_parents[parent]['haplotypes'][0], (dict_parents[parent]['index_haps'][0], "No recombination"), 0)

	# example offspring with 2 crossovers:
	# {'parent1_recomb': (recombinant haplotype, (number of the column in the .hap input that conforms the hap a of this parent, same for hap b), (position of first cut, position of 2nd cut)), (same for parent2)}
	# {'parent1_recomb': ([0, 0, 0, ..., 0, 0, 0], (185, 23), (133624, 206948)), 'parent2_recomb': ([0, 0, 0, ..., 0, 0, 0], (45, 157), [206673, 101316])}

	# return offspring
	return offspring


def mutations_geno(geno_offspring, num_mutations):
	list_positions = [pos for pos in range(len(geno_offspring))]
	for mutation in range(num_mutations):
		# random position for mutation and remove position from list, so it is not chosen again
		mutated_pos = choice(list_positions)
		list_positions.remove(mutated_pos)
		# get locus genotype in that position
		mutated_locus = geno_offspring[mutated_pos]
		# mutate
		if mutated_locus in [0, 2]:
			new_genotype = 1
		else:
			new_genotype = choice([0, 2])
		# apply mutation to phased genotypes
		geno_offspring[mutated_pos] = new_genotype
	return geno_offspring


def mutations_01(phased_geno_offspring, num_mutations):
	list_positions = [pos for pos in range(len(phased_geno_offspring))]
	for mutation in range(num_mutations):
		# random position for mutation and remove position from list, so it is not chosen again
		mutated_pos = choice(list_positions)
		list_positions.remove(mutated_pos)
		# get locus genotype in that position
		mutated_locus = phased_geno_offspring[mutated_pos]
		# mutate
		if mutated_locus is 0:
			new_genotype = 1
		else:
			new_genotype = 0
		# apply mutation to phased genotypes
		phased_geno_offspring[mutated_pos] = new_genotype

	return phased_geno_offspring


def mutations_0123(phased_geno_offspring, num_mutations):
	list_positions = [pos for pos in range(len(phased_geno_offspring))]
	for mutation in range(num_mutations):
		# random position for mutation and remove position from list, so it is not chosen again
		mutated_pos = choice(list_positions)
		list_positions.remove(mutated_pos)
		# get locus genotype in that position
		mutated_locus = phased_geno_offspring[mutated_pos]
		# mutate
		if mutated_locus in [0, 3]:
			new_genotype = choice([1, 2])
		else:
			new_genotype = choice([0, 3])
		# apply mutation to phased genotypes
		phased_geno_offspring[mutated_pos] = new_genotype

	return phased_geno_offspring


def simulated_genotypes_geno(offspring_phased_haplotypes, num_mutations):
	# output genotypes
	parent1_recomb, parent2_recomb = [offspring_phased_haplotypes[parent][0] for parent in offspring_phased_haplotypes.keys()]
	genotypes_offspring = list(map(add, parent1_recomb, parent2_recomb))
	# so that we get the genotype info, which is what we feed to the ANN
	# 0 0 --> 0
	# 0 1 --> 1
	# 1 0 --> 1
	# 1 1 --> 2
	if num_mutations > 0:
		return mutations_geno(genotypes_offspring, num_mutations)
	else:
		return genotypes_offspring


def simulated_genotypes_01(offspring_phased_haplotypes, num_mutations):
	# output the phased genotypes in 2-digit format (poss. digits: 0/1), for feeding to the NN

	phased_geno_offspring = [allele for genotype in zip(offspring_phased_haplotypes[list(offspring_phased_haplotypes.keys())[0]][0],
														offspring_phased_haplotypes[list(offspring_phased_haplotypes.keys())[1]][0]) for allele in genotype]
	if num_mutations > 0:
		return mutations_01(phased_geno_offspring, num_mutations)
	else:
		return phased_geno_offspring


def simulated_genotypes_0123(offspring_phased_haplotypes, num_mutations):
	# output the phased genotypes in single-digit format (possible digits: 0/1/2/3), for feeding to the NN

	parent1_recomb, parent2_recomb = [offspring_phased_haplotypes[parent][0] for parent in offspring_phased_haplotypes.keys()]

	def bin_geno(allele1, allele2):
		# recode haplotypes, read as binary (00,01,10,11), as decimal integers (0,1,2,3)
		return int(f'{allele1}{allele2}', 2)

	phased_geno_offspring = list(map(bin_geno, parent1_recomb, parent2_recomb))

	if num_mutations > 0:
		return mutations_0123(phased_geno_offspring, num_mutations)
	else:
		return phased_geno_offspring


def chunks_per_pop(list_pops, props_pop, num_rcs):

	# normalize so props sum 1
	if sum(props_pop) is not 1:
		sum_props = sum(props_pop)
		props_pop = [prop / sum_props for prop in props_pop]

	chunks_per_pop = dict()
	num_exp_chunks = num_rcs + 1

	# pops with 0 chunks as defined in arguments
	zero_pops = []
	for pop, prop in zip(list_pops, props_pop):
		if prop == 0:
			zero_pops.append(pop)
		chunks_per_pop[pop] = round(prop * num_exp_chunks)

	# check that the observed nº of chunks per pop add up to the expected total nº of chunks (nº of recombinations + 1)
	num_obs_chunks = sum(chunks_per_pop.values())
	if num_obs_chunks != num_exp_chunks:
		if num_obs_chunks > num_exp_chunks:
			while num_obs_chunks > num_exp_chunks:
				random_pop = choice(list_pops)
				while random_pop in zero_pops:
					random_pop = choice(list_pops)
				chunks_per_pop[random_pop] = chunks_per_pop[random_pop] - 1
				num_obs_chunks = sum(chunks_per_pop.values())
		else:
			while num_obs_chunks < num_exp_chunks:
				random_pop = choice(list_pops)
				while random_pop in zero_pops:
					random_pop = choice(list_pops)
				chunks_per_pop[random_pop] = chunks_per_pop[random_pop] + 1
				num_obs_chunks = sum(chunks_per_pop.values())

	return chunks_per_pop


def props_per_pop(list_pops, props_pop, num_pops, num_rcs):
	# 1st randomly select prop of chunks per pop in this iteration
	random_props = []
	for i in range(0, num_pops * 2, 2):
		random_props.append(uniform(props_pop[i], props_pop[i + 1]))

	# 2nd determine num of chunks from each pop in this iteration
	props_pop = chunks_per_pop(list_pops, random_props, num_rcs)

	return props_pop


def ancestral_haps(haplotypes):
	# get random indexes and store random haplotypes
	ancestral_haplotypes = dict()
	for pop in haplotypes:
		length_hap_matrix = len(haplotypes[pop])
		# gets 2 ints from 0 to the number of haplotypes for then picking 2 random parent haplotypes per population
		index0, index1 = 0, 0
		# do not allow repeating same index
		while len({index0, index1}) < 2:
			index0 = randint(0, length_hap_matrix - 1)
			index1 = randint(0, length_hap_matrix - 1)
		# get haplotypes corresponding to indexes
		hap0 = haplotypes[pop][index0]
		hap1 = haplotypes[pop][index1]
		ancestral_haplotypes[pop] = {'parent1': hap0, 'parent2': hap1}

	return ancestral_haplotypes


def crossovers_admix(list_pops, props_pop, ancestral_haplotypes, random_cuts, num_cuts, num_snps):

	# create lists of chunks per hap
	lists_chunks = {'parent1': [], 'parent2': []}
	for hap in ['parent1', 'parent2']:

		# calculate chunks_counter
		if len(list_pops) == len(props_pop):
			chunks_counter = chunks_per_pop(list_pops, props_pop, num_cuts)
		else:
			# randomly select prop of chunks per pop in this iteration, and then determine the num of chunks per pop
			chunks_counter = props_per_pop(list_pops, props_pop, len(list_pops), num_cuts)

		# create list of chunks
		for cut in range(num_cuts + 1):
			if sum(chunks_counter.values()) == 0:
				exit("Error: num_cuts+1 larger than nº of chunks in chunks_counter")
			# select random pop from available
			pop = choice(list(chunks_counter.keys()))
			while chunks_counter[pop] == 0:
				# remove pops with 0 chunks left
				chunks_counter.pop(pop, None)
				pop = choice(list(chunks_counter.keys()))
			# append pop name to list
			lists_chunks[hap].append(pop)
			# subtract 1 to the pop used
			chunks_counter[pop] = chunks_counter[pop] - 1

	# create list for storing the amount of SNPs from each population
	final_prop_pops = dict()
	for pop in list_pops:
		final_prop_pops[pop] = 0

	# actually create the recombinant admixed haplotypes
	dict_parents = dict()
	if len(random_cuts['parent1']) != 0:
		# 1 or more recombination events specified
		for parent in ['parent1', 'parent2']:
			# crossovers
			recombinant = []
			previous_cut_pos = 0
			for i, cut in enumerate(random_cuts[parent]):
				chunk_pop = lists_chunks[parent].pop(0)
				# not first cut:
				if i is not 0:
					recombinant = recombinant + ancestral_haplotypes[chunk_pop][parent][previous_cut_pos:cut]
					final_prop_pops[chunk_pop] += cut - previous_cut_pos
					previous_cut_pos = cut
				# first cut
				else:
					recombinant = recombinant + ancestral_haplotypes[chunk_pop][parent][:cut]
					final_prop_pops[chunk_pop] += cut
					previous_cut_pos = cut
			# last chunk
			chunk_pop = lists_chunks[parent].pop(0)
			recombinant = recombinant + ancestral_haplotypes[chunk_pop][parent][previous_cut_pos:]
			final_prop_pops[chunk_pop] += num_snps - previous_cut_pos
			# store the recombinant chromosome
			dict_parents[parent] = recombinant, ['2nd element for compatibility with simulated_genotypes*()']
	else:
		# -rc 0
		for parent in ['parent1', 'parent2']:
			chunk_pop = lists_chunks[parent].pop(0)
			# store the non-recombinant chromosome
			dict_parents[parent] = ancestral_haplotypes[chunk_pop][parent] , ['2nd element for compatibility with simulated_genotypes*()']
			final_prop_pops[chunk_pop] += num_snps

	# calculate final proportion per population
	for pop in final_prop_pops:
		final_prop_pops[pop] = str(final_prop_pops[pop] / (num_snps * 2))
	return dict_parents, list(final_prop_pops.values())
