import numpy as np

def informativeness_assignment(line, inds_list, pops_dict, num_pops):

	for genotype, pop in zip(line, inds_list):
		try:
			# create dict for relative frequencies
			genotype = sum([int(allele) for allele in genotype.split("|")])
			pops_dict[pop].append(genotype)
		except KeyError:
			# excluded pop
			pass

	rel_freqs_allele1 = []
	rel_freqs_allele2 = []
	transf_rel_freqs_allele1 = []
	transf_rel_freqs_allele2 = []

	for pop in pops_dict:
		# relative freq
		num_inds = len(pops_dict[pop])
		rel_freq_allele1 = sum(pops_dict[pop]) / (num_inds*2)
		rel_freq_allele2 = 1 - rel_freq_allele1
		if rel_freq_allele1 != 0 and rel_freq_allele1 != 1:
			transf_rel_freqs_allele1.append(rel_freq_allele1 * np.log(rel_freq_allele1))
			transf_rel_freqs_allele2.append(rel_freq_allele2 * np.log(rel_freq_allele2))
		# mean freq
		rel_freqs_allele1.append(rel_freq_allele1)
		rel_freqs_allele2.append(rel_freq_allele2)

		# clear list genotypes in pop, otherwise they will be stored for the next snp...
		pops_dict[pop].clear()

	mean_freq_allele1 = sum(rel_freqs_allele1) / num_pops
	if mean_freq_allele1 != 0 and mean_freq_allele1 != 1:
		mean_freq_allele2 = 1 - mean_freq_allele1
		first_comp_allele1 = -mean_freq_allele1 * np.log(mean_freq_allele1)
		first_comp_allele2 = -mean_freq_allele2 * np.log(mean_freq_allele2)
	else:
		# monomorphic for the selected populations
		return 0

	second_comp_allele1 = sum(transf_rel_freqs_allele1) / num_pops
	second_comp_allele2 = sum(transf_rel_freqs_allele2) / num_pops

	allele1 = first_comp_allele1 + second_comp_allele1
	allele2 = first_comp_allele2 + second_comp_allele2

	# return In
	return allele1+allele2
