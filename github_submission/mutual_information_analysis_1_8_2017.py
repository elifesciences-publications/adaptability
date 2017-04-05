import numpy
import sys
import matplotlib.pylab as pt
import matplotlib.cm
import numpy.random
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
import scipy.stats
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import matplotlib

matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.size'] = 10.0
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['lines.markersize'] = 3.5
matplotlib.rcParams['lines.linewidth'] = .5
matplotlib.rcParams['legend.fontsize'] = 8.0
matplotlib.rcParams['axes.linewidth']=.5
matplotlib.rcParams['patch.linewidth']=.5

def permute_within_categories(categories, cat_inds):
	#categories: 1d array where each item has an index indicating which category it belongs to. The category indices need not be consecutive.
	#cat_inds: list of category indices.
	n = len(categories)
	inds = numpy.arange(n) #Original order
	
	permuted_order = numpy.zeros((n,),dtype='int')
	for i in range(len(cat_inds)):
		
		items_in_cat_unpermuted = inds[categories == cat_inds[i]]
		permuted_order[items_in_cat_unpermuted] = numpy.random.permutation(items_in_cat_unpermuted)
	
	return permuted_order

def permute_within_categories_preserve_num_muts(mutation_table, categories, cat_inds):
	#categories: 1d array where each item has an index indicating which category it belongs to. The category indices need not be consecutive.
	#cat_inds: list of category indices.
	
	n = len(categories)
	inds = numpy.arange(n) #Original order
	n_muts = mutation_table.shape[1]
	
	mut_inds = numpy.arange(n_muts)
	
	permuted_mutation_table = numpy.zeros_like(mutation_table)
	for i in range(len(cat_inds)):
		category_indices = inds[categories == cat_inds[i]]
		
		#Construct ordered pair list of which mutations occurred in which clones in this category
		pop_mut_list = []
		for index in category_indices:
			muts = mut_inds[mutation_table[index,:] > .5]
			for mut in muts:
				pop_mut_list.append([index, mut])
		
		#Permute mutations
		n_muts_category = len(pop_mut_list)
		perm = numpy.random.permutation(numpy.arange(n_muts_category))
		
		pop_mut_list = numpy.array(pop_mut_list)
		pop_mut_list_permuted = pop_mut_list
		pop_mut_list_permuted[:,1] = pop_mut_list[perm,1]
		
		#Construct the section of the permuted mutation table for this category
		
		for j in range(len(pop_mut_list_permuted[:,0])):
			mut_loc = pop_mut_list_permuted[j,:]
			permuted_mutation_table[mut_loc[0],mut_loc[1]] = 1
	
	return permuted_mutation_table
		
def calculate_cat_inds(categories):
	
	categories = numpy.array(categories)		
	return numpy.unique(categories)

def calculate_helper_matrix(categories, cat_inds):
	
	#The helper matrix is a utility for quickly summing over specified rows in a table. It is intended to be matrix multiplied by the original mutation table; hence it is n_cats x n_pops
	num_cats = len(cat_inds)
	num_pops = len(categories)
	helper_matrix = numpy.zeros((num_cats,num_pops))
	
	for i in range(num_cats):
		specific_cat_inds = numpy.where(categories == cat_inds[i])
		helper_matrix[i, specific_cat_inds] = 1
	
	return helper_matrix
	
def calculate_entropy_statistic(mutation_table, helper_matrix):
	
	muts_per_gene = numpy.sum(mutation_table, axis = 0)
	
	collapsed_table = numpy.dot(helper_matrix,mutation_table)
	pops_per_category = numpy.dot(helper_matrix,helper_matrix.T)
	#print pops_per_category
	probs = numpy.dot(numpy.linalg.inv(pops_per_category),collapsed_table)
	
	num_genes = mutation_table.shape[1]
	
	entropies = numpy.zeros((num_genes,))
	total_pops = numpy.float(numpy.sum(pops_per_category))
	for i in range(num_genes):
		nonzero_inds = numpy.all([probs[:,i] > 0 , probs[:,i]< 1], axis = 0)
		nonzero_p_hit = probs[:,i][nonzero_inds]
		nonzero_p_no_hit = 1. - nonzero_p_hit
		pops_per_cat_temp = numpy.diag(pops_per_category)[nonzero_inds]
		entropies[i] = numpy.sum(-1*pops_per_cat_temp/total_pops*(nonzero_p_hit*numpy.log2(nonzero_p_hit) + nonzero_p_no_hit*numpy.log2(nonzero_p_no_hit)))
	
	return numpy.sum(entropies)

def calculate_entropy_statistic2(mutation_table, helper_matrix):
	#This function can be used to weight double-hit mutations less than other mutations, since they carry less information.
	#However, for this dataset including the 2-hit mutations with equal weight was equivalently sensitive.
	muts_per_gene = numpy.sum(mutation_table, axis = 0)
	
	collapsed_table = numpy.dot(helper_matrix,mutation_table)
	pops_per_category = numpy.dot(helper_matrix,helper_matrix.T)
	probs = numpy.dot(numpy.linalg.inv(pops_per_category),collapsed_table) #probability that a population in this category got a mutation
	#print probs
	num_genes = mutation_table.shape[1]
	
	entropies = numpy.zeros((num_genes,))
	weight = 1.
	total_pops = numpy.float(numpy.sum(pops_per_category))
	
	for i in range(num_genes):
		if muts_per_gene[i] > 2.1:
			nonzero_inds = numpy.all([probs[:,i] > 0 , probs[:,i]< 1], axis = 0)
			nonzero_p_hit = probs[:,i][nonzero_inds]
			nonzero_p_no_hit = 1. - nonzero_p_hit
			pops_per_cat_temp = numpy.diag(pops_per_category)[nonzero_inds]
			
			entropies[i] = numpy.sum(-1*pops_per_cat_temp/total_pops*(nonzero_p_hit*numpy.log2(nonzero_p_hit) + nonzero_p_no_hit*numpy.log2(nonzero_p_no_hit)))
		else:
			nonzero_inds = numpy.all([probs[:,i] > 0 , probs[:,i]< 1], axis = 0)
			nonzero_p_hit = probs[:,i][nonzero_inds]
			nonzero_p_no_hit = 1. - nonzero_p_hit
			pops_per_cat_temp = numpy.diag(pops_per_category)[nonzero_inds]
			entropies[i] = weight*numpy.sum(-1*pops_per_cat_temp/total_pops*(nonzero_p_hit*numpy.log2(nonzero_p_hit) + nonzero_p_no_hit*numpy.log2(nonzero_p_no_hit)))
	return numpy.sum(entropies)

def calculate_presence_absence_statistic(mutation_table, helper_matrix):
	
	#This is a simple test statistic based on whether or not a particular mutation was or wasn't hit in some category.
	
	collapsed_table = numpy.dot(helper_matrix,mutation_table)
	
	num_zeros = numpy.sum(collapsed_table < .5)
	
	return num_zeros
	
#Read in the list of mutations that fixed in each population. Filter out snps that occur in multiple descendants of the same founder--these were SGV from the passaging of this segregant well.

input_file = 'data/mutation_lists_with_aa_positions_reannotated.txt'

#First loop to find any mutations that are shared among descendants of the same segregant

file = open(input_file,'r')
file_lines = file.readlines()
file.close()

segregant_mut_dict = {}
common_mut_dict = {}
for line in file_lines:
	linelist = line.strip().split('\t')
	if len(linelist) < 1.5:
		#Go to the next clone
		clone_name = linelist[0]
		
		segregant = clone_name.split('_')[0]
		if segregant not in segregant_mut_dict:
			segregant_mut_dict[segregant] = []
		
	else:
		mutation = ('_').join(str(i) for i in linelist)
		if len(linelist) > 5.5:
			if linelist[6] == 'Non':
				if mutation in segregant_mut_dict[segregant]:
					print segregant, mutation
					if segregant in common_mut_dict:
						common_mut_dict[segregant].append(mutation)
					else:
						common_mut_dict[segregant] = [mutation]
				if mutation not in segregant_mut_dict[segregant]:
					
					segregant_mut_dict[segregant].append(mutation)

##Second loop to identify all de novo nonsynonymous mutations (and indels)

gene_dict_by_sample = {}
mutation_dict_by_sample = {}

for line in file_lines:
	
	linelist = line.strip().split('\t')
	
	if len(linelist) < 1.5:
		#Go to the next clone
		clone_name = linelist[0]
		gene_dict_by_sample[clone_name] = []
		mutation_dict_by_sample[clone_name] = []
		local_gene_names = []
		segregant = clone_name.split('_')[0]
		
	else:
		gene_name = linelist[4]
		mutation = ('_').join(str(i) for i in linelist)
		if len(linelist) > 5.5:
			if linelist[6] == 'Non':
				if segregant in common_mut_dict: #There might be shared ancestral snps
				
					if ((gene_name not in local_gene_names) and (len(gene_name) < 6.5) and (mutation not in common_mut_dict[segregant])): #We have not already counted this mutation, it is not an ancestral mutation, and it is not in a dubious ORF
						local_gene_names.append(gene_name)
						gene_dict_by_sample[clone_name].append(gene_name)
						mutation_dict_by_sample[clone_name].append(mutation)
				elif ((gene_name not in local_gene_names) and (len(gene_name) < 6.5)): #We have not already counted this mutation, it is not an ancestral mutation, and it is not in a dubious ORF
					
					local_gene_names.append(gene_name)
					gene_dict_by_sample[clone_name].append(gene_name)
					mutation_dict_by_sample[clone_name].append(mutation)
		


#Determine how many independent times each gene was mutated, and make list of genes for each segregant.
gene_name_counts = []
gene_names = []
samples = sorted(gene_dict_by_sample.keys())

for sample in samples:
	for gene in gene_dict_by_sample[sample]:
		if gene in gene_names:
			index = numpy.where(numpy.array(gene_names) == gene)[0]
			gene_name_counts[index] += 1
		else:
			gene_names.append(gene)
			gene_name_counts.append(1)

gene_name_counts_sc = numpy.zeros_like( numpy.array(gene_name_counts) )
gene_name_counts_ypd = numpy.zeros_like( numpy.array(gene_name_counts) )

for sample in samples:
	env = sample.split('_')[2]
	if env == 'sc':
		for gene in gene_dict_by_sample[sample]:
			index = numpy.where(numpy.array(gene_names) == gene)[0]
			gene_name_counts_sc[index] += 1
	elif env == 'ypd':
		for gene in gene_dict_by_sample[sample]:
			index = numpy.where(numpy.array(gene_names) == gene)[0]
			gene_name_counts_ypd[index] += 1

#print gene_name_counts_sc
#print gene_name_counts_ypd

##Import fitness and founder genotype data


#Import fitness and genotype data
filename1 = 'data/fitness_measurements_with_population_names_12_29_2016.csv'
filename2 = 'data/control_replicate_measurements.csv'
filename3 = 'data/segregant_genotypes_deduplicated_with_header.csv'

segregant_vector = []
init_fits_ypd = []
init_std_errs_ypd = []
init_fits_sc = []
init_std_errs_sc = []

final_fits_ypd_pops_in_ypd = []
segregant_vector_ypd_pops = []
final_fits_sc_pops_in_sc = []
segregant_vector_sc_pops = []

final_fits_sc_pops_in_ypd = []
final_fits_ypd_pops_in_sc = []

file1 = open(filename1,'r')
firstline = 0
for line in file1:
	if firstline < .5:
		firstline += 1
		continue
	linestrs = line.strip().split(';')
	segregant_vector.append(linestrs[0])
	init_fits_ypd.append(float(linestrs[1]))
	init_std_errs_ypd.append(float(linestrs[2]))
	init_fits_sc.append(float(linestrs[3]))
	init_std_errs_sc.append(float(linestrs[4]))
	
	ypd_evolved_pops = linestrs[5].split(',')
	for entry in ypd_evolved_pops:
		segregant_vector_ypd_pops.append(linestrs[0])
		final_fits_ypd_pops_in_ypd.append(float(entry.split()[1]))
		final_fits_ypd_pops_in_sc.append(float(entry.split()[2]))
	sc_evolved_pops = linestrs[6].split(',')
	for entry in sc_evolved_pops:
		segregant_vector_sc_pops.append(linestrs[0])
		final_fits_sc_pops_in_ypd.append(float(entry.split()[1]))
		final_fits_sc_pops_in_sc.append(float(entry.split()[2]))

file1.close()

init_fits_ypd = numpy.array(init_fits_ypd)
init_std_errs_ypd = numpy.array(init_std_errs_ypd)
init_fits_sc = numpy.array(init_fits_sc)
init_std_errs_sc = numpy.array(init_std_errs_sc)

final_fits_ypd_pops_in_ypd = numpy.array(final_fits_ypd_pops_in_ypd)
final_fits_ypd_pops_in_sc = numpy.array(final_fits_ypd_pops_in_sc)
segregant_vector_ypd_pops = numpy.array(segregant_vector_ypd_pops)
segregant_vector_sc_pops = numpy.array(segregant_vector_sc_pops)

final_fits_sc_pops_in_ypd = numpy.array(final_fits_sc_pops_in_ypd)
final_fits_ypd_pops_in_sc = numpy.array(final_fits_ypd_pops_in_sc)


ypd_controls = {}
sc_controls = {}

file2 = open(filename2,'r')

firstline = 0
for line in file2:
	if firstline < .5:
		firstline += 1
		continue
	linestrs = line.strip().split(';')
	ypd_controls[linestrs[0]] = [float(i) for i in linestrs[1].split(',')]
	sc_controls[linestrs[0]] = [float(i) for i in linestrs[2].split(',')]

file2.close()
genotype_mat = []
file3 = open(filename3,'r')
firstline = 0
for line in file3:
	if firstline < .5:
		firstline += 1
		continue
	linelist = line.strip().split(';')
	genotype = [int(i) for i in linelist[1].split(',')]
	genotype_mat.append(genotype)

genotype_mat = numpy.array(genotype_mat)
rm_allele = numpy.array(genotype_mat[:,3777],dtype='Bool')
by_allele = numpy.array(1 - genotype_mat[:,3777],dtype='Bool')
				
####Set up mutation table

gene_names = numpy.array(gene_names)
double_hit_genes = gene_names[numpy.array(gene_name_counts) > 1.5]
#print seg_samples
#print samples
num_double_hit_genes = len(double_hit_genes)

#print double_hit_genes

gene_names_reordered = ['KRE33', 'ENP2', 'BFR2', 'BMS1', 'UTP20', 'RPS8A', 'RPS6A','CRM1', 'ECM16', 'BUD23', 'IRA1', 'IRA2', 'GPB1', 'GPB2', 'PDE2','SIR2', 'SIR3', 'SIR4', 'RXT3', 'NNK1',  'YPK9', 'LTE1', 'SRS2','PAR32','STE11','RRP3','RQC2']

#print set(double_hit_genes) == set(gene_names_reordered)

new_gene_order = []
for i in range(len(gene_names_reordered)):
	index = numpy.where(numpy.array(double_hit_genes) == gene_names_reordered[i])[0][0]
	new_gene_order.append(index)

mutation_table = numpy.zeros((254,num_double_hit_genes))
indel_table = numpy.zeros((254,num_double_hit_genes))
for i in range(len(samples)):
	for j in range(num_double_hit_genes):
		if double_hit_genes[new_gene_order[j]] in gene_dict_by_sample[samples[i]]:
			mutation_table[i,j] = 1
			gene_ind = gene_dict_by_sample[samples[i]].index(double_hit_genes[new_gene_order[j]])
			mutation = mutation_dict_by_sample[samples[i]][gene_ind]
			
			mutation_list = mutation.split('_')
			#print mutation_list
			if (len(mutation_list[3].split(':')) > 1.5 or 'Stop' in mutation_list[-1]): #indels and premature stops
				indel_table[i,j] = 1
				mutation_table[i,j] -= 1

###Determine genotype of sequenced populations

genotype_mat_sequenced_populations = numpy.zeros((len(samples),len(genotype_mat[0,:])))
i = 0
seg_samples = []
for clone in samples:
	
	name_strs = clone.split('_')
	seg = name_strs[0]
	seg_samples.append(seg)
	genotype_index = segregant_vector.index(seg)
	genotype_mat_sequenced_populations[i,:] = genotype_mat[genotype_index,:]
	i += 1

###Set up an indicator variable for the environment
	
env_list = []
for sample in samples:
	env = sample.split('_')[2]
	if env=='sc':
		env_list.append(1)
	elif env=='ypd':
		env_list.append(0)

env_list = numpy.array(env_list, dtype='Bool')

##Set up an indicator variable for the Kre33 allele

kre33_allele = genotype_mat_sequenced_populations[:, 9596]
kre33_allele = numpy.array(kre33_allele, dtype='Bool')

##Set up a variable for the segregant
seg_counter = 0
prev_seg = seg_samples[0]
founder_num_key = []
for i in range(len(samples)):
	seg = seg_samples[i]
	if seg == prev_seg:
		founder_num_key.append(seg_counter)
	else:
		seg_counter += 1
		prev_seg = seg
		founder_num_key.append(seg_counter)
founder_num_key = numpy.array(founder_num_key)

###Determine mutations per gene for 4 categories: Kre33-RM/30C; Kre33-BY/30C; Kre33-RM/37C; Kre33-BY/37C

group4 = numpy.array(env_list*kre33_allele, dtype='Bool')
group3 = numpy.array(env_list*(1 - kre33_allele), dtype='Bool')
group2 = numpy.array((1 - env_list)*kre33_allele, dtype='Bool')
group1 = numpy.array((1 - env_list)*(1 - kre33_allele), dtype='Bool')

counts_grp1_mutations = numpy.sum(mutation_table[group1, :], axis=0)
counts_grp1_indels = numpy.sum(indel_table[group1,:], axis=0)

counts_grp2_mutations = numpy.sum(mutation_table[group2, :], axis=0)
counts_grp2_indels = numpy.sum(indel_table[group2,:], axis=0)

counts_grp3_mutations = numpy.sum(mutation_table[group3, :], axis=0)
counts_grp3_indels = numpy.sum(indel_table[group3,:], axis=0)

counts_grp4_mutations = numpy.sum(mutation_table[group4, :], axis=0)
counts_grp4_indels = numpy.sum(indel_table[group4,:], axis=0)

###Basic counting

num_nonsyn_muts_sc = numpy.sum(gene_name_counts_sc)
num_nonsyn_muts_ypd = numpy.sum(gene_name_counts_ypd)

#print 'num_nonsyn_muts_sc', num_nonsyn_muts_sc, 'num_nonsyn_muts_ypd', num_nonsyn_muts_ypd

num_kre33_ass_muts_sc = numpy.sum( counts_grp1_mutations[0:10] ) + numpy.sum( counts_grp2_mutations[0:10] )
num_kre33_ass_muts_ypd = numpy.sum( counts_grp3_mutations[0:10] ) + numpy.sum( counts_grp4_mutations[0:10] )

frac_kre33_ass_muts_sc = num_kre33_ass_muts_sc/float(num_nonsyn_muts_sc)
frac_kre33_ass_muts_ypd = num_kre33_ass_muts_ypd/float(num_nonsyn_muts_ypd)

#print 'kre33_muts_sc', num_kre33_ass_muts_sc
#print 'kre33_muts_ypd', num_kre33_ass_muts_ypd

#print 'kre33 frac muts sc', frac_kre33_ass_muts_sc
#print 'kre33 frac muts ypd', frac_kre33_ass_muts_ypd

###Basic counting per population

num_pops_kre33_mut = numpy.sum( mutation_table[:, 0] > .5 )
frac_pops_kre33_mut = num_pops_kre33_mut/float(mutation_table.shape[0])

#print num_pops_kre33_mut, frac_pops_kre33_mut



####
iter = 10000
mut_table = mutation_table + indel_table
categories_kre33 = numpy.array(kre33_allele) 

categories_null1 = numpy.zeros((len(categories_kre33),))


cat_inds_kre33 = calculate_cat_inds(categories_kre33)
cat_inds_null1 = calculate_cat_inds(categories_null1)

helper_mat_kre33 = calculate_helper_matrix(categories_kre33, cat_inds_kre33)
helper_mat_null1 = calculate_helper_matrix(categories_null1, cat_inds_null1)
mi_stat_kre33 = calculate_entropy_statistic(mut_table, helper_mat_null1) - calculate_entropy_statistic(mut_table, helper_mat_kre33)

mi_stat_permutations = []
for i in range(iter):
	
	permuted_mut_table = permute_within_categories_preserve_num_muts(mut_table, categories_null1, cat_inds_null1)
	mi_stat_permutations.append(calculate_entropy_statistic(permuted_mut_table, helper_mat_null1) - calculate_entropy_statistic(permuted_mut_table, helper_mat_kre33))

print 'Kre33 allele effect'
print 'p =', numpy.sum(mi_stat_permutations > mi_stat_kre33)/float(iter), '(10,000 permutations)'
print 'mi, treatment', mi_stat_kre33
print 'mi, mean of null', numpy.mean(mi_stat_permutations)
print 'difference', mi_stat_kre33 - numpy.mean(mi_stat_permutations), mi_stat_kre33 - numpy.percentile(mi_stat_permutations, 2.5), mi_stat_kre33 - numpy.percentile(mi_stat_permutations, 97.5)
#Control for the kre33 allele and look for additional effects of the environment

categories_env = env_list*10 + categories_kre33

cat_inds_env = calculate_cat_inds(categories_env)

helper_mat_env = calculate_helper_matrix(categories_env, cat_inds_env)

mi_stat_env = calculate_entropy_statistic(mut_table, helper_mat_kre33) - calculate_entropy_statistic(mut_table, helper_mat_env)

mi_stat_permutations_null2 = []
for i in range(iter):
	permuted_mut_table = permute_within_categories_preserve_num_muts(mut_table, categories_kre33, cat_inds_kre33)
	mi_stat_permutations_null2.append(calculate_entropy_statistic(permuted_mut_table, helper_mat_kre33) - calculate_entropy_statistic(permuted_mut_table, helper_mat_env))

#print permuted_order

print 'Environment effect'
print 'p =', numpy.sum(mi_stat_permutations_null2 > mi_stat_env)/float(iter)
print 'mi, treatment', mi_stat_env
print 'mi, mean of null', numpy.mean(mi_stat_permutations_null2)
print 'difference', mi_stat_env - numpy.mean(mi_stat_permutations_null2), mi_stat_env - numpy.percentile(mi_stat_permutations_null2, 2.5), mi_stat_env - numpy.percentile(mi_stat_permutations_null2, 97.5)
#Control for kre33 and environment and look for additional effects of genotype

categories_genotype = founder_num_key*100 + categories_env #This will sum over the environment groups within a founder class

cat_inds_genotype = calculate_cat_inds(categories_genotype)
helper_mat_genotype = calculate_helper_matrix(categories_genotype, cat_inds_genotype)

mi_stat_genotype = calculate_entropy_statistic(mut_table, helper_mat_env) - calculate_entropy_statistic(mut_table, helper_mat_genotype)

mi_stat_permutations_null4 = []

for i in range(iter):
	permuted_mut_table = permute_within_categories_preserve_num_muts(mut_table, categories_env, cat_inds_env)
	mi_stat_permutations_null4.append(calculate_entropy_statistic(permuted_mut_table, helper_mat_env) - calculate_entropy_statistic(permuted_mut_table, helper_mat_genotype))


print 'Genotype effect'
print 'p =', numpy.sum(mi_stat_permutations_null4 > mi_stat_genotype)/float(iter)
print 'mi, treatment', mi_stat_genotype
print 'mi, mean of null', numpy.mean(mi_stat_permutations_null4)
print 'difference', mi_stat_genotype - numpy.mean(mi_stat_permutations_null4), mi_stat_genotype - numpy.percentile(mi_stat_permutations_null4, 2.5), mi_stat_genotype - numpy.percentile(mi_stat_permutations_null4, 97.5)

#For reference, calculate the overall entropy
hits_per_gene = numpy.sum(mut_table, axis=0)
pi = hits_per_gene/254.
print 'Total entropy'
print numpy.sum(-1*(pi*numpy.log2(pi) + (1.- pi)*numpy.log2(1.-pi)))

f, (ax1, ax2, ax3) = pt.subplots(1, 3, figsize = (16,6))
ax1.set_title('Kre33 allele effect')
ax1.set_ylabel('Probability')
ax1.set_xlabel('Mutual information')
ax1.hist(mi_stat_permutations, color='grey',normed=True,bins=20,edgecolor="none")
ax1.axvline(mi_stat_kre33,0,.2,color='blue')
ymax = ax1.get_ylim()[1]
ax1.plot(mi_stat_kre33,.2*ymax,marker='o',color='b',markersize=5,markeredgewidth=0)


ax2.set_title('Environment effect')
ax2.set_ylabel('Probability')
ax2.set_xlabel('Mutual information')
ax2.hist(mi_stat_permutations_null2, color='grey',normed=True,bins=20,edgecolor="none")

ymax = ax2.get_ylim()[1]
ax2.axvline(mi_stat_env,0,.2,color='blue')
ax2.plot(mi_stat_env,.2*ymax,marker='o',color='b',markersize=5,markeredgewidth=0)

ax3.set_title('Genotype effect')
ax3.set_ylabel('Probability')
ax3.set_xlabel('Mutual information')
ax3.hist(mi_stat_permutations_null4, color='grey',normed=True,bins=20,edgecolor="none")
ax3.axvline(mi_stat_genotype,0,.2,color='blue')

ymax = ax3.get_ylim()[1]
ax3.plot(mi_stat_genotype,.2*ymax,marker='o',color='b',markersize=5,markeredgewidth=0)
pt.savefig('Mutual_information_detection.pdf',bbox_inches='tight')