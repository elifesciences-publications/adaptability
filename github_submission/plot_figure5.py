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

print gene_name_counts_sc
print gene_name_counts_ypd

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

print double_hit_genes

gene_names_reordered = ['KRE33', 'ENP2', 'BFR2', 'BMS1', 'UTP20', 'RPS8A', 'RPS6A','CRM1', 'ECM16', 'BUD23', 'IRA1', 'IRA2', 'GPB1', 'GPB2', 'PDE2','SIR2', 'SIR3', 'SIR4', 'RXT3', 'NNK1',  'YPK9', 'LTE1', 'SRS2','PAR32','STE11','RRP3','RQC2']

print set(double_hit_genes) == set(gene_names_reordered)

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

print counts_grp1_mutations + counts_grp1_indels
print counts_grp2_mutations + counts_grp2_indels
print counts_grp3_mutations + counts_grp3_indels
print counts_grp4_mutations + counts_grp4_indels

print numpy.sum(group1)
print numpy.sum(group2)
print numpy.sum(group3)
print numpy.sum(group4)

###Basic counting

num_nonsyn_muts_sc = numpy.sum(gene_name_counts_sc)
num_nonsyn_muts_ypd = numpy.sum(gene_name_counts_ypd)

num_kre33_ass_muts_sc = numpy.sum( counts_grp1_mutations[0:10] ) + numpy.sum( counts_grp2_mutations[0:10] )
num_kre33_ass_muts_ypd = numpy.sum( counts_grp3_mutations[0:10] ) + numpy.sum( counts_grp4_mutations[0:10] )

frac_kre33_ass_muts_sc = num_kre33_ass_muts_sc/float(num_nonsyn_muts_sc)
frac_kre33_ass_muts_ypd = num_kre33_ass_muts_ypd/float(num_nonsyn_muts_ypd)

print 'kre33_muts_sc', num_kre33_ass_muts_sc
print 'kre33_muts_ypd', num_kre33_ass_muts_ypd

print 'kre33 frac muts sc', frac_kre33_ass_muts_sc
print 'kre33 frac muts ypd', frac_kre33_ass_muts_ypd

###Basic counting per population

num_pops_kre33_mut = numpy.sum( mutation_table[:, 0] > .5 )
frac_pops_kre33_mut = num_pops_kre33_mut/float(mutation_table.shape[0])

print num_pops_kre33_mut, frac_pops_kre33_mut

print numpy.sum( mutation_table )
print numpy.sum( mutation_table[:,0:10] )


####
space = .01
middle_space = .07
buffer = 4*space + middle_space

grp1_width = (1 - buffer)*numpy.sum(group1)/254.
grp2_width = (1 - buffer)*numpy.sum(group2)/254.
grp3_width = (1 - buffer)*numpy.sum(group3)/254.
grp4_width = (1 - buffer)*numpy.sum(group4)/254.

####

bg_alpha = .6
kre_border = num_double_hit_genes + .4 - 10
camp_border = num_double_hit_genes + .4 - 10 - 5
sir_border = num_double_hit_genes + .4 - 10 - 5 - 3
x_max = numpy.max([numpy.max(counts_grp1_mutations + counts_grp1_indels), numpy.max(counts_grp2_mutations + counts_grp2_indels), numpy.max(counts_grp3_mutations + counts_grp3_indels), numpy.max(counts_grp4_mutations + counts_grp4_indels)])



#color1 = [0,158/255.,115/255.]	
#color2 = [230./255.,159/255.,0./255.]
color1 = 'k'
color2 = 'grey'
my_cmap = matplotlib.colors.ListedColormap([[1,1,1],color1,color2],name='my_colormap')

bg_color = 'MediumSlateBlue'

fig = pt.figure(figsize=(8,4))
	
#bounding_rect = [.02, .02, .2, .98]
bounding_rect = [space, .02, grp2_width, .98]
ax = fig.add_axes(bounding_rect,axis_bgcolor= bg_color, alpha= bg_alpha)#, frame_on=False)
#ax.axvline(x_max,color='k',linewidth=2)

for i in range(len(gene_names_reordered)):
	nmut = int(counts_grp2_mutations[::-1][i])
	nindel = int(counts_grp2_indels[::-1][i])
	for j in range(nmut):
		bar1 = ax.barh(left = x_max - j - 1, bottom = i + .5, height=.8, width = .8, color = color1, linewidth=0)
	for j in range(nindel):
		bar2 = ax.barh(left = x_max - nmut - j - 1, bottom = i + .5, height=.8, width = .8, color = color2, linewidth=.5)

ax.set_ylim(.5,num_double_hit_genes + .5)

pt.tick_params(
		axis='x',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom='off',      # ticks along the bottom edge are off
		top='off',		   # ticks along the top edge are on
		labeltop='off',
		direction='in',        
		labelbottom='off',
		length=2)
pt.tick_params(
	axis='y',          # changes apply to the x-axis
	which='both',      # both major and minor ticks are affected
	left='off',      # ticks along the bottom edge are off
	right='off',         # ticks along the top edge are on
	labelleft='off',
	labelright='off')
ax.set_xlim(0,x_max-.2)
ax.axhline(kre_border, color='k')
#for i in numpy.arange(x_max):
#	ax.axvline(i, color=bg_color,alpha=bg_alpha,linewidth=.5)


ax.patch.set_alpha(bg_alpha)
ax.text(5.5, num_double_hit_genes+1, 'Evolved at OT')
ax.axhline(kre_border, color='k',linewidth=.8)
ax.axhline(camp_border, color='k',linewidth=.8)
ax.axhline(sir_border, color='k',linewidth=.8)
ax.legend([bar1, bar2],['missense','nonsense or indel'],loc='lower left')
ax.text(.8, kre_border + 8, '90s preribosomal', rotation=90,fontsize=9)
ax.text(.8, camp_border + 3, 'cAMP', rotation=90,fontsize=9)
ax.text(.8, sir_border + 1.5, 'SIR', rotation=90,fontsize=9)
#ax.text(1, sir_border, 'SIR pathway', rotation=90)

###########

bounding_rect = [grp2_width + grp4_width + 2*space + middle_space, .02, grp1_width, .98]
ax = fig.add_axes(bounding_rect,axis_bgcolor='DarkSlateBlue',alpha=bg_alpha)#,frame_on=False)
#ax.axvline(0,color='k',linewidth=2)

for i in range(len(gene_names_reordered)):
	nmut = int(counts_grp1_mutations[::-1][i])
	nindel = int(counts_grp1_indels[::-1][i])
	for j in range(nmut):
		bar1 = ax.barh(left = j, bottom = i + .5, height=.8, width = .8, color = color1, linewidth=0)
	for j in range(nindel):
		bar2 = ax.barh(left =  nmut + j, bottom = i + .5, height=.8, width = .8, color = color2, linewidth=.5)
		
ax.set_ylim(.5,num_double_hit_genes + .5)

pt.tick_params(
		axis='x',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom='off',      # ticks along the bottom edge are off
		top='off',		   # ticks along the top edge are on
		labeltop='off',
		direction='in',        
		labelbottom='off',
		length=2)
pt.tick_params(
	axis='y',          # changes apply to the x-axis
	which='both',      # both major and minor ticks are affected
	left='off',      # ticks along the bottom edge are off
	right='off',         # ticks along the top edge are on
	labelleft='off',
	labelright='off')
ax.set_xlim(0,x_max)



ax.patch.set_alpha(bg_alpha)
ax.text(5.5, num_double_hit_genes+1, 'Evolved at OT')
ax.axhline(kre_border, color='k',linewidth=.8)
ax.axhline(camp_border, color='k',linewidth=.8)
ax.axhline(sir_border, color='k',linewidth=.8)
########	
bounding_rect = [grp2_width + 2*space, .02, grp4_width, .98]
ax = fig.add_axes(bounding_rect,axis_bgcolor='Tomato',alpha=bg_alpha)#,frame_on=False)
#ax.axvline(x_max,color='k',linewidth=2)

for i in range(len(gene_names_reordered)):
	nmut = int(counts_grp4_mutations[::-1][i])
	nindel = int(counts_grp4_indels[::-1][i])
	for j in range(nmut):
		bar1 = ax.barh(left = x_max - j - 1, bottom = i + .5, height=.8, width = .8, color = color1, linewidth=0)
	for j in range(nindel):
		bar2 = ax.barh(left = x_max - nmut - j - 1, bottom = i + .5, height=.8, width = .8, color = color2, linewidth=.5)
		
ax.set_ylim(.5,num_double_hit_genes + .5)

pt.tick_params(
		axis='x',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom='off',      # ticks along the bottom edge are off
		top='off',		   # ticks along the top edge are on
		labeltop='off',
		direction='in',        
		labelbottom='off',
		length=2)
pt.tick_params(
	axis='y',          # changes apply to the x-axis
	which='both',      # both major and minor ticks are affected
	left='off',      # ticks along the bottom edge are off
	right='off',         # ticks along the top edge are on
	labelleft='off',
	labelright='off')
ax.set_xlim(0,x_max-.2)

ax.patch.set_alpha(bg_alpha)
ax.text(5.5, num_double_hit_genes + 1, 'Evolved at HT')
ax.text(-6, num_double_hit_genes + 2, 'RM KRE33 allele')
ax.axhline(kre_border, color='k')
for i in numpy.arange(num_double_hit_genes):
	gene = gene_names_reordered[num_double_hit_genes-(i+1)]
	ax.text(x_max + 3.25, i+.6, gene, ha='center', fontsize=9)
ax.axhline(kre_border, color='k',linewidth=.8)
ax.axhline(camp_border, color='k',linewidth=.8)
ax.axhline(sir_border, color='k',linewidth=.8)
#######

bounding_rect = [grp2_width+grp4_width+grp1_width+3*space+middle_space, .02, grp1_width, .98]
ax = fig.add_axes(bounding_rect,axis_bgcolor='Brown',alpha=bg_alpha)#,frame_on=False)

#ax.axvline(0,color='k',linewidth=2)
for i in range(len(gene_names_reordered)):
	nmut = int(counts_grp3_mutations[::-1][i])
	nindel = int(counts_grp3_indels[::-1][i])
	for j in range(nmut):
		bar1 = ax.barh(left = j, bottom = i + .5, height=.8, width = .8, color = color1, linewidth=0)
	for j in range(nindel):
		bar2 = ax.barh(left =  nmut + j, bottom = i + .5, height=.8, width = .8, color = color2, linewidth=.5)
ax.set_ylim(.5,num_double_hit_genes + .5)

pt.tick_params(
		axis='x',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom='off',      # ticks along the bottom edge are off
		top='off',		   # ticks along the top edge are on
		labeltop='off',
		direction='in',        
		labelbottom='off',
		length=2)
pt.tick_params(
	axis='y',          # changes apply to the x-axis
	which='both',      # both major and minor ticks are affected
	left='off',      # ticks along the bottom edge are off
	right='off',         # ticks along the top edge are on
	labelleft='off',
	labelright='off')
ax.set_xlim(0,x_max)

ax.set_yticks(numpy.arange(1,num_double_hit_genes + 1,1))

ax.patch.set_alpha(bg_alpha)
ax.text(5.5, num_double_hit_genes+1, 'Evolved at HT')
ax.text(-6, num_double_hit_genes+2, 'BY KRE33 allele')
ax.axhline(kre_border, color='k',linewidth=.8)
ax.axhline(camp_border, color='k',linewidth=.8)
ax.axhline(sir_border, color='k',linewidth=.8)

pt.savefig('mutation_histograms_1_9_2017.pdf',bbox_inches='tight')