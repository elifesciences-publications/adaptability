import numpy
import sys
import matplotlib.pylab as pt
import matplotlib.cm
import numpy.random
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
import scipy.stats

matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.size'] = 10.0
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['lines.markersize'] = 3
matplotlib.rcParams['lines.linewidth'] = 1
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

##Import fitness and founder genotype data

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
clone_vector_ypd_pops = []
final_fits_sc_pops_in_sc = []
segregant_vector_sc_pops = []
clone_vector_sc_pops = []

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
		templist = entry.split()
		segregant_vector_ypd_pops.append(linestrs[0])
		clone_vector_ypd_pops.append(templist[0])
		final_fits_ypd_pops_in_ypd.append(float(templist[1]))
		final_fits_ypd_pops_in_sc.append(float(templist[2]))
		
	sc_evolved_pops = linestrs[6].split(',')
	for entry in sc_evolved_pops:
		templist = entry.split()
		segregant_vector_sc_pops.append(linestrs[0])
		clone_vector_sc_pops.append(templist[0])
		final_fits_sc_pops_in_ypd.append(float(templist[1]))
		final_fits_sc_pops_in_sc.append(float(templist[2]))

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

#Use controls (~8 technical replicates each of 24 final populations) to estimate the error variance on the final fitness measurements in each environment
n_control_pops = 24.
var_sum = 0
n_total_reps = 0
for pop in ypd_controls:
	fit = numpy.mean(ypd_controls[pop])
	var_sum += numpy.sum((ypd_controls[pop] - fit)**2)
	n_total_reps += len(ypd_controls[pop])
	
measurement_error_var_sc = var_sum/float(n_total_reps - n_control_pops)

var_sum = 0
n_total_reps = 0
for pop in sc_controls:
	fit = numpy.mean(sc_controls[pop])
	var_sum += numpy.sum((sc_controls[pop] - fit)**2)
	n_total_reps += len(sc_controls[pop])

measurement_error_var_ypd = var_sum/float(n_total_reps - n_control_pops)

###

#Set up 'helper matrix' utilities to conveniently calculate averages and variances over segregant groups
num_segs = len(segregant_vector)
num_pops_ypd = len(segregant_vector_ypd_pops)
num_pops_sc = len(segregant_vector_sc_pops)

helper_matrix_ypd_pops = numpy.zeros((num_pops_ypd,num_segs))
helper_matrix_sc_pops = numpy.zeros((num_pops_sc,num_segs))

for i in range(num_segs):
	current_seg = segregant_vector[i]
	helper_matrix_ypd_pops[numpy.where(segregant_vector_ypd_pops == current_seg)[0],i] = 1.
	helper_matrix_sc_pops[numpy.where(segregant_vector_sc_pops == current_seg)[0],i] = 1.

pops_per_seg_ypd = numpy.diag(numpy.dot(helper_matrix_ypd_pops.T,helper_matrix_ypd_pops))
pops_per_seg_sc = numpy.diag(numpy.dot(helper_matrix_sc_pops.T,helper_matrix_sc_pops))

rm_allele_pops_sc = numpy.array(numpy.dot(helper_matrix_sc_pops, rm_allele),dtype='Bool')
by_allele_pops_sc = numpy.array(numpy.dot(helper_matrix_sc_pops, by_allele),dtype='Bool')

rm_allele_pops_ypd = numpy.array(numpy.dot(helper_matrix_ypd_pops, rm_allele),dtype='Bool')
by_allele_pops_ypd = numpy.array(numpy.dot(helper_matrix_ypd_pops, by_allele),dtype='Bool')


# #Use the helper matrix to average among populations descended from a particular segregant:

delta_fits_ypd = final_fits_ypd_pops_in_ypd - numpy.dot(helper_matrix_ypd_pops,init_fits_ypd)
delta_fits_sc = final_fits_sc_pops_in_sc - numpy.dot(helper_matrix_sc_pops,init_fits_sc)

delta_fits_ypd_in_sc = final_fits_ypd_pops_in_sc - numpy.dot(helper_matrix_ypd_pops,init_fits_sc)
delta_fits_sc_in_ypd = final_fits_sc_pops_in_ypd - numpy.dot(helper_matrix_sc_pops,init_fits_ypd)

delta_fits_ypd_means = numpy.dot(delta_fits_ypd,helper_matrix_ypd_pops)/pops_per_seg_ypd
delta_fits_sc_means = numpy.dot(delta_fits_sc,helper_matrix_sc_pops)/pops_per_seg_sc
delta_fits_sc_in_ypd_means = numpy.dot(delta_fits_sc_in_ypd,helper_matrix_sc_pops)/pops_per_seg_sc
delta_fits_ypd_in_sc_means = numpy.dot(delta_fits_ypd_in_sc,helper_matrix_ypd_pops)/pops_per_seg_ypd

#Delta fits inherit variance from the initial fitness and final fitness measurements
delta_fits_sc_vars = numpy.dot((delta_fits_sc - numpy.dot(helper_matrix_sc_pops,delta_fits_sc_means))**2, helper_matrix_sc_pops)/(pops_per_seg_sc - 1.) + init_std_errs_sc**2
delta_fits_sc_std_errs = numpy.sqrt(delta_fits_sc_vars/pops_per_seg_sc)

delta_fits_ypd_vars = numpy.dot((delta_fits_ypd - numpy.dot(helper_matrix_ypd_pops,delta_fits_ypd_means))**2, helper_matrix_ypd_pops)/(pops_per_seg_ypd - 1.) + init_std_errs_ypd**2
delta_fits_ypd_std_errs = numpy.sqrt(delta_fits_ypd_vars/pops_per_seg_ypd)

delta_fits_ypd_in_sc_vars = numpy.dot((delta_fits_ypd_in_sc - numpy.dot(helper_matrix_ypd_pops,delta_fits_ypd_in_sc_means))**2, helper_matrix_ypd_pops)/(pops_per_seg_ypd - 1.) + init_std_errs_sc**2
delta_fits_ypd_in_sc_std_errs = numpy.sqrt(delta_fits_ypd_in_sc_vars/pops_per_seg_ypd)

delta_fits_sc_in_ypd_vars = numpy.dot((delta_fits_sc_in_ypd - numpy.dot(helper_matrix_sc_pops,delta_fits_sc_in_ypd_means))**2, helper_matrix_sc_pops)/(pops_per_seg_sc - 1.) + init_std_errs_ypd**2 #- measurement_error_ypd
delta_fits_sc_in_ypd_std_errs = numpy.sqrt(delta_fits_sc_in_ypd_vars/pops_per_seg_sc)
				
####First calculation: number of nonsynonymous, genic mutations vs. fitness in each evolution condition
clone_vector_sc_pops = numpy.array(clone_vector_sc_pops)
clone_vector_ypd_pops = numpy.array(clone_vector_ypd_pops)


clones_sc = [ i + '_' + j for [i,j] in numpy.stack((segregant_vector_sc_pops, clone_vector_sc_pops)).T]
clones_ypd = [i + '_' + j for [i,j] in numpy.stack((segregant_vector_ypd_pops, clone_vector_ypd_pops)).T]

num_muts_fit_array_sc = []
num_muts_fit_array_ypd = []
seg_list_sc_seq = []
seg_list_ypd_seq = []

num_muts_seg_dict_sc = {}
num_muts_seg_dict_ypd = {}

for clone in gene_dict_by_sample:
	name_strs = clone.split('_')
	seg = name_strs[0]
	evol_env = name_strs[2]
	clone_num = name_strs[1]
	seg_index = segregant_vector.index(seg)
	full_name = seg + '_' + clone_num
	
	if seg not in num_muts_seg_dict_sc:
		num_muts_seg_dict_sc[seg] = {}
		num_muts_seg_dict_sc[seg]['init_fit'] = init_fits_sc[seg_index]
		num_muts_seg_dict_sc[seg]['kre'] = rm_allele[seg_index]
		num_muts_seg_dict_sc[seg]['num_muts'] = []
	if seg not in num_muts_seg_dict_ypd:
		num_muts_seg_dict_ypd[seg] = {}
		num_muts_seg_dict_ypd[seg]['init_fit'] = init_fits_ypd[seg_index]
		num_muts_seg_dict_ypd[seg]['kre'] = rm_allele[seg_index]
		num_muts_seg_dict_ypd[seg]['num_muts'] = []
	if (evol_env == 'sc' and full_name in clones_sc):
		index_sc = clones_sc.index(full_name)
		seg_list_sc_seq.append(seg)
		
	elif (evol_env == 'ypd' and full_name in clones_ypd):
		index_ypd = clones_ypd.index(full_name)
		seg_list_ypd_seq.append(seg)

	num_muts = len(gene_dict_by_sample[clone])
	
	if (evol_env == 'sc' and full_name in clones_sc):
		num_muts_seg_dict_sc[seg]['num_muts'].append(num_muts)
		num_muts_fit_array_sc.append([init_fits_sc[seg_index], delta_fits_sc[index_sc], num_muts, init_std_errs_sc[seg_index], rm_allele[seg_index]])
	elif (evol_env == 'ypd' and full_name in clones_ypd):
		num_muts_seg_dict_ypd[seg]['num_muts'].append(num_muts)
		num_muts_fit_array_ypd.append([init_fits_ypd[seg_index], delta_fits_ypd[index_ypd], num_muts, init_std_errs_ypd[seg_index], rm_allele[seg_index]])
	
##To plot just the sequenced segregants:

num_muts_fit_array_sc = numpy.array(num_muts_fit_array_sc)
num_muts_fit_array_ypd = numpy.array(num_muts_fit_array_ypd)

msc, bsc = numpy.polyfit(num_muts_fit_array_sc[:,0], num_muts_fit_array_sc[:,2], 1)
mypd, bypd = numpy.polyfit(num_muts_fit_array_ypd[:,0], num_muts_fit_array_ypd[:,2], 1)

r_sc, p_sc = scipy.stats.pearsonr(num_muts_fit_array_sc[:,0], num_muts_fit_array_sc[:,2])
r_ypd, p_ypd = scipy.stats.pearsonr(num_muts_fit_array_ypd[:,0], num_muts_fit_array_ypd[:,2])

print p_sc, p_ypd
###

fig, (ax1, ax2) = pt.subplots(1,2,figsize=(8,4))

colors = ['brown','Tomato']
for seg in num_muts_seg_dict_sc:
	init_fit = num_muts_seg_dict_sc[seg]['init_fit']
	num_muts = num_muts_seg_dict_sc[seg]['num_muts']
	kre_status = num_muts_seg_dict_sc[seg]['kre']
	color1 = colors[kre_status]
	
	offset = 0
	for num in sorted(num_muts):
		ax2.plot(init_fit, num + offset, 'o', color=color1, alpha=.9, markeredgewidth=0)
		offset += .08
	ax2.plot([init_fit,init_fit], [min(num_muts), max(num_muts) + offset - .08], color=color1,alpha=.9,linewidth=.5)


ax2.set_ylim(-.5,13)

ax2.set_xlim(-.2,.18)


ax2.set_xlabel('Initial fitness, 37 C (%)')
ax2.set_ylabel('Number of nonsynon. muts')

#ax2.set_frame_on(False)
ax2.set_axisbelow(True)
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()
ax2.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are on
    top='off',         # ticks along the top edge are off
    labelbottom='on')
xmin, xmax = ax2.get_xaxis().get_view_interval()
ymin, ymax = ax2.get_yaxis().get_view_interval()
#ax2.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=2))
#ax2.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=2))


ax2.set_xticks(numpy.arange(-.2,.19,.05))
ax2.plot(numpy.arange(-.19,.18,.01), msc*numpy.arange(-.19,.18,.01) + bsc, 'k')
ax2.set_xticklabels(numpy.arange(-20,19,5))
ax2.text(.06,10,'$r^2=$' + str(round(r_sc**2,2)), fontsize=12)
ax2.set_title('Evolved at 37 C')
colors=['DarkSlateBlue','MediumSlateBlue']

for seg in num_muts_seg_dict_ypd:
	init_fit = num_muts_seg_dict_ypd[seg]['init_fit']
	num_muts = num_muts_seg_dict_ypd[seg]['num_muts']
	kre_status = num_muts_seg_dict_ypd[seg]['kre']
	color1 = colors[kre_status]
	
	offset = 0
	for num in sorted(num_muts):
		ax1.plot(init_fit, num + offset, 'o', color=color1, alpha=.9, markeredgewidth=0)
		offset += .08
	ax1.plot([init_fit,init_fit], [min(num_muts), max(num_muts) + offset - .08], color=color1,alpha=.9,linewidth=.5)


ax1.set_ylim(-.25,6)
ax1.set_xlim(-.15,.1)

#ax1.set_frame_on(False)
ax1.set_axisbelow(True)
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
ax1.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are on
    top='off',         # ticks along the top edge are off
    labelbottom='on')
xmin, xmax = ax1.get_xaxis().get_view_interval()
ymin, ymax = ax1.get_yaxis().get_view_interval()
#ax1.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=2))
#ax1.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=2))


ax1.set_xticks(numpy.arange(-.15,.11,.05))
ax1.plot(numpy.arange(-.13,.1,.01), mypd*numpy.arange(-.13,.1,.01) + bypd, 'k')
ax1.set_xticklabels(numpy.arange(-15,11,5))
ax1.set_xlabel('Initial fitness, 30 C (%)')
ax1.set_ylabel('Number of nonsynon. muts')
ax1.set_title('Evolved at 30 C')
ax1.text(.025,4.5,'$r^2=$' + str(round(r_ypd**2,2)), fontsize=12)
pt.savefig('Init_fit_v_num_mut_1_9_2016.pdf',bbox_inches='tight')		