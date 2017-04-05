import matplotlib.pylab as pt
import numpy
import matplotlib.cm as cm
import qtl_detection_adaptability
from calculate_narrow_sense_heritability import narrow_sense_hsq

from calculate_narrow_sense_heritability import narrow_sense_hsq_REML

from calculate_narrow_sense_heritability import rsq
from calculate_narrow_sense_heritability import rsq_linear_model
from calculate_narrow_sense_heritability import broad_sense_Hsq
from calculate_narrow_sense_heritability import broad_sense_Hsq_means

from calculate_narrow_sense_heritability import broad_sense_Hsq_init_fits

import scipy.stats
from matplotlib import colors

import matplotlib

matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.size'] = 10.0
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['lines.markersize'] = 3.5
matplotlib.rcParams['lines.linewidth'] = .7
matplotlib.rcParams['legend.fontsize'] = 8.0
matplotlib.rcParams['axes.linewidth']= .5
matplotlib.rcParams['patch.linewidth']= .5

def calculate_rsq(predictor_mat, delta_fit_vector, helper_matrix):
	
	predictor_mat_expanded = numpy.dot(helper_matrix, predictor_mat)
	beta = numpy.dot(numpy.linalg.inv(numpy.dot(predictor_mat_expanded.T, predictor_mat_expanded)), numpy.dot(predictor_mat_expanded.T, delta_fit_vector))
	rsq = scipy.stats.pearsonr(numpy.dot(predictor_mat_expanded, beta),delta_fit_vector)[0]**2
	
	
	return rsq

def difference_statistic(jacknife_values1, jacknife_values2):
	n_segs = 230
	n_delete = 115
	lb = numpy.percentile(jacknife_values1 - jacknife_values2, 5)
	return lb

def jacknife_sigma(jacknife_values):
	n_segs = 230
	n_delete = 115
	sigma = numpy.sqrt((n_segs - n_delete)/n_delete)*numpy.std(jacknife_values) #5th percentile using jacknife-estimated std deviation for the difference
	return sigma
		
##This script compares different models of delta fitness in the away environment, and calculates narrow-sense heritability.
##Modified 2/7/2017 to calculate statistics for mean segregant phenotypes; estimate h^2 using the REML method; and calculate error bars using a jacknife over segregants.

#Import fitness and genotype data
filename1 = 'data/fitness_measurements_with_population_names_12_29_2016.csv'
filename2 = 'data/control_replicate_measurements.csv'
filename3 = 'data/segregant_genotypes_deduplicated.csv'

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
for line in file3:
	linelist = line.strip().split(';')
	genotype = [int(i) for i in linelist[1].split(',')]
	genotype_mat.append(genotype)

genotype_mat = numpy.array(genotype_mat)
rm_allele = numpy.array(genotype_mat[:,3777],dtype='Bool')
by_allele = numpy.array(1 - genotype_mat[:,3777],dtype='Bool')
#Symmetrize genotype_mat

#genotype_mat = genotype_mat - .5
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

#Delta fits inherit variance from the initial fitness and final fitness measurements, in addition to technical measurement error
delta_fits_sc_vars = numpy.dot((delta_fits_sc - numpy.dot(helper_matrix_sc_pops,delta_fits_sc_means))**2, helper_matrix_sc_pops)/(pops_per_seg_sc - 1.)
delta_fits_sc_std_errs = numpy.sqrt(delta_fits_sc_vars/pops_per_seg_sc)

delta_fits_ypd_vars = numpy.dot((delta_fits_ypd - numpy.dot(helper_matrix_ypd_pops,delta_fits_ypd_means))**2, helper_matrix_ypd_pops)/(pops_per_seg_ypd - 1.)
delta_fits_ypd_std_errs = numpy.sqrt(delta_fits_ypd_vars/pops_per_seg_ypd)

delta_fits_ypd_in_sc_vars = numpy.dot((delta_fits_ypd_in_sc - numpy.dot(helper_matrix_ypd_pops,delta_fits_ypd_in_sc_means))**2, helper_matrix_ypd_pops)/(pops_per_seg_ypd - 1.)
delta_fits_ypd_in_sc_std_errs = numpy.sqrt(delta_fits_ypd_in_sc_vars/pops_per_seg_ypd)

delta_fits_sc_in_ypd_vars = numpy.dot((delta_fits_sc_in_ypd - numpy.dot(helper_matrix_sc_pops,delta_fits_sc_in_ypd_means))**2, helper_matrix_sc_pops)/(pops_per_seg_sc - 1.)
delta_fits_sc_in_ypd_std_errs = numpy.sqrt(delta_fits_sc_in_ypd_vars/pops_per_seg_sc)

##QTL models
n_segs = num_segs

file1 = open('data/initial_fitness_qtl_table_1_18_2017.csv','r')
init_fit_qtls_sc = []
init_fit_qtls_ypd = []
env_counter = 0
for line in file1:
	if line.startswith('#'):
		env_counter += 1
	if (not line.startswith('#') and env_counter < 2.5):
		init_fit_qtls_sc.append( int( line.split(',')[0] ) )
		
	elif not line.startswith('#'):
		init_fit_qtls_ypd.append( int( line.split(',')[0] ) )
		
file1.close()

file1 = open('data/delta_fitness_qtl_table_1_18_2017.csv','r')
delta_fit_qtls_sc = []
delta_fit_qtls_ypd = []
env_counter = 0
for line in file1:
	if line.startswith('#'):
		env_counter += 1
	if (not line.startswith('#') and env_counter < 2.5):
		delta_fit_qtls_sc.append( int( line.split(',')[0] ) )
		
	elif not line.startswith('#'):
		delta_fit_qtls_ypd.append( int( line.split(',')[0] ) )
		
file1.close()

###Note that we're going to test all previously identified loci for significant influence on pleiotropy over and above initial fitness, because loci affecting fitness or adaptability in either environment could be involved.

merged_qtl_list = init_fit_qtls_ypd

for qtl in init_fit_qtls_sc:
	if not numpy.min(numpy.abs(numpy.array(merged_qtl_list) - qtl)) < 20 : ##Already have this locus
		merged_qtl_list.append(qtl)

for qtl in delta_fit_qtls_ypd:
	if not numpy.min( numpy.abs(numpy.array(merged_qtl_list) - qtl )) < 20:
		merged_qtl_list.append(qtl)

for qtl in delta_fit_qtls_sc:
	if not numpy.min( numpy.abs(numpy.array(merged_qtl_list) - qtl )) < 20:
		merged_qtl_list.append(qtl)
		
###Cross environment fitness plus qtl model. Note that here 'sc' refers to sc evolved populations; 'ypd' refers to ypd-evolved populations.

rsq_sc_list = []
rsq_sc_interval_list = []

rsq_ypd_list = []
rsq_ypd_interval_list = []

n_segs = len(pops_per_seg_sc)
X_fit_sc_segs = numpy.concatenate((numpy.ones((n_segs,1)),  init_fits_sc.reshape((n_segs,1)), init_fits_ypd.reshape((n_segs,1))), axis = 1)
X_fit_ypd_segs = numpy.concatenate((numpy.ones((n_segs,1)),  init_fits_ypd.reshape((n_segs,1)), init_fits_sc.reshape((n_segs,1))), axis = 1)

X_qtls_sc_linear = X_fit_sc_segs
X_qtls_ypd_linear = X_fit_ypd_segs
for qtl in merged_qtl_list:
	
	X_qtls_sc_linear = numpy.append(X_qtls_sc_linear, genotype_mat[:, qtl].reshape((n_segs,1)), axis=1)
	X_qtls_ypd_linear = numpy.append(X_qtls_ypd_linear, genotype_mat[:, qtl].reshape((n_segs,1)), axis=1)

print merged_qtl_list
print X_qtls_sc_linear
rsq_all_linear_sc = calculate_rsq(X_qtls_sc_linear, delta_fits_sc_in_ypd, helper_matrix = helper_matrix_sc_pops)
rsq_all_linear_ypd = calculate_rsq(X_qtls_ypd_linear, delta_fits_ypd_in_sc, helper_matrix = helper_matrix_ypd_pops)

print rsq_all_linear_sc
print rsq_all_linear_ypd

X_qtls_sc_pops_linear = numpy.dot( helper_matrix_sc_pops, X_qtls_sc_linear)
X_qtls_ypd_pops_linear = numpy.dot( helper_matrix_ypd_pops, X_qtls_ypd_linear)

beta_sc_linear = numpy.dot( numpy.linalg.inv( numpy.dot(X_qtls_sc_pops_linear.T, X_qtls_sc_pops_linear) ), numpy.dot(X_qtls_sc_pops_linear.T, delta_fits_sc_in_ypd) )
beta_ypd_linear = numpy.dot( numpy.linalg.inv( numpy.dot(X_qtls_ypd_pops_linear.T, X_qtls_ypd_pops_linear) ), numpy.dot(X_qtls_ypd_pops_linear.T, delta_fits_ypd_in_sc) )

print beta_sc_linear
print beta_ypd_linear
####Significance of coefficients for linear qtl model

rsq_sc_linear_list = []
rsq_ypd_linear_list = []

F_sc_linear_list = []
F_ypd_linear_list = []
num_free_parameters = X_qtls_sc_linear.shape[1] + 1
##Calculate rsq with each coefficient left out in turn

for i in range(X_qtls_sc_linear.shape[1]):
	X_qtls_sc_segs_temp = numpy.delete(X_qtls_sc_linear, i, axis=1)
	X_qtls_ypd_segs_temp = numpy.delete(X_qtls_ypd_linear, i, axis=1)
	
	rsq_sc_temp = calculate_rsq(X_qtls_sc_segs_temp, delta_fits_sc_in_ypd, helper_matrix = helper_matrix_sc_pops)
	rsq_ypd_temp = calculate_rsq(X_qtls_ypd_segs_temp, delta_fits_ypd_in_sc, helper_matrix = helper_matrix_ypd_pops)

	rsq_sc_linear_list.append(rsq_sc_temp)
	rsq_ypd_linear_list.append(rsq_ypd_temp)
	
	F_sc_linear_list.append((rsq_all_linear_sc - rsq_sc_temp)/rsq_sc_temp*(num_pops_sc - num_free_parameters))
	F_ypd_linear_list.append((rsq_all_linear_ypd - rsq_ypd_temp)/rsq_ypd_temp*(num_pops_ypd - num_free_parameters))

##Bonferroni correction: to be 95% confident of each coefficient being real, we want them to be significant at a p<.004 level.

alpha = 1. - (.95)**(1./float(num_free_parameters))


p_vals_sc = 1. - scipy.stats.f.cdf(F_sc_linear_list, 1, num_pops_sc - num_free_parameters)
p_vals_ypd = 1. - scipy.stats.f.cdf(F_ypd_linear_list, 1, num_pops_ypd - num_free_parameters)

###Identify significant coefficients

coefficient_inds_sc = numpy.arange(num_free_parameters)[p_vals_sc < alpha]
coefficient_inds_ypd = numpy.arange(num_free_parameters)[p_vals_ypd < alpha]

X_fitness_and_qtls_sc = X_qtls_sc_linear[:, coefficient_inds_sc]
X_fitness_and_qtls_ypd = X_qtls_ypd_linear[:, coefficient_inds_ypd]

print coefficient_inds_sc
print coefficient_inds_ypd
print merged_qtl_list



###

##For plotting, rsq for the means for the combined model

X_fitness_and_qtls_sc = X_qtls_sc_linear[:, coefficient_inds_sc]
X_fitness_and_qtls_ypd = X_qtls_ypd_linear[:, coefficient_inds_ypd]
X_fitness_and_qtls_sc_pops = numpy.dot(helper_matrix_sc_pops, X_fitness_and_qtls_sc)
X_fitness_and_qtls_ypd_pops = numpy.dot(helper_matrix_ypd_pops, X_fitness_and_qtls_ypd)
beta_fitness_and_qtls_sc = numpy.dot( numpy.linalg.inv( numpy.dot(X_fitness_and_qtls_sc_pops.T, X_fitness_and_qtls_sc_pops) ), numpy.dot(X_fitness_and_qtls_sc_pops.T, delta_fits_sc_in_ypd) )
beta_fitness_and_qtls_ypd = numpy.dot( numpy.linalg.inv( numpy.dot(X_fitness_and_qtls_ypd_pops.T, X_fitness_and_qtls_ypd_pops) ), numpy.dot(X_fitness_and_qtls_ypd_pops.T, delta_fits_ypd_in_sc) )

###Print in a table

rsq_sc_nested = []
rsq_ypd_nested = []

F_sc_nested = []
F_ypd_nested = []

for i in range(X_fitness_and_qtls_sc.shape[1]):
	num_free_parameters = i + 1
	rsq = calculate_rsq( X_fitness_and_qtls_sc[:,0:i+1].reshape((n_segs,i+1)), delta_fits_sc_in_ypd, helper_matrix = helper_matrix_sc_pops )
	if i > 1.5:
		F = (rsq - rsq_sc_nested[i-1])/rsq_sc_nested[i-1]*(num_pops_sc - num_free_parameters)
	else:
		F = 'NA'
	#print rsq, rsq2
	rsq_sc_nested.append(rsq)
	F_sc_nested.append(F)

for i in range(X_fitness_and_qtls_ypd.shape[1]):
	num_free_parameters = i + 1
	rsq = calculate_rsq(X_fitness_and_qtls_ypd[:,0:i+1].reshape((n_segs,i+1)), delta_fits_ypd_in_sc, helper_matrix = helper_matrix_ypd_pops )
	if i > 1.5:
		F = (rsq - rsq_ypd_nested[i-1])/rsq_ypd_nested[i-1]*(num_pops_ypd - num_free_parameters)
	else:
		F = 'NA'
	rsq_ypd_nested.append(rsq)
	F_ypd_nested.append(F)

filename = 'data/Pleiotropy_model_parameters_2_19_2017.txt'

# # 
file = open(filename,'w')
file.write('#Trait: Delta fitness in YPD at 30C for populations evolved in SC at 37C' + '\n')
file.write((',').join(('#Term (qtl marker location)','term number','Additional variance explained','Model Coefficient', 'F statistic')) + '\n')
file.write((',').join(('intercept',str(0),'NA', str(beta_fitness_and_qtls_sc[0]), str(F_sc_nested[0]))) + '\n')
file.write((',').join(('initial fitness in YPD at 30C',str(1),str(rsq_sc_nested[1] - rsq_sc_nested[0]), str(beta_fitness_and_qtls_sc[1]), str(F_sc_nested[1]))) + '\n')
counter = 2
for coefficient in coefficient_inds_sc[2:]:
	
	locus = merged_qtl_list[int(coefficient) - 3] #First 3 terms in coefficient inds are the intercept, init fit in SC, init fit in YPD
	file.write((',').join(('linear',str(locus),str(coefficient),str(rsq_ypd_nested[counter] - rsq_ypd_nested[counter-1]), str(beta_fitness_and_qtls_sc[counter]), str(F_sc_nested[counter]) )) + '\n')
	counter += 1
		
# # 
file.write('#Trait: Delta fitness in SC at 37C for populations evolved in YPD at 30C' + '\n')
file.write((',').join(('#Term (qtl marker location)','term number','Additional variance explained','Model Coefficient', 'F statistic')) + '\n')
file.write((',').join(('intercept',str(0),'NA', str(beta_fitness_and_qtls_ypd[0]), str(F_ypd_nested[0]))) + '\n')
file.write((',').join(('initial fitness in SC at 37C',str(1),str(rsq_ypd_nested[1] - rsq_ypd_nested[0]), str(beta_fitness_and_qtls_ypd[1]), str(F_ypd_nested[1]))) + '\n')
counter = 2
for coefficient in coefficient_inds_ypd[2:]:
	
	locus = merged_qtl_list[int(coefficient) - 3]
	file.write((',').join(('linear',str(locus),str(coefficient),str(rsq_ypd_nested[counter] - rsq_ypd_nested[counter-1]), str(beta_fitness_and_qtls_ypd[counter]), str(F_ypd_nested[counter]) )) + '\n')
	counter += 1
file.close()

###

predictor_fit_and_qtls_ypd_unscaled = numpy.dot( X_fitness_and_qtls_ypd, beta_fitness_and_qtls_ypd )
predictor_fit_and_qtls_ypd = (predictor_fit_and_qtls_ypd_unscaled - numpy.mean(predictor_fit_and_qtls_ypd_unscaled))/numpy.std(predictor_fit_and_qtls_ypd_unscaled)
predictor_fit_and_qtls_sc_unscaled = numpy.dot( X_fitness_and_qtls_sc, beta_fitness_and_qtls_sc )
predictor_fit_and_qtls_sc = (predictor_fit_and_qtls_sc_unscaled - numpy.mean(predictor_fit_and_qtls_sc_unscaled))/numpy.std(predictor_fit_and_qtls_sc_unscaled)

rsq_full_model_ypd_for_means = scipy.stats.pearsonr(predictor_fit_and_qtls_ypd, delta_fits_ypd_in_sc_means)[0]**2
rsq_full_model_sc_for_means = scipy.stats.pearsonr(predictor_fit_and_qtls_sc, delta_fits_sc_in_ypd_means)[0]**2


####Calculate heritabilities

n_segs_by = numpy.sum(by_allele)
n_segs_rm = numpy.sum(rm_allele)
n_pops_by_sc = numpy.sum(by_allele_pops_sc)
n_pops_rm_sc = numpy.sum(rm_allele_pops_sc)
n_pops_by_ypd = numpy.sum(by_allele_pops_ypd)
n_pops_rm_ypd = numpy.sum(rm_allele_pops_ypd)

#Set up lists of statistics to be calculated each jacknife round

H2_delta_fit_sc_bs = []
H2_delta_fit_ypd_bs = []
h2_delta_fit_sc_bs = []
h2_delta_fit_ypd_bs = []
r2_delta_fit_sc_fit_model_bs = []
r2_delta_fit_ypd_fit_model_bs = []
r2_delta_fit_sc_qtl_model_bs = []
r2_delta_fit_ypd_qtl_model_bs = []
r2_delta_fit_sc_combined_model_bs = []
r2_delta_fit_ypd_combined_model_bs = []


##Calculate statistics on the whole set of segregants

H2_delta_fit_ypd=broad_sense_Hsq(delta_fits_ypd_in_sc, pops_per_seg_ypd, delta_fits_ypd_in_sc_vars, helper_matrix_ypd_pops)
H2_delta_fit_sc=broad_sense_Hsq(delta_fits_sc_in_ypd, pops_per_seg_sc, delta_fits_sc_in_ypd_vars, helper_matrix_sc_pops)

h2_delta_fit_sc=narrow_sense_hsq_REML(genotype_mat, delta_fits_sc_in_ypd, helper_matrix = helper_matrix_sc_pops, variance_comps_initial_guess=[.4,.2])
h2_delta_fit_ypd=narrow_sense_hsq_REML(genotype_mat, delta_fits_ypd_in_sc, helper_matrix = helper_matrix_ypd_pops, variance_comps_initial_guess=[.4,.2])

r2_delta_fit_sc_fit_model=rsq_linear_model(X_fit_sc_segs, delta_fits_sc_in_ypd, helper_matrix = helper_matrix_sc_pops)
r2_delta_fit_ypd_fit_model=rsq_linear_model(X_fit_ypd_segs, delta_fits_ypd_in_sc, helper_matrix = helper_matrix_ypd_pops)
r2_delta_fit_sc_combined_model=rsq_linear_model(X_fitness_and_qtls_sc, delta_fits_sc_in_ypd, helper_matrix = helper_matrix_sc_pops)
r2_delta_fit_ypd_combined_model=rsq_linear_model(X_fitness_and_qtls_ypd, delta_fits_ypd_in_sc, helper_matrix = helper_matrix_ypd_pops)
	
print H2_delta_fit_sc, H2_delta_fit_ypd, h2_delta_fit_sc, h2_delta_fit_ypd
print r2_delta_fit_sc_fit_model, r2_delta_fit_ypd_fit_model, r2_delta_fit_sc_combined_model, r2_delta_fit_ypd_combined_model

###Jacknife over the 230 segregants

n_iter = 1000

for k in range(n_iter):
	n_delete = 115 #half the segregants are chosen on each round
	bootstrap_inds = numpy.random.choice(numpy.arange(n_segs), n_delete, replace=False) #These are the segregants chosen for this round
	expanded_indices_ypd = [] #These are the indices of the YPD evolved populations that descend from these segregants
	expanded_indices_sc = [] #These are the indices of the SC evolved populations that descend from these segregants
	for bootstrap_ind in bootstrap_inds:
		expanded_indices_ypd.extend(numpy.where(helper_matrix_ypd_pops[:,bootstrap_ind] == 1)[0])
		expanded_indices_sc.extend(numpy.where(helper_matrix_sc_pops[:,bootstrap_ind] == 1)[0])
	
	##Construct lists of phenotypes for the set of segregants/populations for this bootstrap round
	
	init_fits_sc_bs = init_fits_sc[bootstrap_inds]
	#print init_fits_sc_bs
	init_fits_ypd_bs = init_fits_ypd[bootstrap_inds]
	
	
	delta_fits_sc_in_ypd_means_bs = delta_fits_sc_in_ypd_means[bootstrap_inds]
	delta_fits_ypd_in_sc_means_bs = delta_fits_ypd_in_sc_means[bootstrap_inds]
	
	delta_fits_sc_in_ypd_bs = delta_fits_sc_in_ypd[expanded_indices_sc]
	delta_fits_ypd_in_sc_bs = delta_fits_ypd_in_sc[expanded_indices_ypd]
	
	pops_per_seg_ypd_bs = pops_per_seg_ypd[bootstrap_inds]
	pops_per_seg_sc_bs = pops_per_seg_sc[bootstrap_inds]
	
	delta_fits_ypd_in_sc_vars_bs = delta_fits_ypd_in_sc_vars[bootstrap_inds]
	delta_fits_sc_in_ypd_vars_bs = delta_fits_sc_in_ypd_vars[bootstrap_inds]
	##Update the dependent variables for the set of chosen segregants
	
	X_fit_sc_segs_bs = X_fit_sc_segs[bootstrap_inds,:]
	X_fit_ypd_segs_bs = X_fit_ypd_segs[bootstrap_inds,:]
	
	X_fitness_and_qtls_sc_bs = X_fitness_and_qtls_sc[bootstrap_inds,:]
	X_fitness_and_qtls_ypd_bs = X_fitness_and_qtls_ypd[bootstrap_inds,:]
	
	##
	
	##Update genotypes and helper matrices for the segregants in this bootstrap
	
	genotype_mat_bs = genotype_mat[bootstrap_inds, :]
	helper_matrix_sc_pops_bs = helper_matrix_sc_pops[expanded_indices_sc,:][:,bootstrap_inds]
	helper_matrix_ypd_pops_bs = helper_matrix_ypd_pops[expanded_indices_ypd,:][:,bootstrap_inds]
	
	##
	#pt.figure()
	#pt.plot(numpy.dot(helper_matrix_sc_pops_bs, X_init_fit_sc[:,0].reshape((n_segs,1))), delta_fits_sc_bs, 'o')
	
	##Calculate the statistics that we're interested in for this bootstrap round
	
	H2_delta_fit_ypd_bs.append(broad_sense_Hsq(delta_fits_ypd_in_sc_bs, pops_per_seg_ypd_bs, delta_fits_ypd_in_sc_vars_bs, helper_matrix_ypd_pops_bs))
	H2_delta_fit_sc_bs.append(broad_sense_Hsq(delta_fits_sc_in_ypd_bs, pops_per_seg_sc_bs, delta_fits_sc_in_ypd_vars_bs, helper_matrix_sc_pops_bs))
	
	h2_delta_fit_sc_bs.append(narrow_sense_hsq_REML(genotype_mat_bs, delta_fits_sc_in_ypd_bs, helper_matrix = helper_matrix_sc_pops_bs, variance_comps_initial_guess=[.4,.2]))
	h2_delta_fit_ypd_bs.append(narrow_sense_hsq_REML(genotype_mat_bs, delta_fits_ypd_in_sc_bs, helper_matrix =  helper_matrix_ypd_pops_bs, variance_comps_initial_guess=[.4,.2]))
	
	r2_delta_fit_sc_fit_model_bs.append(rsq_linear_model(X_fit_sc_segs_bs, delta_fits_sc_in_ypd_bs, helper_matrix =  helper_matrix_sc_pops_bs))
	r2_delta_fit_ypd_fit_model_bs.append(rsq_linear_model(X_fit_ypd_segs_bs, delta_fits_ypd_in_sc_bs, helper_matrix =  helper_matrix_ypd_pops_bs))
	
	#print rsq_linear_model(X_qtls_delta_bs, delta_fits_ypd_bs, helper_matrix = helper_matrix_ypd_pops_bs)
	r2_delta_fit_sc_combined_model_bs.append(rsq_linear_model(X_fitness_and_qtls_sc_bs, delta_fits_sc_in_ypd_bs, helper_matrix =  helper_matrix_sc_pops_bs))
	r2_delta_fit_ypd_combined_model_bs.append(rsq_linear_model(X_fitness_and_qtls_ypd_bs, delta_fits_ypd_in_sc_bs, helper_matrix =  helper_matrix_ypd_pops_bs))

H2_delta_fit_sc_bs = numpy.array(H2_delta_fit_sc_bs)
H2_delta_fit_ypd_bs = numpy.array(H2_delta_fit_ypd_bs)
h2_delta_fit_sc_bs = numpy.array(h2_delta_fit_sc_bs)
h2_delta_fit_ypd_bs = numpy.array(h2_delta_fit_ypd_bs)
r2_delta_fit_sc_fit_model_bs = numpy.array(r2_delta_fit_sc_fit_model_bs)
r2_delta_fit_ypd_fit_model_bs = numpy.array(r2_delta_fit_ypd_fit_model_bs)
r2_delta_fit_sc_combined_model_bs = numpy.array(r2_delta_fit_sc_combined_model_bs)
r2_delta_fit_ypd_combined_model_bs = numpy.array(r2_delta_fit_ypd_combined_model_bs)

###Calculate differences between statistics over bootstrapped samples

mean_diff_h2_r2b_delta_ypd = numpy.nanmedian(r2_delta_fit_ypd_fit_model_bs - h2_delta_fit_ypd_bs)
lb_diff_h2_r2b_delta_ypd = numpy.nanpercentile(r2_delta_fit_ypd_fit_model_bs - h2_delta_fit_ypd_bs, 5)

mean_diff_h2_r2b_delta_sc = numpy.nanmedian(r2_delta_fit_sc_fit_model_bs - h2_delta_fit_sc_bs)
lb_diff_h2_r2b_delta_sc =numpy.nanpercentile(r2_delta_fit_sc_fit_model_bs - h2_delta_fit_sc_bs, 5)

mean_diff_r2b_r2c_delta_ypd = numpy.nanmedian(r2_delta_fit_ypd_combined_model_bs - r2_delta_fit_ypd_fit_model_bs)
lb_diff_r2b_r2c_delta_ypd = numpy.nanpercentile(r2_delta_fit_ypd_combined_model_bs-r2_delta_fit_ypd_fit_model_bs,5)

mean_diff_r2b_r2c_delta_sc = numpy.nanmedian(r2_delta_fit_sc_combined_model_bs - r2_delta_fit_sc_fit_model_bs)
lb_diff_r2b_r2c_delta_sc = numpy.nanpercentile(r2_delta_fit_sc_combined_model_bs - r2_delta_fit_sc_fit_model_bs,5)

###Record everything in a table

filename_out1 = 'data/heritability_statistics_cross_environments_with_differences_2_19_2017b.csv'
file_out1 = open(filename_out1,'w')
file_out1.write('#H^2' + '\n')
file_out1.write((',').join( ('OT populations, delta fitness at HT',str(H2_delta_fit_ypd), str(numpy.percentile(H2_delta_fit_ypd_bs,2.5)),  str(numpy.percentile(H2_delta_fit_ypd_bs,97.5))) ) + '\n')
file_out1.write((',').join( ('HT populations, delta fitness at OT',str(H2_delta_fit_sc), str(numpy.percentile(H2_delta_fit_sc_bs,2.5)),  str(numpy.percentile(H2_delta_fit_sc_bs,97.5))) ) + '\n')
file_out1.write('#h^2' + '\n')
file_out1.write((',').join( ('OT populations, delta fitness at HT',str(h2_delta_fit_ypd), str(numpy.percentile(h2_delta_fit_ypd_bs,2.5)),  str(numpy.percentile(h2_delta_fit_ypd_bs,97.5))) ) + '\n')
file_out1.write((',').join( ('HT populations, delta fitness at OT',str(h2_delta_fit_sc), str(numpy.percentile(h2_delta_fit_sc_bs,2.5)),  str(numpy.percentile(h2_delta_fit_sc_bs,97.5))) ) + '\n')

file_out1.write('#r^2, fitness models' + '\n')
file_out1.write((',').join( ('OT populations, delta fitness at HT',str(r2_delta_fit_ypd_fit_model), str(numpy.percentile(r2_delta_fit_ypd_fit_model_bs,2.5)),  str(numpy.percentile(r2_delta_fit_ypd_fit_model_bs,97.5))) ) + '\n')
file_out1.write((',').join( ('HT populations, delta fitness at OT',str(r2_delta_fit_sc_fit_model), str(numpy.percentile(r2_delta_fit_sc_fit_model_bs,2.5)),  str(numpy.percentile(r2_delta_fit_sc_fit_model_bs,97.5))) ) + '\n')

file_out1.write('#r^2-h^2, fitness models' + '\n')

file_out1.write((',').join( ('OT populations, delta fitness at HT',str(mean_diff_h2_r2b_delta_ypd), str(lb_diff_h2_r2b_delta_ypd) )) + '\n')
file_out1.write((',').join( ('HT populations, delta fitness at OT',str(mean_diff_h2_r2b_delta_sc), str(lb_diff_h2_r2b_delta_sc) )) + '\n')

file_out1.write('#r^2, combined models' + '\n')
file_out1.write((',').join( ('OT populations, delta fitness at HT',str(r2_delta_fit_ypd_combined_model), str(numpy.percentile(r2_delta_fit_ypd_combined_model_bs,2.5)),  str(numpy.percentile(r2_delta_fit_ypd_combined_model_bs,97.5))) ) + '\n')
file_out1.write((',').join( ('HT populations, delta fitness at OT',str(r2_delta_fit_sc_combined_model), str(numpy.percentile(r2_delta_fit_sc_combined_model_bs,2.5)),  str(numpy.percentile(r2_delta_fit_sc_combined_model_bs,97.5))) ) + '\n')

file_out1.write('#r^2 combined-r^2 fitness' + '\n')

file_out1.write((',').join( ('OT populations, delta fitness at HT',str(mean_diff_r2b_r2c_delta_ypd), str(lb_diff_r2b_r2c_delta_ypd) )) + '\n')
file_out1.write((',').join( ('HT populations, delta fitness at OT',str(mean_diff_r2b_r2c_delta_sc), str(lb_diff_r2b_r2c_delta_sc) )) + '\n')

file_out1.close()

###95% confidence intervals

lp = 2.5
up = 97.5

H2_delta_fit_sc_err = [H2_delta_fit_sc - numpy.percentile(H2_delta_fit_sc_bs, lp), numpy.percentile(H2_delta_fit_sc_bs, up) - H2_delta_fit_sc]
H2_delta_fit_ypd_err = [H2_delta_fit_ypd - numpy.percentile(H2_delta_fit_ypd_bs, lp), numpy.percentile(H2_delta_fit_ypd_bs, up) - H2_delta_fit_ypd]
h2_delta_fit_sc_err = [h2_delta_fit_sc - numpy.nanpercentile(h2_delta_fit_sc_bs, lp), numpy.nanpercentile(h2_delta_fit_sc_bs, up) - h2_delta_fit_sc]
h2_delta_fit_ypd_err = [h2_delta_fit_ypd - numpy.nanpercentile(h2_delta_fit_ypd_bs, lp), numpy.nanpercentile(h2_delta_fit_ypd_bs, up) - h2_delta_fit_ypd]
r2_delta_fit_sc_fit_model_err = [r2_delta_fit_sc_fit_model - numpy.percentile(r2_delta_fit_sc_fit_model_bs, lp), numpy.percentile(r2_delta_fit_sc_fit_model_bs, up) - r2_delta_fit_sc_fit_model]
r2_delta_fit_ypd_fit_model_err = [r2_delta_fit_ypd_fit_model - numpy.percentile(r2_delta_fit_ypd_fit_model_bs, lp), numpy.percentile(r2_delta_fit_ypd_fit_model_bs, up) - r2_delta_fit_ypd_fit_model]

r2_delta_fit_sc_combined_model_err = [r2_delta_fit_sc_combined_model - numpy.percentile(r2_delta_fit_sc_combined_model_bs, lp), numpy.percentile(r2_delta_fit_sc_combined_model_bs, up) - r2_delta_fit_sc_combined_model]
r2_delta_fit_ypd_combined_model_err = [r2_delta_fit_ypd_combined_model - numpy.percentile(r2_delta_fit_ypd_combined_model_bs, lp), numpy.percentile(r2_delta_fit_ypd_combined_model_bs, up) - r2_delta_fit_ypd_combined_model]

model_difference_sc = numpy.nanmedian(r2_delta_fit_sc_combined_model_bs - r2_delta_fit_sc_fit_model_bs)
model_difference_ypd = numpy.nanmedian(r2_delta_fit_ypd_combined_model_bs - r2_delta_fit_ypd_fit_model_bs)
model_difference_sc_err = [model_difference_sc - numpy.percentile(r2_delta_fit_sc_combined_model_bs - r2_delta_fit_sc_fit_model_bs, 5), 0]

model_difference_ypd_err = [model_difference_ypd - numpy.percentile(r2_delta_fit_ypd_combined_model_bs - r2_delta_fit_ypd_fit_model_bs, 5), 0]

#####Plotting

##Colorbrewer2 pastel set 2

color1b = numpy.array((166,206,227),dtype='float')/255.
color2b = numpy.array((31,120,180),dtype='float')/255.
color3b = numpy.array((178,223,138),dtype='float')/255.
color4b = numpy.array((51,160,44),dtype='float')/255.
color5b = numpy.array((251,154,153),dtype='float')/255.

###Colorbrewer2 set 3

color1c = numpy.array((27,158,119),dtype='float')/255.
color2c = numpy.array((217,95,2),dtype='float')/255.
color3c = numpy.array((117,112,179),dtype='float')/255.
color4c = numpy.array((231,41,138),dtype='float')/255.
color5c = numpy.array((102,166,30),dtype='float')/255.

##Colorbrewer set 4
color1d = numpy.array((141,211,199),dtype='float')/255.
color2d = numpy.array((255,255,179),dtype='float')/255.
color3d = numpy.array((190,186,218),dtype='float')/255.
color4d = numpy.array((251,128,114),dtype='float')/255.
color5d = numpy.array((128,177,211),dtype='float')/255.


##Colorbrewer set 5

color1e = numpy.array((166,97,26),dtype='float')/255.
color2e = numpy.array((223,194,125),dtype='float')/255.
color3e = numpy.array((245,245,245),dtype='float')/255.
color4e = numpy.array((128,205,193),dtype='float')/255.
color5e = numpy.array((1,133,113),dtype='float')/255.

fig, ((ax3, ax4),(ax1, ax2)) = pt.subplots(2,2, figsize = (6,6))
#print H2_delta_fit_ypd_err
b1 = ax3.bar(0,H2_delta_fit_ypd, yerr = numpy.array([H2_delta_fit_ypd_err]).T, color=color1e, ecolor='k',capsize=0)#, labels=['$H^2$, delta fit 30 C','$H^2$, delta fit 37 C'] )
b2 = ax3.bar(1,h2_delta_fit_ypd, yerr = numpy.array([h2_delta_fit_ypd_err]).T, color=color2e, ecolor='k',capsize=0)
b4 = ax3.bar(2,r2_delta_fit_ypd_fit_model,yerr = numpy.array([r2_delta_fit_ypd_fit_model_err]).T, color=color4e, ecolor='k',capsize=0)
b5 = ax3.bar(3,r2_delta_fit_ypd_combined_model,yerr = numpy.array([r2_delta_fit_ypd_combined_model_err]).T, color=color5e, ecolor='k',capsize=0)
b6 = ax3.bar(4,model_difference_ypd,yerr = numpy.array([model_difference_ypd_err]).T, color=color3e, ecolor='k',capsize=0)

ax3.legend([b1, b2, b4, b5, b6], ['Broad-sense heritability','Narrow-sense heritability','Fitness model','Combined model','Combined - Fitness'],ncol=1,fontsize=7)
ax3.set_xticks([])
ax3.set_ylim(0,1)
ax3.set_title('OT populations')
#ax1.set_xticklabels(['Fitness gain, OT populations at HT'],fontsize=9)
ax3.set_xlim(-.2, 5.2)
b1 = ax4.bar(0,H2_delta_fit_sc, yerr = numpy.array([H2_delta_fit_sc_err]).T, color=color1e, ecolor='k',capsize=0)#, labels=['$H^2$, delta fit 30 C','$H^2$, delta fit 37 C'] )
b2 = ax4.bar(1,h2_delta_fit_sc, yerr = numpy.array([h2_delta_fit_sc_err]).T, color=color2e, ecolor='k',capsize=0)
b4 = ax4.bar(2,r2_delta_fit_sc_fit_model,yerr = numpy.array([r2_delta_fit_sc_fit_model_err]).T, color=color4e, ecolor='k',capsize=0)
b5 = ax4.bar(3,r2_delta_fit_sc_combined_model,yerr = numpy.array([r2_delta_fit_sc_combined_model_err]).T, color=color5e, ecolor='k',capsize=0)
b6 = ax4.bar(4,model_difference_sc,yerr = numpy.array([model_difference_sc_err]).T, color=color3e, ecolor='k',capsize=0)

ax4.legend([b1, b2, b4, b5, b6], ['Broad-sense heritability','Narrow-sense heritability','Fitness model','Combined model','Combined - Fitness'],ncol=1,fontsize=7)
ax4.set_xticks([])
ax4.set_ylim(0,1)
ax4.set_title('HT populations')
#ax2.set_xticklabels(['Fitness gain, HT populations at OT'])
ax4.set_xlim(-.2, 5.2)

ax1.set_xlim(-2.4,2.5)
ax2.set_xlim(-2.5,2.5)

#ax4.set_yticks(numpy.arange(0,.39,.05))
#ax4.set_yticklabels(numpy.arange(0,39,5),fontsize=8)

ax1.errorbar(-1*predictor_fit_and_qtls_ypd[rm_allele], delta_fits_ypd_in_sc_means[rm_allele], xerr = init_std_errs_sc[rm_allele], yerr = delta_fits_ypd_in_sc_std_errs[rm_allele], fmt='o', color= 'MediumSlateBlue', markeredgewidth=0,capsize=0,alpha=.9,elinewidth=.5)
ax1.errorbar(-1*predictor_fit_and_qtls_ypd[by_allele], delta_fits_ypd_in_sc_means[by_allele], xerr = init_std_errs_sc[by_allele], yerr = delta_fits_ypd_in_sc_std_errs[by_allele], fmt='o', color= 'DarkSlateBlue', markeredgewidth=0,capsize=0,alpha=.9,elinewidth=.5)
ax1.plot(-1*predictor_fit_and_qtls_ypd, predictor_fit_and_qtls_ypd_unscaled, 'k')
ax1.text(1,.25,'$r^2=$' + str(round(rsq_full_model_ypd_for_means, 2)),fontsize=8)
ax1.set_ylim(-.3,.35)
ax1.set_yticks(numpy.arange(-.3,.31,.1))
ax1.set_yticklabels(numpy.arange(-30,31,10))

ax2.errorbar(-1*predictor_fit_and_qtls_sc[rm_allele], delta_fits_sc_in_ypd_means[rm_allele], xerr = init_std_errs_ypd[rm_allele], yerr = delta_fits_sc_in_ypd_std_errs[rm_allele], fmt='o', color= 'Tomato', markeredgewidth=0,capsize=0,alpha=.9,elinewidth=.5)
ax2.errorbar(-1*predictor_fit_and_qtls_sc[by_allele], delta_fits_sc_in_ypd_means[by_allele], xerr = init_std_errs_ypd[by_allele], yerr = delta_fits_sc_in_ypd_std_errs[by_allele], fmt='o', color= 'Brown', markeredgewidth=0,capsize=0,alpha=.9,elinewidth=.5)
ax2.plot(-1*predictor_fit_and_qtls_sc, predictor_fit_and_qtls_sc_unscaled, 'k')
ax2.text(1,.15,'$r^2=$' + str(round(rsq_full_model_sc_for_means, 2)),fontsize=8)
ax2.set_ylim(-.05,.25)
ax2.set_yticks(numpy.arange(0,.21,.1))
ax2.set_yticklabels(numpy.arange(0,21,10))
#ax3.set_ylim(0,.18)
#ax4.set_ylim(0,.38)

#ax1.text(-2.65,.15,'Evolved at 30 C',rotation=90,fontsize=10)
#ax4.text(-2.65,.30,'Evolved at 37 C',rotation=90,fontsize=10)

ax1.set_xlabel('best predictor')
ax1.set_ylabel('Fitness changes at HT (%)')
ax2.set_ylabel('Fitness changes at OT (%)')
ax2.set_xlabel('best predictor')

pt.tight_layout()

pt.savefig('fig9_pleiotropy_2_19_2017b.pdf',bbox_inches='tight')