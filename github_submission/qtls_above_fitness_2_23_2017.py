import matplotlib.pylab as pt
import numpy
import matplotlib.cm as cm
import scipy.stats
from qtl_detection_one_trait import detect_qtls
from qtl_detection_one_trait import detect_qtls_above_fitness
from qtl_detection_one_trait import detect_qtls_with_epistasis
from qtl_detection_one_trait import detect_qtls_with_epistasis2 #This is for detecting qtls with epistasis, above fitness
from qtl_detection_one_trait import calculate_qtl_confidence_intervals
#from calculate_narrow_sense_heritability import rsq_with_bootstrap_linear_model

def rsq_with_bootstrap_nested_model(predictor_mat, delta_fit_vector, helper_matrix):
	#Calculates the difference in rsq value between the model with and without the final column of X, bootstrapping over segregants.
	predictor_mat_expanded = numpy.dot(helper_matrix, predictor_mat)
	beta = numpy.dot(numpy.linalg.inv(numpy.dot(predictor_mat_expanded.T, predictor_mat_expanded)), numpy.dot(predictor_mat_expanded.T, delta_fit_vector))
	rsq = scipy.stats.pearsonr(numpy.dot(predictor_mat_expanded, beta),delta_fit_vector)[0]**2
	
	predictor_mat_inner = numpy.dot(helper_matrix, predictor_mat[:,0:-1])
	beta_inner = numpy.dot(numpy.linalg.inv(numpy.dot(predictor_mat_inner.T, predictor_mat_inner)), numpy.dot(predictor_mat_inner.T, delta_fit_vector))
	rsq_inner = scipy.stats.pearsonr(numpy.dot(predictor_mat_inner, beta_inner),delta_fit_vector)[0]**2
	
	n_segs = len(helper_matrix[0,:])
	rsq_bs_list = []
	for k in range(100):
		bootstrap_inds = numpy.random.choice(numpy.arange(n_segs), size=n_segs)
		expanded_indices = []
		for bootstrap_ind in bootstrap_inds:
			expanded_indices.extend(numpy.where(helper_matrix[:,bootstrap_ind] == 1)[0])
		predictor_mat_expanded_bs = predictor_mat_expanded[expanded_indices,:]
		delta_fit_vector_bs = delta_fit_vector[expanded_indices]
		beta_bs = numpy.dot(numpy.linalg.inv(numpy.dot(predictor_mat_expanded_bs.T, predictor_mat_expanded_bs)), numpy.dot(predictor_mat_expanded_bs.T, delta_fit_vector_bs))
		rsq_bs = scipy.stats.pearsonr(numpy.dot(predictor_mat_expanded_bs, beta_bs),delta_fit_vector_bs)[0]**2
		
		predictor_mat_inner_bs = predictor_mat_inner[expanded_indices,:]
		beta_inner_bs = numpy.dot(numpy.linalg.inv(numpy.dot(predictor_mat_inner_bs.T, predictor_mat_inner_bs)), numpy.dot(predictor_mat_inner_bs.T, delta_fit_vector_bs))
		rsq_inner_bs = scipy.stats.pearsonr(numpy.dot(predictor_mat_inner_bs, beta_inner_bs),delta_fit_vector_bs)[0]**2
		
		rsq_bs_list.append(rsq_bs - rsq_inner_bs)
		
	return rsq - rsq_inner, [numpy.percentile(rsq_bs_list,2.5), numpy.percentile(rsq_bs_list,97.5)]

def calculate_rsq(predictor_mat, delta_fit_vector, helper_matrix):
	
	predictor_mat_expanded = numpy.dot(helper_matrix, predictor_mat)
	beta = numpy.dot(numpy.linalg.inv(numpy.dot(predictor_mat_expanded.T, predictor_mat_expanded)), numpy.dot(predictor_mat_expanded.T, delta_fit_vector))
	rsq = scipy.stats.pearsonr(numpy.dot(predictor_mat_expanded, beta),delta_fit_vector)[0]**2
	
	#check
	# x = numpy.dot(predictor_mat_expanded, beta)
# 	
# 	xavg = numpy.mean(numpy.dot(predictor_mat_expanded, beta))
# 	xstd = numpy.sqrt( numpy.sum( x**2/float(len(x)) ) - xavg**2 )
# 	
# 	y = delta_fit_vector
# 	yavg = numpy.mean(delta_fit_vector)
# 	ystd = numpy.sqrt( numpy.sum( y**2/float(len(y)) ) - yavg**2 )
# 	
# 	r2 = ( numpy.sum(x*y)/float(len(x)) - xavg*yavg )/( xstd*ystd ) ##This is sqrt(rsq)
	
	return rsq
	
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
marker_locations = []
firstline=1
for line in file3:
	if firstline:
		marker_locations = line.strip().split(';')[1].split(',')
		firstline = 0
	else:
		linelist = line.strip().split(';')
		genotype = [int(i) for i in linelist[1].split(',')]
		genotype_mat.append(genotype)

genotype_mat = numpy.array(genotype_mat)
file3.close()

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

rm_allele = numpy.array(genotype_mat[:,3777],dtype='Bool')
by_allele = numpy.array(1 - genotype_mat[:,3777],dtype='Bool')

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
delta_fits_sc_vars = numpy.dot((delta_fits_sc - numpy.dot(helper_matrix_sc_pops,delta_fits_sc_means))**2, helper_matrix_sc_pops)/(pops_per_seg_sc - 1.) + init_std_errs_sc**2 #- measurement_error_var_sc
delta_fits_sc_std_errs = numpy.sqrt(delta_fits_sc_vars/pops_per_seg_sc)

delta_fits_ypd_vars = numpy.dot((delta_fits_ypd - numpy.dot(helper_matrix_ypd_pops,delta_fits_ypd_means))**2, helper_matrix_ypd_pops)/(pops_per_seg_ypd - 1.) + init_std_errs_ypd**2 #- measurement_error_var_ypd
delta_fits_ypd_std_errs = numpy.sqrt(delta_fits_ypd_vars/pops_per_seg_ypd)

delta_fits_ypd_in_sc_vars = numpy.dot((delta_fits_ypd_in_sc - numpy.dot(helper_matrix_ypd_pops,delta_fits_ypd_in_sc_means))**2, helper_matrix_ypd_pops)/(pops_per_seg_ypd - 1.) + init_std_errs_sc**2 #- measurement_error_sc
delta_fits_ypd_in_sc_std_errs = numpy.sqrt(delta_fits_ypd_in_sc_vars/pops_per_seg_ypd)

delta_fits_sc_in_ypd_vars = numpy.dot((delta_fits_sc_in_ypd - numpy.dot(helper_matrix_sc_pops,delta_fits_sc_in_ypd_means))**2, helper_matrix_sc_pops)/(pops_per_seg_sc - 1.) + init_std_errs_ypd**2 #- measurement_error_ypd
delta_fits_sc_in_ypd_std_errs = numpy.sqrt(delta_fits_sc_in_ypd_vars/pops_per_seg_sc)

###Import delta fitness qtls from files

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

merged_qtl_list = list(init_fit_qtls_ypd)
for qtl in init_fit_qtls_sc:
	if min(abs(numpy.array(merged_qtl_list) - qtl)) > 20:
		merged_qtl_list.append(qtl)
		

for qtl in delta_fit_qtls_ypd:
	if min(abs(numpy.array(merged_qtl_list) - qtl)) > 20:
		merged_qtl_list.append(qtl)

for qtl in delta_fit_qtls_sc:
	if min(abs(numpy.array(merged_qtl_list) - qtl)) > 20:
		merged_qtl_list.append(qtl)
						
print merged_qtl_list


###Iterate through the list of qtls; test whether the qtl improves predictions over and above fitness
###This time we will look at significance/variance with respect to mean fitness changes for each segregant
##Fitness only model

X_fit_sc = numpy.concatenate((numpy.ones((num_segs,1)), init_fits_sc.reshape((num_segs,1))), axis = 1)
X_fit_ypd = numpy.concatenate(( numpy.ones((num_segs,1)), init_fits_ypd.reshape((num_segs,1))), axis = 1)

X_fit_sc_pops = numpy.dot(helper_matrix_sc_pops, X_fit_sc)
X_fit_ypd_pops = numpy.dot(helper_matrix_ypd_pops, X_fit_ypd)

beta_sc = numpy.dot(numpy.linalg.inv(numpy.dot(X_fit_sc_pops.T, X_fit_sc_pops)), numpy.dot(X_fit_sc_pops.T, delta_fits_sc))
beta_ypd = numpy.dot(numpy.linalg.inv(numpy.dot(X_fit_ypd_pops.T, X_fit_ypd_pops)), numpy.dot(X_fit_ypd_pops.T, delta_fits_ypd))

rsq_sc_fit = scipy.stats.pearsonr(numpy.dot(X_fit_sc_pops, beta_sc),delta_fits_sc)[0]**2
rsq_ypd_fit = scipy.stats.pearsonr(numpy.dot(X_fit_ypd_pops, beta_ypd),delta_fits_ypd)[0]**2

###

rsq_sc_list = []
rsq_sc_interval_list = []

rsq_ypd_list = []
rsq_ypd_interval_list = []

n_segs = len(pops_per_seg_sc)
X_fit_sc_segs = numpy.concatenate((numpy.ones((n_segs,1)),  init_fits_sc.reshape((n_segs,1))), axis = 1)
X_fit_ypd_segs = numpy.concatenate((numpy.ones((n_segs,1)),  init_fits_ypd.reshape((n_segs,1))), axis = 1)

for qtl in merged_qtl_list:
	
	X_qtls_sc = numpy.append(X_fit_sc_segs, genotype_mat[:, qtl].reshape((n_segs,1)), axis=1)
	
	
	rsq_sc, rsq_sc_interval = rsq_with_bootstrap_nested_model(X_qtls_sc, delta_fits_sc, helper_matrix = helper_matrix_sc_pops)
	
	rsq_sc_list.append(rsq_sc)
	rsq_sc_interval_list.append(rsq_sc_interval)
	
for qtl in merged_qtl_list:	
	X_qtls_ypd = numpy.append(X_fit_ypd_segs, genotype_mat[:, qtl].reshape((n_segs,1)), axis=1)
	rsq_ypd, rsq_ypd_interval = rsq_with_bootstrap_nested_model(X_qtls_ypd, delta_fits_ypd, helper_matrix = helper_matrix_ypd_pops)
	
	rsq_ypd_list.append(rsq_ypd)
	rsq_ypd_interval_list.append(rsq_ypd_interval)
	
print rsq_sc_list
print rsq_ypd_list

X_qtls_sc_linear = X_fit_sc_segs
X_qtls_ypd_linear = X_fit_ypd_segs
for qtl in merged_qtl_list:
	
	X_qtls_sc_linear = numpy.append(X_qtls_sc_linear, genotype_mat[:, qtl].reshape((n_segs,1)), axis=1)
	
for qtl in merged_qtl_list:
	X_qtls_ypd_linear = numpy.append(X_qtls_ypd_linear, genotype_mat[:, qtl].reshape((n_segs,1)), axis=1)

rsq_all_linear_sc = calculate_rsq(X_qtls_sc_linear, delta_fits_sc, helper_matrix =  helper_matrix_sc_pops)
rsq_all_linear_ypd = calculate_rsq(X_qtls_ypd_linear, delta_fits_ypd, helper_matrix =  helper_matrix_ypd_pops)

print rsq_all_linear_sc
print rsq_all_linear_ypd

F_sc_list = []
F_ypd_list = []
num_free_parameters_sc = X_qtls_sc_linear.shape[1] + 1 ##Is this right?
num_free_parameters_ypd = X_qtls_ypd_linear.shape[1] + 1
##Calculate rsq with each coefficient left out in turn

for i in range(X_qtls_sc_linear.shape[1]):
	X_qtls_sc_linear_temp = numpy.delete(X_qtls_sc_linear, i, axis=1)
	
	rsq_sc_temp = calculate_rsq(X_qtls_sc_linear_temp, delta_fits_sc, helper_matrix = helper_matrix_sc_pops)

	rsq_sc_list.append(rsq_sc_temp)
	
	F_sc_list.append((rsq_all_linear_sc - rsq_sc_temp)/rsq_sc_temp*(num_pops_sc - num_free_parameters_sc))

for i in range(X_qtls_ypd_linear.shape[1]):
	X_qtls_ypd_linear_temp = numpy.delete(X_qtls_ypd_linear, i, axis=1)
	rsq_ypd_temp = calculate_rsq(X_qtls_ypd_linear_temp, delta_fits_ypd, helper_matrix = helper_matrix_ypd_pops)
	rsq_ypd_list.append(rsq_ypd_temp)
	F_ypd_list.append((rsq_all_linear_ypd - rsq_ypd_temp)/rsq_ypd_temp*(num_pops_ypd - num_free_parameters_ypd))


##Bonferroni correction: to be 95% confident of each coefficient being real, we want them to be significant at a p<.004 level.

alpha_sc = 1. - (.95)**(1./float(num_free_parameters_sc))
alpha_ypd = 1. - (.95)**(1./float(num_free_parameters_ypd))


p_vals_sc = 1. - scipy.stats.f.cdf(F_sc_list, 1, num_pops_sc - num_free_parameters_sc)
p_vals_ypd = 1. - scipy.stats.f.cdf(F_ypd_list, 1, num_pops_ypd - num_free_parameters_ypd)

###Identify significant coefficients

coefficient_inds_sc = numpy.arange(num_free_parameters_sc)[p_vals_sc < alpha_sc]
coefficient_inds_ypd = numpy.arange(num_free_parameters_ypd)[p_vals_ypd < alpha_ypd]

print coefficient_inds_sc
print coefficient_inds_ypd

###Finally, fit the full model with these coefficients and record everything in a table.
###Reorder the qtl coefficients by F statistic

sig_coeffs_sc = coefficient_inds_sc

sig_coeffs_ypd = coefficient_inds_ypd
X_qtls_sc_linear_restricted = X_qtls_sc_linear[:,sig_coeffs_sc]
X_qtls_ypd_linear_restricted = X_qtls_ypd_linear[:,sig_coeffs_ypd]

X_qtls_sc_pops_restricted = numpy.dot(helper_matrix_sc_pops, X_qtls_sc_linear_restricted)
X_qtls_ypd_pops_restricted = numpy.dot(helper_matrix_ypd_pops, X_qtls_ypd_linear_restricted)


beta_sc = numpy.dot( numpy.linalg.inv( numpy.dot(X_qtls_sc_pops_restricted.T, X_qtls_sc_pops_restricted) ), numpy.dot(X_qtls_sc_pops_restricted.T, delta_fits_sc) )
beta_ypd = numpy.dot( numpy.linalg.inv( numpy.dot(X_qtls_ypd_pops_restricted.T, X_qtls_ypd_pops_restricted) ), numpy.dot(X_qtls_ypd_pops_restricted.T, delta_fits_ypd) )

###For each coefficient, record the additional rsq and the F-statistic associated with adding this coefficient.
###We will add the intercept first, then the fitness term, then the other coefficients in descending order.

rsq_sc_nested = []
rsq_ypd_nested = []

F_sc_nested = []
F_ypd_nested = []

for i in range(X_qtls_sc_pops_restricted.shape[1]):
	num_free_parameters = i + 1
	rsq = calculate_rsq( X_qtls_sc_linear_restricted[:,0:i+1].reshape((n_segs,i+1)), delta_fits_sc, helper_matrix = helper_matrix_sc_pops )
	if i > 1.5:
		F = (rsq - rsq_sc_nested[i-1])/rsq_sc_nested[i-1]*(num_pops_sc - num_free_parameters)
	else:
		F = 'NA'
	#print rsq, rsq2
	rsq_sc_nested.append(rsq)
	F_sc_nested.append(F)

for i in range(X_qtls_ypd_pops_restricted.shape[1]):
	num_free_parameters = i + 1
	rsq = calculate_rsq(X_qtls_ypd_linear_restricted[:,0:i+1].reshape((n_segs,i+1)), delta_fits_ypd, helper_matrix = helper_matrix_ypd_pops)
	if i > 1.5:
		F = (rsq - rsq_ypd_nested[i-1])/rsq_ypd_nested[i-1]*(num_pops_ypd - num_free_parameters)
	else:
		F = 'NA'
	rsq_ypd_nested.append(rsq)
	F_ypd_nested.append(F)


print merged_qtl_list
print rsq_sc_nested

print merged_qtl_list
print rsq_ypd_nested

###Record results in a table

merged_qtl_list = numpy.array(merged_qtl_list)

filename = 'data/Above_fitness_qtl_model_parameters_linear_2_23_2017.txt'
# 
file = open(filename,'w')
file.write('#Trait: Delta fitness in SC at 37C' + '\n')
file.write((',').join(('#Term (qtl marker location)','term number','Additional variance explained','Model Coefficient', 'F statistic')) + '\n')
file.write((',').join(('intercept',str(0),'NA', str(beta_sc[0]), str(F_sc_nested[0]))) + '\n')
file.write((',').join(('initial fitness',str(1),str(rsq_sc_nested[1] - rsq_sc_nested[0]), str(beta_sc[1]), str(F_sc_nested[1]))) + '\n')
counter = 1
for qtl in merged_qtl_list[sig_coeffs_sc[2:]-2]:
	
	file.write((',').join((str(qtl),str(rsq_sc_nested[counter + 1] - rsq_sc_nested[counter]), str(beta_sc[counter + 1]), str(F_sc_nested[counter + 1]))) + '\n')
	counter += 1
# 
file.write('#Trait: Delta fitness in YPD at 30C' + '\n')
file.write((',').join(('#Term (qtl marker location)','term number','Additional variance explained','Model Coefficient', 'F statistic')) + '\n')
file.write((',').join(('intercept',str(0),'NA', str(beta_ypd[0]), str(F_ypd_nested[0]))) + '\n')
file.write((',').join(('initial fitness',str(1),str(rsq_ypd_nested[1] - rsq_ypd_nested[0]), str(beta_ypd[1]), str(F_ypd_nested[1]))) + '\n')

counter = 1
for qtl in merged_qtl_list[sig_coeffs_ypd[2:]-2]:
	
	file.write((',').join((str(qtl),str(rsq_ypd_nested[counter + 1] - rsq_ypd_nested[counter]), str(beta_ypd[counter + 1]), str(F_ypd_nested[counter + 1]))) + '\n')
	counter += 1
file.close()
##

###Epistatic model
kre33_loc = merged_qtl_list[0]
kre33_allele = genotype_mat[:,kre33_loc].reshape((n_segs,1)) - .5
#print kre33_allele
X_qtls_sc_segs = numpy.concatenate((numpy.ones((n_segs,1)),  init_fits_sc.reshape((n_segs,1)), kre33_allele), axis = 1)
X_qtls_ypd_segs = numpy.concatenate((numpy.ones((n_segs,1)),  init_fits_ypd.reshape((n_segs,1)), kre33_allele), axis = 1)

for qtl in merged_qtl_list[1:]:
	
	X_qtls_sc_segs = numpy.concatenate((X_qtls_sc_segs, genotype_mat[:,qtl].reshape((n_segs,1)) - .5, (genotype_mat[:,qtl].reshape((n_segs,1))- .5)*kre33_allele), axis=1)

for qtl in merged_qtl_list[1:]:
	X_qtls_ypd_segs = numpy.concatenate((X_qtls_ypd_segs, genotype_mat[:,qtl].reshape((n_segs,1)) - .5, (genotype_mat[:,qtl].reshape((n_segs,1))- .5)*kre33_allele), axis=1)

rsq_sc_all = calculate_rsq(X_qtls_sc_segs, delta_fits_sc, helper_matrix = helper_matrix_sc_pops)
rsq_ypd_all = calculate_rsq(X_qtls_ypd_segs, delta_fits_ypd, helper_matrix = helper_matrix_ypd_pops)

rsq_sc_epistatic_list = []
rsq_ypd_epistatic_list = []

F_sc_epistatic_list = []
F_ypd_epistatic_list = []
num_free_parameters_sc = X_qtls_sc_segs.shape[1] + 1 ##Is this right?
num_free_parameters_ypd = X_qtls_ypd_segs.shape[1] + 1
##Calculate rsq with each coefficient left out in turn

for i in range(X_qtls_sc_segs.shape[1]):
	X_qtls_sc_segs_temp = numpy.delete(X_qtls_sc_segs, i, axis=1)
	
	rsq_sc_temp = calculate_rsq(X_qtls_sc_segs_temp, delta_fits_sc, helper_matrix = helper_matrix_sc_pops)

	rsq_sc_epistatic_list.append(rsq_sc_temp)
	
	F_sc_epistatic_list.append((rsq_sc_all - rsq_sc_temp)/rsq_sc_temp*(num_pops_sc - num_free_parameters_sc))

for i in range(X_qtls_ypd_segs.shape[1]):
	X_qtls_ypd_segs_temp = numpy.delete(X_qtls_ypd_segs, i, axis=1)
	rsq_ypd_temp = calculate_rsq(X_qtls_ypd_segs_temp, delta_fits_ypd, helper_matrix = helper_matrix_ypd_pops)
	rsq_ypd_epistatic_list.append(rsq_ypd_temp)
	F_ypd_epistatic_list.append((rsq_ypd_all - rsq_ypd_temp)/rsq_ypd_temp*(num_pops_ypd - num_free_parameters_ypd))


##Bonferroni correction: to be 95% confident of each coefficient being real, we want them to be significant at a p<.004 level.

alpha_sc = 1. - (.95)**(1./float(num_free_parameters_sc))
alpha_ypd = 1. - (.95)**(1./float(num_free_parameters_ypd))


p_vals_sc = 1. - scipy.stats.f.cdf(F_sc_epistatic_list, 1, num_pops_sc - num_free_parameters_sc)
p_vals_ypd = 1. - scipy.stats.f.cdf(F_ypd_epistatic_list, 1, num_pops_ypd - num_free_parameters_ypd)

###Identify significant coefficients

coefficient_inds_sc = numpy.arange(num_free_parameters_sc)[p_vals_sc < alpha_sc]
coefficient_inds_ypd = numpy.arange(num_free_parameters_ypd)[p_vals_ypd < alpha_ypd]

print coefficient_inds_sc
print coefficient_inds_ypd

###Finally, fit the full model with these coefficients and record everything in a table.
###Reorder the qtl coefficients by F statistic

sig_coeffs_sc = coefficient_inds_sc

sig_coeffs_ypd = coefficient_inds_ypd
X_qtls_sc_segs_restricted = X_qtls_sc_segs[:,sig_coeffs_sc]
X_qtls_ypd_segs_restricted = X_qtls_ypd_segs[:,sig_coeffs_ypd]

print X_qtls_ypd_segs_restricted.shape

X_qtls_sc_pops = numpy.dot(helper_matrix_sc_pops, X_qtls_sc_segs_restricted)
X_qtls_ypd_pops = numpy.dot(helper_matrix_ypd_pops, X_qtls_ypd_segs_restricted)


beta_sc = numpy.dot( numpy.linalg.inv( numpy.dot(X_qtls_sc_pops.T, X_qtls_sc_pops) ), numpy.dot(X_qtls_sc_pops.T, delta_fits_sc) )
beta_ypd = numpy.dot( numpy.linalg.inv( numpy.dot(X_qtls_ypd_pops.T, X_qtls_ypd_pops) ), numpy.dot(X_qtls_ypd_pops.T, delta_fits_ypd) )

###For each coefficient, record the additional rsq and the F-statistic associated with adding this coefficient.
###We will add the intercept first, then the fitness term, then the other coefficients in descending order.

rsq_sc_nested = []
rsq_ypd_nested = []

F_sc_nested = []
F_ypd_nested = []

for i in range(X_qtls_sc_pops.shape[1]):
	num_free_parameters = i + 1
	rsq = calculate_rsq( X_qtls_sc_segs_restricted[:,0:i+1].reshape((n_segs,i+1)), delta_fits_sc, helper_matrix = helper_matrix_sc_pops )
	if i > 1.5:
		F = (rsq - rsq_sc_nested[i-1])/rsq_sc_nested[i-1]*(num_pops_sc - num_free_parameters)
	else:
		F = 'NA'
	#print rsq, rsq2
	rsq_sc_nested.append(rsq)
	F_sc_nested.append(F)

for i in range(X_qtls_ypd_pops.shape[1]):
	num_free_parameters = i + 1
	rsq = calculate_rsq(X_qtls_ypd_segs_restricted[:,0:i+1].reshape((n_segs,i+1)), delta_fits_ypd, helper_matrix = helper_matrix_ypd_pops)
	if i > 1.5:
		F = (rsq - rsq_ypd_nested[i-1])/rsq_ypd_nested[i-1]*(num_pops_ypd - num_free_parameters)
	else:
		F = 'NA'
	rsq_ypd_nested.append(rsq)
	F_ypd_nested.append(F)


print merged_qtl_list
print rsq_sc_nested

print merged_qtl_list
print rsq_ypd_nested

filename = 'data/Above_fitness_epistatic_qtl_model_parameters_epistatic_2_23_2017.txt'
# # 
file = open(filename,'w')
file.write('#Trait: Delta fitness in SC at 37C' + '\n')
file.write((',').join(('#Term (qtl marker location)','term number','Additional variance explained','Model Coefficient', 'F statistic')) + '\n')
file.write((',').join(('#intercept',str(0),'NA', str(beta_sc[0]), str(F_sc_nested[0]))) + '\n')
file.write((',').join(('#initial fitness',str(1),str(rsq_sc_nested[1] - rsq_sc_nested[0]), str(beta_sc[1]), str(F_sc_nested[1]))) + '\n')
file.write((',').join(('Kre33 (3777)',str(3777),str(2),str(rsq_sc_nested[2] - rsq_sc_nested[1]), str(beta_sc[2]), str(F_sc_nested[2]))) + '\n')
counter = 3
for coefficient in coefficient_inds_sc[3:]:
	
	if coefficient % 2 < .5: #even
		locus = merged_qtl_list[int(coefficient/2.) - 1]
		file.write((',').join(('epistatic with Kre33',str(locus),str(coefficient),str(rsq_sc_nested[counter] - rsq_sc_nested[counter-1]), str(beta_sc[counter]), str(F_sc_nested[counter]) )) + '\n')
	else:
		locus = merged_qtl_list[int((coefficient+1)/2.) - 1]
		file.write((',').join(('linear',str(locus),str(coefficient),str(rsq_sc_nested[counter] - rsq_sc_nested[counter-1]), str(beta_sc[counter]), str(F_sc_nested[counter]) )) + '\n')
	counter += 1
		
# # 
file.write('#Trait: Delta fitness in YPD at 30C' + '\n')
file.write((',').join(('#Term (qtl marker location)','term number','Additional variance explained','Model Coefficient', 'F statistic')) + '\n')
file.write((',').join(('#intercept',str(0),'NA', str(beta_ypd[0]), str(F_ypd_nested[0]))) + '\n')
file.write((',').join(('#initial fitness',str(1),str(rsq_ypd_nested[1] - rsq_ypd_nested[0]), str(beta_ypd[1]), str(F_ypd_nested[1]))) + '\n')
file.write((',').join(('Kre33 (3777)',str(3777),str(2),str(rsq_ypd_nested[2] - rsq_ypd_nested[1]), str(beta_ypd[2]), str(F_ypd_nested[2]))) + '\n')
counter = 3
for coefficient in coefficient_inds_ypd[3:]:
	
	if coefficient % 2 < .5: #even
		locus = merged_qtl_list[int(coefficient/2.) - 1]
		file.write((',').join(('epistatic with Kre33',str(locus),str(coefficient),str(rsq_ypd_nested[counter] - rsq_ypd_nested[counter-1]), str(beta_ypd[counter]), str(F_ypd_nested[counter]) )) + '\n')
	else:
		locus = merged_qtl_list[int((coefficient+1)/2.) - 1]
		file.write((',').join(('linear',str(locus),str(coefficient),str(rsq_ypd_nested[counter] - rsq_ypd_nested[counter-1]), str(beta_ypd[counter]), str(F_ypd_nested[counter]) )) + '\n')
	counter += 1
file.close()