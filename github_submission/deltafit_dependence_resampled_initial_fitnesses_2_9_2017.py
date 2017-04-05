import numpy
import scipy.stats
import matplotlib.pylab as pt

def calculate_rsq(X, y):
	
	beta = numpy.dot(numpy.linalg.inv(numpy.dot(X.T, X)), numpy.dot(X.T,y))
	rsq = scipy.stats.pearsonr( numpy.dot(X, beta), y )[0]**2
	
	return rsq

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
genotype_mat = []
file3 = open(filename3,'r')
for line in file3:
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
n_segs = num_segs
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
delta_fits_sc_vars = numpy.dot((delta_fits_sc - numpy.dot(helper_matrix_sc_pops,delta_fits_sc_means))**2, helper_matrix_sc_pops)/(pops_per_seg_sc - 1.) #- measurement_error_var_sc
delta_fits_sc_std_errs = numpy.sqrt(delta_fits_sc_vars/pops_per_seg_sc)

delta_fits_ypd_vars = numpy.dot((delta_fits_ypd - numpy.dot(helper_matrix_ypd_pops,delta_fits_ypd_means))**2, helper_matrix_ypd_pops)/(pops_per_seg_ypd - 1.) #- measurement_error_var_ypd
delta_fits_ypd_std_errs = numpy.sqrt(delta_fits_ypd_vars/pops_per_seg_ypd)

delta_fits_ypd_in_sc_vars = numpy.dot((delta_fits_ypd_in_sc - numpy.dot(helper_matrix_ypd_pops,delta_fits_ypd_in_sc_means))**2, helper_matrix_ypd_pops)/(pops_per_seg_ypd - 1.) #- measurement_error_sc
delta_fits_ypd_in_sc_std_errs = numpy.sqrt(delta_fits_ypd_in_sc_vars/pops_per_seg_ypd)

delta_fits_sc_in_ypd_vars = numpy.dot((delta_fits_sc_in_ypd - numpy.dot(helper_matrix_sc_pops,delta_fits_sc_in_ypd_means))**2, helper_matrix_sc_pops)/(pops_per_seg_sc - 1.) #- measurement_error_ypd
delta_fits_sc_in_ypd_std_errs = numpy.sqrt(delta_fits_sc_in_ypd_vars/pops_per_seg_sc)

#####Here we are going to test the dependence of delta_fits_sc_in_ypd (fitness gains of populations evolved in sc measured in ypd) on initial fitness in ypd, in the presence of measurement noise.
#We are checking whether initial fitness in ypd is more predictive of delta_fits_sc_in_ypd than initial fitness in sc, even when the errors in these two quantities are uncorrelated.
#Similar for delta_fits_ypd_in_sc

X_sc = numpy.concatenate( (numpy.ones((n_segs, 1)), init_fits_sc.reshape((n_segs, 1))), axis=1)
X_ypd = numpy.concatenate( (numpy.ones((n_segs, 1)), init_fits_ypd.reshape((n_segs, 1))), axis=1)

rsq_sc_home=calculate_rsq(X_sc, delta_fits_sc_means)
rsq_ypd_home=calculate_rsq(X_ypd, delta_fits_ypd_means)
	
n_iter = 1000
init_fits_sc_temp = numpy.zeros_like(init_fits_sc)
init_fits_ypd_temp = numpy.zeros_like(init_fits_ypd)

rsq_resampled_sc_pops = []
rsq_resampled_ypd_pops = []
rsq_sc_resampled_home = []
rsq_ypd_resampled_home = []
for i in range(n_iter):
	
	resampled_errs_sc = []
	resampled_errs_ypd = []
	
	#resample errors
	
	for j in range(num_segs):
		if numpy.isnan(init_std_errs_sc[j]):
			
			resampled_errs_sc.append(numpy.random.normal(0, scale=numpy.nanmean(init_std_errs_sc)))
			resampled_errs_ypd.append(numpy.random.normal(0, scale=init_std_errs_ypd[j]))
		else:
			resampled_errs_sc.append(numpy.random.normal(0, scale=init_std_errs_sc[j]))
			resampled_errs_ypd.append(numpy.random.normal(0, scale=init_std_errs_ypd[j]))
			
	resampled_errs_sc = numpy.array(resampled_errs_sc)
	resampled_errs_ypd = numpy.array(resampled_errs_ypd)
	
	init_fits_sc_temp = init_fits_sc + resampled_errs_sc
	init_fits_ypd_temp = init_fits_ypd + resampled_errs_ypd
	
	#measure correlation between initial fitnesses with resampled errors and delta fitnesses
	
	X_both_fits = numpy.concatenate( (numpy.ones((n_segs, 1)), init_fits_sc_temp.reshape((n_segs, 1)), init_fits_ypd_temp.reshape((n_segs, 1))), axis=1)
	X_sc_temp = numpy.concatenate( (numpy.ones((n_segs, 1)), init_fits_sc_temp.reshape((n_segs, 1))), axis=1)
	X_ypd_temp = numpy.concatenate( (numpy.ones((n_segs, 1)), init_fits_ypd_temp.reshape((n_segs, 1))), axis=1)
	
	rsq_sc_resampled_home.append(calculate_rsq(X_sc_temp, delta_fits_sc_means))
	rsq_ypd_resampled_home.append(calculate_rsq(X_ypd_temp, delta_fits_ypd_means))
	
	rsq_sc_both = calculate_rsq(X_both_fits, delta_fits_sc_in_ypd_means)
	rsq_ypd_both = calculate_rsq(X_both_fits, delta_fits_ypd_in_sc_means)
	rsq_sc_sc = calculate_rsq(X_sc_temp, delta_fits_sc_in_ypd_means)
	rsq_sc_ypd = calculate_rsq(X_ypd_temp, delta_fits_sc_in_ypd_means)
	rsq_ypd_sc = calculate_rsq(X_sc_temp, delta_fits_ypd_in_sc_means)
	rsq_ypd_ypd = calculate_rsq(X_ypd_temp, delta_fits_ypd_in_sc_means)
	
	rsq_resampled_sc_pops.append([rsq_sc_both, rsq_sc_ypd, rsq_sc_sc])
	rsq_resampled_ypd_pops.append([rsq_ypd_both, rsq_ypd_sc, rsq_ypd_ypd])

rsq_resampled_sc_pops = numpy.array(rsq_resampled_sc_pops)
rsq_resampled_ypd_pops = numpy.array(rsq_resampled_ypd_pops)

print rsq_sc_home, numpy.mean(rsq_sc_resampled_home)
print rsq_ypd_home, numpy.mean(rsq_ypd_resampled_home)
print numpy.mean(rsq_sc_home) - numpy.mean(rsq_sc_resampled_home)
print numpy.mean(rsq_ypd_home) - numpy.mean(rsq_ypd_resampled_home)

#print rsq_resampled_sc_pops
##Extract p-values from distributions

p1 = numpy.sum(rsq_resampled_sc_pops[:,0] - rsq_resampled_sc_pops[:,1] > 0)/float(n_iter)
p2 = numpy.sum(rsq_resampled_ypd_pops[:,0] - rsq_resampled_ypd_pops[:,1] > 0)/float(n_iter)

p3 = numpy.sum(rsq_resampled_sc_pops[:,1] - rsq_resampled_sc_pops[:,2] > 0)/float(n_iter)
p4 = numpy.sum(rsq_resampled_ypd_pops[:,1] - rsq_resampled_ypd_pops[:,2] > 0)/float(n_iter)

print 'measurement env better than home env, 30C, 37C', p3, p4

##Likelihood ratio test for nested models (measurement envt model vs. both model)

r1 = -2*num_pops_sc*numpy.log((1.- numpy.mean(rsq_resampled_sc_pops[:,0]))/(1.- numpy.mean(rsq_resampled_sc_pops[:,1])))
r2 = -2*num_pops_ypd*numpy.log((1.-numpy.mean(rsq_resampled_ypd_pops[:,0]))/(1.-numpy.mean(rsq_resampled_ypd_pops[:,1])))

print r1, r2

##Sanity check histograms
pt.figure()
pt.title('Predicting fitness gains of sc pops in ypd')
pt.hist(rsq_resampled_sc_pops[:,0] - rsq_resampled_sc_pops[:,1], bins=50)
pt.xlabel('$r^2$, both fitnesses model - $r^2$, ypd only')
pt.savefig('sc_pops_in_ypd_model_eval.pdf',bbox_inches='tight')
pt.figure()
pt.title('Predicting fitness gains of ypd pops in sc')
pt.hist(rsq_resampled_ypd_pops[:,0] - rsq_resampled_ypd_pops[:,1], bins=50)
pt.xlabel('$r^2$, both fitnesses model - $r^2$, ypd only')
pt.savefig('ypd_pops_in_sc_model_eval.pdf',bbox_inches='tight')
pt.figure()
pt.title('Predicting fitness gains of sc pops in ypd')
pt.hist(rsq_resampled_sc_pops[:,1] - rsq_resampled_sc_pops[:,2], bins=50)
pt.xlabel('$r^2$, ypd only - $r^2$, sc only')
pt.savefig('sc_pops_in_ypd_model_eval2.pdf',bbox_inches='tight')
pt.figure()
pt.title('Predicting fitness gains of ypd pops in sc')
pt.hist(rsq_resampled_ypd_pops[:,1] - rsq_resampled_ypd_pops[:,2], bins=50)
pt.xlabel('$r^2$, sc only - $r^2$, ypd only')
pt.savefig('ypd_pops_in_sc_model_eval2.pdf',bbox_inches='tight')
	

