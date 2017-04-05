import matplotlib.pylab as pt
import numpy
import matplotlib.cm as cm
import scipy.stats
from qtl_detection_one_trait import detect_qtls_one_envt
from qtl_detection_one_trait import detect_qtls_above_fitness
from qtl_detection_one_trait import detect_qtls_with_epistasis
from qtl_detection_one_trait import detect_qtls_with_epistasis2 #This is for detecting qtls with epistasis, above fitness
from qtl_detection_one_trait import calculate_qtl_confidence_intervals_lods

#This file was modified on 1/18/2017 to do the following things: first, to detect QTLs using only informative, not redundant loci.
#Second, to detect qtls separately for the two environments, in addition to jointly, to compare.
#Further modified on 2/6/2017 to detect qtls on the mean trait value rather than including all replicates, for comparison

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


########QTL detection

########We output a table in the format 'marker number', 'chromosome location', 'additional fraction of variance explained', 'confidence intervals'. We will later combine this table with a .gff to add 'genes within confidence intervals'.
########For each trait, we will first calculate the QTL locations iteratively.
########We will then calculate confidence intervals for the QTL locations by bootstrapping over the 230 segregants and measuring the distribution of the location of the QTL peak

n_segs = len(pops_per_seg_ypd)

###Initial fitness qtl detection (linear model)

#Detect QTL locations
qtls_init_sc, beta_sc_init, intervals_sc = detect_qtls_one_envt(genotype_mat, init_fits_sc, helper_matrix = numpy.identity(n_segs), pops_per_seg =numpy.ones((n_segs,)))
qtls_init_ypd, beta_ypd_init, intervals_ypd = detect_qtls_one_envt(genotype_mat, init_fits_ypd, helper_matrix = numpy.identity(n_segs), pops_per_seg =numpy.ones((n_segs,)))

lower_CIs_sc = intervals_sc[:,0]
upper_CIs_sc = intervals_sc[:,1]

lower_CIs_ypd = intervals_ypd[:,0]
upper_CIs_ypd = intervals_ypd[:,1]

#Fraction of variance explained
X_qtls_sc = numpy.ones((n_segs,1))
rsq_sc_list = []

for qtl in qtls_init_sc:
	X_qtls_sc = numpy.append(X_qtls_sc, genotype_mat[:, qtl].reshape((n_segs,1)), axis=1)
	beta_sc = numpy.dot(numpy.linalg.inv(numpy.dot(X_qtls_sc.T, X_qtls_sc)), numpy.dot(X_qtls_sc.T, init_fits_sc))
	rsq_sc_list.append(scipy.stats.pearsonr(numpy.dot(X_qtls_sc, beta_sc),init_fits_sc)[0]**2)
	
rsq_sc_list = numpy.array(rsq_sc_list)
rsq_sc_list_diff = numpy.append([rsq_sc_list[0]],rsq_sc_list[1:]-rsq_sc_list[0:-1],axis=0)

X_qtls_ypd = numpy.ones((n_segs,1))
rsq_ypd_list = []

for qtl in qtls_init_ypd:
	X_qtls_ypd = numpy.append(X_qtls_ypd, genotype_mat[:, qtl].reshape((n_segs,1)), axis=1)
	beta_ypd = numpy.dot(numpy.linalg.inv(numpy.dot(X_qtls_ypd.T, X_qtls_ypd)), numpy.dot(X_qtls_ypd.T, init_fits_ypd))
	rsq_ypd_list.append(scipy.stats.pearsonr(numpy.dot(X_qtls_ypd, beta_ypd),init_fits_ypd)[0]**2)

rsq_ypd_list = numpy.array(rsq_ypd_list)
rsq_ypd_list_diff = numpy.append([rsq_ypd_list[0]],rsq_ypd_list[1:]-rsq_ypd_list[0:-1],axis=0)

#Write everything in a table
filename_out1 = 'data/initial_fitness_qtl_table_2_5_2017.csv'
file_out1 = open(filename_out1,'w')
file_out1.write('#Trait: initial segregant fitness, 37 C' + '\n')
file_out1.write('#' +(',').join(('marker number','location','location lower bound','location upper bound','var. explained YPD 30C', 'var. explained SC 37C','\n')))
for i in range(len(qtls_init_sc)):
	marker = qtls_init_sc[i]
	loc = marker_locations[marker]
	loc_lower = marker_locations[lower_CIs_sc[i]]
	loc_upper = marker_locations[upper_CIs_sc[i]]
	var_sc = rsq_sc_list_diff[i]
	
	file_out1.write((',').join((str(marker), loc, loc_lower, loc_upper, 'NA', str(var_sc),'\n')))
file_out1.write('#Trait: initial segregant fitness, 30 C' + '\n')
for i in range(len(qtls_init_ypd)):
	marker = qtls_init_ypd[i]
	loc = marker_locations[marker]
	loc_lower = marker_locations[lower_CIs_ypd[i]]
	loc_upper = marker_locations[upper_CIs_ypd[i]]
	var_ypd = rsq_ypd_list_diff[i]
	
	file_out1.write((',').join((str(marker), loc, loc_lower, loc_upper, str(var_ypd),'NA','\n')))
	
file_out1.close()

##Delta fitness qtl detection (linear model)

qtls_delta_sc, beta_sc_delta, intervals_sc = detect_qtls_one_envt(genotype_mat, delta_fits_sc_means, helper_matrix = numpy.identity(n_segs), pops_per_seg = numpy.ones((n_segs,)))
qtls_delta_ypd, beta_ypd_delta, intervals_ypd = detect_qtls_one_envt(genotype_mat, delta_fits_ypd_means, helper_matrix = numpy.identity(n_segs), pops_per_seg = numpy.ones((n_segs,)))

lower_CIs_sc = intervals_sc[:,0]
upper_CIs_sc = intervals_sc[:,1]

lower_CIs_ypd = intervals_ypd[:,0]
upper_CIs_ypd = intervals_ypd[:,1]

#Fraction of variance explained
X_qtls_sc = numpy.ones((n_segs,1))
rsq_sc_list = []

for qtl in qtls_delta_sc:
	X_qtls_sc = numpy.append(X_qtls_sc, genotype_mat[:, qtl].reshape((n_segs,1)), axis=1)
	#X_qtls_sc_extended = numpy.dot(helper_matrix_sc_pops, X_qtls_sc)
	beta_sc = numpy.dot(numpy.linalg.inv(numpy.dot(X_qtls_sc.T, X_qtls_sc)), numpy.dot(X_qtls_sc.T, delta_fits_sc_means))
	rsq_sc_list.append(scipy.stats.pearsonr(numpy.dot(X_qtls_sc, beta_sc),delta_fits_sc_means)[0]**2)
	
rsq_sc_list = numpy.array(rsq_sc_list)
rsq_sc_list_diff = numpy.append([rsq_sc_list[0]],rsq_sc_list[1:]-rsq_sc_list[0:-1],axis=0)

X_qtls_ypd = numpy.ones((n_segs,1))
rsq_ypd_list = []

for qtl in qtls_delta_ypd:
	X_qtls_ypd = numpy.append(X_qtls_ypd, genotype_mat[:, qtl].reshape((n_segs,1)), axis=1)
	#X_qtls_ypd_extended = numpy.dot(helper_matrix_ypd_pops, X_qtls_ypd)
	beta_ypd = numpy.dot(numpy.linalg.inv(numpy.dot(X_qtls_ypd.T, X_qtls_ypd)), numpy.dot(X_qtls_ypd.T, delta_fits_ypd_means))
	rsq_ypd_list.append(scipy.stats.pearsonr(numpy.dot(X_qtls_ypd, beta_ypd),delta_fits_ypd_means)[0]**2)
	
rsq_ypd_list = numpy.array(rsq_ypd_list)
rsq_ypd_list_diff = numpy.append([rsq_ypd_list[0]],rsq_ypd_list[1:]-rsq_ypd_list[0:-1],axis=0)

#Write everything in a table
filename_out1 = 'data/delta_fitness_qtl_table_2_5_2017.csv'
file_out1 = open(filename_out1,'w')
file_out1.write('#Trait: mean delta segregant fitness 37C' + '\n')
file_out1.write('#' +(',').join(('marker number','location','location lower bound','location upper bound','var. explained YPD 30C', 'var. explained SC 37C','\n')))
for i in range(len(qtls_delta_sc)):
	marker = qtls_delta_sc[i]
	loc = marker_locations[marker]
	loc_lower = marker_locations[lower_CIs_sc[i]]
	loc_upper = marker_locations[upper_CIs_sc[i]]
	var_sc = rsq_sc_list_diff[i]
	
	file_out1.write((',').join((str(marker), loc, loc_lower, loc_upper, 'NA', str(var_sc),'\n')))
	
file_out1.write('#Trait: mean delta segregant fitness 30C' + '\n')
for i in range(len(qtls_delta_ypd)):
	marker = qtls_delta_ypd[i]
	loc = marker_locations[marker]
	loc_lower = marker_locations[lower_CIs_ypd[i]]
	loc_upper = marker_locations[upper_CIs_ypd[i]]
	var_ypd = rsq_ypd_list_diff[i]
	
	file_out1.write((',').join((str(marker), loc, loc_lower, loc_upper, str(var_ypd),'NA','\n')))
	
file_out1.close()