import matplotlib.pylab as pt
import numpy
import matplotlib.cm as cm
import qtl_detection_adaptability
from calculate_narrow_sense_heritability import narrow_sense_hsq

from calculate_narrow_sense_heritability import narrow_sense_hsq_replicates

from calculate_narrow_sense_heritability import rsq
from calculate_narrow_sense_heritability import rsq_linear_model
from calculate_narrow_sense_heritability import broad_sense_Hsq

from calculate_narrow_sense_heritability import broad_sense_Hsq_init_fits
import scipy.stats
from matplotlib import colors

import matplotlib

matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.size'] = 8.0
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['lines.markersize'] = 3.5
matplotlib.rcParams['lines.linewidth'] = .5
matplotlib.rcParams['axes.linewidth']= .5

##This script compares different models of initial fitness and delta fitness, and calculates narrow-sense heritability.
##The script was updated on 10/26/2016 to bootstrap jointly for all statistics (i.e. to calculate each statistic of interest on each
##resampled set of segregants). This allows us to directly test whether certain models are better than others at particular significance levels,
##and also it's faster than bootstrapping for each statistic separately.

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

delta_fits_sc_in_ypd_vars = numpy.dot((delta_fits_sc_in_ypd - numpy.dot(helper_matrix_sc_pops,delta_fits_sc_in_ypd_means))**2, helper_matrix_sc_pops)/(pops_per_seg_sc - 1.) + init_std_errs_ypd**2
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

##Significant loci above fitness

file1 = open('data/Above_fitness_epistatic_qtl_model_parameters_epistatic_2_23_2017.txt','r')

X_fitness_and_qtls_sc = numpy.concatenate(( numpy.ones((n_segs,1)),  init_fits_sc.reshape((n_segs,1))), axis = 1)
X_fitness_and_qtls_ypd = numpy.concatenate((numpy.ones((n_segs,1)),  init_fits_ypd.reshape((n_segs,1))), axis = 1)

header_blocks = 0
for line in file1:
	if line.startswith('#'):
		header_blocks += 1
	elif (header_blocks > 3.5 and header_blocks < 5.5):
		locus = int(line.split(',')[1])
		if line.split(',')[0] == 'epistatic with Kre33':
			
			X_fitness_and_qtls_sc = numpy.concatenate( (X_fitness_and_qtls_sc, ((genotype_mat[:,3777] - .5)*(genotype_mat[:,locus] - .5)).reshape((n_segs,1))), axis = 1)
		else:
			X_fitness_and_qtls_sc = numpy.concatenate( (X_fitness_and_qtls_sc, (genotype_mat[:,locus] - .5).reshape((n_segs,1))), axis = 1)
	elif header_blocks > 7.5:
		locus = int(line.split(',')[1])
		X_fitness_and_qtls_ypd = numpy.concatenate( (X_fitness_and_qtls_ypd, (genotype_mat[:,locus] - .5).reshape((n_segs,1))), axis = 1)

####

n_segs_by = numpy.sum(by_allele)
n_segs_rm = numpy.sum(rm_allele)
n_pops_by_sc = numpy.sum(by_allele_pops_sc)
n_pops_rm_sc = numpy.sum(rm_allele_pops_sc)
n_pops_by_ypd = numpy.sum(by_allele_pops_ypd)
n_pops_rm_ypd = numpy.sum(rm_allele_pops_ypd)

#Set up lists of statistics to be calculated each bootstrap round

H2_init_fit_sc_bs = []
H2_init_fit_ypd_bs = []
h2_init_fit_sc_bs = []
h2_init_fit_ypd_bs = []
r2_init_fit_sc_bs = []
r2_init_fit_ypd_bs = []

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

#Dependent variables for the actual set of segregants

X_qtls_init_sc = numpy.append(genotype_mat[:,init_fit_qtls_sc], numpy.ones((n_segs,1)), axis = 1)
X_qtls_delta_sc_segs = numpy.append(genotype_mat[:,delta_fit_qtls_sc], numpy.ones((n_segs,1)), axis = 1)

X_qtls_init_ypd = numpy.append(genotype_mat[:,init_fit_qtls_ypd], numpy.ones((n_segs,1)), axis = 1)
X_qtls_delta_ypd_segs = numpy.append(genotype_mat[:,delta_fit_qtls_ypd], numpy.ones((n_segs,1)), axis = 1)

#print X_qtls_delta[0:10,:]

X_init_fit_sc = numpy.append(init_fits_sc.reshape((n_segs,1)), numpy.ones((n_segs,1)), axis = 1)
X_init_fit_ypd = numpy.append(init_fits_ypd.reshape((n_segs,1)), numpy.ones((n_segs,1)), axis = 1)

X_init_fit_sc_pops = numpy.concatenate((numpy.dot(helper_matrix_sc_pops, init_fits_sc).reshape((num_pops_sc,1)), numpy.ones((num_pops_sc,1))),axis=1)
X_init_fit_ypd_pops = numpy.concatenate((numpy.dot(helper_matrix_ypd_pops, init_fits_ypd).reshape((num_pops_ypd,1)), numpy.ones((num_pops_ypd,1))),axis=1)

###

beta_qtl_ypd_delta_means = numpy.dot ( numpy.linalg.inv( numpy.dot( X_qtls_delta_ypd_segs.T, X_qtls_delta_ypd_segs) ), numpy.dot(X_qtls_delta_ypd_segs.T, delta_fits_ypd_means) )
beta_qtl_sc_delta_means = numpy.dot ( numpy.linalg.inv( numpy.dot( X_qtls_delta_sc_segs.T, X_qtls_delta_sc_segs) ), numpy.dot(X_qtls_delta_sc_segs.T, delta_fits_sc_means) )

beta_fitness_and_qtls_ypd = numpy.dot ( numpy.linalg.inv( numpy.dot( X_fitness_and_qtls_ypd.T, X_fitness_and_qtls_ypd) ), numpy.dot(X_fitness_and_qtls_ypd.T, delta_fits_ypd_means) )
beta_fitness_and_qtls_sc = numpy.dot ( numpy.linalg.inv( numpy.dot( X_fitness_and_qtls_sc.T, X_fitness_and_qtls_sc) ), numpy.dot(X_fitness_and_qtls_sc.T, delta_fits_sc_means) )

beta_init_ypd = numpy.dot ( numpy.linalg.inv( numpy.dot( X_init_fit_ypd.T, X_init_fit_ypd) ), numpy.dot(X_init_fit_ypd.T, delta_fits_ypd_means) )
beta_init_sc = numpy.dot ( numpy.linalg.inv( numpy.dot( X_init_fit_sc.T, X_init_fit_sc) ), numpy.dot(X_init_fit_sc.T, delta_fits_sc_means) )

qtl_predictor_linear_ypd_unscaled = numpy.dot(X_qtls_delta_ypd_segs, beta_qtl_ypd_delta_means)
qtl_predictor_linear_ypd = (qtl_predictor_linear_ypd_unscaled - numpy.mean(qtl_predictor_linear_ypd_unscaled))/numpy.std(qtl_predictor_linear_ypd_unscaled)
qtl_predictor_linear_sc_unscaled = numpy.dot(X_qtls_delta_sc_segs, beta_qtl_sc_delta_means)
qtl_predictor_linear_sc = (qtl_predictor_linear_sc_unscaled - numpy.mean(qtl_predictor_linear_sc_unscaled))/numpy.std(qtl_predictor_linear_sc_unscaled)

rsq_qtls_ypd_delta_means = scipy.stats.pearsonr(qtl_predictor_linear_ypd, delta_fits_ypd_means)[0]**2
rsq_qtls_sc_delta_means = scipy.stats.pearsonr(qtl_predictor_linear_sc, delta_fits_sc_means)[0]**2


delta_fits_ypd_init_fit_prediction = numpy.dot( X_init_fit_ypd, beta_init_ypd )
delta_fits_sc_init_fit_prediction = numpy.dot( X_init_fit_sc, beta_init_sc )


rsq_fit_ypd_for_means = scipy.stats.pearsonr( init_fits_ypd, delta_fits_ypd_means )[0]**2
rsq_fit_sc_for_means = scipy.stats.pearsonr( init_fits_sc, delta_fits_sc_means )[0]**2

predictor_fit_and_qtls_ypd_unscaled = numpy.dot( X_fitness_and_qtls_ypd, beta_fitness_and_qtls_ypd )
predictor_fit_and_qtls_ypd = (predictor_fit_and_qtls_ypd_unscaled - numpy.mean(predictor_fit_and_qtls_ypd_unscaled))/numpy.std(predictor_fit_and_qtls_ypd_unscaled)
predictor_fit_and_qtls_sc_unscaled = numpy.dot( X_fitness_and_qtls_sc, beta_fitness_and_qtls_sc )
predictor_fit_and_qtls_sc = (predictor_fit_and_qtls_sc_unscaled - numpy.mean(predictor_fit_and_qtls_sc_unscaled))/numpy.std(predictor_fit_and_qtls_sc_unscaled)

rsq_full_model_ypd_for_means = scipy.stats.pearsonr(predictor_fit_and_qtls_ypd, delta_fits_ypd_means)[0]**2
rsq_full_model_sc_for_means = scipy.stats.pearsonr(predictor_fit_and_qtls_sc, delta_fits_sc_means)[0]**2

fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = pt.subplots(2,3,figsize=(8,4))

ax1.set_title('QTL models')
ax2.set_title('Fitness models')


#ax4.set_title('QTL model, 37 C populations',fontsize=9)
#ax5.set_title('Fitness model, 37 C populations',fontsize=9)

ax1.set_xlabel('best linear predictor',fontsize=8)
ax4.set_xlabel('best linear predictor',fontsize=8)
ax2.set_xlabel('founder fitness, OT',fontsize=8)
ax5.set_xlabel('founder fitness, HT',fontsize=8)
#ax1.set_ylabel('fitness gains (%), 30 C',fontsize=8)
#ax2.set_ylabel('fitness gains (%), 30 C',fontsize=8)
#ax4.set_ylabel('fitness gains (%), 37 C',fontsize=8)
#ax5.set_ylabel('fitness gains (%), 37 C',fontsize=8)
ax1.errorbar(-1*qtl_predictor_linear_ypd, delta_fits_ypd_means, yerr = delta_fits_ypd_std_errs, fmt='o', color= 'MediumSlateBlue', markeredgewidth=0,capsize=0,alpha=.9,elinewidth=.5)
ax1.plot(-1*qtl_predictor_linear_ypd, qtl_predictor_linear_ypd_unscaled, 'k')
ax1.text(1,.15,'$r^2=$' + str(round(rsq_qtls_ypd_delta_means, 2)),fontsize=8)

ax4.errorbar(-1*qtl_predictor_linear_sc, delta_fits_sc_means, yerr = delta_fits_sc_std_errs, fmt='o', color= 'Tomato', markeredgewidth=0,capsize=0,alpha=.9,elinewidth=.5)
ax4.plot(-1*qtl_predictor_linear_sc, qtl_predictor_linear_sc_unscaled, 'k')
ax4.text(1,.3,'$r^2=$' + str(round(rsq_qtls_sc_delta_means, 2)),fontsize=8)
ax1.set_ylim(0,.18)
ax2.errorbar(init_fits_ypd[rm_allele], delta_fits_ypd_means[rm_allele], xerr = init_std_errs_ypd[rm_allele], yerr = delta_fits_ypd_std_errs[rm_allele], fmt='o', color= 'MediumSlateBlue', markeredgewidth=0,capsize=0,alpha=.9,elinewidth=.5)
ax2.errorbar(init_fits_ypd[by_allele], delta_fits_ypd_means[by_allele], xerr = init_std_errs_ypd[by_allele], yerr = delta_fits_ypd_std_errs[by_allele], fmt='o', color= 'DarkSlateBlue', markeredgewidth=0,capsize=0,alpha=.9,elinewidth=.5)
ax2.plot(init_fits_ypd, delta_fits_ypd_init_fit_prediction, 'k')
ax2.text(.05,.15,'$r^2=$' + str(round(rsq_fit_ypd_for_means, 2)),fontsize=8)

ax2.set_xlim(-.17,.13)
ax2.set_ylim(0,.18)
ax1.set_yticks(numpy.arange(0,.19,.02))
ax1.set_yticklabels(numpy.arange(0,19,2),fontsize=8)
ax2.set_yticks(numpy.arange(0,.19,.02))
ax2.set_yticklabels(numpy.arange(0,19,2),fontsize=8)
ax3.set_yticks(numpy.arange(0,.19,.02))
ax3.set_yticklabels(numpy.arange(0,19,2),fontsize=8)

ax1.set_xticks(numpy.arange(-2,2.1,1))
ax1.set_xticklabels(numpy.arange(-2,2.1,1))
ax4.set_xticks(numpy.arange(-2,2.1,1))
ax4.set_xticklabels(numpy.arange(-2,2.1,1))

ax5.set_ylim(0,.38)
ax5.set_ylim(0,.38)
ax5.errorbar(init_fits_sc[rm_allele], delta_fits_sc_means[rm_allele], xerr = init_std_errs_sc[rm_allele], yerr = delta_fits_sc_std_errs[rm_allele], fmt='o', color= 'Tomato', markeredgewidth=0,capsize=0,alpha=.9,elinewidth=.5)
ax5.errorbar(init_fits_sc[by_allele], delta_fits_sc_means[by_allele], xerr = init_std_errs_sc[by_allele], yerr = delta_fits_sc_std_errs[by_allele], fmt='o', color= 'Brown', markeredgewidth=0,capsize=0,alpha=.9,elinewidth=.5)
ax5.plot(init_fits_sc, delta_fits_sc_init_fit_prediction, 'k')
ax5.text(.05,.3,'$r^2=$' + str(round(rsq_fit_sc_for_means, 2)),fontsize=8)

ax1.set_xlim(-2,2.2)
ax4.set_xlim(-2,2.2)
ax3.set_xlim(-2.4,2.5)
ax6.set_xlim(-2.5,2.5)

ax4.set_yticks(numpy.arange(0,.39,.05))
ax4.set_yticklabels(numpy.arange(0,39,5),fontsize=8)
ax5.set_yticks(numpy.arange(0,.39,.05))
ax5.set_yticklabels(numpy.arange(0,39,5),fontsize=8)
ax6.set_yticks(numpy.arange(0,.39,.05))
ax6.set_yticklabels(numpy.arange(0,39,5),fontsize=8)

ax3.errorbar(-1*predictor_fit_and_qtls_ypd[rm_allele], delta_fits_ypd_means[rm_allele], xerr = init_std_errs_ypd[rm_allele], yerr = delta_fits_ypd_std_errs[rm_allele], fmt='o', color= 'MediumSlateBlue', markeredgewidth=0,capsize=0,alpha=.9,elinewidth=.5)
ax3.errorbar(-1*predictor_fit_and_qtls_ypd[by_allele], delta_fits_ypd_means[by_allele], xerr = init_std_errs_ypd[by_allele], yerr = delta_fits_ypd_std_errs[by_allele], fmt='o', color= 'DarkSlateBlue', markeredgewidth=0,capsize=0,alpha=.9,elinewidth=.5)
ax3.plot(-1*predictor_fit_and_qtls_ypd, predictor_fit_and_qtls_ypd_unscaled, 'k')
ax3.text(1,.15,'$r^2=$' + str(round(rsq_full_model_ypd_for_means, 2)),fontsize=8)

ax6.errorbar(-1*predictor_fit_and_qtls_sc[rm_allele], delta_fits_sc_means[rm_allele], xerr = init_std_errs_sc[rm_allele], yerr = delta_fits_sc_std_errs[rm_allele], fmt='o', color= 'Tomato', markeredgewidth=0,capsize=0,alpha=.9,elinewidth=.5)
ax6.errorbar(-1*predictor_fit_and_qtls_sc[by_allele], delta_fits_sc_means[by_allele], xerr = init_std_errs_sc[by_allele], yerr = delta_fits_sc_std_errs[by_allele], fmt='o', color= 'Brown', markeredgewidth=0,capsize=0,alpha=.9,elinewidth=.5)
ax6.plot(-1*predictor_fit_and_qtls_sc, predictor_fit_and_qtls_sc_unscaled, 'k')
ax6.text(1,.3,'$r^2=$' + str(round(rsq_full_model_sc_for_means, 2)),fontsize=8)

ax3.set_ylim(0,.18)
ax6.set_ylim(0,.38)

ax1.text(-3,.15,'Evolved at OT',rotation=90,fontsize=10)
ax4.text(-3,.30,'Evolved at HT',rotation=90,fontsize=10)
ax1.set_ylabel('Fitness gains at OT (%)')
ax4.set_ylabel('Fitness gains at HT (%)')


ax3.set_title('Combined models')
#ax6.set_title('Combined model, 37 C populations',fontsize=9)
ax3.set_xlabel('best predictor',fontsize=8)
#ax3.set_ylabel('fitness gains (%), 30 C',fontsize=8)
ax6.set_xlabel('best predictor')
#ax6.set_ylabel('fitness gains (%), 37 C',fontsize=8)
matplotlib.rcParams['font.sans-serif'] = 'Helvetica'
ax1.text(-2.4, .2, 'A', fontsize=16)
ax4.text(-2.4, .4, 'D', fontsize=16)
ax2.text(-.22, .2, 'B', fontsize=16)
ax3.text(-2.6, .2, 'C', fontsize=16)
ax5.text(-.34, .4, 'E', fontsize=16)
ax6.text(-2.6, .4, 'F', fontsize=16)

pt.tight_layout()
pt.savefig('delta_fitness_predictors_2_5_2017.pdf',bbox_inches='tight')