import matplotlib.pylab as pt
import numpy
import matplotlib.cm as cm

import scipy.stats
from matplotlib import colors

import matplotlib

matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.size'] = 8.0
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['lines.markersize'] = 3.5
matplotlib.rcParams['lines.linewidth'] = .5
matplotlib.rcParams['axes.linewidth']= .5

##This script compares fitness increments for the same segregant in different environments.
##It tests the hypothesis that segregants with low initial fitness in one condition but high initial fitness in the other
##adapt more quickly in the environment where they had low initial fitness.


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

###Find how far above or below the average initial fitness each segregant is in both environments

init_fit_deviations_sc = (init_fits_sc - numpy.mean(init_fits_sc))/numpy.std(init_fits_sc)
init_fit_deviations_ypd = (init_fits_ypd - numpy.mean(init_fits_ypd))/numpy.std(init_fits_ypd)

###Find how far above or below the average fitness increment each segregant was in both environments (note this is calculated on the average fitness gains across pops descended from each segregant in both environments)

delta_fit_deviations_sc = (delta_fits_sc_means - numpy.mean(delta_fits_sc_means))/numpy.std(delta_fits_sc_means)
delta_fit_deviations_ypd = (delta_fits_ypd_means - numpy.mean(delta_fits_ypd_means))/numpy.std(delta_fits_ypd_means)

###Propagate errors

delta_fits_init_errs = numpy.sqrt( init_std_errs_sc**2/numpy.var(init_fits_sc) + init_std_errs_ypd**2/numpy.var(init_fits_ypd) )

delta_fits_deviations_errs = numpy.sqrt( delta_fits_ypd_vars/(pops_per_seg_ypd*numpy.var(delta_fits_ypd_means)) + delta_fits_sc_vars/(pops_per_seg_sc*numpy.var(delta_fits_sc_means)) )

##Calculate r^2 and bootstrap, including errors, to calculate significance

r = scipy.stats.pearsonr(init_fit_deviations_sc - init_fit_deviations_ypd, delta_fit_deviations_sc - delta_fit_deviations_ypd)[0]

init_diffs = init_fit_deviations_sc - init_fit_deviations_ypd
delta_diffs = delta_fit_deviations_sc - delta_fit_deviations_ypd
r_bs_list = []

n_segs = len(init_diffs)

for i in range(10000):
	chosen_indices = numpy.random.randint(0,n_segs,size=n_segs)
	init_diffs_bs = init_diffs[chosen_indices] + numpy.nanmean(delta_fits_init_errs)*numpy.random.randn(n_segs)
	delta_diffs_bs = delta_diffs[chosen_indices] + numpy.nanmean(delta_fits_deviations_errs)*numpy.random.randn(n_segs)
	r_bs = scipy.stats.pearsonr(init_diffs_bs, delta_diffs_bs)[0]
	r_bs_list.append(r_bs)

r_bs_list = numpy.array(r_bs_list)
#print r_bs_list
p_val = numpy.sum(r_bs_list > 0)/float(10000.)

print p_val
print r**2
pt.figure()
pt.hist(r_bs_list)
pt.show()
###Plot differences in these quantities

fig, ax1 = pt.subplots(1,1,figsize=(4,4))

ax1.errorbar(init_fit_deviations_sc - init_fit_deviations_ypd, delta_fit_deviations_sc - delta_fit_deviations_ypd, xerr= delta_fits_init_errs, yerr=delta_fits_deviations_errs, fmt='o', color= 'k', markeredgewidth=0,alpha=.9,capsize=0, elinewidth=.5)
pt.xlabel('normalized init fit, HT - normalized init fit, OT')
pt.ylabel('normalized delta fit, HT - normalized delta fit, OT')
pt.savefig('delta_fitness_seg_by_seg_comparison.pdf',bbox_inches='tight')

##Make some sanity check plots

##Look at this for each Kre33 allele group separately.
##Plot up the normalized initial fitnesses 

fig, (ax1, ax2, ax3) = pt.subplots(1,3,figsize=(12,4))

ax1.errorbar(init_fit_deviations_ypd[rm_allele], init_fit_deviations_sc[rm_allele],  yerr=(init_std_errs_sc**2/numpy.var(init_fits_sc))[rm_allele], xerr=(init_std_errs_ypd**2/numpy.var(init_fits_ypd))[rm_allele], fmt='o', color= 'Grey', markeredgewidth=0,alpha=.9,capsize=0, elinewidth=.5)
ax1.plot([-2,2],[-2,2],'k')
ax1.set_ylabel('Initial fitness at HT (normalized)',fontsize=12)
ax1.set_xlabel('Initial fitness at OT (normalized)',fontsize=12)

ax1.errorbar(init_fit_deviations_ypd[by_allele], init_fit_deviations_sc[by_allele],  yerr=(init_std_errs_sc**2/numpy.var(init_fits_sc))[by_allele], xerr=(init_std_errs_ypd**2/numpy.var(init_fits_ypd))[by_allele], fmt='o', color= 'Black', markeredgewidth=0,alpha=.9,capsize=0, elinewidth=.5)
#ax1.axis('equal')

ax2.errorbar(delta_fit_deviations_ypd[rm_allele], delta_fit_deviations_sc[rm_allele],  yerr=(delta_fits_sc_vars/(pops_per_seg_sc*numpy.var(delta_fits_sc_means)))[rm_allele], xerr=(delta_fits_ypd_vars/(pops_per_seg_ypd*numpy.var(delta_fits_ypd_means)))[rm_allele], fmt='o', color= 'Grey', markeredgewidth=0,alpha=.9,capsize=0, elinewidth=.5)
ax2.errorbar(delta_fit_deviations_ypd[by_allele], delta_fit_deviations_sc[by_allele],  yerr=(delta_fits_sc_vars/(pops_per_seg_sc*numpy.var(delta_fits_sc_means)))[by_allele], xerr=(delta_fits_ypd_vars/(pops_per_seg_ypd*numpy.var(delta_fits_ypd_means)))[by_allele], fmt='o', color= 'Black', markeredgewidth=0,alpha=.9,capsize=0, elinewidth=.5)

ax2.plot([-2,2],[-2,2],'k')
ax2.set_ylabel('Fitness gain at HT (normalized)',fontsize=12)
ax2.set_xlabel('Fitness gain at OT (normalized)', fontsize=12)
#ax2.axis('equal')

ax3.errorbar(init_fit_deviations_sc[rm_allele] - init_fit_deviations_ypd[rm_allele], delta_fit_deviations_sc[rm_allele] - delta_fit_deviations_ypd[rm_allele], xerr= delta_fits_init_errs[rm_allele], yerr=delta_fits_deviations_errs[rm_allele], fmt='o', color= 'Grey', markeredgewidth=0,alpha=.9,capsize=0, elinewidth=.5)
ax3.errorbar(init_fit_deviations_sc[by_allele] - init_fit_deviations_ypd[by_allele], delta_fit_deviations_sc[by_allele] - delta_fit_deviations_ypd[by_allele], xerr= delta_fits_init_errs[by_allele], yerr=delta_fits_deviations_errs[by_allele], fmt='o', color= 'Black', markeredgewidth=0,alpha=.9,capsize=0, elinewidth=.5)
#ax3.axis('equal')
pt.xlabel('Initial fitness at HT - initial fitness at OT',fontsize=12)
pt.ylabel('Fitness gain at HT - fitness gain at OT', fontsize=12)
pt.savefig('delta_fitness_seg_by_seg_comparison2.pdf',bbox_inches='tight')

##Normalize separately for the two Kre33 allele groups

###Find how far above or below the average initial fitness each segregant is in both environments

init_fit_deviations_sc_rm = (init_fits_sc[rm_allele] - numpy.mean(init_fits_sc[rm_allele]))/numpy.std(init_fits_sc[rm_allele])
init_fit_deviations_ypd_rm = (init_fits_ypd[rm_allele] - numpy.mean(init_fits_ypd[rm_allele]))/numpy.std(init_fits_ypd[rm_allele])

init_fit_deviations_sc_by = (init_fits_sc[by_allele] - numpy.mean(init_fits_sc[by_allele]))/numpy.std(init_fits_sc[by_allele])
init_fit_deviations_ypd_by = (init_fits_ypd[by_allele] - numpy.mean(init_fits_ypd[by_allele]))/numpy.std(init_fits_ypd[by_allele])

###Find how far above or below the average fitness increment each segregant was in both environments (note this is calculated on the average fitness gains across pops descended from each segregant in both environments)

delta_fit_deviations_sc_rm = (delta_fits_sc_means[rm_allele] - numpy.mean(delta_fits_sc_means[rm_allele]))/numpy.std(delta_fits_sc_means[rm_allele])
delta_fit_deviations_ypd_rm = (delta_fits_ypd_means[rm_allele] - numpy.mean(delta_fits_ypd_means[rm_allele]))/numpy.std(delta_fits_ypd_means[rm_allele])

delta_fit_deviations_sc_by = (delta_fits_sc_means[by_allele] - numpy.mean(delta_fits_sc_means[by_allele]))/numpy.std(delta_fits_sc_means[by_allele])
delta_fit_deviations_ypd_by = (delta_fits_ypd_means[by_allele] - numpy.mean(delta_fits_ypd_means[by_allele]))/numpy.std(delta_fits_ypd_means[by_allele])

###Propagate errors

delta_fits_init_errs_rm = numpy.sqrt( init_std_errs_sc[rm_allele]**2/numpy.var(init_fits_sc[rm_allele]) + init_std_errs_ypd[rm_allele]**2/numpy.var(init_fits_ypd[rm_allele]) )
delta_fits_deviations_errs_rm = numpy.sqrt( delta_fits_ypd_vars[rm_allele]/(pops_per_seg_ypd[rm_allele]*numpy.var(delta_fits_ypd_means[rm_allele])) + delta_fits_sc_vars[rm_allele]/(pops_per_seg_sc[rm_allele]*numpy.var(delta_fits_sc_means[rm_allele])) )
delta_fits_init_errs_by = numpy.sqrt( init_std_errs_sc[by_allele]**2/numpy.var(init_fits_sc[by_allele]) + init_std_errs_ypd[by_allele]**2/numpy.var(init_fits_ypd[by_allele]) )
delta_fits_deviations_errs_by = numpy.sqrt( delta_fits_ypd_vars[by_allele]/(pops_per_seg_ypd[by_allele]*numpy.var(delta_fits_ypd_means[by_allele])) + delta_fits_sc_vars[by_allele]/(pops_per_seg_sc[by_allele]*numpy.var(delta_fits_sc_means[by_allele])) )

fig, (ax1,ax2,ax3) = pt.subplots(1,3,figsize=(12,4))
colors_rm = (init_fit_deviations_sc_rm - init_fit_deviations_ypd_rm) - numpy.min( init_fit_deviations_sc_rm - init_fit_deviations_ypd_rm )
colors_by = (init_fit_deviations_sc_by - init_fit_deviations_ypd_by) - numpy.min( init_fit_deviations_sc_by - init_fit_deviations_ypd_by )

cnorm_rm= colors_rm/numpy.max((numpy.max(colors_rm),numpy.max(colors_by)))
cnorm_by= colors_by/numpy.max((numpy.max(colors_rm),numpy.max(colors_by)))

ax1.scatter(init_fit_deviations_ypd_rm, init_fit_deviations_sc_rm, c = cnorm_rm, cmap='PRGn', edgecolors='Grey',linewidths=.5)
ax1.plot([-2,2],[-2,2],'k')
ax1.set_ylabel('Initial fitness at HT (normalized)',fontsize=12)
ax1.set_xlabel('Initial fitness at OT (normalized)',fontsize=12)
ax1.scatter(init_fit_deviations_ypd_by, init_fit_deviations_sc_by, c = cnorm_by, cmap='PRGn', edgecolors='Grey',linewidths=.5)
# #ax1.axis('equal')


# 
ax2.scatter(delta_fit_deviations_ypd_rm, delta_fit_deviations_sc_rm, c = cnorm_rm, cmap='PRGn', edgecolors='Grey',linewidths=.5)
ax2.scatter(delta_fit_deviations_ypd_by, delta_fit_deviations_sc_by, c = cnorm_by, cmap='PRGn', edgecolors='Grey',linewidths=.5)
# 
ax2.plot([-2,2],[-2,2],'k')
ax2.set_ylabel('Fitness gain at HT (normalized)',fontsize=12)
ax2.set_xlabel('Fitness gain at OT (normalized)', fontsize=12)
# #ax2.axis('equal')

ax3.errorbar(init_fit_deviations_sc_rm - init_fit_deviations_ypd_rm, delta_fit_deviations_sc_rm - delta_fit_deviations_ypd_rm, xerr= delta_fits_init_errs_rm, yerr=delta_fits_deviations_errs_rm, fmt='o', color= 'Black', markeredgewidth=0,alpha=.9,capsize=0, elinewidth=.5)
ax3.errorbar(init_fit_deviations_sc_by - init_fit_deviations_ypd_by, delta_fit_deviations_sc_by - delta_fit_deviations_ypd_by, xerr= delta_fits_init_errs_by, yerr=delta_fits_deviations_errs_by, fmt='o', color= 'Black', markeredgewidth=0,alpha=.9,capsize=0, elinewidth=.5)
init_fit_diffs_rm_by = numpy.concatenate((init_fit_deviations_sc_rm - init_fit_deviations_ypd_rm, init_fit_deviations_sc_by - init_fit_deviations_ypd_by), axis=0)

delta_fit_diffs_rm_by = numpy.concatenate((delta_fit_deviations_sc_rm - delta_fit_deviations_ypd_rm, delta_fit_deviations_sc_by - delta_fit_deviations_ypd_by), axis=0)
r_rm_by = scipy.stats.pearsonr(init_fit_diffs_rm_by, delta_fit_diffs_rm_by)[0]
#ax3.axis('equal')
pt.xlabel('Initial fitness at HT - initial fitness at OT',fontsize=12)
pt.ylabel('Fitness gain at HT - fitness gain at OT', fontsize=12)
r2_rm_by = numpy.round(r_rm_by**2, 2)
pt.text(.7,3,("").join(("$r^2=$",str(r2_rm_by))),fontsize=12)
pt.savefig('delta_fitness_seg_by_seg_comparison_KRE33_control.pdf',bbox_inches='tight')

fig, (ax1,ax2,ax3) = pt.subplots(1,3,figsize=(12,4))
###Now let's try to make a plot colored by y-x in panel A (i.e. normalized fitness at OT - normalized fitness at HT)
colors = (init_fit_deviations_sc - init_fit_deviations_ypd) - numpy.min( init_fit_deviations_sc - init_fit_deviations_ypd )
cnorm = colors/numpy.max(colors)
print cnorm
ax1.scatter(init_fit_deviations_ypd, init_fit_deviations_sc, c= cnorm, cmap='PRGn', edgecolors='Grey',linewidths=.5)#,capsize=0, elinewidth=.5)fmt='o', yerr=(init_std_errs_sc**2/numpy.var(init_fits_sc)), xerr=(init_std_errs_ypd**2/numpy.var(init_fits_ypd)),
ax1.plot([-2,2],[-2,2],'k')
ax1.set_ylabel('Initial fitness at HT (normalized)',fontsize=12)
ax1.set_xlabel('Initial fitness at OT (normalized)',fontsize=12)

ax2.scatter(delta_fit_deviations_ypd, delta_fit_deviations_sc,  c= cnorm, cmap='PRGn', edgecolors='Grey', linewidths=.5)#,capsize=0, elinewidth=.5)yerr=(delta_fits_sc_vars/(pops_per_seg_sc*numpy.var(delta_fits_sc_means))), xerr=(delta_fits_ypd_vars/(pops_per_seg_ypd*numpy.var(delta_fits_ypd_means))), 

ax2.plot([-2,2],[-2,2],'k')
ax2.set_ylabel('Fitness gain at HT (normalized)',fontsize=12)
ax2.set_xlabel('Fitness gain at OT (normalized)', fontsize=12)

ax3.errorbar(init_fit_deviations_sc - init_fit_deviations_ypd, delta_fit_deviations_sc - delta_fit_deviations_ypd, xerr= delta_fits_init_errs, yerr=delta_fits_deviations_errs, fmt='o', color= 'k', markeredgewidth=0,alpha=.9,capsize=0, elinewidth=.5)
r2 = numpy.round(r**2, 2)
pt.text(.5,2,("").join(("$r^2=$",str(r2))),fontsize=12)
#ax3.axis('equal')
pt.xlabel('Initial fitness at HT - initial fitness at OT',fontsize=12)
pt.ylabel('Fitness gain at HT - fitness gain at OT', fontsize=12)
pt.savefig('delta_fitness_seg_by_seg_comparison_colored.pdf',bbox_inches='tight')