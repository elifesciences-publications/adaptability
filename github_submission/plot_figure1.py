import matplotlib.pylab as pt
import numpy
import matplotlib.cm as cm
import scipy.stats
from matplotlib import colors
import matplotlib

matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.size'] = 10.0
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['lines.markersize'] = 1.7
matplotlib.rcParams['lines.linewidth'] = .3
matplotlib.rcParams['legend.fontsize'] = 8.0
matplotlib.rcParams['legend.markerscale'] = .7
matplotlib.rcParams['legend.frameon'] = False
matplotlib.rcParams['axes.linewidth']= .5
matplotlib.rcParams['patch.linewidth']= .5

#Import fitness and genotype data
filename1 = 'data/fitness_measurements_with_population_names_12_29_2016.csv'
filename2 = 'data/control_replicate_measurements.csv'
filename3 = 'data/segregant_genotypes_12_29_2016.csv'

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
rm_allele = numpy.array(genotype_mat[:,9596],dtype='Bool')
by_allele = numpy.array(1 - genotype_mat[:,9596],dtype='Bool')

##Order by number of RM loci

rm_content = numpy.sum(genotype_mat, axis = 1)
print rm_content[0:10]

rm_order = numpy.argsort(rm_content) ##segregants ordered based on the fraction of their genome that is RM as opposed to BY--from less to more RM

#Symmetrize genotype_mat
genotype_mat = genotype_mat - .5
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
delta_fits_sc_vars = numpy.dot((delta_fits_sc - numpy.dot(helper_matrix_sc_pops,delta_fits_sc_means))**2, helper_matrix_sc_pops)/(pops_per_seg_sc - 1.) + init_std_errs_sc**2 #- measurement_error_var_sc
delta_fits_sc_std_errs = numpy.sqrt(delta_fits_sc_vars/pops_per_seg_sc)

delta_fits_ypd_vars = numpy.dot((delta_fits_ypd - numpy.dot(helper_matrix_ypd_pops,delta_fits_ypd_means))**2, helper_matrix_ypd_pops)/(pops_per_seg_ypd - 1.) + init_std_errs_ypd**2 #- measurement_error_var_ypd
delta_fits_ypd_std_errs = numpy.sqrt(delta_fits_ypd_vars/pops_per_seg_ypd)

delta_fits_ypd_in_sc_vars = numpy.dot((delta_fits_ypd_in_sc - numpy.dot(helper_matrix_ypd_pops,delta_fits_ypd_in_sc_means))**2, helper_matrix_ypd_pops)/(pops_per_seg_ypd - 1.) + init_std_errs_sc**2 #- measurement_error_sc
delta_fits_ypd_in_sc_std_errs = numpy.sqrt(delta_fits_ypd_in_sc_vars/pops_per_seg_ypd)

delta_fits_sc_in_ypd_vars = numpy.dot((delta_fits_sc_in_ypd - numpy.dot(helper_matrix_sc_pops,delta_fits_sc_in_ypd_means))**2, helper_matrix_sc_pops)/(pops_per_seg_sc - 1.) + init_std_errs_ypd**2 #- measurement_error_ypd
delta_fits_sc_in_ypd_std_errs = numpy.sqrt(delta_fits_sc_in_ypd_vars/pops_per_seg_sc)

##Plot the delta fitness for YPD pops and SC pops at home
##dummy figure for legends

fig, ((ax1, ax3),(ax4,ax2)) = pt.subplots(2,2,figsize=(8,4), sharex=True)
for i in range(num_segs):
	seg_idx = rm_order[i]
	pop_idxs = numpy.where(helper_matrix_ypd_pops[:,seg_idx]==1)
	#print pop_idxs
	delta_fits_ypd_these_pops = delta_fits_ypd[pop_idxs]
	#print delta_fits_ypd_these_pops
	bottom = numpy.min(delta_fits_ypd_these_pops)
	top = numpy.max(delta_fits_ypd_these_pops)
	
	ax1.plot([i,i],[bottom,top], color='MediumSlateBlue',alpha=.9)
	for point in delta_fits_ypd_these_pops:
		ax1.plot(i, point, color='MediumSlateBlue', alpha=.9, marker='o', markeredgewidth=0)
	
	pop_idxs_sc = numpy.where(helper_matrix_sc_pops[:,seg_idx]==1)
	delta_fits_sc_these_pops = delta_fits_sc[pop_idxs_sc]
	bottom = numpy.min(delta_fits_sc_these_pops)
	top = numpy.max(delta_fits_sc_these_pops)
	
	ax2.plot([i,i],[bottom,top], color='Tomato',alpha=.9)
	for point in delta_fits_sc_these_pops:
		ax2.plot(i, point, color='Tomato', alpha=.9, marker='o', markeredgewidth=0)
	
	delta_fits_sc_in_ypd_these_pops = delta_fits_sc_in_ypd[pop_idxs_sc]
	
	bottom = numpy.min(delta_fits_sc_in_ypd_these_pops)
	top = numpy.max(delta_fits_sc_in_ypd_these_pops)
	
	ax3.plot([i,i],[bottom,top], color='Tomato',alpha=.9)
	for point in delta_fits_sc_in_ypd_these_pops:
		ax3.plot(i, point, color='Tomato', alpha=.9, marker='o', markeredgewidth=0)
	
	delta_fits_ypd_in_sc_these_pops = delta_fits_ypd_in_sc[pop_idxs]
	
	bottom = numpy.min(delta_fits_ypd_in_sc_these_pops)
	top = numpy.max(delta_fits_ypd_in_sc_these_pops)
	
	ax4.plot([i,i],[bottom,top], color='MediumSlateBlue',alpha=.9)
	for point in delta_fits_ypd_in_sc_these_pops:
		ax4.plot(i, point, color='MediumSlateBlue', alpha=.9, marker='o', markeredgewidth=0)
		
ax1.set_ylim(-.1,.25)
ax1.set_xlim(-2,num_segs+2)
#ax1.set_xlabel('Founder (ordered by % RM)')
ax1.set_ylabel('Fitness gain at OT (%)')
ax1.set_yticks(numpy.arange(-.1,.26,.05))
ax1.set_yticklabels(numpy.arange(-10,26,5))
ax1.set_title('Evolved at OT')
#ax1.legend(dot1,'Evolved at 30 C', numpoints = 1)

ax1.axhline(0,color='k')

ax3.set_ylim(-.1,.25)
ax3.set_xlim(-2,num_segs+2)
#ax1.set_xlabel('Founder (ordered by % RM)')
#ax1.set_ylabel('Fitness gain at 30C (%)')
ax3.set_yticks(numpy.arange(-.1,.26,.05))
ax3.set_yticklabels(numpy.arange(-10,26,5))
#ax3.set_ylabel('Fitness gain at 30C (%)')
ax3.set_title('Evolved at HT')
#ax3.legend(dot2,'Evolved at 37 C', numpoints = 1)

#ax3.legend((dot1,dot2),('Evolved at 30 C','Evolved at 37 C'), numpoints = 1, ncol=2, fontsize=10, loc=(-.2,-.2))
ax3.axhline(0,color='k')

ax2.set_ylim(-.2,.45)
ax2.set_xlim(-2,num_segs+2)
ax2.set_xlabel('Founder index (ordered by % RM)')

ax2.set_yticks(numpy.arange(-.2,.41,.1))
ax2.set_yticklabels(numpy.arange(-20,41,10))

#ax2.set_title('Evolved at 37C')

ax2.axhline(0,color='k')

ax4.set_ylim(-.2,.45)
ax4.set_xlim(-2,num_segs+2)
ax4.set_xlabel('Founder index (ordered by % RM)')
#ax4.set_ylabel('Fitness gain at 37C (%)')
ax4.set_yticks(numpy.arange(-.2,.41,.1))
ax4.set_yticklabels(numpy.arange(-20,41,10))
ax4.set_ylabel('Fitness gain at HT (%)')
#ax4.set_title('Evolved at 30C')
ax4.axhline(0,color='k')
#ax4.legend(dot1,'Evolved at 30 C', numpoints = 1)

pt.tight_layout
pt.savefig('Variation_in_delta_fitness_outcomes_by_founder_four_by_four_12_29_2016.pdf',bbox_inches='tight')