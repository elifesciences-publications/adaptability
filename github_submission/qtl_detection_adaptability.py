import numpy
import numpy.random

from math import log
#import parse_file_2_1_15
import sys

import matplotlib.pylab as pylab

def calculate_lods(genotype_matrix, phenotype_matrix):
    
    # Calculate LODS for the entire genome
    
    # first calculate means
    
    yavgs = phenotype_matrix[:,0]
    y2s = phenotype_matrix[:,1]
    ns = phenotype_matrix[:,2]
    ntot = ns.sum() #4 replicates per segregant (note that this only affects the scale of the peaks, which is relevant for the confidence intervals but not for identifying QTLs)
    
    Yavg = numpy.dot(yavgs,ns)/ntot
    Ystd = ( numpy.dot(y2s,ns)/ntot - Yavg**2 )**0.5
    
    Xavg = numpy.dot(ns, genotype_matrix)/ntot    
    Xstd = numpy.sqrt( numpy.dot(ns, genotype_matrix**2)/ntot - Xavg**2 )
    XYavg = numpy.dot(yavgs*ns, genotype_matrix)/ntot
    
    rs = (XYavg-Xavg*Yavg)/Xstd/Ystd
    
    # should be this
    lods = -ntot*numpy.log(1.0-numpy.square(rs))/2/log(10.0)
    
    # but bloom et al used this instead (about a 6 fold reduction for their data)
    #lods = -len(ns)*numpy.log(1.0-numpy.square(rs))/2/log(10.0)

    return lods

def find_peaks(lods, peak_radius=1000, lod_threshold=2):

    peaks = {chromosome: [[],[]] for chromosome in lods.keys()}
    
    for chromosome in lods.keys():
        
        chromosome_positions = numpy.array(lods[chromosome][0])
        chromosome_lods = numpy.array(lods[chromosome][1])
        
        while chromosome_lods.max() > lod_threshold:
            
            max_idx = chromosome_lods.nanargmax()
            max_position = chromosome_positions[max_idx]
            
            idxs = numpy.searchsorted(chromosome_positions,[max_position-peak_radius,max_position+peak_radius])
            
            peaks[chromosome][0].append(max_position)
            peaks[chromosome][1].append(chromosome_lods.max())
            
            chromosome_lods[idxs[0]:idxs[1]] = 0
    return peaks

def calculate_chromosome_lods(genotypes, phenotypes, chromosome_boundaries):
    
    if chromosome_boundaries == None:
        chromosome_boundaries = parse_file.get_chromosome_boundaries(loci)

    lods = calculate_lods(genotypes, phenotypes)
    
    chromosome_peak_lods = []
    chromosome_peak_lod_idxs = []
    for chromosome_start,chromosome_end in chromosome_boundaries:
        chromosome_peak_lods.append(numpy.nanmax(lods[chromosome_start:chromosome_end]))
        chromosome_peak_lod_idxs.append(numpy.nanargmax(lods[chromosome_start:chromosome_end])+chromosome_start)
    
    return numpy.array(chromosome_peak_lods), numpy.array(chromosome_peak_lod_idxs)

def calculate_lod_threshold(bootstrapped_peak_lods, observed_peak_lods, FDR=0.05):
    
    print "FDR =", FDR
    
    num_bootstraps = len(bootstrapped_peak_lods)
    print num_bootstraps
    all_lods = []
    bootstrapped_lods = []
    
    all_lods.extend(observed_peak_lods)
    for i in xrange(0,num_bootstraps):
        all_lods.extend(bootstrapped_peak_lods[i])
        bootstrapped_lods.extend(bootstrapped_peak_lods[i])
    
    
    bootstrapped_lods = numpy.array(bootstrapped_lods)
    
    all_lods.sort()
    
    for lam in reversed(all_lods):
        false_positives = (bootstrapped_lods >= lam).sum()*1.0/(num_bootstraps*1.0)
        
        observed = (observed_peak_lods >= lam).sum()+1e-012
        observed_FDR = false_positives/observed
        print lam, observed_FDR
        if observed_FDR >= FDR:
            return lam
            
    return all_bootstrapped_lods[0]       
    
if __name__=='__main__':
    import pylab
    import regression
    import parse_file_2_1_15
	
    ej_data, loci, environments = parse_file_2_1_15.get_ej_data()
    chromosome_boundaries = parse_file_2_1_15.get_chromosome_boundaries(loci)
    genotype_matrix, initial_fitnesses, final_fitnesses, segregant_vector = parse_file_2_1_15.get_genotype_phenotype_matrices_ej(ej_data,environment='SC')
	#genotype_matrix_collapsed, initfit_collapsed, finalfit_collapsed = parse_file_2_1_15.get_collapsed_genotype_init_final_fitnesses(genotype_matrix, initial_fitnesses, final_fitnesses, segregant_vector)
    
    beta_current, F_current =  regression.ordinary_linear_regression(final_fitnesses,initial_fitnesses)
    final_fitnesses = numpy.vstack((final_fitnesses,final_fitnesses**2,numpy.ones((len(final_fitnesses),)))).T
    dependent_variables = initial_fitnesses
    residuals_for_detection = numpy.copy(final_fitnesses)
    residuals_for_detection[:,0] = residuals_for_detection[:,0] - F_current(dependent_variables)
    residuals_for_detection[:,1] = residuals_for_detection[:,1] - 2*residuals_for_detection[:,0]*F_current(dependent_variables) + F_current(dependent_variables)**2
	
    genotype_matrix_collapsed, residuals_collapsed = parse_file_2_1_15.get_collapsed_genotype_residual_matrices(genotype_matrix, residuals_for_detection, segregant_vector)
    pylab.figure(1)
    pylab.xlabel('Position on chromosome',fontsize=12)
    pylab.ylabel('LOD',fontsize=12)
    
    #chromosome_peak_lods, chromosome_peak_lod_idxs = qtl_detection_adaptability.calculate_chromosome_lods(genotype_matrix_collapsed, residuals_collapsed, chromosome_boundaries)
    
    lods = calculate_lods(genotype_matrix_collapsed, residuals_collapsed)
    locations = numpy.array([location for chromosome,location in loci])
    
    permuted_phenotype_matrix = numpy.random.permutation(residuals_collapsed)
    permuted_lods = calculate_lods(genotype_matrix_collapsed, permuted_phenotype_matrix)
    
    chromosome_peak_lods, chromosome_peak_lod_idxs = calculate_chromosome_lods(genotype_matrix_collapsed, residuals_collapsed, chromosome_boundaries)
    
    bootstrapped_lods = []
    for i in xrange(0,1000):
        #print i
        #print "Generating permutation..."
        permuted_phenotype_matrix = numpy.random.permutation(residuals_collapsed)
        #numpy.random.shuffle(permuted_phenotype_matrix)
        
        #print "Calculating LODS..."
        bootstrapped_lods.append(calculate_chromosome_lods(genotype_matrix_collapsed, permuted_phenotype_matrix, chromosome_boundaries)[0])
        #print bootstrapped_lods[-1].shape, bootstrapped_lods[-1]
        #print "Done"   
     
    lam = calculate_lod_threshold(bootstrapped_lods, chromosome_peak_lods)
    print lam
    
    pylab.figure(1)
    for chromosome_start,chromosome_end in chromosome_boundaries:
        pylab.plot(locations[chromosome_start:chromosome_end], lods[chromosome_start:chromosome_end])    
    pylab.plot([0,1600000],[lam,lam],'k:')
    
    pylab.figure(2)
    for chromosome_start,chromosome_end in chromosome_boundaries:
        pylab.plot(locations[chromosome_start:chromosome_end], permuted_lods[chromosome_start:chromosome_end])    
    pylab.plot([0,1600000],[lam,lam],'k:')
    
    pylab.show()
    
        