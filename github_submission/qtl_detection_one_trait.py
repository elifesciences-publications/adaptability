import qtl_detection_adaptability
import numpy
import matplotlib.pylab as pt
import regression

##Updated 11-2-2016 to include a second epistasis model, 'detect_qtls_with_epistasis2'
##Updated 12-21-2016 to calculate confidence intervals based on LOD drop-off during QTL detection
##Updated 1-18-2016 to include a function for finding QTLs on the two environments separately

def detect_qtls(genotype_mat, phenotype_list_sc,  phenotype_list_ypd, helper_matrix_sc, helper_matrix_ypd, pops_per_seg_sc, pops_per_seg_ypd):
	
	#Initialize residuals as phenotypes; format is [<phenotype>,<phenotype**2>,pops_per_seg]
	
	n_segs = len(pops_per_seg_sc)
	
	phenotype_means_sc = numpy.dot(phenotype_list_sc, helper_matrix_sc)/pops_per_seg_sc
	phenotype_second_moments_sc = numpy.dot(phenotype_list_sc**2, helper_matrix_sc)/pops_per_seg_sc
	
	phenotypes_sc = numpy.append(phenotype_means_sc.reshape((n_segs,1)), phenotype_second_moments_sc.reshape((n_segs,1)), axis=1)
	phenotypes_sc = numpy.append(phenotypes_sc, pops_per_seg_sc.reshape((n_segs,1)), axis = 1)
	
	phenotype_means_ypd = numpy.dot(phenotype_list_ypd, helper_matrix_ypd)/pops_per_seg_ypd
	phenotype_second_moments_ypd = numpy.dot(phenotype_list_ypd**2, helper_matrix_ypd)/pops_per_seg_ypd
	
	phenotypes_ypd = numpy.append(phenotype_means_ypd.reshape((n_segs,1)), phenotype_second_moments_ypd.reshape((n_segs,1)), axis=1)
	phenotypes_ypd = numpy.append(phenotypes_ypd, pops_per_seg_ypd.reshape((n_segs,1)), axis = 1)
	
	expanded_genotype_mat_sc = numpy.dot(helper_matrix_sc, genotype_mat)
	expanded_genotype_mat_ypd = numpy.dot(helper_matrix_ypd, genotype_mat)
	
	lods_sc = qtl_detection_adaptability.calculate_lods(genotype_mat, phenotypes_sc)
	lods_ypd = qtl_detection_adaptability.calculate_lods(genotype_mat, phenotypes_ypd)
	lods = lods_sc + lods_ypd
	
	top_lod = numpy.nanmax(lods)
	top_lod_idx = numpy.nanargmax(lods)
	
	##Confidence intervals around this peak
	
	intervals = []
	
	relative_height_lb = top_lod
	relative_height_ub = top_lod
	lb_index = top_lod_idx
	ub_index = top_lod_idx
	
	consecutive_low_scores = 0
	
	while consecutive_low_scores < 40:
		
		lb_index -= 1
		relative_height_lb = lods[lb_index]
		if relative_height_lb < top_lod - 1.5:
			consecutive_low_scores += 1
		else:
			consecutive_low_scores = 0
		if consecutive_low_scores == 1:
			first_consecutive_low_idx = lb_index
			
	consecutive_low_scores = 0
		
	while consecutive_low_scores < 40:
	
		ub_index += 1
		relative_height_ub = lods[ub_index]	
		if relative_height_ub < top_lod - 1.5:
			consecutive_low_scores += 1
		else:
			consecutive_low_scores = 0
		if consecutive_low_scores == 1:
			first_consecutive_high_idx = ub_index
			
	bootstrapped_lods = []
	n_iter = 1000
	for i in xrange(0,n_iter):
		#print i
		#print "Generating permutation..."
		permutation = numpy.random.permutation(numpy.arange(n_segs))
		
		permuted_phenotype_matrix_ypd = phenotypes_ypd[permutation,:]
		permuted_phenotype_matrix_sc = phenotypes_sc[permutation,:]
		
		permuted_lods_ypd = qtl_detection_adaptability.calculate_lods(genotype_mat, permuted_phenotype_matrix_ypd)
		permuted_lods_sc = qtl_detection_adaptability.calculate_lods(genotype_mat, permuted_phenotype_matrix_sc)
		permuted_lods = permuted_lods_sc + permuted_lods_ypd
		#print permuted_lods
		bootstrapped_lods.append(numpy.nanmax(permuted_lods))
	
	#print bootstrapped_lods	
	sig_threshold = numpy.sort(bootstrapped_lods)[::-1][49]
	#sig_threshold = numpy.sort(bootstrapped_lods)[::-1][99]

	# pt.plot(lods)
#  	pt.plot(permuted_lods,'g')
#  	pt.axhline(sig_threshold, 0, 1, 'k')
#  	pt.show()
	print 'sig_threshold =', sig_threshold

	#print chromosome_peak_lod_idxs

	if top_lod > sig_threshold:
		new_QTLs = [top_lod_idx]
		intervals.append([first_consecutive_low_idx, first_consecutive_high_idx])
		print intervals	
	current_QTLs = new_QTLs
	all_QTLs_found = False
	while all_QTLs_found ==False:
		print current_QTLs
		#Fit a linear model using the current QTL list--or a nonlinear model
		
		qtl_matrix_sc = expanded_genotype_mat_sc[:,current_QTLs]
		qtl_matrix_ypd = expanded_genotype_mat_ypd[:,current_QTLs]
		
		beta_sc, betanorm_sc, F_sc = regression.ordinary_linear_regression(phenotype_list_sc,qtl_matrix_sc)
		beta_ypd, betanorm_ypd, F_ypd = regression.ordinary_linear_regression(phenotype_list_ypd,qtl_matrix_ypd)
		
		residuals_sc = phenotype_list_sc - F_sc(numpy.dot(beta_sc,qtl_matrix_sc.T))
		residuals_ypd = phenotype_list_ypd - F_ypd(numpy.dot(beta_ypd,qtl_matrix_ypd.T))
		
		residuals_means_sc = numpy.dot(residuals_sc, helper_matrix_sc)/pops_per_seg_sc
		residuals_second_moments_sc = numpy.dot(residuals_sc**2, helper_matrix_sc)/pops_per_seg_sc
		
		phenotypes_new_sc = numpy.append(residuals_means_sc.reshape((n_segs,1)),residuals_second_moments_sc.reshape((n_segs,1)), axis=1)
		phenotypes_new_sc = numpy.append(phenotypes_new_sc, pops_per_seg_sc.reshape((n_segs,1)), axis = 1)
		
		residuals_means_ypd = numpy.dot(residuals_ypd, helper_matrix_ypd)/pops_per_seg_ypd
		residuals_second_moments_ypd = numpy.dot(residuals_ypd**2, helper_matrix_ypd)/pops_per_seg_ypd
		
		phenotypes_new_ypd = numpy.append(residuals_means_ypd.reshape((n_segs,1)), residuals_second_moments_ypd.reshape((n_segs,1)), axis=1)
		phenotypes_new_ypd = numpy.append(phenotypes_new_ypd, pops_per_seg_ypd.reshape((n_segs,1)), axis = 1)
		
		lods_sc = qtl_detection_adaptability.calculate_lods(genotype_mat, phenotypes_new_sc)
		lods_ypd = qtl_detection_adaptability.calculate_lods(genotype_mat, phenotypes_new_ypd)
		lods = lods_sc + lods_ypd

		top_lod = numpy.nanmax(lods)
		top_lod_idx = numpy.nanargmax(lods)
		
		relative_height_lb = top_lod
		relative_height_ub = top_lod
		lb_index = top_lod_idx
		ub_index = top_lod_idx
		
		consecutive_low_scores = 0
	
		while consecutive_low_scores < 20:
		
			lb_index -= 1
			relative_height_lb = lods[lb_index]
			
			if relative_height_lb < top_lod - 2:
				consecutive_low_scores += 1
			else:
				consecutive_low_scores = 0
			if consecutive_low_scores == 1:
				first_consecutive_low_idx = lb_index
				
		consecutive_low_scores = 0
		
		while consecutive_low_scores < 20:
	
			ub_index += 1
			relative_height_ub = lods[ub_index]	
			if relative_height_ub < top_lod - 2:
				consecutive_low_scores += 1
			else:
				consecutive_low_scores = 0
			if consecutive_low_scores == 1:
				first_consecutive_high_idx = ub_index
				
		bootstrapped_lods = []
		n_iter = 1000
		for i in xrange(0,n_iter):
			
			permutation = numpy.random.permutation(numpy.arange(n_segs))
			
			permuted_phenotype_matrix_ypd = phenotypes_new_ypd[permutation,:]
			permuted_phenotype_matrix_sc = phenotypes_new_sc[permutation,:]
			

			permuted_lods_ypd = qtl_detection_adaptability.calculate_lods(genotype_mat, permuted_phenotype_matrix_ypd)
			permuted_lods_sc = qtl_detection_adaptability.calculate_lods(genotype_mat, permuted_phenotype_matrix_sc)
			permuted_lods = permuted_lods_sc + permuted_lods_ypd
			
			bootstrapped_lods.append(numpy.nanmax(permuted_lods))
		
		
		sig_threshold = numpy.sort(bootstrapped_lods)[-50] #p < .05
		#sig_threshold = numpy.sort(bootstrapped_lods)[-100] #p < .05



		print 'sig_threshold =', sig_threshold
		#pt.plot(lods)
 		#pt.plot(permuted_lods,'g')
 		#pt.axhline(sig_threshold, 0, 1, 'k')
 		#pt.show()
		if top_lod > sig_threshold:
			current_QTLs.append(top_lod_idx)
			intervals.append([first_consecutive_low_idx, first_consecutive_high_idx])
		else:
			print 'all_QTLs_found'
			all_QTLs_found = True
					
	return current_QTLs, beta_ypd*betanorm_ypd, beta_sc*betanorm_sc, numpy.array(intervals)

def detect_qtls_one_envt(genotype_mat, phenotype_list, helper_matrix, pops_per_seg):
	
	#Initialize residuals as phenotypes; format is [<phenotype>,<phenotype**2>,pops_per_seg]
	
	n_segs = len(pops_per_seg)
	
	phenotype_means = numpy.dot(phenotype_list, helper_matrix)/pops_per_seg
	phenotype_second_moments = numpy.dot(phenotype_list**2, helper_matrix)/pops_per_seg
	
	phenotypes = numpy.append(phenotype_means.reshape((n_segs,1)), phenotype_second_moments.reshape((n_segs,1)), axis=1)
	phenotypes = numpy.append(phenotypes, pops_per_seg.reshape((n_segs,1)), axis = 1)
	
	expanded_genotype_mat = numpy.dot(helper_matrix, genotype_mat)
	
	
	lods = qtl_detection_adaptability.calculate_lods(genotype_mat, phenotypes)
	
	top_lod = numpy.nanmax(lods)
	top_lod_idx = numpy.nanargmax(lods)
	
	##Confidence intervals around this peak
	
	intervals = []
	
	relative_height_lb = top_lod
	relative_height_ub = top_lod
	lb_index = top_lod_idx
	ub_index = top_lod_idx
	
	consecutive_lowores = 0
	
	while consecutive_lowores < 20:
		
		lb_index -= 1
		relative_height_lb = lods[lb_index]
		if relative_height_lb < top_lod - 1.5:
			consecutive_lowores += 1
		else:
			consecutive_lowores = 0
		if consecutive_lowores == 1:
			first_consecutive_low_idx = lb_index
			
	consecutive_lowores = 0
		
	while consecutive_lowores < 20:
	
		ub_index += 1
		relative_height_ub = lods[ub_index]	
		if relative_height_ub < top_lod - 1.5:
			consecutive_lowores += 1
		else:
			consecutive_lowores = 0
		if consecutive_lowores == 1:
			first_consecutive_high_idx = ub_index
			
	bootstrapped_lods = []
	n_iter = 1000
	for i in xrange(0,n_iter):
		#print i
		#print "Generating permutation..."
		permutation = numpy.random.permutation(numpy.arange(n_segs))
		
		
		permuted_phenotype_matrix = phenotypes[permutation,:]
		
		permuted_lods = qtl_detection_adaptability.calculate_lods(genotype_mat, permuted_phenotype_matrix)
		
		#print permuted_lods
		bootstrapped_lods.append(numpy.nanmax(permuted_lods))
	
	#print bootstrapped_lods	
	sig_threshold = numpy.sort(bootstrapped_lods)[::-1][49]
	#sig_threshold = numpy.sort(bootstrapped_lods)[::-1][99]

	# pt.plot(lods)
#  	pt.plot(permuted_lods,'g')
#  	pt.axhline(sig_threshold, 0, 1, 'k')
#  	pt.show()
	print 'sig_threshold =', sig_threshold

	#print chromosome_peak_lod_idxs

	if top_lod > sig_threshold:
		new_QTLs = [top_lod_idx]
		intervals.append([first_consecutive_low_idx, first_consecutive_high_idx])
		print intervals	
	current_QTLs = new_QTLs
	all_QTLs_found = False
	while all_QTLs_found ==False:
		print current_QTLs
		#Fit a linear model using the current QTL list--or a nonlinear model
		
		qtl_matrix = expanded_genotype_mat[:,current_QTLs]
		
		
		beta, betanorm, F = regression.ordinary_linear_regression(phenotype_list,qtl_matrix)
		
		
		residuals = phenotype_list - F(numpy.dot(beta,qtl_matrix.T))
		
		
		residuals_means = numpy.dot(residuals, helper_matrix)/pops_per_seg
		residuals_second_moments = numpy.dot(residuals**2, helper_matrix)/pops_per_seg
		
		phenotypes_new = numpy.append(residuals_means.reshape((n_segs,1)),residuals_second_moments.reshape((n_segs,1)), axis=1)
		phenotypes_new = numpy.append(phenotypes_new, pops_per_seg.reshape((n_segs,1)), axis = 1)
		
		lods = qtl_detection_adaptability.calculate_lods(genotype_mat, phenotypes_new)

		top_lod = numpy.nanmax(lods)
		top_lod_idx = numpy.nanargmax(lods)
		
		relative_height_lb = top_lod
		relative_height_ub = top_lod
		lb_index = top_lod_idx
		ub_index = top_lod_idx
		
		consecutive_lowores = 0
	
		while consecutive_lowores < 20:
		
			lb_index -= 1
			relative_height_lb = lods[lb_index]
			
			if relative_height_lb < top_lod - 1.5:
				consecutive_lowores += 1
			else:
				consecutive_lowores = 0
			if consecutive_lowores == 1:
				first_consecutive_low_idx = lb_index
				
		consecutive_lowores = 0
		
		while consecutive_lowores < 20:
	
			ub_index += 1
			relative_height_ub = lods[ub_index]	
			if relative_height_ub < top_lod - 1.5:
				consecutive_lowores += 1
			else:
				consecutive_lowores = 0
			if consecutive_lowores == 1:
				first_consecutive_high_idx = ub_index
				
		bootstrapped_lods = []
		n_iter = 1000
		for i in xrange(0,n_iter):
			
			permutation = numpy.random.permutation(numpy.arange(n_segs))
			
			
			permuted_phenotype_matrix = phenotypes_new[permutation,:]
			

			
			permuted_lods = qtl_detection_adaptability.calculate_lods(genotype_mat, permuted_phenotype_matrix)
			
			bootstrapped_lods.append(numpy.nanmax(permuted_lods))
		
		
		sig_threshold = numpy.sort(bootstrapped_lods)[-50] #p < .05
		#sig_threshold = numpy.sort(bootstrapped_lods)[-100] #p < .05



		print 'sig_threshold =', sig_threshold
		#pt.plot(lods)
 		#pt.plot(permuted_lods,'g')
 		#pt.axhline(sig_threshold, 0, 1, 'k')
 		#pt.show()
		if top_lod > sig_threshold:
			current_QTLs.append(top_lod_idx)
			intervals.append([first_consecutive_low_idx, first_consecutive_high_idx])
		else:
			print 'all_QTLs_found'
			all_QTLs_found = True
					
	return current_QTLs, beta*betanorm, numpy.array(intervals)
	

def calculate_qtl_confidence_intervals_lods(qtl_locs, genotype_mat, phenotype_sc, phenotype_ypd, helper_matrix_sc=numpy.identity(229), helper_matrix_ypd=numpy.identity(229), pops_per_seg_sc=numpy.ones((229,)), pops_per_seg_ypd=numpy.ones((229,))):
	
	#This function takes an arbitrary number of phenotypes (columns of phenotype_mat) and assumes qtls have been detected on them jointly
	#evol_env_vector records which environment populations with a given phenotype evolved, if applicable; 1=sc at 37 C, 0=ypd at 30 C.
	#Confidence intervals are calculated based on the location at which the LOD score (log-likelihood) falls to half its maximum value.
	#Initialize residuals as phenotypes; format is [<phenotype>,<phenotype**2>,pops_per_seg]
	
	n_segs = len(pops_per_seg_sc)
	#n_phenotypes = len(evol_env_vector)
	n_loci = genotype_mat.shape[1]
	
	lod_idxs = []
	
	intervals = []
	real_qtl_locs = []
	
	##Set up phenotype matrixes
	
	phenotype_means_sc = numpy.dot(phenotype_sc, helper_matrix_sc)/pops_per_seg_sc
	phenotype_second_moments_sc = numpy.dot(phenotype_sc**2, helper_matrix_sc)/pops_per_seg_sc
	phenotypes_sc = numpy.append(phenotype_means_sc.reshape((n_segs,1)), phenotype_second_moments_sc.reshape((n_segs,1)), axis=1)
	phenotypes_sc = numpy.append(phenotypes_sc, pops_per_seg_sc.reshape((n_segs,1)), axis = 1)
		
	phenotype_means_ypd = numpy.dot(phenotype_ypd, helper_matrix_ypd)/pops_per_seg_ypd
	phenotype_second_moments_ypd = numpy.dot(phenotype_ypd**2, helper_matrix_ypd)/pops_per_seg_ypd
	phenotypes_ypd = numpy.append(phenotype_means_ypd.reshape((n_segs,1)), phenotype_second_moments_ypd.reshape((n_segs,1)), axis=1)
	phenotypes_ypd = numpy.append(phenotypes_ypd, pops_per_seg_ypd.reshape((n_segs,1)), axis = 1)
	
	lods = numpy.zeros((n_loci,))
		
	lods += qtl_detection_adaptability.calculate_lods(genotype_mat, phenotypes_sc)
			
	lods += qtl_detection_adaptability.calculate_lods(genotype_mat, phenotypes_ypd)
	
	for qtl in qtl_locs:
		
		if qtl > 20:
			real_peak = qtl - 20 + numpy.nanargmax(lods[qtl-20:qtl+20])
		else:
			real_peak = numpy.nanargmax(lods[0:qtl+20])
		
		peak_height = lods[real_peak]
		relative_height_lb = peak_height
		relative_height_ub = peak_height
		lb_index = real_peak
		ub_index = real_peak
		
		print real_peak
		print peak_height
		
		while relative_height_lb > .5*peak_height:
			lb_index -= 1
			relative_height_lb = lods[lb_index]
		while relative_height_ub > .5*peak_height:
			ub_index += 1
			relative_height_ub = lods[ub_index]
		
		intervals.append([lb_index, ub_index])	
		real_qtl_locs.append(real_peak)
		
	return real_qtl_locs, numpy.array(intervals)
			
def detect_qtls_above_fitness(genotype_mat, phenotype_list_sc,  phenotype_list_ypd, initfit_list_sc, initfit_list_ypd, helper_matrix_sc, helper_matrix_ypd, pops_per_seg_sc, pops_per_seg_ypd):
	
	#Initialize residuals as phenotypes; format is [<phenotype>,<phenotype**2>,pops_per_seg]
	
	n_segs = len(pops_per_seg_sc)
	n_pops_sc = sum(pops_per_seg_sc)
	n_pops_ypd = sum(pops_per_seg_ypd)
	
	expanded_genotype_mat_sc = numpy.dot(helper_matrix_sc, genotype_mat)
	expanded_genotype_mat_ypd = numpy.dot(helper_matrix_ypd, genotype_mat)
	
	#symmetrize the genotype matrix
	genotype_mat = 1./2.*(genotype_mat - (1 - genotype_mat))
	expanded_genotype_mat_sc = 1./2.*(expanded_genotype_mat_sc - (1 - expanded_genotype_mat_sc))
	expanded_genotype_mat_ypd = 1./2.*(expanded_genotype_mat_ypd - (1 - expanded_genotype_mat_ypd))
	
	#Initial dependent variables are initial fitnesses
	X_sc = numpy.dot(helper_matrix_sc, initfit_list_sc).reshape((n_pops_sc,1))
	X_ypd = numpy.dot(helper_matrix_ypd, initfit_list_ypd).reshape((n_pops_ypd,1))

	current_QTLs = []
	all_QTLs_found = False
	while all_QTLs_found ==False:
		print current_QTLs
		#Fit a linear model using the current QTL list--or a nonlinear model
		
		qtl_matrix_sc = numpy.append(X_sc, expanded_genotype_mat_sc[:,current_QTLs], axis = 1)
		qtl_matrix_ypd = numpy.append(X_ypd, expanded_genotype_mat_ypd[:,current_QTLs], axis = 1)
		
		beta_sc, betanorm_sc, F_sc = regression.ordinary_linear_regression(phenotype_list_sc,qtl_matrix_sc)
		beta_ypd, betanorm_ypd, F_ypd = regression.ordinary_linear_regression(phenotype_list_ypd,qtl_matrix_ypd)
		
		residuals_sc = phenotype_list_sc - F_sc(numpy.dot(beta_sc,qtl_matrix_sc.T))
		residuals_ypd = phenotype_list_ypd - F_ypd(numpy.dot(beta_ypd,qtl_matrix_ypd.T))
		
		residuals_means_sc = numpy.dot(residuals_sc, helper_matrix_sc)/pops_per_seg_sc
		residuals_second_moments_sc = numpy.dot(residuals_sc**2, helper_matrix_sc)/pops_per_seg_sc
		
		phenotypes_new_sc = numpy.append(residuals_means_sc.reshape((n_segs,1)),residuals_second_moments_sc.reshape((n_segs,1)), axis=1)
		phenotypes_new_sc = numpy.append(phenotypes_new_sc, pops_per_seg_sc.reshape((n_segs,1)), axis = 1)
		
		residuals_means_ypd = numpy.dot(residuals_ypd, helper_matrix_ypd)/pops_per_seg_ypd
		residuals_second_moments_ypd = numpy.dot(residuals_ypd**2, helper_matrix_ypd)/pops_per_seg_ypd
		
		phenotypes_new_ypd = numpy.append(residuals_means_ypd.reshape((n_segs,1)), residuals_second_moments_ypd.reshape((n_segs,1)), axis=1)
		phenotypes_new_ypd = numpy.append(phenotypes_new_ypd, pops_per_seg_ypd.reshape((n_segs,1)), axis = 1)
		
		lods_sc = qtl_detection_adaptability.calculate_lods(genotype_mat, phenotypes_new_sc)
		lods_ypd = qtl_detection_adaptability.calculate_lods(genotype_mat, phenotypes_new_ypd)
		lods = lods_sc + lods_ypd
		#pt.plot(lods)
		
		top_lod = numpy.nanmax(lods)
		top_lod_idx = numpy.nanargmax(lods)
		#print top_lod
		bootstrapped_lods = []
		n_iter = 1000
		for i in xrange(0,n_iter):
			
			permutation = numpy.random.permutation(numpy.arange(n_segs))
			
			permuted_phenotype_matrix_ypd = phenotypes_new_ypd[permutation,:]
			permuted_phenotype_matrix_sc = phenotypes_new_sc[permutation,:]
			

			permuted_lods_ypd = qtl_detection_adaptability.calculate_lods(genotype_mat, permuted_phenotype_matrix_ypd)
			permuted_lods_sc = qtl_detection_adaptability.calculate_lods(genotype_mat, permuted_phenotype_matrix_sc)
			permuted_lods = permuted_lods_sc + permuted_lods_ypd
			
			bootstrapped_lods.append(numpy.nanmax(permuted_lods))
		
		
		sig_threshold = numpy.sort(bootstrapped_lods)[-49]

		print 'sig_threshold =', sig_threshold
		#print numpy.sort(bootstrapped_lods)
		#pt.plot(permuted_lods,'g')
		#pt.axhline(sig_threshold,0,1,'k')
		#pt.show()
		
		if top_lod > sig_threshold:
			current_QTLs.append(top_lod_idx)
			
		else:
			print 'all_QTLs_found'
			all_QTLs_found = True
					
	return current_QTLs, beta_ypd*betanorm_ypd, beta_sc*betanorm_sc
		
def detect_qtls_with_epistasis(genotype_mat, phenotype_list_sc,  phenotype_list_ypd, helper_matrix_sc, helper_matrix_ypd, pops_per_seg_sc, pops_per_seg_ypd):
	
	#Initialize residuals as phenotypes; format is [<phenotype>,<phenotype**2>,pops_per_seg]
	
	n_segs = len(pops_per_seg_sc)
	n_pops_sc = sum(pops_per_seg_sc)
	n_pops_ypd = sum(pops_per_seg_ypd)
	
	expanded_genotype_mat_sc = numpy.dot(helper_matrix_sc, genotype_mat)
	expanded_genotype_mat_ypd = numpy.dot(helper_matrix_ypd, genotype_mat)
	
	#symmetrize the genotype matrix
	genotype_mat = 1./2.*(genotype_mat - (1 - genotype_mat))
	expanded_genotype_mat_sc = 1./2.*(expanded_genotype_mat_sc - (1 - expanded_genotype_mat_sc))
	expanded_genotype_mat_ypd = 1./2.*(expanded_genotype_mat_ypd - (1 - expanded_genotype_mat_ypd))
	
	kre33_loc = 9596
	kre_genotypes = genotype_mat[:,kre33_loc]
	kre_genotypes_sc = expanded_genotype_mat_sc[:,kre33_loc]
	kre_genotypes_ypd = expanded_genotype_mat_ypd[:,kre33_loc]
	current_main_effect_QTLs = []#new_QTLs
	current_epistatic_QTLs = []
	all_QTLs_found = False
	while all_QTLs_found ==False:
		print current_main_effect_QTLs
		print current_epistatic_QTLs
		#Fit a linear model using the current QTL list--or a nonlinear model
		
		coefficient_matrix_sc = kre_genotypes_sc.reshape((n_pops_sc,1))
		
		if len(current_main_effect_QTLs) > .5:
			coefficient_matrix_sc = numpy.append(coefficient_matrix_sc, expanded_genotype_mat_sc[:,current_main_effect_QTLs], axis=1)
		if len(current_epistatic_QTLs) > .5:
			coefficient_matrix_sc = numpy.append(coefficient_matrix_sc, kre_genotypes_sc.reshape((n_pops_sc,1))*expanded_genotype_mat_sc[:,current_epistatic_QTLs], axis=1)
		
		coefficient_matrix_ypd = kre_genotypes_ypd.reshape((n_pops_ypd,1))
		if len(current_main_effect_QTLs) > .5:
			coefficient_matrix_ypd = numpy.append(coefficient_matrix_ypd, expanded_genotype_mat_ypd[:,current_main_effect_QTLs], axis=1)
		if len(current_epistatic_QTLs) > .5:
			coefficient_matrix_ypd = numpy.append(coefficient_matrix_ypd, kre_genotypes_ypd.reshape((n_pops_ypd,1))*expanded_genotype_mat_ypd[:,current_epistatic_QTLs], axis=1)
			
		beta_sc, betanorm_sc, F_sc = regression.ordinary_linear_regression(phenotype_list_sc,coefficient_matrix_sc)
		beta_ypd, betanorm_ypd, F_ypd = regression.ordinary_linear_regression(phenotype_list_ypd,coefficient_matrix_ypd)
		
		residuals_sc = phenotype_list_sc - F_sc(numpy.dot(beta_sc,coefficient_matrix_sc.T))
		residuals_ypd = phenotype_list_ypd - F_ypd(numpy.dot(beta_ypd,coefficient_matrix_ypd.T))
		
		residuals_means_sc = numpy.dot(residuals_sc, helper_matrix_sc)/pops_per_seg_sc
		residuals_second_moments_sc = numpy.dot(residuals_sc**2, helper_matrix_sc)/pops_per_seg_sc
		
		phenotypes_new_sc = numpy.append(residuals_means_sc.reshape((n_segs,1)),residuals_second_moments_sc.reshape((n_segs,1)), axis=1)
		phenotypes_new_sc = numpy.append(phenotypes_new_sc, pops_per_seg_sc.reshape((n_segs,1)), axis = 1)
		
		residuals_means_ypd = numpy.dot(residuals_ypd, helper_matrix_ypd)/pops_per_seg_ypd
		residuals_second_moments_ypd = numpy.dot(residuals_ypd**2, helper_matrix_ypd)/pops_per_seg_ypd
		
		phenotypes_new_ypd = numpy.append(residuals_means_ypd.reshape((n_segs,1)), residuals_second_moments_ypd.reshape((n_segs,1)), axis=1)
		phenotypes_new_ypd = numpy.append(phenotypes_new_ypd, pops_per_seg_ypd.reshape((n_segs,1)), axis = 1)
		
		#print phenotypes_new_sc
		
		##Calculate lods for new main-effect loci
		lods_sc = qtl_detection_adaptability.calculate_lods(genotype_mat, phenotypes_new_sc)
		lods_ypd = qtl_detection_adaptability.calculate_lods(genotype_mat, phenotypes_new_ypd)
		lods = lods_sc + lods_ypd
		# pt.figure()
		# pt.plot(lods)
		# pt.show()
		top_lod = numpy.nanmax(lods)
		top_lod_idx = numpy.nanargmax(lods)
		print top_lod
		##Calculate potential epistatic effects of loci already in the model
		
		if len(current_main_effect_QTLs) > .5:
			genotype_mat_interactions = kre_genotypes.reshape((n_segs,1))*genotype_mat[:,current_main_effect_QTLs]
			#print genotype_mat_interactions
			lods_sc_ints = qtl_detection_adaptability.calculate_lods(genotype_mat_interactions, phenotypes_new_sc)
			lods_ypd_ints = qtl_detection_adaptability.calculate_lods(genotype_mat_interactions, phenotypes_new_ypd)
			lods_interactions = lods_sc_ints + lods_ypd_ints
		
			top_lod_int = numpy.nanmax(lods_interactions)
			top_lod_int_idx = current_main_effect_QTLs[numpy.nanargmax(lods_interactions)]
			print top_lod_int
		else:
			top_lod_int = 0
		
		
		bootstrapped_lods = []
		n_iter = 1000
		for i in xrange(0,n_iter):
			
			permutation = numpy.random.permutation(numpy.arange(n_segs))
			
			permuted_phenotype_matrix_ypd = phenotypes_new_ypd[permutation,:]
			permuted_phenotype_matrix_sc = phenotypes_new_sc[permutation,:]
			

			permuted_lods_ypd = qtl_detection_adaptability.calculate_lods(genotype_mat, permuted_phenotype_matrix_ypd)
			permuted_lods_sc = qtl_detection_adaptability.calculate_lods(genotype_mat, permuted_phenotype_matrix_sc)
			permuted_lods = permuted_lods_sc + permuted_lods_ypd
			if len(current_main_effect_QTLs) > .5:
				permuted_lods_sc_ints = qtl_detection_adaptability.calculate_lods(genotype_mat_interactions, permuted_phenotype_matrix_sc)
				permuted_lods_ypd_ints = qtl_detection_adaptability.calculate_lods(genotype_mat_interactions, permuted_phenotype_matrix_ypd)
				permuted_lods_interactions = permuted_lods_sc_ints + permuted_lods_ypd_ints
			
				all_permuted_lods = numpy.append(permuted_lods, permuted_lods_interactions)
			else:
				all_permuted_lods = permuted_lods
				
			bootstrapped_lods.append(numpy.nanmax(all_permuted_lods))
		
		
		sig_threshold = numpy.sort(bootstrapped_lods)[-49]

		print 'sig_threshold =', sig_threshold
		
		if (top_lod > sig_threshold or top_lod_int > sig_threshold):
		
			if top_lod > top_lod_int:
			
				current_main_effect_QTLs.append(top_lod_idx)
				
			elif top_lod_int > top_lod:
			
				current_epistatic_QTLs.append(top_lod_int_idx)
		else:
			print 'all_QTLs_found'
			all_QTLs_found = True
					
	return current_main_effect_QTLs, current_epistatic_QTLs, beta_ypd*betanorm_ypd, beta_sc*betanorm_sc

def detect_qtls_with_epistasis2(genotype_mat, phenotype_list_sc,  initfit_list_sc, phenotype_list_ypd, initfit_list_ypd, helper_matrix_sc, helper_matrix_ypd, pops_per_seg_sc, pops_per_seg_ypd):
	
	#Initialize residuals as phenotypes; format is [<phenotype>,<phenotype**2>,pops_per_seg]
	
	n_segs = len(pops_per_seg_sc)
	n_pops_sc = sum(pops_per_seg_sc)
	n_pops_ypd = sum(pops_per_seg_ypd)
	
	expanded_genotype_mat_sc = numpy.dot(helper_matrix_sc, genotype_mat)
	expanded_genotype_mat_ypd = numpy.dot(helper_matrix_ypd, genotype_mat)
	
	#symmetrize the genotype matrix
	genotype_mat = 1./2.*(genotype_mat - (1 - genotype_mat))
	expanded_genotype_mat_sc = expanded_genotype_mat_sc - .5
	expanded_genotype_mat_ypd = expanded_genotype_mat_ypd - .5
	#expanded_genotype_mat_sc = 1./2.*(expanded_genotype_mat_sc - (1 - expanded_genotype_mat_sc))
	#expanded_genotype_mat_ypd = 1./2.*(expanded_genotype_mat_ypd - (1 - expanded_genotype_mat_ypd))
	
	kre33_loc = 9596
	kre_genotypes = genotype_mat[:,kre33_loc]
	kre_genotypes_sc = expanded_genotype_mat_sc[:,kre33_loc]
	kre_genotypes_ypd = expanded_genotype_mat_ypd[:,kre33_loc]
	current_QTLs = [] #new_QTLs
	
	#At each step we are going to fit the model: delta_X = a + bX + c*kre_genotypes + sum_i=1^n_qtls d_i1*kre_genotypes*g_i + d_i2*(1-kre_genotypes)*g_i
	#At the final step, we will fit the full model and determine if all the coefficients are significant.
	##Initialize: fit delta_X = a + bX + c*kre_genotypes
	
	X_sc = numpy.concatenate((numpy.dot(helper_matrix_sc, initfit_list_sc).reshape((n_pops_sc,1)), kre_genotypes_sc.reshape((n_pops_sc,1)), numpy.ones((n_pops_sc,1))), axis = 1)
	X_ypd = numpy.concatenate((numpy.dot(helper_matrix_ypd, initfit_list_ypd).reshape((n_pops_ypd,1)), kre_genotypes_ypd.reshape((n_pops_ypd,1)),numpy.ones((n_pops_ypd,1))), axis = 1)
	
	all_QTLs_found = False
	while all_QTLs_found ==False:
		
		#If this is not the first iteration, add the (potentially epistatic) qtls to the model
		if len(current_QTLs) > .5:
			qtl_mat_sc = expanded_genotype_mat_sc[:, current_QTLs]
			qtl_mat_ypd = expanded_genotype_mat_ypd[:, current_QTLs]
			
			X_sc_temp = numpy.concatenate((X_sc, qtl_mat_sc, qtl_mat_sc*kre_genotypes_sc.reshape((n_pops_sc,1))), axis=1)
			X_ypd_temp = numpy.concatenate((X_ypd, qtl_mat_ypd, qtl_mat_ypd*kre_genotypes_ypd.reshape((n_pops_ypd,1))), axis=1)
			#print X_sc_temp.shape
		else:
			X_sc_temp = X_sc
			X_ypd_temp = X_ypd
		
		#Calculate residuals:
		beta_sc = numpy.dot(numpy.linalg.inv(numpy.dot(X_sc_temp.T, X_sc_temp)), numpy.dot(X_sc_temp.T, phenotype_list_sc))
		residuals_sc = phenotype_list_sc - numpy.dot(X_sc_temp, beta_sc) #check dot product direction
		
		beta_ypd = numpy.dot(numpy.linalg.inv(numpy.dot(X_ypd_temp.T, X_ypd_temp)), numpy.dot(X_ypd_temp.T, phenotype_list_ypd))
		residuals_ypd = phenotype_list_ypd - numpy.dot(X_ypd_temp, beta_ypd) 
		
		residuals_means_sc = numpy.dot(residuals_sc, helper_matrix_sc)/pops_per_seg_sc
		residuals_second_moments_sc = numpy.dot(residuals_sc**2, helper_matrix_sc)/pops_per_seg_sc
		
		residual_mat_sc = numpy.concatenate((residuals_means_sc.reshape((n_segs,1)), residuals_second_moments_sc.reshape((n_segs,1)), pops_per_seg_sc.reshape((n_segs,1))), axis=1)
		
		residuals_means_ypd = numpy.dot(residuals_ypd, helper_matrix_ypd)/pops_per_seg_ypd
		residuals_second_moments_ypd = numpy.dot(residuals_ypd**2, helper_matrix_ypd)/pops_per_seg_ypd
		
		residual_mat_ypd = numpy.concatenate((residuals_means_ypd.reshape((n_segs,1)), residuals_second_moments_ypd.reshape((n_segs,1)), pops_per_seg_ypd.reshape((n_segs,1))), axis=1)
		
		##Calculate lods for new loci
		lods_sc = qtl_detection_adaptability.calculate_lods(genotype_mat, residual_mat_sc)
		lods_ypd = qtl_detection_adaptability.calculate_lods(genotype_mat, residual_mat_ypd)
		lods = lods_sc + lods_ypd
		
		top_lod = numpy.nanmax(lods)
		top_lod_idx = numpy.nanargmax(lods)
		
		##Bootstrap over segregants
		bootstrapped_lods = []
		n_iter = 1000
		for i in xrange(0,n_iter):
			
			permutation = numpy.random.permutation(numpy.arange(n_segs))
			
			permuted_phenotype_matrix_ypd = residual_mat_ypd[permutation,:]
			permuted_phenotype_matrix_sc = residual_mat_sc[permutation,:]
			

			permuted_lods_ypd = qtl_detection_adaptability.calculate_lods(genotype_mat, permuted_phenotype_matrix_ypd)
			permuted_lods_sc = qtl_detection_adaptability.calculate_lods(genotype_mat, permuted_phenotype_matrix_sc)
			permuted_lods = permuted_lods_sc + permuted_lods_ypd
			
			bootstrapped_lods.append(numpy.nanmax(permuted_lods))
		
		
		sig_threshold = numpy.sort(bootstrapped_lods)[-49]

		print 'sig_threshold =', sig_threshold
		#print numpy.sort(bootstrapped_lods)
		#pt.plot(permuted_lods,'g')
		#pt.plot(lods,'b')
		#pt.axhline(sig_threshold,0,1,'k')
		#pt.show()
		
		if top_lod > sig_threshold:
			current_QTLs.append(top_lod_idx)
			
		else:
			print 'all_QTLs_found'
			all_QTLs_found = True
		
	return current_QTLs, beta_sc, beta_ypd