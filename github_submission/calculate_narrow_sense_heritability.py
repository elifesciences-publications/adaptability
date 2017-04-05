import numpy
import scipy.stats
import matplotlib.pylab as pt
import statsmodels.api as sm
#import scipy.optimize.minimize

def narrow_sense_hsq(genotype_mat, phenotype_list):
	##Calculates narrow-sense heritability for a given phenotype (length n_segs), with genotypes specified by genotype_mat (n_segs x n_loci).
	##This is one iteration; functions with built-in bootstrapping are found below.
	##This is a method-of-moments based calculation, currently deprecated in favor of REML (below)
	##Calculate relatedness matrix
	#genotype_mat = .5 + genotype_mat
	n_loci = numpy.float(genotype_mat.shape[1])
	n_segs = len(phenotype_list)
	
	rm_allele_frequencies = numpy.mean(genotype_mat, axis = 0)
	by_allele_frequencies = 1 - rm_allele_frequencies
	
	fixed_loci = numpy.arange(n_loci)[(rm_allele_frequencies < .1) + (by_allele_frequencies < .1)]
	genotype_mat = numpy.delete(genotype_mat, fixed_loci, axis = 1)
	rm_allele_frequencies = numpy.delete(rm_allele_frequencies, fixed_loci)
	by_allele_frequencies = numpy.delete(by_allele_frequencies, fixed_loci)
	n_loci -= len(fixed_loci)
	freq_scaled_genotype_mat = -1*by_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*genotype_mat + rm_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*(1 - genotype_mat)
	
	
	relatedness_mat = numpy.dot(freq_scaled_genotype_mat, freq_scaled_genotype_mat.T)/numpy.float(n_loci)
	phenotype_pairwise_difs = numpy.ones((n_segs,n_segs),dtype='float')
	correct_inds = numpy.where(relatedness_mat < .9) #This just excludes the i=j case, or the comparison between a segregant and itself.
	
	relatedness_list = relatedness_mat[correct_inds] - numpy.mean(relatedness_mat[correct_inds])
	for i in range(n_segs):
		for j in range(n_segs):
			phenotype_pairwise_difs[i,j] = (phenotype_list[i] - phenotype_list[j])**2

	# #Center y data
	phenotype_pairwise_difs_list = phenotype_pairwise_difs[correct_inds] - numpy.mean(phenotype_pairwise_difs[correct_inds])

	#Calculate h^2
	#print n_loci, n_segs
	betahat = 1./numpy.sum(relatedness_list*relatedness_list)*numpy.sum(relatedness_list*phenotype_pairwise_difs_list)
	#(binned_relatedness_list, bin_locs, bin_nums) = scipy.stats.binned_statistic(relatedness_list, phenotype_pairwise_difs_list, bins=100)
	#pt.scatter(bin_locs[:-1], binned_relatedness_list)
 	#pt.plot(numpy.arange(-1,1,.01),betahat*numpy.arange(-1,1,.01),'k')
# 	pt.figure()
# 	pt.scatter(relatedness_list, phenotype_pairwise_difs_list)
# 	pt.plot(numpy.arange(-1,1,.01),betahat*numpy.arange(-1,1,.01),'k')
 	#pt.show()
	#print betahat, numpy.var(phenotype_list)
	
	h2 = betahat/-2./numpy.var(phenotype_list)
	
	return h2

def narrow_sense_hsq_REML(genotype_mat, phenotype_list, helper_matrix, variance_comps_initial_guess):
	
	##The basic model we are trying to fit is Y = Zu + e, where Z is the genotype matrix, and u and e are normal random variables.
	##In the case where there are replicates, we are fitting Y_i = \sum Z_ik u_k + H_ij d_j + e_i, where Y_i is the phenotype on a population; d_j is a normal random variable associated with each genotype (this is the heritable part of the variance that is not explained by a linear model of loci), H_ij is the helper matrix (n_pops x n_segs) which is 1 if pop j is descended from segregant i, 0 otherwise; and e_i is a random variable for each population.
	##
	##This function calculates h^2 by maximizing the log-likelihood associated with the model KY = K(Zu + e), where e is N(0,sigma_e^2), u is N(0,sigma_u^2), sigma_u^2 + sigma_e^2 = var(Y), and K is a matrix such that KX = 0, where X are the fixed effects (here a vector of ones).
	##And Z is the genotype matrix, scaled to have mean 0 and variance 1. Z is n_pops x n_loci; u has n_loci entries; Y has n_pops entries.
	
	n_loci = numpy.float(genotype_mat.shape[1])
	
	centered_phenotypes = (numpy.array(phenotype_list) - numpy.mean(phenotype_list))/numpy.std(phenotype_list) #For ease, we are taking out the global mean at the start. Note that this fixes one coefficient and thus biases our variance estimates down slightly, but it should be a very small correction (i.e. 230/229)
	
	expanded_genotype_mat = numpy.dot(helper_matrix,genotype_mat)
	n_pops = expanded_genotype_mat.shape[0]
	
	rm_allele_frequencies = numpy.mean(expanded_genotype_mat, axis = 0)
	by_allele_frequencies = 1. - rm_allele_frequencies
	
	fixed_loci = numpy.arange(n_loci)[(rm_allele_frequencies < .02) + (by_allele_frequencies < .02)]
	genotype_mat = numpy.delete(genotype_mat, fixed_loci, axis = 1)
	
	rm_allele_frequencies = numpy.delete(rm_allele_frequencies, fixed_loci)
	by_allele_frequencies = numpy.delete(by_allele_frequencies, fixed_loci)
	
	n_loci -= len(fixed_loci)
	
	expanded_genotype_mat = numpy.dot(helper_matrix,genotype_mat)
	#Scale the genotype matrix so that sigma_u^2 is the additive variance
	
	freq_scaled_genotype_mat = -1*by_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*expanded_genotype_mat + rm_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*(1 - expanded_genotype_mat)
	
	G = freq_scaled_genotype_mat
	
	relatedness_mat = numpy.dot( G, G.T )/numpy.float(n_loci)
	
	##The other covariance matrix contains the replicate structure:
	
	H = numpy.dot(helper_matrix, helper_matrix.T)
	
	##The matrix K is (n-r)xn, such that KX = 0, where X are the fixed effects:
	
	K = numpy.identity(n_pops - 1) - numpy.ones((n_pops - 1, n_pops -1))*1/float(n_pops)
	K = numpy.concatenate( (K, -1*numpy.ones((n_pops - 1, 1))*1./float(n_pops)), axis=1)
	
	##REML: Ky ~ N(0, KVK^T)
	
	Ky = numpy.dot(K, centered_phenotypes).reshape((n_pops-1,1))
	
	n_var_comps = len(variance_comps_initial_guess)
	
	V_mats = [G/numpy.sqrt(n_loci), helper_matrix]
	V_list = []
	
	for i in range(n_var_comps):
		V_list.append( numpy.dot( numpy.dot(K, V_mats[i]), numpy.dot( V_mats[i].T, K.T ) ) - numpy.dot( K, K.T) ) #This is a list that will be useful for calculating derivatives; the second term comes from the fact that we are constraining the variance components to sum to 1
		
	converged = 0
	var_comps = variance_comps_initial_guess
	
	n_steps = 0
	while (converged < .5 and n_steps < 100):
		V_tot = numpy.zeros((n_pops-1,n_pops-1),dtype='float')
		for i in range(n_var_comps):
			
			V_tot += numpy.dot( numpy.dot(K, V_mats[i]), numpy.dot( V_mats[i].T, K.T ) )*var_comps[i]
	
		V_tot += numpy.dot( K, K.T )*( 1 - numpy.sum(var_comps) ) ## we are fixing variance components to sum to 1 during estimation
	
		P = numpy.linalg.inv( V_tot )
		
		ll_deriv = ll_derivative(P, V_list, Ky)
		ll_hess = ll_hessian(P, V_list, Ky)
		
		#Newton's method for direction to descend
		
		
		if ll_hess.size < 1.5: ##second derivative is scalar; optimizing over one variable
			
			var_step = ll_deriv/ll_hess[0]
		else:
			try:
				var_step = numpy.dot( numpy.linalg.inv(ll_hess), ll_deriv )
			except: ##This should be rare, but if the hessian is singular, reset the variance components to [.1,.5]
				print 'Error: singular Hessian'
				var_step = numpy.array([.1,.5]) - var_comps
				
			
		
		v_comps_temp = var_comps + var_step
		
		if numpy.sum(v_comps_temp) > 1:
			var_step_amended = var_step/numpy.sum(v_comps_temp) - 10**(-4)
		if numpy.any(v_comps_temp < 0):
			ind = numpy.where(v_comps_temp < 0)
			var_step_amended = var_step
			for i in ind:
				var_step_amended[i] = var_step[i] - v_comps_temp[i] + 10**(-4)
		else:
			var_step_amended = var_step
			
		var_comps += var_step_amended
		
		if max( abs(var_step) ) < 10**(-4):
			converged = 1
		
		n_steps += 1
	
	
	h2_estimate = var_comps[0]
	
	return h2_estimate

def ll_derivative(P, V_list, Ky):
	
	#P is the inverse variance-covariance matrix of the model; in this case, P = (K(\sum_i V_i sigma_i^2)K^T)^-1
	#Where V_i is a list of matrices, each npops x npops, which contain the random-effects error structures, and sigma_i^2 are the variance components.
	#V_list is a list of the covariance matrices, V_list[i] = KV_iK^T.
	#Returns a list of the partial derivatives of the log-likelihood with respect to the variance components.
	
	
	ll_deriv = numpy.zeros((len(V_list),),dtype='float')
	
	
	for i in range(len(V_list)):
		ll_deriv[i] =  -.5*numpy.trace( numpy.dot(P, V_list[i]) ) + .5*numpy.dot( numpy.dot( numpy.dot( numpy.dot( Ky.T, P), V_list[i] ), P), Ky)
	
	return ll_deriv
	
def ll_hessian(P, V_list, Ky):
	
	Hessian = numpy.zeros((len(V_list),len(V_list)), dtype='float')
	
	for i in range(len(V_list)):
		for j in range(i+1):
			
			Hessian[i,j] = .5*numpy.dot(numpy.dot( numpy.dot( numpy.dot( numpy.dot( numpy.dot(Ky.T, P), V_list[i]), P), V_list[j]), P), Ky)
			#Hessian[i,j] = -.5*numpy.trace( numpy.dot( numpy.dot( numpy.dot( P, V_list[i] ), P ), V_list[j] ) ) + numpy.dot( numpy.dot( numpy.dot( numpy.dot( numpy.dot(Ky.T, P), V_list[i]), V_list[j]), P), Ky)
			Hessian[j,i] = Hessian[i,j]
	
	return Hessian

def narrow_sense_hsq_ML(genotype_mat, phenotype_list, helper_matrix):
	
	##This function calculates h^2 by maximizing the log-likelihood associated with the model Y = Zu + e, where e is N(0,sigma_e^2), and u is N(0,sigma_u^2, with respect to sigma_u^2 and sigma_e^2
	##And Z is the genotype matrix, scaled to have mean 0 and variance 1. Z is n_pops x n_loci; u has n_loci entries; Y has n_pops entries.
	
	n_loci = numpy.float(genotype_mat.shape[1])
	
	centered_phenotypes = (numpy.array(phenotype_list) - numpy.mean(phenotype_list))/numpy.std(phenotype_list) #For ease, we are taking out the global mean at the start. Note that this fixes one coefficient and thus biases our variance estimates down slightly, but it should be a very small correction (i.e. 230/229)
	
	expanded_genotype_mat = numpy.dot(helper_matrix,genotype_mat)
	n_pops = expanded_genotype_mat.shape[0]
	
	rm_allele_frequencies = numpy.mean(expanded_genotype_mat, axis = 0)
	by_allele_frequencies = 1 - rm_allele_frequencies
	
	fixed_loci = numpy.arange(n_loci)[(rm_allele_frequencies < .1) + (by_allele_frequencies < .1)]
	genotype_mat = numpy.delete(genotype_mat, fixed_loci, axis = 1)
	
	rm_allele_frequencies = numpy.delete(rm_allele_frequencies, fixed_loci)
	by_allele_frequencies = numpy.delete(by_allele_frequencies, fixed_loci)
	
	n_loci -= len(fixed_loci)
	
	expanded_genotype_mat = numpy.dot(helper_matrix,genotype_mat)
	#Scale the genotype matrix so that sigma_u^2 is the additive variance
	
	freq_scaled_genotype_mat = -1*by_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*expanded_genotype_mat + rm_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*(1 - expanded_genotype_mat)
	#freq_scaled_genotype_mat = 2*expanded_genotype_mat - 1
	
	G = freq_scaled_genotype_mat
	relatedness_mat = numpy.dot( G, G.T )/numpy.float(n_loci)
	
	##Initialize guesses for sigma2_mu:
	
	#initial_params = .5 ##We could do this more intelligently, using H^2
	
	#Use scipy.minimize to minimize the negative log-likelihood with respect to sigma2_mu and sigma2_e
	
	#min_obj = scipy.optimize.minimize(loglikelihood_random_effects_model, x0 = initial_params, args = (freq_scaled_genotype_mat, phenotype_list), bounds=(0,1))
	
	#sigma2_u = min_obj.x
	
	##Because the function is one-dimensional, let's just sweep across a range of parameters and try to find the minimum by hand
	
	step = .01
	
	n_steps_u = int(.9/step)
	n_steps_e = int(.49/step)
	
	sigma2_u_vec = numpy.arange(.1,1,step)*numpy.var(centered_phenotypes)
	sigma2_e_vec = numpy.arange(.01,.5,step)*numpy.var(centered_phenotypes)
	
	negloglikelihood_array = numpy.zeros((n_steps_u,n_steps_e),dtype='float')
	negloglikelihood_vec = numpy.zeros((n_steps_u,))
	negloglikelihood_vec2 = numpy.zeros((n_steps_u,))
	##The matrix K is (n-r)xn, such that KX = 0, where X are the fixed effects:
	
	K = numpy.identity(n_pops - 1) - numpy.ones((n_pops - 1, n_pops -1))*1/float(n_pops)
	K = numpy.concatenate( (K, -1*numpy.ones((n_pops - 1, 1))*1./float(n_pops)), axis=1)
	
	for i in range(n_steps_u):
		for j in range(n_steps_e):
			sigma2_u = sigma2_u_vec[i]
			sigma2_e = sigma2_e_vec[j]
			#sigma2_e = 1 - sigma2_u
			#negloglikelihood_vec[i] = negloglikelihood_REML(sigma2_u, sigma2_e, freq_scaled_genotype_mat, centered_phenotypes, K, relatedness_mat)
			#negloglikelihood_vec2[i] = negloglikelihood_REML(sigma2_u, sigma2_e, freq_scaled_genotype_mat, centered_phenotypes, K, numpy.dot(helper_matrix,helper_matrix.T))
			#The function calculates -1*log-likelihood
			negloglikelihood_array[i,j] = negloglikelihood_REML(sigma2_u, sigma2_e, freq_scaled_genotype_mat, centered_phenotypes, K, relatedness_mat)
			#negloglikelihood_array[i,j] = negloglikelihood_REML2(sigma2_u, sigma2_e, centered_phenotypes, relatedness_mat)

	#pt.plot(negloglikelihood_vec)
	CS = pt.contour(numpy.log(negloglikelihood_array),20)
	pt.clabel(CS, inline=1, fontsize=10)
	
	#pt.plot(negloglikelihood_vec2,'r')
	#CS = pt.contour(numpy.log(negloglikelihood_array),20)
	#pt.clabel(CS, inline=1, fontsize=10)
	pt.show()
	min_index = numpy.nanargmin( negloglikelihood_array )
	min_index_e = min_index % n_steps_e
	min_index_u = numpy.floor( min_index/float(n_steps_e) )
	
	##If constraining variance components to sum to 1 during optimization:
	
	#min_index_u = numpy.nanargmin( negloglikelihood_vec )
	#h2_estimate = sigma2_u_estimate/numpy.var(centered_phenotypes)
	
	##Otherwise:
	sigma2_u_estimate = sigma2_u_vec[min_index_u]
	sigma2_e_estimate = sigma2_e_vec[min_index_e]
	
	#e_estimate = sigma2_e_estimate/numpy.var(centered_phenotypes)
	h2_estimate = sigma2_u_estimate/(sigma2_u_estimate + sigma2_e_estimate) #This is what Josh Bloom does... does this make more or less sense than constraining the variance components to sum to one in the model?
	#sigma2_e_estimate = sigma2_e_vec[min_index_e]
	
	#min_index_u2 = numpy.nanargmin( negloglikelihood_vec2 )
	#H2_estimate = sigma2_u_vec[min_index_u2]/numpy.var(centered_phenotypes)
	
	ll_min = negloglikelihood_vec[min_index_u]
	
	h2_lb = numpy.min( sigma2_u_vec[ negloglikelihood_vec - ll_min < numpy.log(10) ]) #2 log-likelihood units; note we calculated log-likelihood*2
	h2_ub = numpy.max( sigma2_u_vec[ negloglikelihood_vec - ll_min < numpy.log(10) ])
	
	print h2_estimate, h2_lb, h2_ub
	#print H2_estimate
	
	return h2_estimate #, h2_lb, h2_ub
	
def negloglikelihood_random_effects_model(sigma2_mu, sigma2_e, freq_scaled_genotype_mat, phenotype_list):
	
	##This function calculates -2 times the log-likelihood associated with the model Y = Zu + e, where e is N(0,sigma_e^2), and u is N(0,sigma_u^2
	##And Z is the genotype matrix, scaled to have mean 0, with variance 1/sqrt(n_loci) down each column. Z is n_pops x n_loci; u has n_loci entries; Y has n_pops entries.
	
	#V: variance-covariance matrix of Y. Under this model, this is V = ZZ^T sigma2_mu + I sigma2_eps
	
	n_pops, n_loci = freq_scaled_genotype_mat.shape
	G = freq_scaled_genotype_mat
	relatedness_mat = numpy.dot( G, G.T )
	
	V = relatedness_mat*sigma2_mu + numpy.identity(n_pops)*sigma2_e
	Vinv = numpy.linalg.inv(V)
	Vdet = numpy.linalg.det(V)
	
	if Vdet == 0:
		Vdet = 1
		
	#Log likelihood times -1, without the initial scale factor that is constant in sigma2_mu and sigma2_e
	
	negloglikelihood = .5*(numpy.log( Vdet ) + numpy.dot( numpy.dot(phenotype_list.T, Vinv), phenotype_list ))
	
	return negloglikelihood
	
def negloglikelihood_REML(sigma2_mu, sigma2_e, freq_scaled_genotype_mat, phenotype_list, K, relatedness_mat):
	
	##This function calculates -1 times the log-likelihood associated with the model Y = \beta + Zu + e, where e is N(0,sigma_e^2), and u is N(0,sigma_u^2
	##And Z is the genotype matrix, scaled to have mean 0 and variance 1. Z is n_pops x n_loci; u has n_loci entries; Y has n_pops entries.
	##In REML, we are calculating the likelihood of KY, where Y is the observed data and K is an (n-r) x n matrix of full rank, chosen such that KX = 0, where X is the matrix of fixed effects.
	##Here, we are only fixing the overall mean, so X is an nx1 vector of ones.
	#V: variance-covariance matrix of Y. Under this model, this is V = ZZ^T sigma2_mu + I sigma2_eps
	
	n_pops = relatedness_mat.shape[0]
	
	V = relatedness_mat*sigma2_mu + numpy.identity(n_pops)*sigma2_e
	c = 1.
	#c = 1/(2*(float(n_pops) - 1))
	K_scaled = K/c
	
	
	KVKT = numpy.dot( numpy.dot( K_scaled, V ), K_scaled.T )
	
	KVKTinv = numpy.linalg.inv(KVKT)*(c**2)
	
	KVKTdet = numpy.linalg.det(KVKT)
	
	#if KVKTdet ==0:
	#	KVKTdet = 1
	
	
	logKVKTdet = 2*(n_pops - 1) * numpy.log( c ) + numpy.log( KVKTdet )
	
	#Log likelihood times -1, without the initial scale factor that is constant in sigma2_mu and sigma2_e
	
	phenotype_list = phenotype_list.reshape((n_pops,1))
	negloglikelihood = .5*(logKVKTdet + numpy.dot( numpy.dot( numpy.dot(K, phenotype_list).T, KVKTinv), numpy.dot(K, phenotype_list) ))
	#print logKVKTdet
	#print numpy.dot( numpy.dot( numpy.dot(K, phenotype_list).T, KVKTinv), numpy.dot(K, phenotype_list) )
	return negloglikelihood

def negloglikelihood_REML(sigma2_mu, sigma2_e1, sigma2_e2, phenotype_list, K, relatedness_mat, error_structure):
	
	##This function calculates -1 times the log-likelihood associated with the model Y = \beta + Zu + e, where e is N(0,sigma_e^2), and u is N(0,sigma_u^2
	##And Z is the genotype matrix, scaled to have mean 0 and variance 1. Z is n_pops x n_loci; u has n_loci entries; Y has n_pops entries.
	##In REML, we are calculating the likelihood of KY, where Y is the observed data and K is an (n-r) x n matrix of full rank, chosen such that KX = 0, where X is the matrix of fixed effects.
	##Here, we are only fixing the overall mean, so X is an nx1 vector of ones.
	#V: variance-covariance matrix of Y. Under this model, this is V = ZZ^T sigma2_mu + I sigma2_eps
	
	n_pops = relatedness_mat.shape[0]
	
	V = relatedness_mat*sigma2_mu + error_structure*sigma2_e1 + numpy.identity(n_pops)*sigma2_e2
	c = 1.
	#c = 1/(2*(float(n_pops) - 1))
	K_scaled = K/c
	
	
	KVKT = numpy.dot( numpy.dot( K_scaled, V ), K_scaled.T )
	
	KVKTinv = numpy.linalg.inv(KVKT)*(c**2)
	
	KVKTdet = numpy.linalg.det(KVKT)
	
	#if KVKTdet ==0:
	#	KVKTdet = 1
	
	
	logKVKTdet = 2*(n_pops - 1) * numpy.log( c ) + numpy.log( KVKTdet )
	
	#Log likelihood times -1, without the initial scale factor that is constant in sigma2_mu and sigma2_e
	
	phenotype_list = phenotype_list.reshape((n_pops,1))
	negloglikelihood = .5*(logKVKTdet + numpy.dot( numpy.dot( numpy.dot(K, phenotype_list).T, KVKTinv), numpy.dot(K, phenotype_list) ))
	#print logKVKTdet
	#print numpy.dot( numpy.dot( numpy.dot(K, phenotype_list).T, KVKTinv), numpy.dot(K, phenotype_list) )
	return negloglikelihood

def negloglikelihood_REML2(sigma2_mu, sigma2_e, phenotype_list, relatedness_mat):
	
	##This function calculates -2 times the log-likelihood associated with the model Y = \beta + Zu + e, where e is N(0,sigma_e^2), and u is N(0,sigma_u^2
	##And Z is the genotype matrix, scaled to have mean 0 and variance 1. Z is n_pops x n_loci; u has n_loci entries; Y has n_pops entries.
	##In REML, we are calculating the likelihood of KY, where Y is the observed data and K is an (n-r) x n matrix of full rank, chosen such that KX = 0, where X is the matrix of fixed effects.
	##The variance-covariance matrix is Pinv, where P = Vinv - Vinv X (X^T Vinv X).inv X^T Vinv where X is the matrix of fixed effects and V is the variance-covariance matrix of Y.
	#V: variance-covariance matrix of Y. Under this model, this is V = ZZ^T sigma2_mu + I sigma2_eps
	
	n_pops = relatedness_mat.shape[0]
	X = numpy.ones((n_pops,1),dtype='float')
	
	V = relatedness_mat*sigma2_mu + numpy.identity(n_pops)*sigma2_e
	Vinv = numpy.linalg.inv(V)
	XTVinvXinv = numpy.linalg.inv( numpy.dot(numpy.dot(X.T, Vinv), X) )
	
	P = Vinv - numpy.dot( numpy.dot( numpy.dot( numpy.dot( Vinv, X ), XTVinvXinv), X.T), Vinv )
	#print 'round'
	#print numpy.linalg.det(P)
	#print numpy.dot( numpy.dot( phenotype_list.T, P), phenotype_list )
	ll = .5*numpy.log( numpy.linalg.det(P) ) - .5*numpy.dot( numpy.dot( phenotype_list.T, P), phenotype_list )
	#print ll
	return -1*ll
	
def narrow_sense_hsq_window(genotype_mat, phenotype_list, window):
	##Calculates narrow-sense heritability for a given phenotype (length n_segs), with genotypes specified by genotype_mat (n_segs x n_loci).
	##This is one iteration; functions with built-in bootstrapping are found below.
	##Calculate relatedness matrix
	#genotype_mat = .5 + genotype_mat
	n_loci = numpy.float(genotype_mat.shape[1])
	n_segs = len(phenotype_list)
	
	rm_allele_frequencies = numpy.mean(genotype_mat, axis = 0)
	by_allele_frequencies = 1 - rm_allele_frequencies
	
	fixed_loci = numpy.arange(n_loci)[(rm_allele_frequencies < .1) + (by_allele_frequencies < .1)]
	genotype_mat = numpy.delete(genotype_mat, fixed_loci, axis = 1)
	rm_allele_frequencies = numpy.delete(rm_allele_frequencies, fixed_loci)
	by_allele_frequencies = numpy.delete(by_allele_frequencies, fixed_loci)
	n_loci -= len(fixed_loci)
	freq_scaled_genotype_mat = -1*by_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*genotype_mat + rm_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*(1 - genotype_mat)
	
	
	relatedness_mat = numpy.dot(freq_scaled_genotype_mat, freq_scaled_genotype_mat.T)/numpy.float(n_loci)
	phenotype_pairwise_difs = numpy.ones((n_segs,n_segs),dtype='float')
	correct_inds = numpy.where(numpy.abs(relatedness_mat) < window) #This just excludes the i=j case, or the comparison between a segregant and itself.
	
	relatedness_list = relatedness_mat[correct_inds] - numpy.mean(relatedness_mat[correct_inds])
	for i in range(n_segs):
		for j in range(n_segs):
			phenotype_pairwise_difs[i,j] = (phenotype_list[i] - phenotype_list[j])**2

	# #Center y data
	phenotype_pairwise_difs_list = phenotype_pairwise_difs[correct_inds] - numpy.mean(phenotype_pairwise_difs[correct_inds])

	#Calculate h^2
	#print n_loci, n_segs
	betahat = 1./numpy.sum(relatedness_list*relatedness_list)*numpy.sum(relatedness_list*phenotype_pairwise_difs_list)
	#(binned_relatedness_list, bin_locs, bin_nums) = scipy.stats.binned_statistic(relatedness_list, phenotype_pairwise_difs_list, bins=100)
	#pt.scatter(bin_locs[:-1], binned_relatedness_list)
 	#pt.plot(numpy.arange(-1,1,.01),betahat*numpy.arange(-1,1,.01),'k')
# 	pt.figure()
# 	pt.scatter(relatedness_list, phenotype_pairwise_difs_list)
# 	pt.plot(numpy.arange(-1,1,.01),betahat*numpy.arange(-1,1,.01),'k')
 	#pt.show()
	#print betahat, numpy.var(phenotype_list)
	
	h2 = betahat/-2./numpy.var(phenotype_list)
	
	return h2

def narrow_sense_hsq_bootstrap(genotype_mat, phenotype_list):
	##Narrow-sense heritability calculations

	# #Calculate relatedness matrix
	#genotype_mat = .5 + genotype_mat
	n_loci = numpy.float(genotype_mat.shape[1])
	n_segs = len(phenotype_list)
	
	rm_allele_frequencies = numpy.mean(genotype_mat, axis = 0)
	by_allele_frequencies = 1 - rm_allele_frequencies
	
	fixed_loci = numpy.arange(n_loci)[(rm_allele_frequencies < .1) + (by_allele_frequencies < .1)]
	genotype_mat = numpy.delete(genotype_mat, fixed_loci, axis = 1)
	rm_allele_frequencies = numpy.delete(rm_allele_frequencies, fixed_loci)
	by_allele_frequencies = numpy.delete(by_allele_frequencies, fixed_loci)
	n_loci -= len(fixed_loci)
	freq_scaled_genotype_mat = -1*by_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*genotype_mat + rm_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*(1 - genotype_mat)
	
	
	relatedness_mat = numpy.dot(freq_scaled_genotype_mat, freq_scaled_genotype_mat.T)/numpy.float(n_loci)
	phenotype_pairwise_difs = numpy.ones((n_segs,n_segs),dtype='float')
	correct_inds = numpy.where(relatedness_mat < .9)
	#lower_tri_indices = numpy.tril_indices(n_segs, 0)
	#lower_tri_indices = numpy.tril_indices(n_segs, -1)
	relatedness_list = relatedness_mat[correct_inds] - numpy.mean(relatedness_mat[correct_inds])
	for i in range(n_segs):
		for j in range(n_segs):
			phenotype_pairwise_difs[i,j] = (phenotype_list[i] - phenotype_list[j])**2

	# #Center y data
	phenotype_pairwise_difs_list = phenotype_pairwise_difs[correct_inds] - numpy.mean(phenotype_pairwise_difs[correct_inds])

	#Calculate h^2
	#print n_loci, n_segs
	betahat = 1./numpy.sum(relatedness_list*relatedness_list)*numpy.sum(relatedness_list*phenotype_pairwise_difs_list)
	#(binned_relatedness_list, bin_locs, bin_nums) = scipy.stats.binned_statistic(relatedness_list, phenotype_pairwise_difs_list, bins=100)
	#pt.scatter(bin_locs[:-1], binned_relatedness_list)
 	#pt.plot(numpy.arange(-1,1,.01),betahat*numpy.arange(-1,1,.01),'k')
# 	pt.figure()
# 	pt.scatter(relatedness_list, phenotype_pairwise_difs_list)
# 	pt.plot(numpy.arange(-1,1,.01),betahat*numpy.arange(-1,1,.01),'k')
 	#pt.show()
	#print betahat, numpy.var(phenotype_list)
	
	h2 = betahat/-2./numpy.var(phenotype_list)
	# #Jacknife (leave one out) to calculate confidence intervals

	h2_bs = []

	for k in range(1000):
		bootstrap_inds = numpy.random.choice(numpy.arange(n_segs), size=n_segs)
		phenotype_pairwise_difs_bs = phenotype_pairwise_difs[bootstrap_inds,:][:,bootstrap_inds]
		relatedness_bs = relatedness_mat[bootstrap_inds,:][:,bootstrap_inds]
		phenotype_list_bs = phenotype_list[bootstrap_inds]
		correct_inds = numpy.where(relatedness_bs < .9)
		#re-center y data
		phenotype_pairwise_difs_list_bs = phenotype_pairwise_difs_bs[correct_inds] - numpy.mean(phenotype_pairwise_difs_bs[correct_inds])
		relatedness_list_bs = relatedness_bs[correct_inds] - numpy.mean(relatedness_bs[correct_inds])
		# #((X^TX)^-1)XTy
		betahat = 1./numpy.sum(relatedness_list_bs*relatedness_list_bs)*numpy.sum(relatedness_list_bs*phenotype_pairwise_difs_list_bs)
		
		h2_bs.append(betahat/-2./numpy.var(phenotype_list_bs))
		
	
	#return h2, [numpy.percentile(h2_bs,2.5), numpy.percentile(h2_bs,97.5)]
	return h2, [numpy.percentile(h2_bs,25), numpy.percentile(h2_bs,75)]

def narrow_sense_hsq_bootstrap_ninetyfive(genotype_mat, phenotype_list):
	##Narrow-sense heritability calculations

	# #Calculate relatedness matrix
	#genotype_mat = .5 + genotype_mat
	n_loci = numpy.float(genotype_mat.shape[1])
	n_segs = len(phenotype_list)
	
	rm_allele_frequencies = numpy.mean(genotype_mat, axis = 0)
	by_allele_frequencies = 1 - rm_allele_frequencies
	
	fixed_loci = numpy.arange(n_loci)[(rm_allele_frequencies < .1) + (by_allele_frequencies < .1)]
	genotype_mat = numpy.delete(genotype_mat, fixed_loci, axis = 1)
	rm_allele_frequencies = numpy.delete(rm_allele_frequencies, fixed_loci)
	by_allele_frequencies = numpy.delete(by_allele_frequencies, fixed_loci)
	n_loci -= len(fixed_loci)
	freq_scaled_genotype_mat = -1*by_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*genotype_mat + rm_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*(1 - genotype_mat)
	#freq_scaled_genotype_mat = -1*genotype_mat + (1 - genotype_mat)
	#freq_scaled_genotype_mat = 2*genotype_mat
	
	relatedness_mat = numpy.dot(freq_scaled_genotype_mat, freq_scaled_genotype_mat.T)/numpy.float(n_loci)
	phenotype_pairwise_difs = numpy.ones((n_segs,n_segs),dtype='float')
	correct_inds = numpy.where(relatedness_mat < .9)
	#lower_tri_indices = numpy.tril_indices(n_segs, 0)
	#lower_tri_indices = numpy.tril_indices(n_segs, -1)
	relatedness_list = relatedness_mat[correct_inds] - numpy.mean(relatedness_mat[correct_inds])
	for i in range(n_segs):
		for j in range(n_segs):
			phenotype_pairwise_difs[i,j] = (phenotype_list[i] - phenotype_list[j])**2

	# #Center y data
	phenotype_pairwise_difs_list = phenotype_pairwise_difs[correct_inds] - numpy.mean(phenotype_pairwise_difs[correct_inds])

	#Calculate h^2
	#print n_loci, n_segs
	betahat = 1./numpy.sum(relatedness_list*relatedness_list)*numpy.sum(relatedness_list*phenotype_pairwise_difs_list)
	#(binned_relatedness_list, bin_locs, bin_nums) = scipy.stats.binned_statistic(relatedness_list, phenotype_pairwise_difs_list, bins=100)
	#pt.scatter(bin_locs[:-1], binned_relatedness_list)
 	#pt.plot(numpy.arange(-1,1,.01),betahat*numpy.arange(-1,1,.01),'k')
# 	pt.figure()
# 	pt.scatter(relatedness_list, phenotype_pairwise_difs_list)
# 	pt.plot(numpy.arange(-1,1,.01),betahat*numpy.arange(-1,1,.01),'k')
 	#pt.show()
	#print betahat, numpy.var(phenotype_list)
	
	h2 = betahat/-2./numpy.var(phenotype_list)
	# #Jacknife (leave one out) to calculate confidence intervals

	h2_bs = []

	for k in range(1000):
		bootstrap_inds = numpy.random.choice(numpy.arange(n_segs), size=n_segs)
		phenotype_pairwise_difs_bs = phenotype_pairwise_difs[bootstrap_inds,:][:,bootstrap_inds]
		relatedness_bs = relatedness_mat[bootstrap_inds,:][:,bootstrap_inds]
		phenotype_list_bs = phenotype_list[bootstrap_inds]
		correct_inds = numpy.where(relatedness_bs < .9)
		#re-center y data
		phenotype_pairwise_difs_list_bs = phenotype_pairwise_difs_bs[correct_inds] - numpy.mean(phenotype_pairwise_difs_bs[correct_inds])
		relatedness_list_bs = relatedness_bs[correct_inds] - numpy.mean(relatedness_bs[correct_inds])
		# #((X^TX)^-1)XTy
		betahat = 1./numpy.sum(relatedness_list_bs*relatedness_list_bs)*numpy.sum(relatedness_list_bs*phenotype_pairwise_difs_list_bs)
		
		h2_bs.append(betahat/-2./numpy.var(phenotype_list_bs))
		
	
	return h2, [numpy.percentile(h2_bs,2.5), numpy.percentile(h2_bs,97.5)]
	#return h2, [numpy.percentile(h2_bs,25), numpy.percentile(h2_bs,75)]

def narrow_sense_hsq_replicates(genotype_mat, phenotype_list, helper_matrix):
	##Narrow-sense heritability calculations

	# #Calculate relatedness matrix
	#genotype_mat = .5 + genotype_mat
	n_loci = numpy.float(genotype_mat.shape[1])
	n_segs = numpy.int(genotype_mat.shape[0])
	n_pops = len(phenotype_list)
	
	expanded_genotype_mat = numpy.dot(helper_matrix,genotype_mat)
	
	rm_allele_frequencies = numpy.mean(expanded_genotype_mat, axis = 0)
	by_allele_frequencies = 1 - rm_allele_frequencies
	
	fixed_loci = numpy.arange(n_loci)[(rm_allele_frequencies < .1) + (by_allele_frequencies < .1)]
	expanded_genotype_mat = numpy.delete(expanded_genotype_mat, fixed_loci, axis = 1)
	rm_allele_frequencies = numpy.delete(rm_allele_frequencies, fixed_loci)
	by_allele_frequencies = numpy.delete(by_allele_frequencies, fixed_loci)
	n_loci -= len(fixed_loci)
	freq_scaled_genotype_mat = -1*by_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*expanded_genotype_mat + rm_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*(1 - expanded_genotype_mat)
	#freq_scaled_genotype_mat = -1*expanded_genotype_mat + (1 - expanded_genotype_mat)
	#freq_scaled_genotype_mat = 2*expanded_genotype_mat
	#print freq_scaled_genotype_mat
	relatedness_mat = numpy.dot(freq_scaled_genotype_mat, freq_scaled_genotype_mat.T)/numpy.float(n_loci)
	
	correct_inds = numpy.where(relatedness_mat < .9)
	phenotype_pairwise_difs = numpy.ones((n_pops,n_pops),dtype='float')
	
	relatedness_list = relatedness_mat[correct_inds] - numpy.mean(relatedness_mat[correct_inds],dtype='float')
	
	for i in range(n_pops):
		for j in range(n_pops):
			phenotype_pairwise_difs[i,j] = (phenotype_list[i] - phenotype_list[j])**2

	# #Center y data

	phenotype_pairwise_difs_list = phenotype_pairwise_difs[correct_inds] - numpy.mean(phenotype_pairwise_difs[correct_inds])
	
	#Calculate h^2

	betahat = 1./numpy.sum(relatedness_list*relatedness_list)*numpy.sum(relatedness_list*phenotype_pairwise_difs_list)
	# (binned_relatedness_list, bin_locs, bin_nums) = scipy.stats.binned_statistic(relatedness_list, phenotype_pairwise_difs_list, bins=100)
	# pt.scatter(bin_locs[:-1], binned_relatedness_list)
	# pt.plot(numpy.arange(-1,1,.01),betahat*numpy.arange(-1,1,.01),'k')
	
	# pt.show()
	h2 = betahat/-2./numpy.var(phenotype_list)
	
	return h2

def narrow_sense_hsq_replicates_window(genotype_mat, phenotype_list, helper_matrix, window):
	##Narrow-sense heritability calculations

	# #Calculate relatedness matrix
	#genotype_mat = .5 + genotype_mat
	n_loci = numpy.float(genotype_mat.shape[1])
	n_segs = numpy.int(genotype_mat.shape[0])
	n_pops = len(phenotype_list)
	
	expanded_genotype_mat = numpy.dot(helper_matrix,genotype_mat)
	
	rm_allele_frequencies = numpy.mean(expanded_genotype_mat, axis = 0)
	by_allele_frequencies = 1 - rm_allele_frequencies
	
	fixed_loci = numpy.arange(n_loci)[(rm_allele_frequencies < .1) + (by_allele_frequencies < .1)]
	expanded_genotype_mat = numpy.delete(expanded_genotype_mat, fixed_loci, axis = 1)
	rm_allele_frequencies = numpy.delete(rm_allele_frequencies, fixed_loci)
	by_allele_frequencies = numpy.delete(by_allele_frequencies, fixed_loci)
	n_loci -= len(fixed_loci)
	freq_scaled_genotype_mat = -1*by_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*expanded_genotype_mat + rm_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*(1 - expanded_genotype_mat)
	#freq_scaled_genotype_mat = -1*expanded_genotype_mat + (1 - expanded_genotype_mat)
	#freq_scaled_genotype_mat = 2*expanded_genotype_mat
	#print freq_scaled_genotype_mat
	relatedness_mat = numpy.dot(freq_scaled_genotype_mat, freq_scaled_genotype_mat.T)/numpy.float(n_loci)
	
	correct_inds = numpy.where(numpy.abs(relatedness_mat) < window)
	phenotype_pairwise_difs = numpy.ones((n_pops,n_pops),dtype='float')
	
	relatedness_list = relatedness_mat[correct_inds] - numpy.mean(relatedness_mat[correct_inds],dtype='float')
	
	for i in range(n_pops):
		for j in range(n_pops):
			phenotype_pairwise_difs[i,j] = (phenotype_list[i] - phenotype_list[j])**2

	# #Center y data

	phenotype_pairwise_difs_list = phenotype_pairwise_difs[correct_inds] - numpy.mean(phenotype_pairwise_difs[correct_inds])
	
	#Calculate h^2

	betahat = 1./numpy.sum(relatedness_list*relatedness_list)*numpy.sum(relatedness_list*phenotype_pairwise_difs_list)
	# (binned_relatedness_list, bin_locs, bin_nums) = scipy.stats.binned_statistic(relatedness_list, phenotype_pairwise_difs_list, bins=100)
	# pt.scatter(bin_locs[:-1], binned_relatedness_list)
	# pt.plot(numpy.arange(-1,1,.01),betahat*numpy.arange(-1,1,.01),'k')
	
	# pt.show()
	h2 = betahat/-2./numpy.var(phenotype_list)
	
	return h2
	
def narrow_sense_hsq_replicates_bootstrap(genotype_mat, phenotype_list, helper_matrix):
	##Narrow-sense heritability calculations

	# #Calculate relatedness matrix
	#genotype_mat = .5 + genotype_mat
	n_loci = numpy.float(genotype_mat.shape[1])
	n_segs = numpy.int(genotype_mat.shape[0])
	n_pops = len(phenotype_list)
	
	expanded_genotype_mat = numpy.dot(helper_matrix,genotype_mat)
	
	rm_allele_frequencies = numpy.mean(expanded_genotype_mat, axis = 0)
	by_allele_frequencies = 1 - rm_allele_frequencies
	
	fixed_loci = numpy.arange(n_loci)[(rm_allele_frequencies < .1) + (by_allele_frequencies < .1)]
	expanded_genotype_mat = numpy.delete(expanded_genotype_mat, fixed_loci, axis = 1)
	rm_allele_frequencies = numpy.delete(rm_allele_frequencies, fixed_loci)
	by_allele_frequencies = numpy.delete(by_allele_frequencies, fixed_loci)
	n_loci -= len(fixed_loci)
	freq_scaled_genotype_mat = -1*by_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*expanded_genotype_mat + rm_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*(1 - expanded_genotype_mat)
	#freq_scaled_genotype_mat = -1*expanded_genotype_mat + (1 - expanded_genotype_mat)
	#freq_scaled_genotype_mat = 2*expanded_genotype_mat
	#print freq_scaled_genotype_mat
	relatedness_mat = numpy.dot(freq_scaled_genotype_mat, freq_scaled_genotype_mat.T)/numpy.float(n_loci)
	
	correct_inds = numpy.where(relatedness_mat < .9)
	phenotype_pairwise_difs = numpy.ones((n_pops,n_pops),dtype='float')
	
	relatedness_list = relatedness_mat[correct_inds] - numpy.mean(relatedness_mat[correct_inds],dtype='float')
	
	for i in range(n_pops):
		for j in range(n_pops):
			phenotype_pairwise_difs[i,j] = (phenotype_list[i] - phenotype_list[j])**2

	# #Center y data

	phenotype_pairwise_difs_list = phenotype_pairwise_difs[correct_inds] - numpy.mean(phenotype_pairwise_difs[correct_inds])
	
	#Calculate h^2

	betahat = 1./numpy.sum(relatedness_list*relatedness_list)*numpy.sum(relatedness_list*phenotype_pairwise_difs_list)
	# (binned_relatedness_list, bin_locs, bin_nums) = scipy.stats.binned_statistic(relatedness_list, phenotype_pairwise_difs_list, bins=100)
	# pt.scatter(bin_locs[:-1], binned_relatedness_list)
	# pt.plot(numpy.arange(-1,1,.01),betahat*numpy.arange(-1,1,.01),'k')
	
	# pt.show()
	h2 = betahat/-2./numpy.var(phenotype_list)
	# #Bootstrap to calculate confidence intervals

	h2_bs = []

	for k in range(1000):
		bootstrap_inds = numpy.random.choice(numpy.arange(n_segs), size=n_segs)
		expanded_indices = []
		for bootstrap_ind in bootstrap_inds:
			expanded_indices.extend(numpy.where(helper_matrix[:,bootstrap_ind] == 1)[0])
		
		phenotype_pairwise_difs_bs = phenotype_pairwise_difs[expanded_indices, :][:, expanded_indices]
		#print phenotype_pairwise_difs_bs.shape
		relatedness_bs = relatedness_mat[expanded_indices, :][:, expanded_indices]
		phenotype_list_bs = phenotype_list[expanded_indices]
		correct_inds_bs = numpy.where(relatedness_bs < .9)
		#re-center y data
		phenotype_pairwise_difs_list_bs = (phenotype_pairwise_difs_bs[correct_inds_bs] - numpy.mean(phenotype_pairwise_difs_bs[correct_inds_bs])).flatten()
		relatedness_list_bs = (relatedness_bs[correct_inds_bs] - numpy.mean(relatedness_bs[correct_inds_bs],dtype='float')).flatten()
		# #((X^TX)^-1)XTy
		betahat = 1./numpy.sum(relatedness_list_bs*relatedness_list_bs)*numpy.sum(relatedness_list_bs*phenotype_pairwise_difs_list_bs)
		
		
		h2_bs.append(betahat/-2./numpy.var(phenotype_list_bs))
		
	
	#return h2, [numpy.percentile(h2_bs,2.5), numpy.percentile(h2_bs,97.5)]
	return h2, [numpy.percentile(h2_bs,25), numpy.percentile(h2_bs,75)]

def narrow_sense_hsq_replicates_bootstrap_ninetyfive(genotype_mat, phenotype_list, helper_matrix):
	##Narrow-sense heritability calculations

	# #Calculate relatedness matrix
	#genotype_mat = .5 + genotype_mat
	n_loci = numpy.float(genotype_mat.shape[1])
	n_segs = numpy.int(genotype_mat.shape[0])
	n_pops = len(phenotype_list)
	
	expanded_genotype_mat = numpy.dot(helper_matrix,genotype_mat)
	
	rm_allele_frequencies = numpy.mean(expanded_genotype_mat, axis = 0)
	by_allele_frequencies = 1 - rm_allele_frequencies
	
	fixed_loci = numpy.arange(n_loci)[(rm_allele_frequencies < .1) + (by_allele_frequencies < .1)]
	expanded_genotype_mat = numpy.delete(expanded_genotype_mat, fixed_loci, axis = 1)
	rm_allele_frequencies = numpy.delete(rm_allele_frequencies, fixed_loci)
	by_allele_frequencies = numpy.delete(by_allele_frequencies, fixed_loci)
	n_loci -= len(fixed_loci)
	freq_scaled_genotype_mat = -1*by_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*expanded_genotype_mat + rm_allele_frequencies/numpy.sqrt(by_allele_frequencies*rm_allele_frequencies)*(1 - expanded_genotype_mat)
	#freq_scaled_genotype_mat = -1*expanded_genotype_mat + (1 - expanded_genotype_mat)
	#freq_scaled_genotype_mat = 2*expanded_genotype_mat
	#print freq_scaled_genotype_mat
	relatedness_mat = numpy.dot(freq_scaled_genotype_mat, freq_scaled_genotype_mat.T)/numpy.float(n_loci)
	
	correct_inds = numpy.where(relatedness_mat < .9)
	phenotype_pairwise_difs = numpy.ones((n_pops,n_pops),dtype='float')
	
	relatedness_list = relatedness_mat[correct_inds] - numpy.mean(relatedness_mat[correct_inds],dtype='float')
	
	for i in range(n_pops):
		for j in range(n_pops):
			phenotype_pairwise_difs[i,j] = (phenotype_list[i] - phenotype_list[j])**2

	# #Center y data

	phenotype_pairwise_difs_list = phenotype_pairwise_difs[correct_inds] - numpy.mean(phenotype_pairwise_difs[correct_inds])
	
	#Calculate h^2

	betahat = 1./numpy.sum(relatedness_list*relatedness_list)*numpy.sum(relatedness_list*phenotype_pairwise_difs_list)
	# (binned_relatedness_list, bin_locs, bin_nums) = scipy.stats.binned_statistic(relatedness_list, phenotype_pairwise_difs_list, bins=100)
	# pt.scatter(bin_locs[:-1], binned_relatedness_list)
	# pt.plot(numpy.arange(-1,1,.01),betahat*numpy.arange(-1,1,.01),'k')
	
	# pt.show()
	h2 = betahat/-2./numpy.var(phenotype_list)
	# #Bootstrap to calculate confidence intervals

	h2_bs = []

	for k in range(1000):
		bootstrap_inds = numpy.random.choice(numpy.arange(n_segs), size=n_segs)
		expanded_indices = []
		for bootstrap_ind in bootstrap_inds:
			expanded_indices.extend(numpy.where(helper_matrix[:,bootstrap_ind] == 1)[0])
		
		phenotype_pairwise_difs_bs = phenotype_pairwise_difs[expanded_indices, :][:, expanded_indices]
		#print phenotype_pairwise_difs_bs.shape
		relatedness_bs = relatedness_mat[expanded_indices, :][:, expanded_indices]
		phenotype_list_bs = phenotype_list[expanded_indices]
		correct_inds_bs = numpy.where(relatedness_bs < .9)
		#re-center y data
		phenotype_pairwise_difs_list_bs = (phenotype_pairwise_difs_bs[correct_inds_bs] - numpy.mean(phenotype_pairwise_difs_bs[correct_inds_bs])).flatten()
		relatedness_list_bs = (relatedness_bs[correct_inds_bs] - numpy.mean(relatedness_bs[correct_inds_bs],dtype='float')).flatten()
		# #((X^TX)^-1)XTy
		betahat = 1./numpy.sum(relatedness_list_bs*relatedness_list_bs)*numpy.sum(relatedness_list_bs*phenotype_pairwise_difs_list_bs)
		
		
		h2_bs.append(betahat/-2./numpy.var(phenotype_list_bs))
		
	
	return h2, [numpy.percentile(h2_bs,2.5), numpy.percentile(h2_bs,97.5)]
	#return h2, [numpy.percentile(h2_bs,25), numpy.percentile(h2_bs,75)]


def rsq(init_fit_vector, delta_fit_vector, helper_matrix):

	rsq = scipy.stats.pearsonr(numpy.dot(helper_matrix,init_fit_vector),delta_fit_vector)[0]**2
	
	return rsq
	
def rsq_with_bootstrap(init_fit_vector, delta_fit_vector, helper_matrix):

	rsq = scipy.stats.pearsonr(numpy.dot(helper_matrix,init_fit_vector),delta_fit_vector)[0]**2
	
	n_segs = len(init_fit_vector)
	rsq_bs = []
	for k in range(1000):
		bootstrap_inds = numpy.random.choice(numpy.arange(n_segs), size=n_segs)
		expanded_indices = []
		for bootstrap_ind in bootstrap_inds:
			expanded_indices.extend(numpy.where(helper_matrix[:,bootstrap_ind] == 1)[0])
		init_fit_vector_bs = init_fit_vector[bootstrap_inds]
		delta_fit_vector_bs = delta_fit_vector[expanded_indices]
		helper_matrix_bs = helper_matrix[expanded_indices,:][:,bootstrap_inds]
		rsq_bs.append(scipy.stats.pearsonr(numpy.dot(helper_matrix_bs,init_fit_vector_bs),delta_fit_vector_bs)[0]**2)
	
	#return rsq, [numpy.percentile(rsq_bs,2.5), numpy.percentile(rsq_bs,97.5)]
	return rsq, [numpy.percentile(rsq_bs,25), numpy.percentile(rsq_bs,75)]

def rsq_linear_model(predictor_mat, delta_fit_vector, helper_matrix):
	predictor_mat_expanded = numpy.dot(helper_matrix, predictor_mat)
	#print predictor_mat_expanded.shape
	beta = numpy.dot(numpy.linalg.inv(numpy.dot(predictor_mat_expanded.T, predictor_mat_expanded)), numpy.dot(predictor_mat_expanded.T, delta_fit_vector))
	rsq = scipy.stats.pearsonr(numpy.dot(predictor_mat_expanded, beta),delta_fit_vector)[0]**2
	
	return rsq
	
def rsq_with_bootstrap_linear_model(predictor_mat, delta_fit_vector, helper_matrix):
	predictor_mat_expanded = numpy.dot(helper_matrix, predictor_mat)
	beta = numpy.dot(numpy.linalg.inv(numpy.dot(predictor_mat_expanded.T, predictor_mat_expanded)), numpy.dot(predictor_mat_expanded.T, delta_fit_vector))
	rsq = scipy.stats.pearsonr(numpy.dot(predictor_mat_expanded, beta),delta_fit_vector)[0]**2
	
	n_segs = len(helper_matrix[0,:])
	rsq_bs = []
	for k in range(100):
		bootstrap_inds = numpy.random.choice(numpy.arange(n_segs), size=n_segs)
		expanded_indices = []
		for bootstrap_ind in bootstrap_inds:
			expanded_indices.extend(numpy.where(helper_matrix[:,bootstrap_ind] == 1)[0])
		predictor_mat_expanded_bs = predictor_mat_expanded[expanded_indices,:]
		delta_fit_vector_bs = delta_fit_vector[expanded_indices]
		beta_bs = numpy.dot(numpy.linalg.inv(numpy.dot(predictor_mat_expanded_bs.T, predictor_mat_expanded_bs)), numpy.dot(predictor_mat_expanded_bs.T, delta_fit_vector_bs))
		rsq_bs.append(scipy.stats.pearsonr(numpy.dot(predictor_mat_expanded_bs, beta_bs),delta_fit_vector_bs)[0]**2)
	
	#return rsq, [numpy.percentile(rsq_bs,2.5), numpy.percentile(rsq_bs,97.5)]
	return rsq, [numpy.percentile(rsq_bs,25), numpy.percentile(rsq_bs,75)]

def rsq_with_bootstrap_linear_model_ninetyfive(predictor_mat, delta_fit_vector, helper_matrix):
	predictor_mat_expanded = numpy.dot(helper_matrix, predictor_mat)
	beta = numpy.dot(numpy.linalg.inv(numpy.dot(predictor_mat_expanded.T, predictor_mat_expanded)), numpy.dot(predictor_mat_expanded.T, delta_fit_vector))
	rsq = scipy.stats.pearsonr(numpy.dot(predictor_mat_expanded, beta),delta_fit_vector)[0]**2
	
	n_segs = len(helper_matrix[0,:])
	rsq_bs = []
	for k in range(1000):
		bootstrap_inds = numpy.random.choice(numpy.arange(n_segs), size=n_segs)
		expanded_indices = []
		for bootstrap_ind in bootstrap_inds:
			expanded_indices.extend(numpy.where(helper_matrix[:,bootstrap_ind] == 1)[0])
		predictor_mat_expanded_bs = predictor_mat_expanded[expanded_indices,:]
		delta_fit_vector_bs = delta_fit_vector[expanded_indices]
		beta_bs = numpy.dot(numpy.linalg.inv(numpy.dot(predictor_mat_expanded_bs.T, predictor_mat_expanded_bs)), numpy.dot(predictor_mat_expanded_bs.T, delta_fit_vector_bs))
		rsq_bs.append(scipy.stats.pearsonr(numpy.dot(predictor_mat_expanded_bs, beta_bs),delta_fit_vector_bs)[0]**2)
	
	return rsq, [numpy.percentile(rsq_bs,2.5), numpy.percentile(rsq_bs,97.5)]
	#return rsq, [numpy.percentile(rsq_bs,25), numpy.percentile(rsq_bs,75)]
	
def broad_sense_Hsq(phenotypes, pops_per_seg, non_heritable_variances, helper_matrix):
	
	non_heritable = 1./numpy.sum(pops_per_seg)*numpy.nansum(pops_per_seg*non_heritable_variances)
	Hsq = (1 - non_heritable/numpy.var(phenotypes))
	
	return Hsq

def broad_sense_Hsq_means(phenotypes, pops_per_seg, non_heritable_variances, helper_matrix):
	
	non_heritable = 1./numpy.sum(pops_per_seg)*numpy.nansum(non_heritable_variances)
	Hsq = (1 - non_heritable/numpy.var(phenotypes))
	
	return Hsq
	
def broad_sense_Hsq_with_bootstrap(phenotypes, pops_per_seg, non_heritable_variances, helper_matrix):
	
	non_heritable = 1./numpy.sum(pops_per_seg)*numpy.nansum(pops_per_seg*non_heritable_variances)
	Hsq = (1 - non_heritable/numpy.var(phenotypes))
	
	n_segs = len(pops_per_seg)
	Hsq_bs = []
	for k in range(1000):
		bootstrap_inds = numpy.random.choice(numpy.arange(n_segs), size=n_segs)
		expanded_indices = []
		for bootstrap_ind in bootstrap_inds:
			expanded_indices.extend(numpy.where(helper_matrix[:,bootstrap_ind] == 1)[0])
		phenotypes_bs = phenotypes[expanded_indices]
		non_heritable_variances_bs = non_heritable_variances[bootstrap_inds]
		pops_per_seg_bs = pops_per_seg[bootstrap_inds]
		non_heritable = 1./numpy.sum(pops_per_seg_bs)*numpy.nansum(pops_per_seg_bs*non_heritable_variances_bs)
		Hsq_bs.append(1 - non_heritable/numpy.var(phenotypes_bs))
	
	#return Hsq, [numpy.percentile(Hsq_bs,2.5),numpy.percentile(Hsq_bs,97.5)]
	return Hsq, [numpy.percentile(Hsq_bs,25),numpy.percentile(Hsq_bs,75)]
	
def broad_sense_Hsq_with_bootstrap_ninetyfive(phenotypes, pops_per_seg, non_heritable_variances, helper_matrix):
	
	non_heritable = 1./numpy.sum(pops_per_seg)*numpy.nansum(pops_per_seg*non_heritable_variances)
	Hsq = (1 - non_heritable/numpy.var(phenotypes))
	
	n_segs = len(pops_per_seg)
	Hsq_bs = []
	for k in range(1000):
		bootstrap_inds = numpy.random.choice(numpy.arange(n_segs), size=n_segs)
		expanded_indices = []
		for bootstrap_ind in bootstrap_inds:
			expanded_indices.extend(numpy.where(helper_matrix[:,bootstrap_ind] == 1)[0])
		phenotypes_bs = phenotypes[expanded_indices]
		non_heritable_variances_bs = non_heritable_variances[bootstrap_inds]
		pops_per_seg_bs = pops_per_seg[bootstrap_inds]
		non_heritable = 1./numpy.sum(pops_per_seg_bs)*numpy.nansum(pops_per_seg_bs*non_heritable_variances_bs)
		Hsq_bs.append(1 - non_heritable/numpy.var(phenotypes_bs))
	
	return Hsq, [numpy.percentile(Hsq_bs,2.5),numpy.percentile(Hsq_bs,97.5)]
	#return Hsq, [numpy.percentile(Hsq_bs,25),numpy.percentile(Hsq_bs,75)]

def broad_sense_Hsq_init_fits(init_fits, std_errs):
	num_segs = numpy.sum(1 - numpy.isnan(std_errs))
	H2 = 1 - 1./num_segs*numpy.nansum(std_errs**2)/numpy.var(init_fits)
	
	H2_bs = []
	for k in range(1000):
		bootstrap_inds = numpy.random.choice(numpy.arange(num_segs), size=num_segs)
		init_fits_bs = init_fits[bootstrap_inds]
		std_errs_bs = std_errs[bootstrap_inds]
		H2_bs.append( 1 - 1./num_segs*numpy.nansum(std_errs_bs**2)/numpy.var(init_fits_bs) )
	
	#return H2, [numpy.percentile(H2_bs,2.5),numpy.percentile(H2_bs,97.5)]
	return H2, [numpy.percentile(H2_bs,25),numpy.percentile(H2_bs,75)]