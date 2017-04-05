import sys
import numpy
from scipy.linalg import lstsq

from scipy.interpolate import UnivariateSpline
from numpy.random import normal
import matplotlib.pylab as pt

class KernelSmoother:

    def __init__(self, U, Y, W, dh=0.05):
        
        self.h =  dh*(U.max()-U.min())
        self.U = U
        self.Y = Y
        self.W = W*1.0/W.sum()
        self.Wtot = W.sum()
        self.weighted_Y = self.Y*self.W
        
    def __call__(self, u):
    
        u_residuals = (numpy.tile(u,(len(self.U),1)).T-numpy.tile(self.U,(len(u),1)))/self.h
        log_weights = -numpy.square(u_residuals)/2
        log_weights = (log_weights - numpy.tile(log_weights.max(axis=1),(len(self.U),1)).T)
        weights = numpy.tile(self.W,(len(u),1))*numpy.exp(log_weights)
        weights = weights / ( numpy.tile(weights.sum(axis=1), (len(self.U),1)).T )
        
        weighted_Y = numpy.dot(weights,self.Y)
        gu = weighted_Y 
        
        #u_residuals = (numpy.tile(u,(len(self.U),1)).T-numpy.tile(self.U,(len(u),1)))/self.h
        #log_weights = -numpy.square(u_residuals)/2
        #log_weights = (log_weights - numpy.tile(log_weights.max(axis=1),(len(self.U),1)).T)
        #weights = numpy.tile(self.W,(len(u),1))*numpy.exp(log_weights)
        
        #gu = (numpy.tile(self.Y,(len(u),1))*weights).sum(axis=1)/weights.sum(axis=1)
        
        return gu
        
    def derivative(self, u):
        
        u_residuals = (numpy.tile(u,(len(self.U),1)).T-numpy.tile(self.U,(len(u),1)))/self.h
        log_weights = -numpy.square(u_residuals)/2
        log_weights = (log_weights - numpy.tile(log_weights.max(axis=1),(len(self.U),1)).T)
        weights = numpy.tile(self.W,(len(u),1))*numpy.exp(log_weights)
        
        gu = (numpy.tile(self.Y,(len(u),1))*weights).sum(axis=1)/weights.sum(axis=1)
        
        dgu =  ((numpy.tile(self.Y,(len(u),1))-numpy.tile(gu,(len(self.U),1)).T)*(-1.0*u_residuals)*weights).sum(axis=1)/weights.sum(axis=1)/self.h
        
        return dgu
        
class KernelRegressionSmoother:

    def __init__(self, U, Y, W, dh=0.05):
        
        self.h =  dh*(U.max()-U.min())
        self.U = U
        self.Y = Y
        self.W = W*1.0/W.sum()
        
        
    def __call__(self, u):
        
        u_residuals = (numpy.tile(u,(len(self.U),1)).T-numpy.tile(self.U,(len(u),1)))/self.h
        log_weights = -numpy.square(u_residuals)/2
        log_weights = (log_weights - numpy.tile(log_weights.max(axis=1),(len(self.U),1)).T)
        weights = numpy.tile(self.W,(len(u),1))*numpy.exp(log_weights)
        weights = weights / ( numpy.tile(weights.sum(axis=1), (len(self.U),1)).T )
        
        weighted_Y = numpy.dot(weights,self.Y)
        weighted_U = numpy.dot(weights,self.U)
        weighted_YU = numpy.dot(weights, self.Y*self.U)
        weighted_U2 = numpy.dot(weights, self.U*self.U)
        
        ms = ( weighted_YU - weighted_Y*weighted_U ) / (weighted_U2 - weighted_U * weighted_U )
        bs = weighted_Y - ms * weighted_U 
        
        gu = bs + ms*u
        return gu
        
    def derivative(self, u):
        
        u_residuals = (numpy.tile(u,(len(self.U),1)).T-numpy.tile(self.U,(len(u),1)))/self.h
        log_weights = -numpy.square(u_residuals)/2
        log_weights = (log_weights - numpy.tile(log_weights.max(axis=1),(len(self.U),1)).T)
        weights = numpy.tile(self.W,(len(u),1))*numpy.exp(log_weights)
        weights = weights / ( numpy.tile(weights.sum(axis=1), (len(self.U),1)).T )
        
        weighted_Y = numpy.dot(weights,self.Y)
        weighted_U = numpy.dot(weights,self.U)
        weighted_YU = numpy.dot(weights, self.Y*self.U)
        weighted_U2 = numpy.dot(weights, self.U*self.U)
        
        ms = ( weighted_YU - weighted_Y*weighted_U ) / (weighted_U2 - weighted_U * weighted_U )
        bs = weighted_Y - ms * weighted_U
        
        gu = bs + ms*u
        
        dweights = -1*u_residuals*weights/self.h
        dweights = dweights - weights * numpy.tile( dweights.sum(axis=1), (len(self.U),1)).T
        
        dweighted_Y = numpy.dot(dweights,self.Y)
        dweighted_U = numpy.dot(dweights,self.U)
        dweighted_YU = numpy.dot(dweights,self.Y*self.U)
        dweighted_U2 = numpy.dot(dweights,self.U*self.U)
        
        dms = ( dweighted_YU - weighted_U*dweighted_Y - weighted_Y*dweighted_U - ms*(dweighted_U2-2*weighted_U*dweighted_U) ) / (weighted_U2-weighted_U*weighted_U)
        
        dgu = dweighted_Y + ms*(1-dweighted_U) + (u-weighted_U)*dms
        
        #dgu = (ms + 2 * weighted_U*dweighted_U - dweighted_U - dweighted_U2 + (u-weighted_U)*(dweighted_YU-weighted_U*dweighted_Y-weighted_Y*dweighted_U)/(weighted_U2-weighted_U*weighted_U) )
        
        return dgu

    
def ordinary_linear_regression(Y,X,W=None):
    # Y = (n,1) vector/matrix
    # X = (n,k) matrix
    # W = (n,1) vector/matrix of weights (not normalized)
    
    
    if W==None:
        W=numpy.ones_like(Y)*1.0
    else:
        W=W*1.0
    if len(X.shape) > 1.5:    
    	X_int = numpy.ones((X.shape[0],X.shape[1]+1))
    	X_int[:,1:] = X
    else:
    	X_int = numpy.ones((X.shape[0],2))
    	X_int[:,1] = X
    
    A = numpy.dot(X_int.T,numpy.dot(numpy.diag(W),X_int))
    b = numpy.dot(Y*W,X_int)
    beta = lstsq(A,b)[0]
    print beta
    beta0 = beta[0]
    beta = beta[1:]
    betanorm = numpy.square(beta).sum()**0.5
    beta = beta/betanorm
    
    return beta, betanorm, lambda x: betanorm*x+beta0
       
def generalized_linear_regression(Y,X,W=None,varM=None,type='kernel_regression',num_total_iters=10,num_beta_iters=10):
    
    if W==None:
        W=numpy.ones_like(Y)*1.0
    else:
        W=W*1.0
    
    if varM == None:
        # set varM to be 5% of variance
        varM = 0.05*numpy.var(Y)
    
    # To get an initial estimate of the coefficients,
    # we first estimate beta when when F(u) = u+b
    beta,betanorm,dummy_F = ordinary_linear_regression(Y,X,W)
    
    # iterate over estimating F(u) and beta
    for current_total_iters in xrange(0,num_total_iters):
        
        U = numpy.dot(X,beta)
        Urange = (U.max()-U.min())
        U += normal(0,Urange*1e-04,size=len(U))
        
        # estimate global trait function
        if type=='spline':
            sorted_U, sorted_Y = (list(x) for x in zip(*sorted(zip(U, Y))))
    
            sorted_U = numpy.array(sorted_U)
            sorted_Y = numpy.array(sorted_Y)
    
            smoothing_factor = 0.43*len(sorted_Y)*(sorted_Y.std()**2)
            g = UnivariateSpline(sorted_U,sorted_Y,s=smoothing_factor)
            gprime = g.derivative()
        elif type=='kernel_regression':
            g = KernelRegressionSmoother(U,Y,W,dh=0.05)
            gprime = lambda x: g.derivative(x)   
        else:   
            g = KernelSmoother(U,Y,W)
            gprime = lambda x: g.derivative(x)
    
        residuals = Y-g(U)
        gprimes = gprime(U)
        # estimate variance weight
    
        varA = 0*((numpy.square(residuals)*numpy.square(gprimes)).sum()-varM*numpy.square(gprimes).sum())/(numpy.square(gprimes)*numpy.square(gprimes)).sum()
    
        # estimate coefficients 
        betahat = beta
        for current_beta_iters in xrange(0,num_beta_iters):
            Uhat = numpy.dot(X,betahat)
            residuals = Y-g(Uhat)
            
            gprimes = gprime(Uhat)
            gprimesquareds = numpy.square(gprimes)
            variances = varM+varA*gprimesquareds
            A = numpy.dot(X.T,numpy.dot(numpy.diag(W*gprimesquareds/variances), X))
            b = numpy.dot(residuals*W*gprimes/variances, X)
            
            #gprimes = numpy.tile(gprime(Uhat),(X.shape[1],1)).T
            #Xhat = numpy.multiply(X,gprimes)    
            #A = numpy.dot(Xhat.T,numpy.dot(numpy.diag(W),Xhat))
            #b = numpy.dot(residuals*W, Xhat) 
            
            deltabeta = lstsq(A,b)[0]
            betahat += deltabeta
    
        beta = betahat/numpy.square(betahat).sum()**0.5

    return beta, numpy.square(betahat).sum()**0.5, g

def mixed_regression(Y,X_lin,X_nonlin,W=None,varM=None,type='kernel_regression',num_total_iters=10,num_beta_iters=10):
    #This fits a function of the form Y = beta*X_lin + F(beta2*X_nonlin)
    print X_lin.shape, X_nonlin.shape
    X_augmented = numpy.hstack((X_lin[:,None],X_nonlin))
    print X_augmented.shape
    if W==None:
        W=numpy.ones_like(Y)*1.0
    else:
        W=W*1.0
    
    if varM == None:
        # set varM to be 5% of variance
        varM = 0.05*numpy.var(Y)
        
    # To get an initial estimate of the coefficients,
    # we first estimate beta when when F(u) = u+b
    beta, betanorm, dummy_F = ordinary_linear_regression(Y,X_augmented,W)
    beta = betanorm*beta
    
    # iterate over estimating F(u) and beta
    for current_total_iters in xrange(0,num_total_iters):
        #Note that the trait, U, which is the argument of the nonlinear function g, is now only a function of X_nonlin (loci after the first)
        U = numpy.dot(X_nonlin,beta[1:])
        
        Urange = (U.max()-U.min())
        U += normal(0,Urange*1e-04,size=len(U))
        #print (beta[0]*X_lin).shape, Y.shape
        residuals = Y - beta[0]*X_lin #Note that the intercept will end up encoded in g_nonlin
        #print U.shape, residuals.shape
        # estimate global trait function
        if type=='spline':
            sorted_U, sorted_Y = (list(x) for x in zip(*sorted(zip(U, residuals))))
    
            sorted_U = numpy.array(sorted_U)
            sorted_Y = numpy.array(sorted_Y)
    
            smoothing_factor = 0.43*len(sorted_Y)*(sorted_Y.std()**2)
            g_nonlin = UnivariateSpline(sorted_U,sorted_Y,s=smoothing_factor)
            gprime = g.derivative()
        elif type=='kernel_regression':
            g_nonlin = KernelRegressionSmoother(U,residuals,W,dh=0.1)
            # pt.plot(U,residuals,'o')
            # pt.hold('on')
            
            # xvec = numpy.arange(-.1,.1,.01)
            # pt.plot(xvec,g_nonlin(xvec),'k')
            # pt.show()
            gprime = lambda x: g_nonlin.derivative(x)   
        else:   
            g_nonlin = KernelSmoother(U,residuals,W) #The nonlinear function
            gprime = lambda x: g_nonlin.derivative(x)
    
        
        # estimate variance weight
        varA = 0
        #varA = 0*((numpy.square(residuals2)*numpy.square(gprimes)).sum()-varM*numpy.square(gprimes).sum())/(numpy.square(gprimes)*numpy.square(gprimes)).sum()
    
        # estimate coefficients 
        betahat = beta
        for current_beta_iters in xrange(0,num_beta_iters):
            Uhat_nonlin = numpy.dot(X_nonlin,betahat[1:]) #nonlinear argument
            residuals2 = Y-betahat[0]*X_lin-g_nonlin(Uhat_nonlin)
            gprimes = gprime(Uhat_nonlin)
            gprimesquareds = numpy.square(gprimes)
            variances = varM+varA*gprimesquareds
            X_lin_scaled = X_lin/gprimes#
            
            X_augmented = numpy.hstack((X_lin_scaled[:,None],X_nonlin))
			
            num_traits = X_augmented.shape[1]
            A = numpy.dot(X_augmented.T,numpy.dot(numpy.diag(W*gprimesquareds/variances),X_augmented))
            b = numpy.dot(residuals2*W*gprimes/variances, X_augmented)
            deltabeta = lstsq(A,b)[0]
            betahat += deltabeta
        betanorm = 1.
        beta = betahat
        #beta = betahat/numpy.square(betahat).sum()**0.5

    return beta, numpy.square(betahat).sum()**0.5, g_nonlin

def mixed_regression_sergey_model(Y,X_lin,X_nonlin,W=None,varM=None,type='kernel_regression',num_total_iters=10,num_beta_iters=10):
    #This fits a function of the form Y = beta*X_lin + F(beta2*X_nonlin)
	#Where F takes the 'Sergey' form: F(U) = 
    #print X_lin.shape, X_nonlin.shape
    X_augmented = numpy.hstack((X_lin[:,None],X_nonlin))
    #print X_augmented.shape
    if W==None:
        W=numpy.ones_like(Y)*1.0
    else:
        W=W*1.0
    
    if varM == None:
        # set varM to be 5% of variance
        varM = 0.05*numpy.var(Y)
        
    # To get an initial estimate of the coefficients,
    # we first estimate beta when when F(u) = u+b
    beta, betanorm, dummy_F = ordinary_linear_regression(Y,X_augmented,W)
    beta = betanorm*beta
    
    # iterate over estimating beta
	
    
	#Note that the trait, U, which is the argument of the nonlinear function g, is now only a function of X_nonlin (loci after the first)
    U = numpy.dot(X_nonlin,beta[1:])
	
    Urange = (U.max()-U.min())
    U += normal(0,Urange*1e-04,size=len(U))
    #print (beta[0]*X_lin).shape, Y.shape
    residuals = Y - beta[0]*X_lin #Note that the intercept will end up encoded in g_nonlin
    #print U.shape, residuals.shape
    # estimate global trait function
    X_c = .177
    g_nonlin = lambda x: X_c*(1.-numpy.exp(-1.*x/X_c))
    gprime = lambda x: numpy.exp(-1.*x/X_c)


    # estimate variance weight
    varA = 0
    #varA = 0*((numpy.square(residuals2)*numpy.square(gprimes)).sum()-varM*numpy.square(gprimes).sum())/(numpy.square(gprimes)*numpy.square(gprimes)).sum()

    # estimate coefficients 
    betahat = beta
    for current_beta_iters in xrange(0,num_beta_iters):
		Uhat_nonlin = numpy.dot(X_nonlin,betahat[1:]) #nonlinear argument
		residuals2 = Y-betahat[0]*X_lin-g_nonlin(Uhat_nonlin)
		gprimes = gprime(Uhat_nonlin)
		gprimesquareds = numpy.square(gprimes)
		variances = varM+varA*gprimesquareds
		X_lin_scaled = X_lin/gprimes#
		
		X_augmented = numpy.hstack((X_lin_scaled[:,None],X_nonlin))
		
		num_traits = X_augmented.shape[1]
		A = numpy.dot(X_augmented.T,numpy.dot(numpy.diag(W*gprimesquareds/variances),X_augmented))
		b = numpy.dot(residuals2*W*gprimes/variances, X_augmented)
		deltabeta = lstsq(A,b)[0]
		betahat += deltabeta
    betanorm = 1.
    beta = betahat
	#beta = betahat/numpy.square(betahat).sum()**0.5

    return beta, numpy.square(betahat).sum()**0.5, g_nonlin
    
if __name__=='__main__':
    x = numpy.zeros((1000,1))
    x[:,0] = numpy.linspace(-1,1,1000)
    y = normal(2*numpy.linspace(-1,1,1000),1)
    
    beta,F = generalized_linear_regression(y,x)
    print beta

    import pylab
    
    y_pred = F(numpy.dot(beta,x.T))
    pylab.plot(y_pred,y,'k.')
    pylab.show()
