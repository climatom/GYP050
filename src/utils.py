# Python code 
import numpy as np
import scipy

# Function to comptute the probabilities of values, x,
# Given shape and scale parameters
def weibull_pdf(x,shape,loc,scale):
    p=shape/scale*np.power((loc-x)/scale,(shape-1))*np.exp(-np.power((loc-x)/scale,shape))
    return p

def gumbel_pdf(x,loc,scale):
	x=(x-loc)/scale
	p=1/scale*np.exp(x)*np.exp(-np.exp(x))
	return p

def norm_pdf(xs):
	p=1/(np.sqrt(2*np.pi))*np.exp(-np.power(xs,2)/2)
	return p

def norm_cdf(xs):
	c=1/2.*(1+scipy.special.erf(xs/np.sqrt(2)))
	return c

def skew_norm_pdf(x,loc,scale,skew):
	xs=(x-loc)/scale
	p=2/scale*norm_pdf(xs)*norm_cdf(skew*xs)	
	return p

def ns_skew_norm_cdf(x,scale,skew,alpha,beta,cov):
	loc=alpha+beta*cov
	cdf=norm_cdf((x-loc)/scale) -2*scipy.special.owens_t((x-loc)/scale,skew)
	return cdf

def ns_skew_norm_pdf(x,scale,skew,alpha,beta,cov):
	loc=alpha+beta*cov
	xs=(x-loc)/scale
	p=2/scale*norm_pdf(xs)*norm_cdf(skew*xs)	
	return p

def gumbel_cdf(x,loc,scale):
	x=(x-loc)/scale
	cdf=1-np.exp(-np.exp(x))
	return cdf

def gumbel_nll(x0,x):
	loc=x0[0]
	scale=x0[1]
	nll=-np.log(np.prod(gumbel_pdf(x,loc,scale)))
	return nll

def norm_skew_nll(x0,x):
	loc=x0[0]
	scale=x0[1]
	skew=x0[2]
	nll=-np.log(np.prod(skew_norm_pdf(x,loc,scale,skew)))
	return nll

def ns_norm_skew_nll(x0,x,cov):
	scale=x0[0]
	skew=x0[1]
	alpha=x0[2]
	beta=x0[3]
	nll=-np.log(np.prod(ns_skew_norm_pdf(x,scale,skew,alpha,beta,cov)))
	return nll


def optimize_gumbel(x,loc,scale):
	out=scipy.optimize.minimize(fun=gumbel_nll,
								x0=[loc,scale],
								args=(x),
								bounds=((None,None),(0,None)),
								method='Nelder-Mead')
	return out


def optimize_skew_norm(x,loc,scale,skew):
	out=scipy.optimize.minimize(fun=norm_skew_nll,
								x0=[loc,scale,skew],
								args=(x),
								bounds=((None,None),(None,None)),
								method='Nelder-Mead')
	return out


def ns_optimize_skew_norm(x,scale,skew,alpha,beta,cov):
	out=scipy.optimize.minimize(fun=ns_norm_skew_nll,
								x0=[scale,skew,alpha,beta],
								args=(x,cov),
								method='Nelder-Mead')
	return out

