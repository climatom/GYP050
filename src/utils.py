# Python code 
import numpy as np
import scipy
import numba as nb
from numba import prange

# Assorted stats 
def gev_pdf(x,shape,loc,scale):
    # Support is: -inf < x < (μ − σ / ξ)
    x=np.atleast_1d(x)
    p=np.ones(len(x))*1e-15
    t=np.ones(len(x))*1e-15
    idx=x<(loc-scale/shape)
    t[idx]=np.power((1+shape*(x[idx]-loc)/scale),-1/shape)
    p[idx]=1/scale*np.power(t[idx],shape+1)*np.exp(-t[idx])
    return p

def gev_cdf(x,shape,loc,scale):
    x=np.atleast_1d(x)
    c=np.ones(len(x))
    t=np.ones(len(x))
    idx=x<(loc-scale/shape)
    t[idx]=np.power((1+shape*(x[idx]-loc)/scale),-1/shape)
    c[idx]=np.exp(-t[idx])
    return c


def gev_nll(x0,x):
        shape=x0[0]
        loc=x0[1]
        scale=x0[2]
        p=gev_pdf(x,shape,loc,scale)
        nll=-np.log(np.prod(p[~np.isnan(p)]))
        return nll


def optimize_gev(x,shape,loc,scale):
       out=scipy.optimize.minimize(fun=gev_nll,
                  x0=[shape,loc,scale],
                  args=(x),
                  #bounds=((-0.4,0.4),(None,None),(None,None)),
                  method='Nelder-Mead')
       return out

def ns_gev_pdf(x,shape,scale,alpha,beta,cov):
    # Support is: -inf < x < (μ − σ / ξ)
    # note that cov (and therefore loc) is/are  scalars
    loc=alpha+beta*cov
    x=np.atleast_1d(x)
    #print(cov,x)
    p=np.ones(len(x))*1e-50
    t=np.ones(len(x))*1e-50
    idx=x<(loc-scale/shape)
    t[idx]=np.power((1+shape*((x[idx]-loc)/scale)),-1/shape)
    p[idx]=1/scale*np.power(t[idx],shape+1)*np.exp(-t[idx])
    return p

def ns_gev_cdf(x,shape,scale,alpha,beta,cov):
    loc=alpha+beta*cov
    x=np.atleast_1d(x)
    c=np.zeros(len(x))
    t=np.zeros(len(x)) 
    idx=x<(loc-scale/shape)
    t[idx]=np.power((1+shape*((x[idx]-loc)/scale)),-1/shape)
    c=np.exp(-t)
    c[~idx]=1.
    return c


def ns_gev_nll(x0,x,cov):
        shape=x0[0]
        scale=x0[1]
        alpha=x0[2]
        beta=x0[3]
        loc=alpha+beta*cov
        x=np.atleast_1d(x)
        loc=np.atleast_1d(loc)
        p=np.ones(len(x))*1e-50
        t=np.ones(len(x))*1e-50
        idx=x<(loc-scale/shape)
        t[idx]=np.power((1+shape*((x[idx]-loc[idx])/scale)),-1/shape)
        p[idx]=1/scale*np.power(t[idx],shape+1)*np.exp(-t[idx])
        nll=-np.log(np.max([np.prod(p),1e-50]))
        return nll


def ns_optimize_gev(x,shape,scale,alpha,beta,cov):
       out=scipy.optimize.minimize(fun=ns_gev_nll,
                  x0=[shape,scale,alpha,beta],
                  args=(x,cov),
                  #bounds=((-0.4,0.4),(None,None),(None,None)),
                  method='Nelder-Mead')
       return out



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





@nb.njit(fastmath=True)
def _satVP(t):
    """Saturation specific humidity from temp (k)"""
    # t=np.atleast_1d(t)
    # esat=np.zeros(t.shape)
    t0=273.16
    a1_w=611.21; a3_w=17.502; a4_w=32.19
    a1_i=611.21; a3_i=22.587; a4_i=-0.7
    
    # Sat Vp according to Teten's formula
    if t>273.15:
        esat=a1_w*np.exp(a3_w*(t-t0)/(t-a4_w))
    else:
        esat=a1_i*np.exp(a3_i*(t-t0)/(t-a4_i))   
        
    return esat

@nb.njit(fastmath=True,parallel=True)
def gen_future_conditions(dsst,dt,drh,sst,t,q,p,nt,nlev,nr,nc):
    sst_out=sst+dsst
    t_out=t+dt
    rh_out=np.zeros((nt,nlev,nr,nc))*np.nan
    q_out=np.zeros((nt,nlev,nr,nc))*np.nan
    for _t in prange(nt):
       for _l in prange(nlev):
         for _r in prange(nr):
           for _c in prange(nc): 

             svp=_satVP(t[_t,_l,_r,_c])/100. # hPa
             sq=0.622*svp/p[_l]
             rh=q[_t,_l,_r,_c]/sq;
             rh_out[_t,_l,_r,_c]=rh+drh  

             if rh_out[_t,_l,_r,_c] >1: rh_out[_t,_l,_r,_c]=1.
             if rh_out[_t,_l,_r,_c] <0: rh_out[_t,_l,_r,_c]=0.01

             t_out[_t,_l,_r,_c]=t[_t,_l,_r,_c]+dt      
             svp=_satVP(t_out[_t,_l,_r,_c])/100. # hPa
             sq=0.622*svp/p[_l]

             q_out[_t,_l,_r,_c]=sq*rh_out[_t,_l,_r,_c]

    # Returns all ts in K
    return sst_out,t_out,rh_out,q_out

    

