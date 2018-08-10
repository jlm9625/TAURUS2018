import numpy as np
from   scipy.optimize import curve_fit
import pdb
import mpfit
import matplotlib.pyplot as plt
from   mpfitexpr import mpfitexpr

def fit_func(x, a0, a1, a2, a3, a4, a5):
	z = (x - a1)/a2
	y = a0 * np.exp(-z**2/2.0) + a3 + a4*x + a5*x**2
	return y

def fit_funcmp(p,x=None,y=None,err=None,fjac=None,model=False):
    z     = (x-p[1])/p[2]
    modl  = p[0]*np.exp(-z**2/2.0) + p[3] + p[4]*x + p[5]*x**2
    resid = (y-modl)/err
    if model == False:    return [0,resid]
    if model == True: return modl
    
    
def gaussfit(xdata,ydata,p0):
	pars, cov = curve_fit(fit_func,xdata,ydata,p0)
	thecurve  = fit_func(xdata,pars[0],pars[1],pars[2],pars[3],pars[4],pars[5])
	return pars , cov,thecurve

def simple(p,x=None,y=None,err=None,fjac=None,model=False):
    modl  =p[0] + p[1]*x + p[2]*x**2
    resid = (y-modl)/err
    if model == False:    return [0,resid]
    if model == True: return modl
  
  
def gaussian_2d(p,x=None,y=None,z=None,err=None,fjac=None,model=False):
    ##P is the parameters of the model
    ##p=[normalization, mean1,mean2, stdev1,stdev2,background_const,rotation angle]  the starting point for the fit
    ##x and y are the 2d grids of independent variable points 


    ##first rotate the x/ygrids to new x/y
    if p[5] > 0.0000001: ##I think this is right....
        x = x*np.cos(np.pi/180*p[6]) - y*np.sin(np.pi/180*p[6])
        y = x*np.sin(np.pi/180*p[6]) - y*np.cos(np.pi/180*p[6])
        
    themodel = np.exp(-(x-p[1])**2/2/p[3]**2 -(y-p[2])**2/2/p[4]**2)*p[0]+p[7]
    if model == True: return themodel
    if model == False: return [0,(z-themodel)/err]
            
        
def gaussfit_mp_gridless(zgrid,error,p0,norotation=False):
    ##construct the dependent variables
    x,y = np.meshgrid(np.arange(zgrid.shape[0],dtype=float),np.arange(zgrid[1],dtype=float))
    ##zgrid are the dependent variable measurements to fit to, error is the uncertainty on zgrid.
    ##p0=[normalization, mean1,mean2, stdev1,stdev2,background_const,rotation angle]  the starting point for the fit
    ##set norotation = True to fix the rotation to zero.
    fa      = {'x':xgrid, 'y':ygrid, 'z':zgrid, 'err':error}
    parinfo = [{'fixed':0, 'limited':[0,0], 'limits':[0.,0.]} for dddd in range(p0.shape[0])]
    
    ##fix rotation if asked to
    if norotation == True:
        parinfo[6]['fixed'] = 1
        p0[6]=0.0000
    
    thisfit = mpfit.mpfit(gaussian_2d,p0,functkw=fa,quiet=True,parinfo=parinfo)
    fitpars = thisfit.params
    fitcov  = thisfit.covar
    
    themodel = gaussian_2d(fitpars,x=xgrid,y=ygrid,z=zgrid,err=error,model=True)
        
    return fitpars,fitcov,themodel       
  
def gaussfit2d_mp(xgrid,ygrid,zgrid,error,p0,norotation=False):
    ##xgrid/ygrid are the 2D grid of independent variables on which co calculate the model
    ##zgrid are the dependent variable measurements to fit to, error is the uncertainty on zgrid.
    ##p0=[normalization, mean1,mean2, stdev1,stdev2,background_const,rotation angle]  the starting point for the fit
    ##set norotation = True to fix the rotation to zero.
    fa      = {'x':xgrid, 'y':ygrid, 'z':zgrid, 'err':error}
    parinfo = [{'fixed':0, 'limited':[0,0], 'limits':[0.,0.]} for dddd in range(p0.shape[0])]
    
    if norotation == True:
        parinfo[6]['fixed'] = 1
        p0[6]=0.0000
    
    thisfit = mpfit.mpfit(gaussian_2d,p0,functkw=fa,quiet=True,parinfo=parinfo)
    fitpars = thisfit.params
    fitcov  = thisfit.covar
    
    themodel = gaussian_2d(fitpars,x=xgrid,y=ygrid,z=zgrid,err=error,model=True)
        
    return fitpars,fitcov,themodel
    
    
    
  
def gaussfit_mp(xdata,ydata,error,p0,nocurve=False):
    fa      = {'x':xdata, 'y':ydata, 'err':error}
    parinfo = [{'fixed':0, 'limited':[0,0], 'limits':[0.,0.]} for dddd in range(p0.shape[0])]
    parinfo[2]['limited'] = [1,0]
    parinfo[2]['limits']  = [0.0,100.]
    
    if nocurve==True:
        parinfo[4]['fixed'] = 1
        parinfo[5]['fixed'] = 1
        p0[4] = 0.0
        p0[5] = 0.0
    
    thisfit = mpfit.mpfit(fit_funcmp,p0,functkw=fa,quiet=True,parinfo=parinfo)
    fitpars = thisfit.params
    fitcov  = thisfit.covar
    
    thecurve = fit_funcmp(fitpars,x=xdata,y=ydata,err=error,model=True)
    #pdb.set_trace()

    return fitpars,fitcov,thecurve