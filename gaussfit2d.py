import numpy as np
import pdb

##intented for fitting PSF images, sop amplitude is decoupled from the sigma, and only one sigma is used
def gauss2d_resid(p,xgrid=None,ygrid=None,im=None,sig_im=None,fjac=None,model=False):
    ##make the model
    ##params are p=[mux,muy,sig,ampl]
    themodel = np.exp(-(xgrid-p[0])**2/2.0/p[2] -(ygrid-p[1])**2/2.0/p[2])*p[3]
    resid    = (themodel-im)/sig_im
    if model == True: 
        return resid,themodel
    else:
        return [0,resid.flatten()]
        
        
def gaussfit2d(im,sig_im,pstart):
    xgrid,ygrid = np.meshgrid(np.arange(im.shape[0]),np.arange(im.shape[1]))
   # pdb.set_trace()
    import mpfit
    fa  = {'xgrid':xgrid, 'ygrid':ygrid,'im':im,'sig_im':sig_im}
    

    thisfit = mpfit.mpfit(gauss2d_resid,pstart,functkw=fa,quiet=True)
    
    bestresid,bestmod = gauss2d_resid(thisfit.params,xgrid=xgrid,ygrid=ygrid,im=im,sig_im=sig_im,model=True) 
    return thisfit,bestresid,bestmod
