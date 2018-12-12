import numpy as np
import glob,os,sys,pdb,time,pickle
import pandas as pd
from readcol import readcol


'''
script to read in the Baraffe/BTSETTL models
for solar/gaia-bands
'''
modfile = 'isochrones/baraffe_2015_solar_gaia.txt'


##i guess go line-by-line??

mfl = np.array(open(modfile,'rb').readlines())

hinx = np.where(mfl == mfl[5])[0]
dinx = np.where(mfl == mfl[6])[0]
ainx = hinx-2
ages = ()
model_array = ()
for i in range(len(ainx)):
    ages += (float(mfl[ainx[i]].split(' ')[-1].split('\n')[0])*1000.0,)
    startline = hinx[i]+2
    if i< len(ainx)-1:    
        endline   = hinx[i+1]-2
    else: endline = len(mfl)-1

    thischunk = readcol(modfile,skipline=startline,skipafter=endline)
    ##now format this into a structure similar to the PARSEC models
    ##define a recarray with the right columns but it will have only the columns we care about.
    #headstring = 'Zini Age Mini  Mass   logL    logTe  logg  label   McoreTP C_O  period0 period1 pmode  Mloss  tau1m   X   Y   Xc  Xn  Xo  Cexcess  Z mbolmag  Gmag    G_BPmag  G_RPmag  B_Tmag  V_Tmag  Jmag    Hmag    Ksmag'
    thisarray = np.recarray(len(thischunk[:,0]),dtype=[('Mini',float),('Mass',float),('Zini',float),('logTe',float),('logg',float),('Gmag',float),('G_BPmag',float),('G_RPmag',float)])
    thisarray.Mini = thischunk[:,0]*1.0
    thisarray.Mass = thischunk[:,0]*1.0
    thisarray.Zini[:] = 0.0152
    thisarray.logTe = np.log10(thischunk[:,1])
    thisarray.logg  = thischunk[:,3]*1.0
    thisarray.Gmag  =  thischunk[:,-3]*1.0
    thisarray.G_RPmag  =  thischunk[:,-1]*1.0
    thisarray.G_BPmag  =  thischunk[:,-2]*1.0
    ##now pandas it into a dataframe
    model_array += (pd.DataFrame.from_records(thisarray),)
    print('Model age ' + str(float(mfl[ainx[i]].split(' ')[-1].split('\n')[0])*1000.0))
    
pickle.dump((ages,model_array),open(modfile.split('.')[0]+'_stacked.pkl','wb'))
print('Model saved to '+modfile.split('.')[0]+'_stacked.pkl')
print 'Done'