import numpy as np
import matplotlib.pyplot as plt
import os,sys,pdb,pickle,glob
import pandas as pd
from gcirc import gcirc
from astroquery.vizier import Vizier
import astropy.units as u
import astropy.coordinates as coord

Vizier.columns=['all']

##cross match an input file with 2mass using queryvizier. This is going to be slow I assume.

##The basline Taurus datafile as a first example:
datafile = 'datastorage/gaia_dr2_tauaur.csv'

##Process:
    ##See if anything within 5 arcsec in Gaia, 
    ##if significant source within 2 arcsec, see if 2mass id's this
    ##if not, and if source not super bright, correct for this with JHK mag error term of some sort.
    ##if nothing in Gaia within 5 arcsec all good.
    

##read the data
data = pd.read_csv(datafile,sep=',')

##all the stuff we want to get from Vizier
J     = np.zeros(len(data))
H     = np.zeros(len(data))
K     = np.zeros(len(data))
sig_J = np.zeros(len(data))+1000.0
sig_H = np.zeros(len(data))+1000.0
sig_K = np.zeros(len(data))+1000.0
qj    = np.zeros(len(data),dtype=str)
qh    = np.zeros(len(data),dtype=str)
qk    = np.zeros(len(data),dtype=str)
bj    = np.zeros(len(data),dtype=str)
bh    = np.zeros(len(data),dtype=str)
bk    = np.zeros(len(data),dtype=str)
flag = np.zeros(len(data),dtype=int)

start,J,H,K,sig_J,sig_H,sig_K,qj,qh,qk,bj,bh,bk,flag = pickle.load(open(datafile.split('.')[0]+'_X2MASS_progress.pkl','rb'))

pdb.set_trace()
for i in range(start,len(data)):
    queryname = 'Gaia DR2 ' + str(data.source_id.values[i]) ##The official object name out of Gaia    
    try:
        #stuff0 =  Vizier.query_region(queryname,radius="0d0m5s",catalog='I/345/gaia2')
        stuff1 =  Vizier.query_region(queryname,radius="0d0m3s",catalog='II/246/out')
    except: ##try at least twice before quitting
        print('Excepted on first query??!!??')
        #stuff0 =  Vizier.query_region(queryname,radius="0d0m5s",catalog='I/345/gaia2')
        stuff1 =  Vizier.query_region(queryname,radius="0d0m3s",catalog='II/246/out')
    #if len(stuff0) == 0: 
    #    print('Cant find Gaia DR2 entry, thats bad...')
    #    pdb.set_trace()
        
    if (len(stuff1) > 1): pdb.set_trace()
    if (len(stuff1) == 1): ##if something in 2mass and assuming Gaia source is found
        cl = np.argmin(stuff1[0]['_r'])
        qj[i] = stuff1[0]['Qflg'][cl][0]
        qh[i] = stuff1[0]['Qflg'][cl][1]
        qk[i] = stuff1[0]['Qflg'][cl][2]
        bj[i] = stuff1[0]['Bflg'][cl][0]
        bh[i] = stuff1[0]['Bflg'][cl][1]
        bk[i] = stuff1[0]['Bflg'][cl][2]
        J[i] = stuff1[0]['Jmag'][cl]*1.0
        H[i] = stuff1[0]['Hmag'][cl]*1.0  
        K[i] = stuff1[0]['Kmag'][cl]*1.0  

        ##only write new errors if the flags are reasonable
        if (qj[i] != 'U') & (qj[i] != 'X'): sig_J[i] = stuff1[0]['e_Jmag'][cl]*1.0
        if (qh[i] != 'U') & (qh[i] != 'X'): sig_H[i] = stuff1[0]['e_Hmag'][cl]*1.0
        if (qk[i] != 'U') & (qk[i] != 'X'): sig_K[i] = stuff1[0]['e_Kmag'][cl]*1.0       
        
    else: flag[i] = -1 ##flag for no source in 2mass points source cat.
    
    print('Up to ' + str(i+1) + ' out of ' + str(len(data)))
    if np.mod(i,1000) == 0: pickle.dump((i,J,H,K,sig_J,sig_H,sig_K,qj,qh,qk,bj,bh,bk,flag),open(datafile.split('.')[0]+'_X2MASS_progress.pkl','wb'))
        


##output to the same place as the data was read from, with slightly different filename
pickle.dump((J,H,K,sig_J,sig_H,sig_K,qj,qh,qk,bj,bh,bk,flag),open(datafile.split('.')[0]+'_X2MASS.pkl','wb'))

#pdb.set_trace()
print 'Done'
