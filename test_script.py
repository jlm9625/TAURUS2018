import taurus_functions as tau
import numpy as np
import matplotlib.pyplot  as plt
import matplotlib as mpl
import pdb,os,sys,pickle
import pandas as pd

##here's the data, read it in
datadir = 'datastorage/'
dfile = 'GDR2_TAUAUR_cleaned.csv'
data = pd.read_csv(datadir + dfile)

clean_dfile = tau.clean_gaia_data(data,'dummy.csv')
print(clean_dfile)
##can comment the above out once you have run them and use the following line
clean_dfile = 'datastorage/dummy_20180822-12:04:02.csv'

##now read the clean file
cdata = pd.read_csv(clean_dfile,sep=',') 

##now compute things like absolute mags and errors:
G,BP,RP,sig_G,sig_BP,sig_RP,sourceID = tau.produce_data_gaia(cdata)

##this should be a flipped HR diagram 
fig0,ax0 = plt.subplots()
ax0.plot(BP-RP,G,'k,')


##load in a stacked model for checking, its parsec Gaia + 2mass solar 
#mod_age,mod_stack = pickle.load(open('isochrones/parsec1.2s_z0.0152_GDR2_2MASS_stacked.pkl','rb'))##this should be fine I think

##generate samples without metallicity for a solar metallicity model stack
modfile = 'isochrones/parsec1.2s_z0.0152_GDR2_2MASS_stacked.pkl'
#sample_2d = tau.sample_generate_2d(nsamples=100000,regen=True,use_mini=True,outfilename = 'isochrones/tmp1_dummy.pkl',readfile='isochrones/tmp1_dummy.pkl',modelfile=modfile,agelims=[1.0,3000.0],masslims=[0.1,1.0])
##try loading it instead of regening
sample_2d = tau.sample_generate_2d(nsamples=100000,regen=False,outfilename = None,readfile='isochrones/fullsamp_mini.pkl',modelfile=None,agelims=[1.0,3000.0],masslims=[0.1,1.0])
 
 ##Above code all works!!

##A quick plot of the samples
fig,ax = plt.subplots()
binbin = np.where(sample_2d[:,5] > 0)[0]
ax.plot(sample_2d[:,3]-sample_2d[:,4],sample_2d[:,2],'k,')
ax.plot(sample_2d[binbin,3]-sample_2d[binbin,4],sample_2d[binbin,2],'r,')
plt.show() 


##now probabilities for one star, note the option to show the user histograms if switched on.
##Run test on HP Tau G3
ttt = np.where(sourceID == 145213192171160064)[0] ##HPTau G3 
prob_data_model,age,mass,sig_age,sig_mass = tau.probability_calculation_singlestar(sample_2d,G[ttt[0]],BP[ttt[0]],RP[ttt[0]],sig_G[ttt[0]],sig_BP[ttt[0]],sig_RP[ttt[0]],sourceID[ttt[0]],run_2mass=False,showhist=True)

##how about running things on a few stars ##HP Tau G3 and FF Tau
torun = np.where((sourceID == 145213192171160064) | (sourceID == 145213875069914496) )[0]
age2,mass2,sig_age2,sig_mass2 = tau.probability_calculation_all(sample_2d,G[torun],BP[torun],RP[torun],sig_G[torun],sig_BP[torun],sig_RP[torun],sourceID[torun],run_2mass=False)

##Thats all for now

pdb.set_trace()
print 'Im done'