import taurus_functions as tau
import numpy as np
import matplotlib.pyplot  as plt
import matplotlib as mpl
import pdb,os,sys,pickle
import pandas as pd
##A script to make input data files for running the pipeline on the supercomputer. So that we don't waste 
##supercomputer time on boring stuff like organizing the data or computing magnitudes.



rundata = True
if rundata==True:
    ##Sco-Cen chunks:
    data_dir = '../BAFGKM/gaiadr2_scripts/'
    ##US
    #dfile = 'gaia_dr2_usco_inputformat_20180516pair.pkl'
    ##UCL
    #dfile = 'gaia_dr2_ucl_inputformat_20180511pair.pkl'
    ##LCC
    #dfile = 'gaia_dr2_lcc_inputformat_20180511pair.pkl'
#    data0 = pd.read_pickle(data_dir + dfile)
#    data  = pd.DataFrame.from_records(data0)

    
    dfile = 'gaia_dr2_100pc.csv'
    data_dir = '/Volumes/UT1/'
    data = pd.read_csv(data_dir + dfile)

    #pdb.set_trace()
    clean_dfile = tau.clean_gaia_data(data,dfile.split('.')[0] +'_ageprobcleaned.pkl',noclean=False,arenou_clean=True)
    print(clean_dfile)
    #pdb.set_trace()

    cdata = pd.read_csv(clean_dfile,sep=',') 
    dataname = clean_dfile.split('.')[0] + '_inputformat.pkl'
    ##now remove everything but the essentials:
    G,BP,RP,sig_G,sig_BP,sig_RP,sourceID = tau.produce_data_gaia(cdata,writetofile=dataname)
    pdb.set_trace()
###generate a model sample??
generate_model = True
if generate_model == True:
    modfile,outfilename   = 'isochrones/parsec1.2s_z0.0152_GDR2_2MASS_stacked.pkl','Solar_fullrange_2M.pkl' ##SOLAR METALLICITY PARSEC
    modfile,outfilename   = 'isochrones/baraffe_2015_solar_gaia_stacked.pkl','BHAC15_solar_fullsample.pkl'##SOLAR METALLICITY BHAC15
    sample_2d = tau.sample_generate_2d(nsamples=2000000,regen=True,use_mini=True,outfilename = outfilename,modelfile=modfile,agelims=[1.0,3000.0],masslims=[0.1,1.4])
#sample_2d = tau.sample_generate_2d(nsamples=100000,regen=False,outfilename = None,readfile='isochrones/Solar_fullrange_2M.pkl',modelfile=None,agelims=[1.0,3000.0],masslims=[0.1,1.0])

pdb.set_trace()
print 'Done'