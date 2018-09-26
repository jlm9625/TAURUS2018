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
ages = np.zeros(len(ainx))
pdb.set_trace()
for i in range(len(ainx)):
    ages[i] = float(mfl[ainx[i]].split(' ')[-1].split('\n')[0])*1000.0
    startline = hinx[i]+2
    if i< len(ainx)-1:    
        endline   = hinx[i+1]-2
    else: endline = len(mfl)-1

    thischunk = readcol(modfile,skipline=startline,skipafter=endline)
    ##now format this into a structure similar to the PARSEC models
    



pdb.set_trace()