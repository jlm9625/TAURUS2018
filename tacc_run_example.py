import numpy as np
import os,glob,sys,pdb,time,pickle
import taurus_functions as tau
import pandas as pd

##MPI STUFF
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE

#set up the communicator basics
comm = MPI.COMM_WORLD
rank = comm.Get_rank() ##each core's ID number
size = comm.Get_size() ##the number of cores running
ncor = size
jobID = str(sys.argv[1]) 

##Here is a script that will run on a supercoputer. 
##You can also run it on your own computer (which has multiple cores) 
##from the terminal like this:
#mpiexec -n 2 python tacc_run_example.py JOBNAME
##the "2" means 2 cores, and JOBNAME can be anything you like without spaces in it.
##currently the script just runs 2 stars 


##here's the data, read it in. Assumes all data prep was done before on a local machine.
datadir     = 'datastorage/'
clean_dfile = 'datastorage/GDR2_TAUAUR_cleaned.csv'
cdata       = pd.read_csv(clean_dfile,sep=',') 

##for testing, just run two stars that are Taurus members:
torun = np.where((cdata.source_id.values == 145213192171160064) | (cdata.source_id.values== 145213875069914496) )[0]
cdata =cdata.ix[torun]

#pdb.set_trace()

##Make a new folder for the results. If it doesn't exists already, cancel everything rather than overwrite stuff.
##output with starting timestamp so all runs are distinguishable
datestamp    = time.strftime('%Y%m%d-%H:%M:%S', time.localtime(time.time())) 
outdir = 'outputs/' + jobID + '_'+datestamp+'/'
if rank == 0:
    if os.path.exists(outdir): 
        print 'Failed Due to Pre-existing output!!' ##SHOULD NEVER HAPPEN!!
        comm.Abort()
    else: os.makedirs(outdir)


##now each core has to figure out which data entries it is responsible for:
nstars = len(cdata) ##the number of input stars to run
s   = nstars/ncor
nx  = (s+1)*ncor - nstars
ny  = nstars - s*ncor
if rank <= nx-1: 
    mynumber = s
    myrange = [(rank)*s,(rank+1)*s]
if rank > nx-1: 
    mynumber = s+1
    myrange = [nx*s+(rank-nx)*mynumber,nx*s+(rank-nx+1)*mynumber]

##now it's clear how many test are running on each processor, tell the user
if rank == 0: print 'Using ' + str(nx) + ' cores of size ' + str(s) + ' and ' +str(ny) + ' cores of size ' + str(s+1)

##cut out each core's data:
cdata = cdata.iloc[myrange[0]:myrange[1]]

##now compute abs mags and so on. Each core just works on the data it's responsible for
G,BP,RP,sig_G,sig_BP,sig_RP,sourceID = tau.produce_data_gaia(cdata)

##read in the model file Again, assuming it's pre-generated for us
sample_2d = tau.sample_generate_2d(regen=False,outfilename = None,readfile='isochrones/fullsamp_mini.pkl',modelfile=None,agelims=[1.0,3000.0],masslims=[0.1,1.0])

age,mass,sig_age,sig_mass = tau.probability_calculation_all(sample_2d,G,BP,RP,sig_G,sig_BP,sig_RP,sourceID,run_2mass=False)


##now output everthing correctly
##make a unique filename and put it in the unique folder we already created/checked.
#Just store the source ID's and the results to save space and time.
pickle.dump((sourceID,age,mass,sig_age,sig_mass),open(outdir + 'result_'+str(rank)+'_'+datestamp+'.pkl','wb'))

pdb.set_trace()

print 'Done'



