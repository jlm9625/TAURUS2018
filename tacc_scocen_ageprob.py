import numpy as np
import os,glob,sys,pdb,time,pickle
import taurus_functions as tau
import pandas as pd

##MPI STUFF
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE

datadir     = 'datastorage/'
isodir      = 'isochrones/'
#set up the communicator basics
comm = MPI.COMM_WORLD
rank = comm.Get_rank() ##each core's ID number
size = comm.Get_size() ##the number of cores running
ncor = size
jobID = str(sys.argv[1]) 
subg  = str(sys.argv[2])

##Make a new folder for the results. If it doesn't exists already, cancel everything rather than overwrite stuff.
##output with starting timestamp so all runs are distinguishable
datestamp    = time.strftime('%Y%m%d-%H:%M', time.localtime(time.time())) 
outdir = 'outputs/' + jobID + '_'+datestamp+'/'
if rank == 0:
    if os.path.exists(outdir): 
        print 'Failed Due to Pre-existing output!!' ##SHOULD NEVER HAPPEN!!
        comm.Abort()
    else: os.makedirs(outdir)


input_file='fail' 
if subg == 'US':input_file = 'gaia_dr2_usco_inputformat_20180516pair_ageprobcleaned_20180919-12:58:56_inputformat.pkl'
if subg == 'UCL':input_file = 'gaia_dr2_ucl_inputformat_20180511pair_ageprobcleaned_20180919-13:51:29_inputformat.pkl'
if subg == 'LCC':input_file = 'gaia_dr2_lcc_inputformat_20180511pair_ageprobcleaned_20180919-14:08:34_inputformat.pkl'
if input_file == 'fail':
    print('thats not a valid scocen subgroup!?! ABORTING')
    comm.Abort()

##Here is a script that will run on a supercoputer. 
##You can also run it on your own computer (which has multiple cores) 
##from the terminal like this:
#mpiexec -n 2 python tacc_run_example.py JOBNAME
##the "2" means 2 cores, and JOBNAME can be anything you like without spaces in it.
##currently the script just runs 2 stars 


##here's the data, read it in. Assumes all data prep was done before on a local machine.

G,BP,RP,sig_G,sig_BP,sig_RP,sourceID = pickle.load(open(datadir+input_file,'rb'))

##for testing on specific stars
# torun = np.where((sourceID == 6057558691062793856) | (sourceID== 6237142264484167296) )[0]
# G = G[torun]
# RP=RP[torun]
# BP=BP[torun]
# sig_G = sig_G[torun]
# sig_RP=sig_RP[torun]
# sig_BP=sig_BP[torun]
# sourceID=sourceID[torun]
# pdb.set_trace()




##now each core has to figure out which data entries it is responsible for:
nstars = len(sourceID) ##the number of input stars to run
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
G = G[myrange[0]:myrange[1]]
RP=RP[myrange[0]:myrange[1]]
BP=BP[myrange[0]:myrange[1]]
sig_G = sig_G[myrange[0]:myrange[1]]
sig_RP=sig_RP[myrange[0]:myrange[1]]
sig_BP=sig_BP[myrange[0]:myrange[1]]
sourceID=sourceID[myrange[0]:myrange[1]]


##read in the model file Again, assuming it's pre-generated for us
sample_2d = tau.sample_generate_2d(regen=False,outfilename=None,readfile='isochrones/Solar_fullrange_2M.pkl')

age,mass,sig_age,sig_mass = tau.probability_calculation_all(sample_2d,G,BP,RP,sig_G,sig_BP,sig_RP,sourceID,run_2mass=False)


##now output everthing correctly
##make a unique filename and put it in the unique folder we already created/checked.
#Just store the source ID's and the results to save space and time.
#also store the range of data points used
pickle.dump((sourceID,age,mass,sig_age,sig_mass,myrange),open(outdir + 'result_'+str(rank)+'_'+datestamp+'.pkl','wb'))

print 'Done'



