import numpy as np
import pickle 
import matplotlib.pyplot as plt
import glob,pdb,os,sys
import pandas as pd
from readcol import readcol

##US chunk
# datadir= 'outputs/US180919_20180919-16:13/'
# datafile = 'datastorage/gaia_dr2_usco_inputformat_20180516pair_ageprobcleaned_20180919-12:58:56_inputformat.pkl'
# outname= 'US_ageprob.pkl'
# 
# datadir= 'outputs/UCL180919_20180919-22:32/'
# datafile= 'datastorage/gaia_dr2_ucl_inputformat_20180511pair_ageprobcleaned_20180919-13:51:29_inputformat.pkl'
# outname='UCL_ageprob.pkl'
# 
# datadir= 'outputs/LCC180919_20180919-22:49/'
# datafile='datastorage/gaia_dr2_lcc_inputformat_20180511pair_ageprobcleaned_20180919-14:08:34_inputformat.pkl'
# outname = 'LCC_ageprob.pkl'

##100 pc chunk
datadir = 'outputs/tp100bhac_20180928-15:59/'
datafile = 'datastorage/gaia_dr2_100pc_ageprobcleaned_20180928-14:17:59_arenouclean_inputformat.pkl'
outname = 'tp100bhac_20180928-15:59_arenouclean_ageprob.pkl'
pairgaia = 'datastorage/gaia_dr2_100pc_ageprobcleaned_20180928-14:17:59_arenouclean.csv'


#gaia_dr2_lcc_inputformat_20180511pair_ageprobcleaned_20180919-14:08:34_inputformat.pkl'


files = np.array(glob.glob(datadir+'*pkl'))
fnum = np.array([ff.split('/')[-1].split('_')[1] for ff in files]).astype(int)
order =np.argsort(fnum)
fnum  =fnum[order]
files =files[order]

G,BP,RP,sig_G,sig_BP,sig_RP,ID = pickle.load(open(datafile,'rb'))
aa = readcol('isochrones/bhac_20myr.txt',asRecArray=True )
#pdb.set_trace()
for i in range(len(files)):
    print(str(i))
    if i == 0: 
        sourceID,age,mass,sig_age,sig_mass,myrange=pickle.load(open(files[i],'rb'))
    else:
        sourceID1,age1,mass1,sig_age1,sig_mass1,myrange1=pickle.load(open(files[i],'rb'))
        sourceID= np.concatenate((sourceID,sourceID1))
        age = np.concatenate((age,age1))
        mass = np.concatenate((mass,mass1))
        sig_mass = np.concatenate((sig_mass,sig_mass1))
        sig_age = np.concatenate((sig_age,sig_age1))
        myrange  =np.concatenate((myrange,myrange1))


#pdb.set_trace()
if np.sum(sourceID-ID) != 0:
    pdb.set_trace()



print ('reading pandas data table of gaia cat ')
data = pd.read_csv(pairgaia)




qwe= np.where((age>0)&(age<40) & (G>4.5) )[0]
from gal_xyz import gal_xyz

x,y,z = gal_xyz(data.ra.values,data.dec.values,data.parallax.values,radec=True,plx=True)
asd = np.where((x[qwe] < -35) & (y[qwe] < -64))[0]
pdb.set_trace()
##visual check HRD at age 20
plt.plot(BP-RP,G,',k')
plt.plot(BP[qwe]-RP[qwe],G[qwe],'.b')
plt.plot(aa.G_BP-aa.G_RP,aa.G,'b')
plt.ylim([15,-2])
plt.show()




pickle.dump((ID,age,mass,sig_age,sig_mass),open(outname,'wb'))

#qwe= np.where(age < 20)[0]








#pdb.set_trace()
print 'Done'