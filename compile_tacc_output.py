import numpy as np
import pickle 
import matplotlib.pyplot as plt
import glob,pdb,os,sys

##US chunk
datadir= 'outputs/US180919_20180919-16:13/'
datafile = 'datastorage/gaia_dr2_usco_inputformat_20180516pair_ageprobcleaned_20180919-12:58:56_inputformat.pkl'

#gaia_dr2_lcc_inputformat_20180511pair_ageprobcleaned_20180919-14:08:34_inputformat.pkl'


files = np.array(glob.glob(datadir+'*pkl'))
fnum = np.array([ff.split('/')[-1].split('_')[1] for ff in files]).astype(int)
order =np.argsort(fnum)
fnum  =fnum[order]
files =files[order]

G,BP,RP,sig_G,sig_BP,sig_RP,ID = pickle.load(open(datafile,'rb'))


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


qwe= np.where(age < 20)[0]
plt.plot(BP-RP,G,',k')
plt.plot(BP[qwe]-RP[qwe],G[qwe],'.b')
plt.ylim([15,05])
plt.show()







pdb.set_trace()
print 'Done'