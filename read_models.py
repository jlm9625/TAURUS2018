import numpy as np
import os,sys,glob,pdb,pickle
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
sys.path.append('pyutils')
#from readcol import readcol

def readfile(filename):
	output = pd.read_csv(filename, sep=',')
	return output

from numpy.lib.recfunctions import stack_arrays


headstring = 'Zini Age Mini  Mass   logL    logTe  logg  label   McoreTP C_O  period0 period1 pmode  Mloss  tau1m   X   Y   Xc  Xn  Xo  Cexcess  Z mbolmag  Gmag    G_BPmag  G_RPmag  B_Tmag  V_Tmag  Jmag    Hmag    Ksmag'
headarr = headstring.split(' ')
clean = np.where(np.array(headarr) != '')[0]
headarr = np.array(headarr)[clean]
#pdb.set_trace()
files= glob.glob('isochrones/parsec1.2*txt')
print(files)
for i in range(len(files)):
    modfile=open(files[i],'r')
    linenumber=0
    done = 1
    ages =()
    model_array = ()
    while True:
        while True:
            thisline=modfile.readline()	
            linenumber += 1
            if len(thisline) == 0: 
                #pdb.set_trace()
                break
                print ('here')

            #if done == -1: pdb.set_trace()
            if thisline[0] != '#':
                startline=linenumber*1
                break
        #pdb.set_trace()

        while True: 
            thisline = modfile.readline()
            linenumber +=1
            if len(thisline) == 0: 
                numlines=1
                break
            #pdb.set_trace()
            if thisline[0] == '#':
                numlines = linenumber-startline
                break
        
        if numlines < 10: break
     #  
    #  pdb.set_trace()
        thisiso = pd.read_table(files[i],skiprows=startline-1,nrows=numlines,skipinitialspace=True,sep='\s+',names=headarr.astype(list),header=None)
        ages += (thisiso.Age[0]/1.0e6,)
        model_array += (thisiso.copy(),)
        print(linenumber,numlines)
       # pdb.set_trace()

    outname = files[i][0:len(files[i])-4] + '_stacked.pkl'
    pickle.dump((ages,model_array),open(outname,'wb'))
	
# pdb.set_trace()
# for i in range(len(linenums)):
# 	print('up to ' + str(i+1) + ' out of ' + str(len(linenums)))
# 	finishline = linenums[i]
# 	if i > 0:
# 		startline = linenums[i-1]
# 	else: startline =0
# 	
# 	
# 	thisiso = pd.read_csv('parsec_isochrones.txt',skiprows=startline-1,nrows=finishline-startline-2,skipinitialspace=True,delimiter=' ').to_records()
# 	ages[i] = thisiso.Age[0]/1.0e6
# 	#if i == 0: 
# 	model_array += (thisiso.copy(),)
# 	#pdb.set_trace()
# 	#else:
# 		#model_array = stack_arrays((model_array,thisiso),asrecarray=True,usemask=False)
# 	#	pdb.set_trace()
# 	#	model_array = np.array([model_array,thisiso])
# #		model_array = np.concatenate((model_array,np.array([thisiso])),axis=1)
# 
# 
# #example		stack_arrays((a, b), asrecarray=True, usemask=False)
# 	
# 	##remove the age columns, stack with previous array.
# 	
# 	
# 	#pdb.set_trace()
# #model_array = np.array([model_array])
# pickle.dump((ages,model_array),open('stacked_parsec_models.pkl','wb'))
# 
# ##reading it back int
# model_ages,model_array=pickle.load(open('stacked_parsec_models.pkl','rb'))
# 
# 
# ##1. Generate a random population of stars with ages/mass
# #used numpy.random.uniform(
# random_age = np.random.uniform(1.,3000,size=1000)
# random_mass = np.random.uniform(0.1,1.0,size=1000)
# ##initalize a storage array: model_stars = np.zeros((npoints,5),dtype=float)
# ##bound age on 1 million to 3 billion in age
# ##and 0.1 to 1 solar masses in mass
# #loop here
# ##find the nearest two ages 
# ##low_indx = np.where(model_ages <= random_age[i])[0][-1]
# ##high_indx = np.where(model_ages > random_age[i])[0][0]
# 
# ##then interpolated based on the random_mass[i] for both low/high models (numpy.interp)
# ##agelow_gmag = np.interp(random_mass[i],model_array[low_indx].Mass,model_array[low_indx].Gmag)
# ##agehigh_gmag = np.interp(random_mass[i],model_array[high_indx].Mass,model_array[high_indx].Gmag)
# ##model_Gmag = np.interp(random_age[i],[model_ages[low_indx],model_ages[high_indx]],[agelow_gmag,agehigh_gmag]) 
# 
# ##then interpolated the two outputs of that to the right age.
# 
# 
# 
# 
# #store things uniformly
# #model_stars[i,0] = random_mass[i]
# #model_start[i,1] = random_age[i]
# #model_stars[i,3] = model_Gmag*1.0
# #same for BP, RP
# 
# 
# 
# 
# 
# pdb.set_trace()
# 
# dfile = 'GDR2_TAUAUR_cleaned.csv'
# data = readfile(dfile)
# 
# ##basics:
# ##distance = 1000.0/parallax(in mas)
# 
# ##1. calculate errors on phot_bp_mean_mag and phot_rp_mean_mag, and phot_g_mean_mag
# ##check out phot_x_mean_flux and phot_x_mean_flux_error. and google how to calculate mag from flux.
# ##propagation of uncertainty on wiki and look and the log10 example and apparent magnitude on wiki for the equation.
# gmag_efrac = data.phot_g_mean_flux.values/data.phot_g_mean_flux_error.values
# bmag_efrac = data.phot_bp_mean_flux.values/data.phot_bp_mean_flux_error.values
# 
# plt.hist(bmag_efrac,bins=300)
# plt.show()
# 
# 
# ##absolute mag from apparent and parallax (google distance modulus)
# abs_g = data.phot_g_mean_mag.values+5-5*np.log10(1000.0/data.parallax.values)
# pdb.set_trace()
# ##2. Error on absolute mag from error on G (apparent) and error on parallax.
# 
# 
# 
# ##3. clean out bad points (or even low quality points)
# ##read the gaia paper, and look at things like parallax/parallax_error, and 
# ##phot_x_mean_flux/phot_x_mean_flux_error
# 
# 
# 
# ##simple HRD
# plt.plot(data.phot_bp_mean_mag.values-data.phot_rp_mean_mag.values,abs_g,',b')
# plt.ylim([15,-2])
# plt.savefig('coolplot.pdf')
# plt.show()
# pdb.set_trace()


print('Finished')
