##Standard stuff
import numpy as np
import os,sys,glob,pdb,pickle,math,time
import pandas as pd

##Plotting
import matplotlib.pyplot as plt
import matplotlib as mpl


'''
All the completed code goes here as functions for import by other scripts and notebooks.

ACR!! to load all this into a jupyter notebook, 


'''


def read_pickle_file(file):
    '''
    basic function wrapping pandas read pickle    
    again not really nessecary
    '''
    pickle_data = pd.read_pickle(file)
    return pickle_data
    
def readfile(filename):
    '''
    function wrapping pandas read_csv, good for gaia dr2 chunks
    these are not super nessecary
    '''
    output = pd.read_csv(filename, sep=',')
    return output
    
def clean_gaia_data(data,outfilename,outdir='datastorage/',noclean=False,arenou_clean=False):
    '''
    input a pandas table for Gaia data and clean it with some limits
    also input a filename string for the output.
    In future this needs to take in an additional table of 2MASS J/H/K magnitudes.
    '''
    if (noclean == False) & (arenou_clean==False):
        badrow = np.where((np.isnan(data.parallax.values) == True) & (np.isnan(data.phot_g_mean_mag.values) == True) & (np.isnan(data.phot_rp_mean_mag.values) == True) & (np.isnan(data.phot_bp_mean_mag.values) == True))[0]
        if len(badrow) == len(data):
            print('Everything is NaN!?!?!?!?')
            print('Not outputting anything')
            return -1
        print('Removing ' + str(len(badrow)) + ' entries with NaNs in key data')

        data = data.drop(data.index[badrow])

        ##Jordyn this is where this magnitude cuts can be placed. Note that
        ##these are searching in the negative sense.
        badrow = np.where((data.phot_bp_mean_flux_over_error.values) < 5)[0] 
        if len(badrow) == len(data):
            print('Your limit cuts ate all the data! Not outputting anything')
            return -1
        print('Removing ' + str(len(badrow)) + ' entries with bad photometry/parallax')
        data = data.drop(data.index[badrow])
    
    if (arenou_clean==True):
        ##first kill nans
        badrow = np.where((np.isnan(data.parallax.values) == True) & (np.isnan(data.phot_g_mean_mag.values) == True) & (np.isnan(data.phot_rp_mean_mag.values) == True) & (np.isnan(data.phot_bp_mean_mag.values) == True))[0]
        if len(badrow) == len(data):
            print('Everything is NaN!?!?!?!?')
            print('Not outputting anything')
            return -1
        print('Removing ' + str(len(badrow)) + ' entries with NaNs in key data')

        data = data.drop(data.index[badrow])
        
        ##now apply cleaning used by gaia team in lingegren 2018 and arenou 2018 to select "astrometrically clean subsets"
        u_param = np.sqrt(data.astrometric_chi2_al.values/(data.astrometric_n_good_obs_al.values - 5))
        b_r = data.phot_bp_mean_mag.values-data.phot_rp_mean_mag.values
        ef = data.phot_bp_rp_excess_factor.values
        g = data.phot_g_mean_mag.values
        dmod = +5-5*np.log10(1000.0/data.parallax.values)
        gexp = np.exp(-0.2*(g-19.5))
        qwe = np.where(gexp < 1.0)[0]
        gexp[qwe] = 1.0
        badrow = np.where((u_param >= 1.2*gexp) | (ef <=1.+0.015*b_r**2) | (ef >= 1.3+0.06*b_r**2))[0]
        goodrow = np.where((u_param < 1.2*gexp) & (ef >1.+0.015*b_r**2) & (ef < 1.3+0.06*b_r**2))[0]
        gr2 = np.where(u_param < 1.2*gexp)[0]
       # pdb.set_trace()
        data = data.drop(data.index[badrow])
        print('Removing ' + str(len(badrow)) +' that fail the lindegren/arenou cuts')
        
        
    if (noclean == True) & (arenou_clean==False):
        print('not cleaning anything from the input data')
        ##just replace nans with 0 and a huge error
        #actually, nans should just propagate, producing nans in the output right? THat's probably fine for now
        ##badrow = np.where((np.isnan(data.parallax.values) == True) | (np.isnan(data.phot_g_mean_mag.values) == True) | (np.isnan(data.phot_rp_mean_mag.values) == True) | (np.isnan(data.phot_bp_mean_mag.values) == True))[0]

    ##output the results
    ##output with starting timestamp for uniquness
    datestamp    = time.strftime('%Y%m%d-%H:%M:%S', time.localtime(time.time()))
    if arenou_clean == True: datestamp += '_arenouclean'
    usefname = outfilename.split('.')[0] + '_'+datestamp + '.'+outfilename.split('.')[1]
    
    ##Check for overwrite, should never happen with timestamping but who knows...
    if os.path.isfile(outdir+usefname) == True:
        print('Pick a new filename please then try again')
        print('Not outputting anything')
        return -1
        
    ##ok all good lets write the output
    data.to_csv(outdir+usefname,index=False,sep=',')
    
    ##if everything is fine give the user the final filename of the clean data.
    return outdir+usefname

def produce_data_gaia(data,writetofile=False):
    '''
    calculate mags, errors, and all those things for a 
    set of data. Assumes you have already cleaned the data for bad entries in 
    previous steps.  
    
    assumes input is in pandas format.
    

    '''
    
    ##Gaia flux errors, seems like these arent needed??
    #phot_bp_mean_flux_error = data.phot_bp_mean_flux.values *( 1/(data.phot_bp_mean_flux_over_error.values))
    #phot_rp_mean_flux_error = data.phot_rp_mean_flux.values *( 1/(data.phot_rp_mean_flux_over_error.values))
    #phot_g_mean_flux_error  = data.phot_g_mean_flux.values *( 1/(data.phot_g_mean_flux_over_error.values))
    
    ##Gaia mag errors
    mag_bp_error = np.abs(-2.5*(data.phot_bp_mean_flux_error.values/(data.phot_bp_mean_flux.values*np.log(10))))
    mag_rp_error = np.abs(-2.5*(data.phot_rp_mean_flux_error.values/(data.phot_rp_mean_flux.values*np.log(10))))
    mag_g_error  = np.abs(-2.5*(data.phot_g_mean_flux_error.values/(data.phot_g_mean_flux.values*np.log(10))))
    ##Gaia absolute mags 
    abs_mag_g    = data.phot_g_mean_mag.values+5-5*np.log10(1000.0/data.parallax.values)
    abs_mag_bp   = data.phot_bp_mean_mag.values+5-5*np.log10(1000.0/data.parallax.values)
    abs_mag_rp   = data.phot_rp_mean_mag.values+5-5*np.log10(1000.0/data.parallax.values)
    
    ##These shouldnt be needed??
    #parallax_error_np_array= np.array([data_3.parallax_error.values])
    #parallax_np_array = np.array([data_3.parallax.values])
    #mag_g_error_np_array= np.array([mag_g_error])
    #mag_bp_error_np_array= np.array([mag_bp_error])
    #mag_rp_error_np_array= np.array([mag_rp_error])
    
    ##Gaia Absolute mag errors 
    abs_mag_g_error  = np.sqrt((mag_g_error)**2 + np.abs(5*(data.parallax_error.values/data.parallax.values)/(np.log(10)))**2)
    abs_mag_bp_error = np.sqrt((mag_bp_error)**2 + np.abs(5*(data.parallax_error.values/data.parallax.values)/(np.log(10)))**2)
    abs_mag_rp_error = np.sqrt((mag_rp_error)**2 + np.abs(5*(data.parallax_error.values/data.parallax.values)/(np.log(10)))**2)

    if writetofile != False:
        pickle.dump((abs_mag_g,abs_mag_bp,abs_mag_rp,abs_mag_g_error,abs_mag_bp_error,abs_mag_rp_error,data.source_id.values),open(writetofile,'wb'))
        print('wrote the input variable to ' + writetofile)
    return abs_mag_g,abs_mag_bp,abs_mag_rp,abs_mag_g_error,abs_mag_bp_error,abs_mag_rp_error,data.source_id.values
    

def sample_generate_2d(nsamples=2000000, regen=False, outfilename='isochrones/tmp.pkl',readfile=None, modelfile=None,agelims=[1.0,3000.0],masslims=[0.1,1.0],use_mini=True):
    '''
    Function to generate samples for the 2D version of the model (where metallicity isn't included)
    
    ACR!!!!: Dont change this function, make a new one below called e.g. sample_generate_3d for the metallicity verion and copy this code over to it
    '''
   
   
    if regen == False:
        print ('Attempting to read previously generated samples')
        if readfile == None: 
            print('Give me a read file if you dont want to regen!!')
            return -1
        if os.path.isfile(readfile) == False:
            print('The read file is non-existant!!')
            print('You tried to read: ' + readfile + ' please try another ')
            return -1
        f = read_pickle_file(readfile)
#        f = pickle.load(open(readfile, 'rb'))
        return f
    
    if regen == True:
        print ('Generating Samples')
        if modelfile == None:
            print('Please give me a model file if you want to generate samples')
            return -1
        if os.path.isfile(modelfile)==False:
            print('The modelfile is non-existant!!')
            print('You tried to read: ' + modelfile + ' please try another ')
            return -1
            
        model_ages,model_array=pickle.load(open(modelfile,'rb'))
        model_stars = np.zeros((nsamples,6),dtype=float)
        
        random_age = 10**(np.random.uniform(np.log10(agelims[0]), np.log10(agelims[1]),size=nsamples))
        
        
        
        random_mass = np.random.uniform(masslims[0],masslims[1],size=nsamples)
        random_binprob = np.random.uniform(0.0,1.0,size=nsamples)
        binary_flag = np.zeros(len(random_mass)) -1
       # pdb.set_trace()
        for i in range(len(random_age)):
            low_indx = np.where(model_ages <= random_age[i])[0][-1]
            high_indx = np.where(model_ages > random_age[i])[0][0]
            
            lowmodel  = model_array[low_indx]
            highmodel = model_array[high_indx]
        
            ## ACR!!: Probably a good idea to sample on Mini (initial mass) rather 
            ##than Mass (current mass), this will eliminate interpolation problems. 
            if use_mini ==True:
                lmass = lowmodel.Mini
                hmass = highmodel.Mini
            else:
                lmass = lowmodel.Mass
                hmass = highmodel.Mass        
        
        
            low_gmag = np.interp(random_mass[i],lmass,lowmodel.Gmag)
            low_bpmag = np.interp(random_mass[i],lmass,lowmodel.G_BPmag)
            low_rpmag = np.interp(random_mass[i],lmass,lowmodel.G_RPmag)
            high_gmag = np.interp(random_mass[i],hmass,highmodel.Gmag)
            high_bpmag = np.interp(random_mass[i],hmass,highmodel.G_BPmag)
            high_rpmag = np.interp(random_mass[i],hmass,highmodel.G_RPmag)
            
            if random_binprob[i] < 0.44:
                random_companion_mass = np.random.uniform(np.max([lmass[0],hmass[0]]),random_mass[i],size=1)[0]
                binary_flag[i] = random_companion_mass/random_mass[i]
                ##now do the above interpolation on the grid but for the companion random mass
                
                ##this gives you a set of low/high G/R/B mags foe the companion
               
                
               # pdb.set_trace()
                low_gmag_comp  = np.interp(random_companion_mass,lmass,lowmodel.Gmag)
                high_gmag_comp = np.interp(random_companion_mass,hmass,highmodel.Gmag)
                low_bpmag_comp  = np.interp(random_companion_mass,lmass,lowmodel.G_BPmag)
                high_bpmag_comp = np.interp(random_companion_mass,hmass,highmodel.G_BPmag)
                low_rpmag_comp  = np.interp(random_companion_mass,lmass,lowmodel.G_RPmag)
                high_rpmag_comp = np.interp(random_companion_mass,hmass,highmodel.G_RPmag)
                
                low_gmag  = -2.5*np.log10(10**(-0.4*low_gmag) + 10**(-0.4*low_gmag_comp))
                high_gmag = -2.5*np.log10(10**(-0.4*high_gmag) + 10**(-0.4*high_gmag_comp))
                low_bpmag  = -2.5*np.log10(10**(-0.4*low_bpmag) + 10**(-0.4*low_bpmag_comp))
                high_bpmag = -2.5*np.log10(10**(-0.4*high_bpmag) + 10**(-0.4*high_bpmag_comp))
                low_rpmag  = -2.5*np.log10(10**(-0.4*low_rpmag) + 10**(-0.4*low_rpmag_comp))
                high_rpmag = -2.5*np.log10(10**(-0.4*high_rpmag) + 10**(-0.4*high_rpmag_comp))
                
                ##add them to the corresponding values for the primary
                
            
            
            
            model_stars[i,0] =(random_age[i])
            model_stars[i,1] =(random_mass[i])
            model_stars[i,2] = np.interp(random_age[i],[model_ages[low_indx],model_ages[high_indx]],[low_gmag,high_gmag])
            model_stars[i,3] = np.interp(random_age[i],[model_ages[low_indx],model_ages[high_indx]],[low_bpmag,high_bpmag])
            model_stars[i,4] = np.interp(random_age[i],[model_ages[low_indx],model_ages[high_indx]],[low_rpmag,high_rpmag])
            model_stars[i,5] = binary_flag[i]
          
            if np.mod(i+1,10000) == 0: print('up to ' + str(i+1) + ' out of ' + str(nsamples))
            
           
        pickle.dump(model_stars, open(outfilename,'wb'))
        print('outputting samples to ' + outfilename + ' and returning the samples to user')
        return model_stars
    
 
#I think I'm most confused where/how we are readin in the stacked files. I'm unsure how to read in multiple files and have it pull from all of those, thats a little confusing to me. I did make the stacked files and they're labeled by metallicity in the isochrones folder in taurus_project
def sample_generate_3d(nsamples=2000000, regen=False, outfilename='isochrones_metallicity/tmp.pkl',readfile=None, modelfile=None,agelims=[1.0,3000.0],masslims=[0.1,1.0],use_mini=True):
    '''
    Function to generate samples for the 3D version of the model (where metallicity is included)
    '''
   
   
    if regen == False:
        print ('Attempting to read previously generated samples')
        if readfile == None: 
            print('Give me a read file if you dont want to regen!!')
            return -1
        if os.path.isfile(readfile) == False:
            print('The read file is non-existant!!')
            print('You tried to read: ' + readfile + ' please try another ')
            return -1
        f = read_pickle_file(readfile)
#        f = pickle.load(open(readfile, 'rb'))
        return f
    
    if regen == True:
        print ('Generating Samples')
        if modelfile == None:
            print('Please give me a model file if you want to generate samples')
            return -1
        if os.path.isfile(modelfile)==False:
            print('The modelfile is non-existant!!')
            print('You tried to read: ' + modelfile + ' please try another ')
            return -1
            
        model_ages,model_array=pickle.load(open(modelfile,'rb'))
        model_stars = np.zeros((nsamples,9),dtype=float)
        
        random_age = 10**(np.random.uniform(np.log10(agelims[0]), np.log10(agelims[1]),size=nsamples))
        
        random_metallicity = np.random.uniform(-0.3, 0.3, size=nsamples)#actual z ranges from 0.00761804595 to 0.03032798718
        random_mass = np.random.uniform(masslims[0],masslims[1],size=nsamples)
        random_binprob = np.random.uniform(0.0,1.0,size=nsamples)
        random_tripprob = np.random.uniform(0.0,1.0,size=nsamples)
        binary_flag = np.zeros(len(random_mass)) -1
        triple_flag_1 = np.zeros(len(random_mass)) -1 #Not sure if this is the best way to do this
        triple_flag_2 = np.zeros(len(random_mass)) -1 #Not sure if this is the best way to do this
       # pdb.set_trace()
        for i in range(len(random_age)):
            low_indx = np.where(model_ages <= random_age[i])[0][-1]
            high_indx = np.where(model_ages > random_age[i])[0][0]
            
            lowmodel  = model_array[low_indx]
            highmodel = model_array[high_indx]
        
            ## ACR!!: Probably a good idea to sample on Mini (initial mass) rather 
            ##than Mass (current mass), this will eliminate interpolation problems. 
            if use_mini ==True:
                lmass = lowmodel.Mini
                hmass = highmodel.Mini
            else:
                lmass = lowmodel.Mass
                hmass = highmodel.Mass        
        
        
            low_gmag = np.interp(random_mass[i],lmass,lowmodel.Gmag)
            low_bpmag = np.interp(random_mass[i],lmass,lowmodel.G_BPmag)
            low_rpmag = np.interp(random_mass[i],lmass,lowmodel.G_RPmag)
            high_gmag = np.interp(random_mass[i],hmass,highmodel.Gmag)
            high_bpmag = np.interp(random_mass[i],hmass,highmodel.G_BPmag)
            high_rpmag = np.interp(random_mass[i],hmass,highmodel.G_RPmag)
            
            if random_binprob[i] < 0.44:
                random_companion_mass = np.random.uniform(np.max([lmass[0],hmass[0]]),random_mass[i],size=1)[0]
                binary_flag[i] = random_companion_mass/random_mass[i]
                ##now do the above interpolation on the grid but for the companion random mass
                
                ##this gives you a set of low/high G/R/B mags foe the companion
               
         
               # pdb.set_trace()
                low_gmag_comp  = np.interp(random_companion_mass,lmass,lowmodel.Gmag)
                high_gmag_comp = np.interp(random_companion_mass,hmass,highmodel.Gmag)
                low_bpmag_comp  = np.interp(random_companion_mass,lmass,lowmodel.G_BPmag)
                high_bpmag_comp = np.interp(random_companion_mass,hmass,highmodel.G_BPmag)
                low_rpmag_comp  = np.interp(random_companion_mass,lmass,lowmodel.G_RPmag)
                high_rpmag_comp = np.interp(random_companion_mass,hmass,highmodel.G_RPmag)
                
                low_gmag  = -2.5*np.log10(10**(-0.4*low_gmag) + 10**(-0.4*low_gmag_comp))
                high_gmag = -2.5*np.log10(10**(-0.4*high_gmag) + 10**(-0.4*high_gmag_comp))
                low_bpmag  = -2.5*np.log10(10**(-0.4*low_bpmag) + 10**(-0.4*low_bpmag_comp))
                high_bpmag = -2.5*np.log10(10**(-0.4*high_bpmag) + 10**(-0.4*high_bpmag_comp))
                low_rpmag  = -2.5*np.log10(10**(-0.4*low_rpmag) + 10**(-0.4*low_rpmag_comp))
                high_rpmag = -2.5*np.log10(10**(-0.4*high_rpmag) + 10**(-0.4*high_rpmag_comp))
                
                ##add them to the corresponding values for the primary
            if random_tripprob[i]: ##< #I don't have the Raghavan paper with this info, so if you could sent that to me I can add in                                     the value here
                random_companion_1_mass = np.random.uniform(np.max([lmass[0],hmass[0]]),random_mass[i],size=1)[0]
                random_companion_2_mass = np.random.uniform(np.max([lmass[0],hmass[0]]),random_mass[i],size=1)[0]
                triple_flag_1[i] = random_companion_1_mass/random_mass[i] #Not sure the best way/how to do this; I'm assuming we want to just divide the mass of the original target object by 3 and create 3 new objects?
                triple_flag_2[i] =random_companion_2_mass/random_mass[i] #Not sure the best way/how to do this; I'm assuming we want to just divide the mass of the original target object by 3 and create 3 new objects?
                '''low_gmag_comp  = np.interp(random_companion_mass,lmass,lowmodel.Gmag)
                high_gmag_comp = np.interp(random_companion_mass,hmass,highmodel.Gmag)
                low_bpmag_comp  = np.interp(random_companion_mass,lmass,lowmodel.G_BPmag)
                high_bpmag_comp = np.interp(random_companion_mass,hmass,highmodel.G_BPmag)
                low_rpmag_comp  = np.interp(random_companion_mass,lmass,lowmodel.G_RPmag)
                high_rpmag_comp = np.interp(random_companion_mass,hmass,highmodel.G_RPmag)
                
                low_gmag  = -2.5*np.log10(10**(-0.4*low_gmag) + 10**(-0.4*low_gmag_comp))
                high_gmag = -2.5*np.log10(10**(-0.4*high_gmag) + 10**(-0.4*high_gmag_comp))
                low_bpmag  = -2.5*np.log10(10**(-0.4*low_bpmag) + 10**(-0.4*low_bpmag_comp))
                high_bpmag = -2.5*np.log10(10**(-0.4*high_bpmag) + 10**(-0.4*high_bpmag_comp))
                low_rpmag  = -2.5*np.log10(10**(-0.4*low_rpmag) + 10**(-0.4*low_rpmag_comp))
                high_rpmag = -2.5*np.log10(10**(-0.4*high_rpmag) + 10**(-0.4*high_rpmag_comp))''' #again not sure the best way to do this so we can talk about if first
            
            
            
            
            
            model_stars[i,0] =(random_age[i])
            model_stars[i,1] =(random_mass[i])
            model_stars[i,2] = np.interp(random_age[i],[model_ages[low_indx],model_ages[high_indx]],[low_gmag,high_gmag])
            model_stars[i,3] = np.interp(random_age[i],[model_ages[low_indx],model_ages[high_indx]],[low_bpmag,high_bpmag])
            model_stars[i,4] = np.interp(random_age[i],[model_ages[low_indx],model_ages[high_indx]],[low_rpmag,high_rpmag])
            model_stars[i,5] = binary_flag[i]
            model_stars[i,6] = (random_metallicity[i])
            model_stars[i,7] = triple_flag_1[i] #Again, not sure if this is right
            model_stars[i,8] = triple_flag_2[i] #Again, not sure if this is right
          
            if np.mod(i+1,10000) == 0: print('up to ' + str(i+1) + ' out of ' + str(nsamples))
            
           
        pickle.dump(model_stars, open(outfilename,'wb'))
        print('outputting samples to ' + outfilename + ' and returning the samples to user')
        return model_stars

    
    
        
def probability_calculation_singlestar(random_samples,G,BP,RP,sig_G,sig_BP,sig_RP,sourceID,run_2mass=False,showhist=False,histsave=''):
    '''
    Calculate the probability of the data given the star model population 
    Assumes G, BP, RP, sig_G, sig_BP, sig_RP are all floats (not arrays)
    Can be added to to include more magnitudes
       
    '''
    random_age = random_samples[:,0]
    random_mass = random_samples[:,1]
    interp_gmag = random_samples[:,2]
    interp_bpmag = random_samples[:,3]
    interp_rpmag = random_samples[:,4]
    binaryflag  = random_samples[:,5]
    
    prob_data_modgmag = 1/np.sqrt(2*np.pi*sig_G**2)*np.exp(-(interp_gmag - G)**2/2.0/sig_G**2)
    prob_data_modrpmag = 1/np.sqrt(2*np.pi*sig_RP**2)*np.exp(-(interp_rpmag - RP)**2/2.0/sig_RP**2)
    prob_data_modbpmag = 1/np.sqrt(2*np.pi*sig_BP**2)*np.exp(-(interp_bpmag - BP)**2/2.0/sig_BP**2)
    priors =  (1/((1/(-1.35))*((1.0**(-1.35))-(0.1)**((-1.35)))))*(random_mass)**(-2.35)
    
    ##we've missed a prior here: our log-flat age sampling builds in a bias to smaller ages:
    ##to fixe this we need to correct for this by doing a variable transform from Log10X to X
    
    dydx = np.log10(np.exp(1.0))*1/random_age
    prob_age = dydx*1.0 ##f(x) = f(y)dy/dx with f(y) = 1.0 (uniform)
    
    prob_data_mod = prob_data_modgmag*prob_data_modgmag*prob_data_modbpmag*priors/prob_age
    
    ##ACR!!: Put 2mass JHK work here in future!!
    if run_2mass == True:
        print(' 2MASS IR photometry not implemented yet !!!!! figure that out first then try me!')
       

    ##calculate confidence ranges:
    
    if (np.sum(prob_data_mod) == 0):
        age,mass,sig_age,sig_mass = -1000.,-1000.,-1000.,-1000.0
    else:
        age       = np.average(np.log10(random_age), weights = prob_data_mod)
        sig_age   = np.sqrt(np.average((np.log10(random_samples[:,0])-age)**2, weights=prob_data_mod))
        mass      = np.average(random_mass, weights = prob_data_mod)
        sig_mass  = np.sqrt(np.average((random_mass-mass)**2, weights=prob_data_mod))
        ##unlog stuff
        sig_age  = np.sqrt(sig_age**2/np.log(10.0)**2/age**2)
        age = 10**age

##plot histogram to the user??
    if showhist == True: show_me_histograms(prob_data_mod,random_age,random_mass,savename=histsave,ar=[np.max([age-50*sig_age,1.0]),np.min([3000.0,age+50*sig_age])],mr=[mass-20*sig_mass,mass+20*sig_mass])
    
        
    ##output the probability to the caller, also the best numbers
    return prob_data_mod,age,mass,sig_age,sig_mass

    
def show_me_histograms(prob,model_age,model_mass,savename = '',ar=[1.0,3000.],mr=[0.1,1.0]):
    '''
    Make Jordyn's histogram plots for age/mass.
    '''
    agehist, agebin_edges   = np.histogram(np.log10(model_age), weights = prob, bins = 200,density=True, normed=True)
    masshist, massbin_edges = np.histogram(model_mass, weights = prob, bins = 100,density=True, normed=True)
    
    fig,axarr = plt.subplots(2)
    axarr[0].plot(10**(agebin_edges[0:len(agehist)]),agehist,'b-')
    axarr[0].set_xlabel('Age (Myr)')
    axarr[0].set_xlim(ar)
    axarr[1].plot(massbin_edges[0:len(masshist)],masshist,'b-')
    axarr[1].set_xlabel(r'Mass (M$_\mathrm{\odot}$)')
    axarr[1].set_xlim(mr)
    plt.tight_layout()
    if savename != '': plt.savefig(savename)
    plt.show()
    
def probability_calculation_all(random_samples,G,BP,RP,sig_G,sig_BP,sig_RP,sourceID,run_2mass=False,logfile=None):   
    '''
    mag and error inputs are now arrays
    This just loops probability_calculation_singlestar for each stars
    '''
    
    if (len(G) != len(BP)) | (len(G) != len(RP)) | (len(G) != len(sig_G)) | (len(G) != len(sig_RP)) | (len(G) != len(sig_BP)) : 
        print('Input photometry/error arrays are not the same size!!!!!!')
        if logfile != None:
            logfile.write('pcalc: inputs wrong sizes, this cant happen so.... \n') 
            logfile.flush()
        return -1
    
    age_result  = np.zeros(len(G))
    mass_result = np.zeros(len(G))
    age_error   = np.zeros(len(G))
    mass_error  = np.zeros(len(G))
    tstart = time.time()
    for i in range(len(G)):
        if logfile != None:
            logfile.write('pcalc: ' + str(i) + ' of ' + str(len(G)) + ' ' + str(sourceID[i]) + ' \n')
            logfile.flush()
        prob_model_given_data_star,expage,expmass,sigage,sigmass = probability_calculation_singlestar(random_samples,G[i],BP[[i]],RP[[i]],sig_G[[i]],sig_BP[[i]],sig_RP[[i]],sourceID[[i]],run_2mass=run_2mass,showhist=False)
        age_result[i]  = expage*1.0
        mass_result[i] = expmass*1.0
        age_error[i]   = sigage*1.0
        mass_error[i]  = sigmass*1.0
    
        if np.mod(i+1,100) == 0:
            tneed = (time.time()-tstart)/(i+1)*(len(G)-i-1)
            print('Up to ' + str(i+1) + ' out of ' + str(len(G)) + ' stars')
            print('Will finish at ' + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(tneed+time.time())))
    if logfile != None:
        logfile.write('finished pcalc, returning \n')
        logfile.flush()
    return age_result,mass_result,age_error,mass_error


 
  