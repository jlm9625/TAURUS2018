##function to calculate galactic XYZ positions given l,b, (or ra,dec) 
##and parallax (or distance)
import numpy as np
from ACRastro.glactc import glactc
#import time
##import pdb
def gal_xyz(ll,bb,dd,radec=False,plx=True,reverse=False):
	#tt1 = time.clock()
	if reverse == False:
	##convert to galactic coords if needed
		l=ll
		b=bb
		if radec == True:
			l,b = (glactc(ll,bb,2000))
			#l=l*180/np.pi
			#b=b*180/np.pi
	##pdb.set_trace
	##convert to distance if needed
		dist = dd
		if plx == True:
			dist = 1000.0/dd
	##Note the negative sign for x makes it positive away from the galactic centre, to match our 
	## use of U as positive to the galactic anticentre
		x = -dist*np.cos(np.pi/180.0*b)*np.cos(np.pi/180.0*l)
		y = dist*np.cos(np.pi/180.0*b)*np.sin(np.pi/180.0*l)
		z = dist*np.sin(np.pi/180.0*b)	
		return x,y,z
	
	if reverse == True:
		##in the case user supplies x,y,z(pc) and want's gall,galb,plx/dist
		x = ll
		y = bb
		z = dd
		dist = np.sqrt(x**2 + y**2 + z**2)
		##again negative sign in front of x is because x is positive away from galactic centre
		ll   = np.arctan2(y,-x)
		bb   = np.arcsin(z/dist)
		ll   = np.mod(ll+2*np.pi,2*np.pi)
		##print time.clock()-tt1
		if plx == True:
			return ll,bb,1000.0/dist
		if plx == False:
			return ll,bb,dist