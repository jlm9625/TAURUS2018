import numpy as np
def gcirc(ra1,dc1,ra2,dc2,format=2):
	d2r    = np.pi/180.0
	as2r   = np.pi/648000.0
	h2r    = np.pi/12.0
	sin = np.sin
	cos = np.cos
	if format==0:
		rarad1 = ra1
		rarad2 = ra2
		dcrad1 = dc1
		dcrad2 = dc2
	if format==1:
		rarad1 = ra1*h2r
		rarad2 = ra2*h2r
		dcrad1 = dc1*d2r
		dcrad2 = dc2*d2r
	if format==2:
		rarad1 = ra1*d2r
		rarad2 = ra2*d2r
		dcrad1 = dc1*d2r
		dcrad2 = dc2*d2r
	if ((format!=1) & (format!=2) & (format!=0)):
		print "Choose a correct format for gcirc intput/output next time!!!"
		return -1
    
	radif = np.mod(2*np.pi+np.absolute(rarad2-rarad1),2*np.pi)

	cosdis = sin(dcrad1)*sin(dcrad2) + cos(dcrad1)*cos(dcrad2)*cos(radif)
	big    = np.where(np.absolute(cosdis) > 1.0)[0]
	if big.size!=0:
		cosdis[big] = cosdis[big]/np.absolute(cosdis[big])
	dis    = np.arccos(cosdis)
	if format != 0:
		dis = dis/as2r
	return dis