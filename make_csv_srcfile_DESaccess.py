import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

fullcatfile='/Volumes/TARDIS/Work/askap/fullfields/srcfind/9351/AS101_Continuum_Component_Catalogue_9351_78.csv'
fullcat = np.genfromtxt(fullcatfile, delimiter=',',dtype='str')
source_list='/Users/emma/OneDrive/share/EMUzoo/alpha_test/cx_below5000_150.txt'
sources=np.genfromtxt(source_list,dtype='str')

headers=fullcat[0,:]
component_names=fullcat[1:,7]
ra_deg=np.asarray(fullcat[1:,10],dtype=float)
dec_deg=np.asarray(fullcat[1:,11],dtype=float)
stacked=np.stack([ra_deg,dec_deg])
nsrc=len(sources)
ncat=len(component_names)
print(stacked.shape)

radeclist=[]


for i in range(0,nsrc):
	#print(sources[i])
	loc=np.argwhere(component_names==sources[i])[0][0]
	radeclist.append(stacked[:,loc])

out=np.asarray(radeclist)


fname='SB_9351_DES_alpha_cx_below5000_150.csv'
np.savetxt(fname,out,delimiter=',',fmt=['%.6f','%.6f'],header='RA,DEC',comments='')

#use the following if you have a very long source list
'''
maxsrc=2500
nfiles=nsrc/maxsrc

for i in range(0,int(np.floor(nfiles))):
	nmin=i*maxsrc
	nmax=(i+1)*maxsrc
	out=np.transpose(stacked[:,nmin:nmax])
	fname='SB_9351_DES'+str(nmin)+'-'+str(nmax)+'.csv'
	np.savetxt(fname,out,delimiter=',',fmt=['%.6f','%.6f'],header='RA,DEC',comments='')

nmin=maxsrc*int(np.floor(nfiles))
nmax=nsrc
out=np.transpose(stacked[:,nmin:nmax])
print(out.shape)
fname='SB_9351_DES_alpha'+str(nmin)+'-'+str(nmax)+'.csv'
np.savetxt(fname,out,delimiter=',',fmt=['%.6f','%.6f'],header='RA,DEC',comments='')
'''