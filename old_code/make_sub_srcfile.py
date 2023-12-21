import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

fullcatfile='/Volumes/TARDIS/Work/askap/fullfields/srcfind/9351/AS101_Continuum_Component_Catalogue_9351_78.csv'
fullcat = np.genfromtxt(fullcatfile, delimiter=',',dtype='str')
source_list='/Volumes/TARDIS/Work/EMUzoo/alpha_test/cx_above5000_150.txt'
sources=np.genfromtxt(source_list,dtype='str')

headers=fullcat[0,:]
data=fullcat[1:,:]
nsrc=len(sources)
component_names=fullcat[1:,7]

fname='/Volumes/TARDIS/Work/EMUzoo/alpha_test/cat_above5000.txt'

out=[]

for i in range(0,nsrc):
	#print(sources[i])
	loc=np.argwhere(component_names==sources[i])[0][0]
	out.append(data[loc,:])

np.savetxt(fname,out,delimiter=' ',fmt='%s',comments='')
