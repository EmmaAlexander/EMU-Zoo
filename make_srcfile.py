import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

fullcatfile='/Volumes/TARDIS/Work/askap/fullfields/srcfind/9351/AS101_Continuum_Island_Catalogue_9351_79.csv'
fullcat = np.genfromtxt(fullcatfile, delimiter=',',dtype='str')

headers=fullcat[0,:]

component_names=fullcat[1:,6]
coords=SkyCoord(fullcat[1:,10], fullcat[1:,11], unit=u.degree,frame='fk5')
local_rms=fullcat[1:,20]

outputlist=[]
outputlist.append('#srcname coords local_rms')

nsrc=len(component_names)
print(nsrc)
for i in range(0,nsrc):
	#print("Source {}/{}".format(i,nsrc))
	print(i)
	coord=coords[i]
	srcname=component_names[i]
	coordstrhms=coord.to_string('hmsdms').replace(' ',',')
	
	output=srcname+' '+coordstrhms+' '+str(local_rms[i])
	outputlist.append(output)


outputlist=np.asarray(outputlist)
np.savetxt('SB_9351_islandcat.txt',outputlist,fmt='%s')

