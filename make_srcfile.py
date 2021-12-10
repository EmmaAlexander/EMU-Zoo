import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

bmaj=12.0348
bmin=10.712
beamsize=np.pi*0.25*bmaj*bmin #arcsec
fullcatfile='/Volumes/TARDIS/Work/askap/fullfields/srcfind/9351/AS101_Continuum_Component_Catalogue_9351_78.csv'
fullcat = np.genfromtxt(fullcatfile, delimiter=',',dtype='str')

print(fullcat.shape)
keptsrcs=[]

for i in range(1,fullcat.shape[0]):
	component_name=fullcat[i,7]
	maj_axis=float(fullcat[i,19])
	min_axis=float(fullcat[i,20])

	if maj_axis >=1.5*bmaj and min_axis >= 1.5*bmin:
		keptsrcs.append(fullcat[i,:])

keptsrcsarr=np.asarray(keptsrcs)
print(keptsrcsarr.shape)

component_names=keptsrcsarr[:,7]
coords=SkyCoord(keptsrcsarr[:,10], keptsrcsarr[:,11], unit=u.degree,frame='fk5')

doubles=[]
triples=[]

simples=[]

for i in range(0,keptsrcsarr.shape[0]):
	coord=coords[i]
	sep=coord.separation(coords)/u.arcsec
	neighbour=np.where(sep<45,1,np.nan)
	object1=keptsrcsarr[i,7]
	neighbour_objects=component_names[np.isfinite(neighbour)]
	if len(neighbour_objects) ==2:
		#picked up more than just itself
		doubles.append(neighbour_objects)
	elif len(neighbour_objects) ==3:
		triples.append(neighbour_objects)
	elif len(neighbour_objects) ==1:
		#just itself
		coordstr=coord.to_string('hmsdms').replace(' ',',')
		output=object1+' '+coordstr+' 9351'
		simples.append(output)
	else:
		print("ooh maybe a quadruple source??")

doubles=np.unique(np.asarray(doubles))
triples=np.unique(np.asarray(triples))

simples=np.asarray(simples)

print(doubles)
print(triples)
print(simples.shape)

np.savetxt('SB_9351_viewcat.txt',simples,fmt='%s')

