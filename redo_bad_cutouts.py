#----------Setup----------------
#imports
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs import utils
import numpy as np 
from astroquery.skyview import SkyView
import os
import astropy.units as u
from astropy.coordinates import SkyCoord
import warnings
warnings.filterwarnings('ignore')
import glob, os, sys
import img_scale
from matplotlib.patches import Rectangle
import gc
import argparse


gc.enable()



#read in the file with the sources
sourcefile = 'cat_above5000.txt'



SB=9351
directory ='/Volumes/TARDIS/Work/EMUzoo/alpha_test/WISE_DSS_cutouts/'
filedirectory='/Volumes/TARDIS/Work/EMUzoo/alpha_test/'
sources=np.loadtxt(filedirectory+sourcefile,dtype='str')

#specify if you want to overwrite existing images or not
overwrite=False

def fitsopen(filepath):
	#open a fits file and return the data and header
	hdu=fits.open(filepath)
	header=hdu[0].header
	data=hdu[0].data
	hdu.close()
	return(data,header)

def main():
	sources=np.loadtxt(filedirectory+sourcefile,dtype='str')
	baddes=np.loadtxt(filedirectory+'baddes.txt' ,dtype='str')
	badwise=np.loadtxt(filedirectory+'badwise.txt' ,delimiter=' ',dtype='str')
	#badwise=['J203623-551941', 'J202827-571807', 'J202112-563655']
	nsrc=sources.shape[0]

	survey_pix=np.loadtxt('survey_pix.txt',dtype='str')
	commonpixscale=False
	pixscale=None #arcseconds

	imsize=12*60. #arcminutes

	#	specify if you want to overwrite existing images or not
	overwrite=False
	bleh=['J203455-565256', 'J203018-534949', 'J204553-565854','J205720-565940', 'J205319-552733', 'J204124-535646']

	for i in range (0,nsrc):
		src=sources[i,6]
		if src in bleh:
			surveys= ['DSS','WISE 3.4']
		#if src in badwise and src not in baddes:
			#surveys= ['WISE 3.4']
		#elif src in baddes and src not in badwise:
			#surveys= ['DSS']
		#elif src in baddes and src in badwise:
			#surveys= ['WISE 3.4','DSS']
		else:
			surveys=[]

		if len(surveys)!=0:

			print("downloading for src {}".format(src))
			#get source coords and download image
			#coordsstr=sources[i,1]
			ra_deg=sources[i,9]
			dec_deg=sources[i,10]
			coordstr=ra_deg+dec_deg
			coords = SkyCoord(coordstr,unit="deg")  

			EMUSB=float(sources[i,36])

			#loop over surveys
			for j in range (0,len(surveys)):
				survey=surveys[j]
				survey_nospace=survey.replace(" ","_")

				if commonpixscale==False:
					#get the pixel scale of the survey
					#otherwise will use the value already specified in the cell above
					idx=np.where(survey_pix==survey_nospace)
					try:
						pixscale=float(survey_pix[idx[0][0],1])
					except:
						print("Error: no pixscale specified in file for {}".format(survey_nospace))
						break

				# make the image name
				imagename=src+'_'+survey_nospace+'.fits'

				n_pix=int(imsize/pixscale)

				#only download if file doesn't already exist 
				if (os.path.isfile(directory+imagename)==False) or (overwrite==True):
					try:
						url=SkyView.get_image_list(position=coords,survey=[survey],coordinates='J2000',projection='Sin',resolver='NED',pixels=n_pix,height=imsize*u.arcsec,width=imsize*u.arcsec)
						os.system('curl -o '+directory+imagename+' '+url[0])
						print("Downloaded {} image for source {}".format(survey,src))
					except:
						print("Survey {} not found in SkyView".format(survey))
				else:
					print("Not downloading {} image for source {}".format(survey,src)) 
		else:
			print("no issue with this source")

if __name__ == "__main__":
	main()









