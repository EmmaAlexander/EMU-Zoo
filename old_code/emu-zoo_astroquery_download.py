#----------Setup----------------
#imports
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np 
from astroquery.skyview import SkyView
import os
import astropy.units as u
from astropy.coordinates import SkyCoord
import warnings
warnings.filterwarnings('ignore')


#read in the file with the sources
#file format should be (1) Name (2) Centre coordinates (3) EMU SB (not needed for this code)
sourcefile = 'cat_above5000.txt'
directory ='/Volumes/TARDIS/Work/EMUzoo/alpha_test/'
sources=np.loadtxt(directory+sourcefile,dtype='str')

baddes=np.loadtxt(directory+'baddes.txt',dtype='str')
#badwise=np.loadtxt(directory+'badwise.txt',dtype='str')


# list of survey images you wish to download 
# for names see https://astroquery.readthedocs.io/en/latest/skyview/skyview.html
# or use SkyView.list_surveys()
# change these in and out as you wish
# note that if a survey is not available for the coordinates you specify, then it will download an empty fits
surveys= ['WISE 3.4','DSS']#,'SUMSS 843 MHz','2MASS-J','GLEAM 170-231 MHz']

# annoyingly, SkyView doesn't register the native pixel size of surveys, so need to specify
# load in the file with these data from https://skyview.gsfc.nasa.gov/current/cgi/survey.pl
survey_pix=np.loadtxt('survey_pix.txt',dtype='str')

# alternatively, specify a common reolution you would like all the surveys to be sampled at
# beware resampling errors... 
commonpixscale=False
pixscale=None #arcseconds

imsize=12*60. #arcminutes

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
	# begin finding and downloading data
	# loop over sources in source list
	for i in range (0,sources.shape[0]):
		
		# extract source name from source list
		src=sources[i,6]
		print(src)
		if src in baddes or src in badwise:
			print("downloading for src {}".format(src))
			#get source coords and download image
			#coordsstr=sources[i,1]
			ra_deg=sources[i,9]
			dec_deg=sources[i,10]
			coordstr=ra_deg+dec_deg
			print(ra_deg)
			print(dec_deg)
			print("ohhh")
			coords = SkyCoord(coordstr.replace(',', ''))  
			print(coords) 

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