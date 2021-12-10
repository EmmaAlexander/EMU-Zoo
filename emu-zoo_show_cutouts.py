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
sourcefile = 'SB_9351_viewcat.txt'
directory ='/Users/emma/Dropbox/Public/EMUZoo/'
sources=np.loadtxt(sourcefile,dtype='str')

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'cmu serif'
plt.rcParams.update({'lines.linewidth':0.3})
plt.rc('font', size=10)          # controls default text sizes

# list of survey images you wish to download 
# for names see https://astroquery.readthedocs.io/en/latest/skyview/skyview.html
# or use SkyView.list_surveys()
# change these in and out as you wish
# note that if a survey is not available for the coordinates you specify, then it will download an empty fits
surveys= ['WISE 3.4','DSS','SUMSS 843 MHz','2MASS-J','GLEAM 170-231 MHz']

# annoyingly, SkyView doesn't register the native pixel size of surveys, so need to specify
# load in the file with these data from https://skyview.gsfc.nasa.gov/current/cgi/survey.pl
survey_pix=np.loadtxt('survey_pix.txt',dtype='str')

# alternatively, specify a common reolution you would like all the surveys to be sampled at
# beware resampling errors... 
commonpixscale=False
pixscale=None #arcseconds

imsize=3*60. #arcminutes

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
		src=sources[i,0]

		#get source coords and download image
		coordsstr=sources[i,1]
		coords = SkyCoord(coordsstr.replace(',', ''))  
		print(coords) 

		EMUSB=float(sources[i,2])

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

		print("Making .png images for source") 

		# background greyscale
		greyscale_im=directory+src+'_DSS.fits'
		greyscale_label='DSS'
		greyscale_data,greyscale_header=fitsopen(greyscale_im)
		wcs=WCS(greyscale_header)

		#radio contours
		radio_im=directory+src+'_EMU.fits'
		radio_label='EMU'
		radio_data,radio_header=fitsopen(radio_im)
		radio_wcs=WCS(radio_header)

		#WISE contours
		wise_im=directory+src+'_WISE_3.4.fits'
		wise_label='WISE 3.4'
		wise_data,wise_header=fitsopen(wise_im)
		wise_wcs=WCS(wise_header)

		centre=int(radio_data.shape[0]/2.)
		rvmax=np.nanmax(radio_data[centre-20:centre+20,centre-20:centre+20])

		contourexps=np.arange(start=0,stop=32,step=0.5)
		contourmults=np.power(2,contourexps)
		basecont=0.00012 #roughly 3x emu rms noise
		radio_contours = [basecont * i for i in contourmults]
		infernocmap=plt.cm.inferno
		infernocmap.set_bad('grey',1)

		fig=plt.figure(dpi=400,figsize=(9,3))
		ax1=plt.subplot(131,projection=wcs,adjustable='box')
		ax1.imshow(radio_data,origin='lower',transform=ax1.get_transform(radio_wcs),cmap=infernocmap,vmax=rvmax)
		ax1.contour(radio_data, transform=ax1.get_transform(radio_wcs), colors='green',levels=radio_contours,linewidths=0.2)		
		ax1.set_xlabel(' ')
		ax1.set_ylabel(' ')

		ax1.set_title('EMU')
		ax1.tick_params(
			axis='x',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom=True,      
			top=True,         
			labelbottom=False) 
		ax1.tick_params(
			axis='y',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom=True,      
			top=True,         
			left=True,
			right=True,
			labelbottom=False) 

		ax1.set_xlim(0,greyscale_data.shape[0])
		ax1.set_ylim(0,greyscale_data.shape[1])


		ax2=plt.subplot(132,projection=wcs,adjustable='box')
		ax2.imshow(greyscale_data,origin='lower',cmap=plt.cm.gray_r)
		ax2.contour(radio_data, transform=ax2.get_transform(radio_wcs), colors='red',levels=radio_contours,linewidths=0.2)


		ax2.set_xlabel(' ')
		ax2.set_ylabel(' ')

		ax2.tick_params(
			axis='x',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom=True,      
			top=True,         
			labelbottom=False) 
		ax2.tick_params(
			axis='y',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom=True,      
			top=True,         
			left=True,
			right=True,
			labelbottom=False) 

		ax2.set_title('EMU & DSS')

		ax2.set_xlim(0,greyscale_data.shape[0])
		ax2.set_ylim(0,greyscale_data.shape[1])


		ax3=plt.subplot(133,projection=wcs,adjustable='box')
		#ax3.imshow(greyscale_data,origin='lower',cmap=plt.cm.gray,vmax=0.5*np.nanmax(greyscale_data))
		ax3.imshow(wise_data,transform=ax3.get_transform(wise_wcs),origin='lower',cmap=plt.cm.gray_r)
		ax3.contour(radio_data, transform=ax3.get_transform(radio_wcs), colors='red',levels=radio_contours,linewidths=0.2)

		ax3.set_xlabel(' ')
		ax3.set_ylabel(' ')

		ax3.tick_params(
			axis='x',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom=True,      
			top=True,         
			labelbottom=False) 
		ax3.tick_params(
			axis='y',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom=True,      
			top=True,         
			left=True,
			right=True,
			labelbottom=False) 
		ax3.set_xlim(0,greyscale_data.shape[0])
		ax3.set_ylim(0,greyscale_data.shape[1])

		ax3.set_title('EMU & WISE 3.4')
		plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.9, wspace=0.02, hspace=0.02)
		plt.savefig(directory+src+'_overlay.png',dpi=400)
if __name__ == "__main__":
	main()
