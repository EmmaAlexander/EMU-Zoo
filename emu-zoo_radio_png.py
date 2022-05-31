#----------Setup----------------
#imports
import matplotlib.pyplot as plt
import matplotlib.colors as colors

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
sourcefile = 'SB_9351_fullcat.txt'
directory ='/Users/emma/OneDrive/share/EMUzoo/new_cutouts/fitsfiles/batch10/'
sources=np.loadtxt(sourcefile,dtype='str')
nsrc=sources.shape[0]

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'cmu serif'
plt.rcParams.update({'lines.linewidth':0.3})
plt.rc('font', size=10)          # controls default text sizes

imsize=6*60. #arcminutes

#specify if you want to overwrite existing images or not
overwrite=False

def fitsopen(filepath):
	#open a fits file and return the data and header
	hdu=fits.open(filepath)
	header=hdu[0].header
	data=hdu[0].data
	hdu.close()
	return(data,header)

def vlim(array,percent):
    vlim_val=np.percentile(array[np.isfinite(array)],percent)
    return vlim_val

def main():
	# begin finding and downloading data
	# loop over sources in source list
	for i in range (22500,nsrc):
	
		# extract source name from source list
		src=sources[i,0]
		#get source coords and download image
		coordsstr=sources[i,1]
		coords = SkyCoord(coordsstr.replace(',', ''))  

		rms=float(sources[i,2])

		print("Making .png images for source {}".format(i)) 

		#radio contours
		radio_im=directory+src+'_EMU.fits'
		radio_label='EMU'
		radio_data,radio_header=fitsopen(radio_im)
		radio_wcs=WCS(radio_header)

		centrex=int(radio_data.shape[0]/2.)
		centrey=int(radio_data.shape[1]/2.)

		xmin=centrex-75
		xmax=centrex+75
		ymin=centrey-75
		ymax=centrey+75
		#rvmax=np.nanmax(radio_data[centre-20:centre+20,centre-20:centre+20])

		contourexps=np.arange(start=0,stop=32,step=0.5)
		contourmults=np.power(2,contourexps)
		basecont=(3*rms)/1000
		radio_contours = [basecont * i for i in contourmults]
		infernocmap=plt.cm.inferno
		infernocmap.set_bad('grey',1)

		dpi=72
		npix=150

		figsize=(npix/131)*(npix/dpi)


		fig=plt.figure(dpi=dpi,figsize=(figsize,figsize))
		ax1=plt.subplot(111,adjustable='box')
		ax1.imshow(radio_data[xmin:xmax,ymin:ymax],origin='lower',cmap=infernocmap,vmax=vlim(radio_data[int(xmin+1):int(xmax-1),int(ymin+1):int(ymax-1)],99.99))
		#ax1.contour(radio_data, colors='green',levels=radio_contours,linewidths=0.2)		
		ax1.set_xlabel(' ')
		ax1.set_ylabel(' ')
		ax1.tick_params(
			axis='x',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom=False,      
			top=False,         
			labelbottom=False) 
		ax1.tick_params(
			axis='y',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom=False,      
			top=False,         
			left=False,
			right=False,
			labelbottom=False) 
		plt.axis('off')
		plt.tight_layout()
		plt.savefig('/Users/emma/OneDrive/share/EMUzoo/new_cutouts/pngfiles/batch10/'+src+'_5arcmin.png',dpi=dpi,transparent=True,bbox_inches='tight')

if __name__ == "__main__":
	main()