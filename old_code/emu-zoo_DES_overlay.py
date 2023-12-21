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
import glob, os, sys
import img_scale
from matplotlib.patches import Rectangle



#read in the file with the sources
#file format should be (1) Name (2) Centre coordinates (3) EMU SB (not needed for this code)
sourcefile = 'SB_9351_fullcat.txt'
complexitiesfile='batch_1to10_fn_and_complexity_FITS_v2.csv'



SB=9351
directory ='/Users/emma/OneDrive/share/EMUzoo/SB9551/'
filedirectory='/Volumes/TARDIS/Work/EMUzoo/SB9351/'
sources=np.loadtxt(filedirectory+sourcefile,dtype='str')
dpi=250

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'cmu serif'
plt.rc('font', size=0.5)  
plt.rcParams.update({'lines.linewidth':0.5})
plt.rc('font', size=10)          # controls default text sizes
infernocmap=plt.cm.inferno
infernocmap.set_bad('black',1)
redcmap=plt.cm.Reds_r
redcmap.set_bad('black',1)
greycmap=plt.cm.Greys
greycmap.set_bad('white',1)

#specify if you want to overwrite existing images or not
overwrite=False

failedsources=[]

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

	#get the coordintes of the WISE images to know which to load and use
	WISEfiles=glob.glob(filedirectory+'WISE_fits/*.fits')
	WISEcoords=[]
	for file in WISEfiles:
		filename=file.split('/')[-1] #get the actual file name
		ra=float(filename[0:4])/10
		dec=-1*float(filename[5:8])/10
		WISEcoords.append(SkyCoord(ra,dec,unit=u.degree,frame='fk5'))

	WISEcoords=SkyCoord(np.asarray(WISEcoords))	

	complexities=np.loadtxt(filedirectory+complexitiesfile,delimiter=',',dtype='str')

	for j in range (0,complexities.shape[0]):#sources.shape[0]):2500

		filename=complexities[j][0]
		#if 'batch1/' in filename:
		if 'J204923-563857' in filename:

			src=filename[30:44]
			complexity=complexities[j,1]
			print(complexity)

			listpos=np.where(sources[:,0]==src)[0][0]
			i=listpos
	
			# extract source name from source list
			#src=sources[i,0]

			try:

				print(sources[i,:])

				#get source coords 
				coordsstr=sources[i,1]
				ra,dec=coordsstr.split(',')
				coords = SkyCoord(coordsstr.replace(',', ''),frame='fk5')  

				EMUrms=float(sources[i,2]) #mjy

				print("Making .png images for source {}".format(src)) 
				
				racode=ra.split('.')[0].replace('h', '').replace('m', '')
				deccode=dec.split('.')[0].replace('d', '').replace('m', '')
				DESstr=filedirectory+'DES/DESJ{}*{}*'.format(racode,deccode)
				
				#load the DES data
				R,Rhead=fitsopen(glob.glob(DESstr+'_six_i.fits')[0])
				G,Ghead=fitsopen(glob.glob(DESstr+'_six_r.fits')[0])
				B,Bhead=fitsopen(glob.glob(DESstr+'_six_g.fits')[0])
				des_wcs=WCS(Rhead)

				img = np.zeros((R.shape[0], R.shape[1], 3), dtype=float)
				img[:,:,0] = img_scale.log(R, scale_min=5, scale_max=2000)
				img[:,:,1] = img_scale.log(G, scale_min=5, scale_max=2000)
				img[:,:,2] = img_scale.log(B, scale_min=5, scale_max=2000)


				#radio contours
				radio_im=filedirectory+'cutouts/{}_SB{}.fits'.format(src,str(SB))
				print(radio_im)
				radio_label='EMU'
				radio_data,radio_header=fitsopen(radio_im)
				radio_wcs=WCS(radio_header)

				print(radio_data.shape)

				#WISE contours
				# find out which WISE image to load
				WISEsep=coords.separation(WISEcoords)/u.degree
				wise_im=WISEfiles[np.argmin(WISEsep)]
				#wise_im=directory+src+'_WISE_3.4.fits'
				wise_label='WISE 3.4'
				wise_data,wise_header=fitsopen(wise_im)
				wise_wcs=WCS(wise_header)

				print(np.nanmin(wise_data))

				contourexps=np.arange(start=0,stop=32,step=2)
				contourmults=np.power(2,contourexps)
				basecont=3*EMUrms/1000
				print(basecont)
				radio_contours = [basecont * i for i in contourmults]

				plt.figure(dpi= dpi,figsize=(6,6))
				# plot the greyscale
				ax1=plt.subplot(221,projection=radio_wcs)			
				ax1.imshow(radio_data,origin='lower',cmap=infernocmap,norm=colors.LogNorm(vmin=basecont/3, vmax=np.nanmax(radio_data)/2))
				ax1.contour(radio_data,levels=radio_contours,colors=['green'])

				ax1.set_xlim(90,271)
				ax1.set_ylim(90,271)

				ax2=plt.subplot(222,projection=radio_wcs)			
				ax2.imshow(radio_data,origin='lower',cmap=infernocmap,norm=colors.LogNorm(vmin=basecont/3, vmax=np.nanmax(radio_data)/2))
				ax2.imshow(wise_data,transform=ax2.get_transform(wise_wcs),origin='lower',cmap=greycmap,norm=colors.LogNorm(vmin=2.5,vmax=100))
				ax2.contour(radio_data,levels=radio_contours,colors=['green'])

				ax2.set_xlim(90,271)
				ax2.set_ylim(90,271)

				ax3=plt.subplot(223,projection=radio_wcs)
				ax3.imshow(radio_data,origin='lower',cmap=infernocmap,norm=colors.LogNorm(vmin=basecont/3, vmax=np.nanmax(radio_data)/2))
				ax3.contour(radio_data,levels=radio_contours,colors=['green'])

				ax3.set_xlim(0,361)
				ax3.set_ylim(0,361)

				ax4=plt.subplot(224,projection=radio_wcs)
				ax4.imshow(radio_data,origin='lower',cmap=infernocmap,norm=colors.LogNorm(vmin=basecont/3, vmax=np.nanmax(radio_data)/2))
				ax4.imshow(R,origin='lower',cmap=greycmap,transform=ax4.get_transform(des_wcs),vmax=200)
				#ax4.imshow(img,origin='lower',transform=ax4.get_transform(des_wcs))
				ax4.contour(radio_data,levels=radio_contours,colors=['green'])

				#ax4.set_xlim(135,226)
				#ax4.set_ylim(135,226)

				ax4.set_xlim(90,271)
				ax4.set_ylim(90,271)

				ax1.axis('off')
				ax2.axis('off')
				ax3.axis('off')
				ax4.axis('off')
				
				arcmin3=Rectangle(xy=(180-45,180-45),width=90,height=90,edgecolor='white',fc='none',alpha=1,lw=0.5)
				ax3.add_patch(arcmin3)

				arcmin6=Rectangle(xy=(180-90,180-90),width=180,height=180,edgecolor='white',fc='none',alpha=1,lw=0.5)
				ax3.add_patch(arcmin6)

				ax1.annotate(xy=(0.02,0.96),xycoords='figure fraction',text="Radio (6'x6')",c='white',ha="left")  
				ax2.annotate(xy=(0.51,0.96),xycoords='figure fraction',text="WISE + radio (6'x6')",c='black',ha="left")  
				ax3.annotate(xy=(0.02,0.47),xycoords='figure fraction',text="Radio (12'x12')",c='white',ha="left")  
				ax4.annotate(xy=(0.51,0.47),xycoords='figure fraction',text="DES-i + radio (6'x6')",c='black',ha="left")  
				
				ax3.annotate(xy=(0.145,0.3505),xycoords='figure fraction',text="6'x6'",c='white',ha="left")  
				ax3.annotate(xy=(0.205,0.295),xycoords='figure fraction',text="3'x3'",c='white',ha="left")  


				plt.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.99, wspace=0.01, hspace=0.01)



				plt.savefig(filedirectory+'png/{}_cplx{}_wbg_6x6des.png'.format(src,str(complexity)),dpi=dpi,bbox_inches='tight',transparent=True)

			except:
				failedsources.append(src)
				#sys.exit()

	print("Sources failed:")
	print(failedsources)

if __name__ == "__main__":
	main()