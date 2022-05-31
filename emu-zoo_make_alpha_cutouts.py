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
sourcefile = 'cat_below5000.txt'



SB=9351
directory ='/Volumes/TARDIS/Work/EMUzoo/alpha_test/'
filedirectory='/Volumes/TARDIS/Work/EMUzoo/alpha_test/'
sources=np.loadtxt(filedirectory+sourcefile,dtype='str')
dpi=684

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'cmu serif'
plt.rc('font', size=0.5)  
plt.rcParams.update({'lines.linewidth':0.5})
plt.rc('font', size=10)          # controls default text sizes
infernocmap=plt.cm.inferno
infernocmap.set_bad('black',1)
magmacmap=plt.cm.magma
magmacmap.set_bad('black',1)
redcmap=plt.cm.Reds_r
redcmap.set_bad('black',1)
orangecmap=plt.cm.Oranges_r
orangecmap.set_bad('black',1)
greycmap=plt.cm.Greys_r
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

def main(args):
	# begin finding and downloading data
	# loop over sources in source list
	start=args.start[0]
	end=args.end[0]

	print(start)
	print(end)

	#get the coordintes of the WISE images to know which to load and use
	WISEfiles=glob.glob('/Volumes/TARDIS/Work/EMUzoo/SB9351_prealpha/WISE_fits/*.fits')
	WISEcoords=[]
	for file in WISEfiles:
		filename=file.split('/')[-1] #get the actual file name
		ra=float(filename[0:4])/10
		dec=-1*float(filename[5:8])/10
		WISEcoords.append(SkyCoord(ra,dec,unit=u.degree,frame='fk5'))
	WISEcoords=SkyCoord(np.asarray(WISEcoords))	

	sources=np.loadtxt(filedirectory+sourcefile,dtype='str')
	baddes=np.loadtxt(filedirectory+'baddes.txt' ,dtype='str')
	badwise=np.loadtxt(filedirectory+'badwise.txt' ,delimiter=' ',dtype='str')
	nodes=['J202120-572034', 'J202154-561930', 'J202343-560746', 'J203047-562751', 'J202707-571728', 'J202827-571807', 'J202214-563129', 'J203022-563134', 'J202959-560557', 'J202112-563655','J202826-563818', 'J201838-572542', 'J201921-565155']
	nsrc=sources.shape[0]
	bleh=['J203455-565256', 'J203018-534949', 'J204553-565854','J205720-565940', 'J205319-552733', 'J204124-535646']



	for i in range (start,end):
			print(i)
			gc.collect()
			#check file
			src=sources[i,6]
			if src =='J202827-571807':
			#if 1+1==2:
			#if src  in badwise:
				#try:
				#if 1+1==2:
					ra_deg=sources[i,9]
					dec_deg=sources[i,10]
					#get source coords 
					coordstr=ra_deg+dec_deg
					coords = SkyCoord(coordstr,frame='fk5',unit=u.degree)  
					ra=coords.ra.hms
					dec=coords.dec.dms

					EMUrms=float(sources[i,36]) #mjy
					
					print("Making .png images for source {}".format(src)) 
					
					#racode=ra.split('.')[0].replace('h', '').replace('m', '')
					#deccode=dec.split('.')[0].replace('d', '').replace('m', '')

					racode="{:02d}".format(int(ra[0]))+"{:02d}".format(int(ra[1]))+"{:02d}".format(int(ra[2]))
					deccode="{:02d}".format(-1*int(dec[0]))+"{:02d}".format(-1*int(dec[1]))+"{:02d}".format(-1*int(dec[2]))
					DESstr=filedirectory+'DES/DESJ{}*{}*'.format(racode,deccode)
					#print(DESstr)
					#print(ra_deg,',',dec_deg)
					
					#load the DES data
					if src not in nodes:
						R,Rhead=fitsopen(glob.glob(DESstr+'_i.fits')[0])
						G,Ghead=fitsopen(glob.glob(DESstr+'_r.fits')[0])
						B,Bhead=fitsopen(glob.glob(DESstr+'_g.fits')[0])
						des_wcs=WCS(Rhead)



						img = np.zeros((R.shape[0], R.shape[1], 3), dtype=float)
						img[:,:,0] = img_scale.log(R, scale_min=5, scale_max=5000)
						img[:,:,1] = img_scale.log(G, scale_min=5, scale_max=5000)
						img[:,:,2] = img_scale.log(B, scale_min=5, scale_max=5000)

					#radio contours
					radio_im=filedirectory+'radio/{}_SB{}.fits'.format(src,str(SB))
					radio_label='EMU'
					radio_data,radio_header=fitsopen(radio_im)
					radio_wcs=WCS(radio_header)

					xp,yp=utils.skycoord_to_pixel(coords,radio_wcs)
					print(xp,yp)

					#WISE contours
					# find out which WISE image to load
					
					if src  not in badwise:
						WISEsep=coords.separation(WISEcoords)/u.degree
						wise_im=WISEfiles[np.argmin(WISEsep)]
						#
					else:
						#wise_label='WISE 3.4'
						wise_im=directory+'/WISE_DSS_cutouts/'+src+'_WISE_3.4.fits'

					wise_data,wise_header=fitsopen(wise_im)
					wise_wcs=WCS(wise_header)
					
					#wise_im=directory+'WISE_DSS_cutouts/{}_WISE_3.4.fits'.format(src)
					#wise_data,wise_header=fitsopen(wise_im)
					#wise_wcs=WCS(wise_header)

					#DSS_im=directory+'WISE_DSS_cutouts/{}_DSS.fits'.format(src)
					#DSS_data,DSS_header=fitsopen(DSS_im)
					#DSS_wcs=WCS(DSS_header)
					#print(np.nanmin(DSS_data))
					#print(np.nanmax(DSS_data))

					contourexps=np.arange(start=0,stop=32,step=2)
					contourmults=np.power(2,contourexps)
					basecont=3*EMUrms/1000
					radio_contours = [basecont * i for i in contourmults]

					plt.figure(dpi= dpi,figsize=(6,6))
					# plot the greyscale
					ax1=plt.subplot(131,projection=radio_wcs)		
					ax1.set_facecolor('k')		
					ax1.imshow(radio_data,origin='lower',cmap=magmacmap,norm=colors.LogNorm(vmin=basecont/3, vmax=np.nanmax(radio_data)/2))
					ax1.contour(radio_data,levels=radio_contours,colors=['green'])

					ax1.set_xlim(90,271)
					ax1.set_ylim(90,271)

					ax2=plt.subplot(132,projection=radio_wcs)
					ax2.set_facecolor('k')				
					#ax2.imshow(radio_data,origin='lower',cmap=infernocmap,norm=colors.LogNorm(vmin=basecont/3, vmax=np.nanmax(radio_data)/2))
					ax2.imshow(wise_data,transform=ax2.get_transform(wise_wcs),origin='lower',cmap=orangecmap,norm=colors.LogNorm(vmin=2.5,vmax=150))
					ax2.contour(radio_data,levels=radio_contours,colors=['green'])
					ax2.set_xlim(90,271)
					ax2.set_ylim(90,271)

					ax3=plt.subplot(133,projection=radio_wcs)
					ax3.set_facecolor('k')	
					#ax3.imshow(radio_data,origin='lower',cmap=infernocmap,norm=colors.LogNorm(vmin=basecont/3, vmax=np.nanmax(radio_data)/2))
					#ax3.imshow(DSS_data,origin='lower',transform=ax3.get_transform(DSS_wcs),cmap=greycmap,vmin=3000,vmax=17500)
					if src not in nodes:
						ax3.imshow(img,origin='lower',transform=ax3.get_transform(des_wcs))
					ax3.contour(radio_data,levels=radio_contours,colors=['green'])

					ax3.set_xlim(90,271)
					ax3.set_ylim(90,271)

					ax1.axis('off')
					ax2.axis('off')
					ax3.axis('off')

					#ax1.annotate(xy=(0.02,0.96),xycoords='figure fraction',text="Radio",c='white',ha="left")  
					#ax2.annotate(xy=(0.51,0.96),xycoords='figure fraction',text="IR",c='black',ha="left")  
					#ax3.annotate(xy=(0.02,0.47),xycoords='figure fraction',text="Optical",c='white',ha="left")  

					plt.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.99, wspace=0.01, hspace=0.01)
					plt.savefig(filedirectory+'png/{}_alpha_basicview.png'.format(src),dpi=dpi,bbox_inches='tight',pad_inches=0.0,transparent=True)

					### individual 

					plt.figure(dpi= dpi,figsize=(2,2))
					# plot the greyscale
					ax1=plt.subplot(111,projection=radio_wcs)
					ax1.set_facecolor('k')				
					ax1.imshow(radio_data,origin='lower',cmap=magmacmap,norm=colors.LogNorm(vmin=basecont/3, vmax=np.nanmax(radio_data)/2))
					ax1.contour(radio_data,levels=radio_contours,colors=['green'])
					ax1.set_xlim(90,271)
					ax1.set_ylim(90,271)
					ax1.axis('off')
					plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
					plt.savefig(filedirectory+'png/{}_alpha_6x6_radio.png'.format(src),dpi=dpi,bbox_inches='tight',pad_inches=0.0,transparent=True)

					plt.figure(dpi= dpi,figsize=(2,2))
					# plot the greyscale
					ax1=plt.subplot(111,projection=radio_wcs)	
					ax1.set_facecolor('k')			
					ax1.imshow(radio_data,origin='lower',cmap=infernocmap,norm=colors.LogNorm(vmin=basecont/3, vmax=np.nanmax(radio_data)/2))
					ax1.imshow(wise_data,transform=ax1.get_transform(wise_wcs),origin='lower',cmap=orangecmap,norm=colors.LogNorm(vmin=2.5,vmax=150))
					ax1.contour(radio_data,levels=radio_contours,colors=['green'])
					ax1.set_xlim(90,271)
					ax1.set_ylim(90,271)
					ax1.axis('off')
					plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
					plt.savefig(filedirectory+'png/{}_alpha_6x6_IR.png'.format(src),dpi=dpi,bbox_inches='tight',pad_inches=0.0,transparent=True)

					plt.figure(dpi= dpi,figsize=(2,2))
					# plot the greyscale
					ax1=plt.subplot(111,projection=radio_wcs)	
					ax1.set_facecolor('k')	
					#ax1.imshow(DSS_data,origin='lower',transform=ax1.get_transform(DSS_wcs),cmap=greycmap,vmin=3000,vmax=17500)
					if src not in nodes:		
						ax1.imshow(img,origin='lower',transform=ax1.get_transform(des_wcs))
					ax1.contour(radio_data,levels=radio_contours,colors=['green'])
					ax1.set_xlim(90,271)
					ax1.set_ylim(90,271)
					ax1.axis('off')
					plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
					plt.savefig(filedirectory+'png/{}_alpha_6x6_optical.png'.format(src),dpi=dpi,bbox_inches='tight',pad_inches=0.0,transparent=True)

					plt.figure(dpi= dpi,figsize=(2,2))
					# plot the greyscale
					ax1=plt.subplot(111,projection=radio_wcs)
					ax1.set_facecolor('k')				
					ax1.imshow(radio_data,origin='lower',cmap=infernocmap,norm=colors.LogNorm(vmin=basecont/3, vmax=np.nanmax(radio_data)/2))
					ax1.imshow(wise_data,transform=ax1.get_transform(wise_wcs),origin='lower',cmap=orangecmap,norm=colors.LogNorm(vmin=2.5,vmax=150))
					ax1.contour(radio_data,levels=radio_contours,colors=['green'])
					ax1.set_xlim(135,226)
					ax1.set_ylim(135,226)
					ax1.axis('off')
					plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
					plt.savefig(filedirectory+'png/{}_alpha_3x3_IR.png'.format(src),dpi=dpi,bbox_inches='tight',pad_inches=0.0,transparent=True)

					plt.figure(dpi= dpi,figsize=(2,2))
					# plot the greyscale
					ax1=plt.subplot(111,projection=radio_wcs)	
					ax1.set_facecolor('k')			
					ax1.imshow(radio_data,origin='lower',cmap=infernocmap,norm=colors.LogNorm(vmin=basecont/3, vmax=np.nanmax(radio_data)/2))
					ax1.imshow(wise_data,transform=ax1.get_transform(wise_wcs),origin='lower',cmap=orangecmap,norm=colors.LogNorm(vmin=2.5,vmax=150))
					ax1.contour(radio_data,levels=radio_contours,colors=['green'])
					ax1.set_xlim(0,361)
					ax1.set_ylim(0,361)
					ax1.axis('off')
					plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
					plt.savefig(filedirectory+'png/{}_alpha_12x12_IR.png'.format(src),dpi=dpi,bbox_inches='tight',pad_inches=0.0,transparent=True)

					### crosshair 
					
					plt.figure(dpi= dpi,figsize=(2,2))
					# plot the greyscale
					ax1=plt.subplot(111,projection=radio_wcs)
					ax1.set_facecolor('k')				
					ax1.imshow(radio_data,origin='lower',cmap=magmacmap,norm=colors.LogNorm(vmin=basecont/3, vmax=np.nanmax(radio_data)/2))
					ax1.contour(radio_data,levels=radio_contours,colors=['green'])
					ax1.set_xlim(90,271)
					ax1.set_ylim(90,271)
					ax1.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
					ax1.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
					ax1.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
					ax1.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')
					ax1.axis('off')
					plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
					plt.savefig(filedirectory+'png/{}_alpha_6x6_radio_cross.png'.format(src),dpi=dpi,bbox_inches='tight',pad_inches=0.0,transparent=True)

					plt.figure(dpi= dpi,figsize=(2,2))
					# plot the greyscale
					ax1=plt.subplot(111,projection=radio_wcs)	
					ax1.set_facecolor('k')			
					ax1.imshow(radio_data,origin='lower',cmap=infernocmap,norm=colors.LogNorm(vmin=basecont/3, vmax=np.nanmax(radio_data)/2))
					ax1.imshow(wise_data,transform=ax1.get_transform(wise_wcs),origin='lower',cmap=orangecmap,norm=colors.LogNorm(vmin=2.5,vmax=150))
					ax1.contour(radio_data,levels=radio_contours,colors=['green'])
					ax1.set_xlim(90,271)
					ax1.set_ylim(90,271)
					ax1.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
					ax1.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
					ax1.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
					ax1.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')
					ax1.axis('off')
					plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
					plt.savefig(filedirectory+'png/{}_alpha_6x6_IR_cross.png'.format(src),dpi=dpi,bbox_inches='tight',pad_inches=0.0,transparent=True)

					plt.figure(dpi= dpi,figsize=(2,2))
					# plot the greyscale
					ax1=plt.subplot(111,projection=radio_wcs)
					ax1.set_facecolor('k')
					#ax1.imshow(DSS_data,origin='lower',transform=ax1.get_transform(DSS_wcs),cmap=greycmap,vmin=3000,vmax=17500)
					if src not in nodes:			
						ax1.imshow(img,origin='lower',transform=ax1.get_transform(des_wcs))
					ax1.contour(radio_data,levels=radio_contours,colors=['green'])
					ax1.set_xlim(90,271)
					ax1.set_ylim(90,271)
					ax1.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
					ax1.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
					ax1.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
					ax1.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')
					ax1.axis('off')
					plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
					plt.savefig(filedirectory+'png/{}_alpha_6x6_optical_cross.png'.format(src),dpi=dpi,bbox_inches='tight',pad_inches=0.0,transparent=True)

					plt.figure(dpi= dpi,figsize=(2,2))
					# plot the greyscale
					ax1=plt.subplot(111,projection=radio_wcs)
					ax1.set_facecolor('k')			
					ax1.imshow(radio_data,origin='lower',cmap=infernocmap,norm=colors.LogNorm(vmin=basecont/3, vmax=np.nanmax(radio_data)/2))
					ax1.imshow(wise_data,transform=ax1.get_transform(wise_wcs),origin='lower',cmap=orangecmap,norm=colors.LogNorm(vmin=2.5,vmax=150))
					ax1.contour(radio_data,levels=radio_contours,colors=['green'])
					ax1.set_xlim(135,226)
					ax1.set_ylim(135,226)
					ax1.axis('off')
					ax1.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
					ax1.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
					ax1.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
					ax1.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')
					plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
					plt.savefig(filedirectory+'png/{}_alpha_3x3_IR_cross.png'.format(src),dpi=dpi,bbox_inches='tight',pad_inches=0.0,transparent=True)

					plt.figure(dpi= dpi,figsize=(2,2))
					# plot the greyscale
					ax1=plt.subplot(111,projection=radio_wcs)	
					ax1.set_facecolor('k')			
					ax1.imshow(radio_data,origin='lower',cmap=infernocmap,norm=colors.LogNorm(vmin=basecont/3, vmax=np.nanmax(radio_data)/2))
					ax1.imshow(wise_data,transform=ax1.get_transform(wise_wcs),origin='lower',cmap=orangecmap,norm=colors.LogNorm(vmin=2.5,vmax=150))
					ax1.contour(radio_data,levels=radio_contours,colors=['green'])
					ax1.set_xlim(0,361)
					ax1.set_ylim(0,361)
					ax1.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
					ax1.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
					ax1.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
					ax1.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')

					ax1.axis('off')
					plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
					plt.savefig(filedirectory+'png/{}_alpha_12x12_IR_cross.png'.format(src),dpi=dpi,bbox_inches='tight',pad_inches=0.0,transparent=True)
					
					#plt.clear()
					plt.close()
					plt.clf()


					if src not in nodes:
						del R
						del G
						del B
						del img
						del Rhead
						del Bhead
						del Ghead
					del radio_im
					del wise_im

					gc.collect()
					
				#except:
					#failedsources.append(src)
					#sys.exit()
			#else:
				#print("nop")
				#print("Already exists for src {}".format(src))

	print("Sources failed:")
	print(failedsources)

if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument('start', metavar='s', type=int, nargs='+')
	ap.add_argument('end', metavar='e', type=int, nargs='+')
	args = ap.parse_args()
	main(args)