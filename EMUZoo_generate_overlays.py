# EMUZoo cutout generation code
# Author: Emma Alexander (emma.alexander@gmail.com)
# https://github.com/EmmaAlexander/EMU-Zoo

# NOTE: currently limited to EMU Pilot I survey fields area. Will needs tweaks for main survey data. 
# Also note: copied over from a jupyter notebook. Should hopefully still work.
# Original code (and data) can be found on AusSRC:
# scripts at /mnt/shared/home/ealexander/
# data at /mnt/shared/des

#-----------------------------------------------------------
#imports
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs import utils
import numpy as np 
from astroquery.skyview import SkyView
import astropy.units as u
from astropy.coordinates import SkyCoord
import glob, os, sys
from astropy.nddata import Cutout2D
import math
import cmasher as cmr
from astropy.visualization import lupton_rgb
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
from matplotlib.patches import Ellipse

#-----------------------------------------------------------
def ashinh_scale(array,zeropoint=0,scale=1):
    scaled=(array-zeropoint)*np.arcsinh(array*scale)/(array*np.arcsinh(scale))
    return scaled

def percentile(array,percent):
    val=np.percentile(array[np.isfinite(array)],percent)
    return val
#-----------------------------------------------------------

# Hard coded variables (data locations etc)

dataloc='/mnt/shared/des/askap_data/' #askap data location
WISEfiles=glob.glob('/mnt/shared/des/WISEtiles/*.fits') #WISE raw data tiles location

SB='9351'

#hard coded image names due to inconsistent naming conventions but could be streamlined
image='image.i.SB'+SB+'.cont.taylor.0.restored.fits'

#note that there are two catalogue types: island and component
#cat=dataloc+'catalogues/AS101_Continuum_Island_Catalogue_'+SB+'.csv'
cat=dataloc+'catalogues/AS101_Continuum_Component_Catalogue_9351_78.csv'
island=False #change to true if using island catalogue

# folder for cutout data output if you want to save radio fits files
# not currently implemented due to data storage issues
radiooutloc='/mnt/shared/des/radio_cutouts/SB'+SB
if os.path.isdir(radiooutloc) ==False:
    os.system('mkdir '+radiooutloc+SB)

overlayloc='/mnt/shared/des/radio_cutouts/SB'+SB
if os.path.isdir(overlayloc) ==False:
    os.system('mkdir '+overlayloc+SB)
overlay_suffix='.png'
#you can add to the above if you need different versions 
#e.g. change to '_v2.png'

arcmins=12. #maximum size of cutout in armins
# disclaimer: things may break if you change this...
# number of pixel dimensions
npix_edge=int(15*arcmins)
deg_edge=arcmins*0.025/3.

#-----------------------------------------------------------
# set up some plot parameters
# note that many of these are no longer used but kept for posterity

dpi=300
plt.rc('font', size=0.5)  
plt.rcParams.update({'lines.linewidth':0.8})
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
viridis=plt.cm.viridis
viridis.set_bad('black',1)
gist_heat=plt.cm.gist_heat
gist_heat.set_bad('black',1)
greencmap=plt.cm.Greens_r
greencmap.set_bad('black',1)

twl_blue = cmr.get_sub_cmap('twilight_shifted', 0, 0.5)
twl_red = cmr.get_sub_cmap('twilight', 0.5, 1)
twl_red.set_bad(twl_red(0),1)
twl_red_r = cmr.get_sub_cmap('twilight_shifted', 0.5, 1)
twl_blue_r = cmr.get_sub_cmap('twilight',0,0.5)

#-----------------------------------------------------------

#read in full ASKAP image 
hdu=fits.open(dataloc+image)
image=hdu[0].data.squeeze()
header=hdu[0].header
wcs= WCS(hdu[0].header).celestial
hdu.close()
#some useful parameters that don't get used here
pixscale=header['CDELT2']
bmaj_pix=header['BMAJ']/header['CDELT2']
bmin_pix=header['BMIN']/header['CDELT2']
bpa=header['BPA']

# read in catalogue
catalogue=np.genfromtxt(cat,dtype='str',delimiter=',')
headers=catalogue[0,:]
data=catalogue[1:,:]
#below lines for if you'd like to check catalogue headers
#for i in range(0,45):
    #print(i, headers[i],data[0,i])

#sort data by major axis (optional, but data_sorted is array used going forward
#note that 19 needs changing to 31 for an island catalogue because they have different columns
#TODO: use header name rather than column index (pandas...?)
data_sorted = data[data[:,19].argsort()[::-1]] #componenet cat, major axis

#-----------------------------------------------------------

#get the coordintes of the WISE images to know which to load and use
WISEcoords=[]
for file in WISEfiles:
    filename=file.split('/')[-1] #get the actual file name
    ra=float(filename[0:4])/10
    dec=-1*float(filename[5:8])/10
    WISEcoords.append(SkyCoord(ra,dec,unit=u.degree,frame='fk5'))
WISEcoords=SkyCoord(np.asarray(WISEcoords))

WISEtiles=np.genfromtxt('WISE_tiles.txt',dtype='str') #this file should be in the same location as the code
DEStiles=np.genfromtxt('/DES_tiles_dims.txt',dtype='str') #this file should be in the same location as the code

#-----------------------------------------------------------

#MAIN BIT

#loop over sources in the list

for i in range(0,len(data_sorted)):
    if i % 100 == 0:
        print(i)
        #or use a better tracking method idk

    if island==True:
    	#get values from catalogue file
    	# not all of these are used but I have left them in
        src = data_sorted[i,6]    
        n_components = data_sorted[i,7]
        ra_hms_cont = data_sorted[i,8]
        dec_dms_cont = data_sorted[i,9]
        ra_deg_cont = data_sorted[i,10]
        dec_deg_cont = data_sorted[i,11]
        maj_axis = float(data_sorted[i,13])
        min_axis = float(data_sorted[i,14])
        pos_ang = float(data_sorted[i,15])
        flux_int = float(data_sorted[i,16])
        flux_int_err = float(data_sorted[i,16])
        flux_peak = float(data_sorted[i,18])
        background_noise = float(data_sorted[i,20])
        n_pix  = float(data_sorted[i,30])
        solid_angle = float(data_sorted[i,31])
        x_cen = float(data_sorted[i,35])
        y_cen = float(data_sorted[i,36])
    else:
    	#component catalogue being used
    
        src=data_sorted[i,7]
        ra_hms_cont	=data_sorted[i,8]
        dec_dms_cont=data_sorted[i,9]
        ra_deg_cont	=data_sorted[i,10]
        dec_deg_cont=data_sorted[i,11]      
        background_noise =float(data_sorted[i,37])    

        #re-iterating todo: change to header names rather than index (likely to break)	

    filename=overlayloc+src+overlay_suffix
    filename_cross=overlayloc+src+'_cross_'+overlay_suffix

    # check to see if file already exists for this source
    # aka you shouldn't be able to overwrite existing files unless you tweak this!
    if os.path.isfile(filename) == False:
    	#print(i,src)
        coords=SkyCoord(ra_deg_cont,dec_deg_cont,frame='fk5',unit=u.degree)   

        ra_max=float((coords.ra/u.degree)+deg_edge) 
        ra_min=float((coords.ra/u.degree)-deg_edge)
        dec_max=float((coords.dec/u.degree)+deg_edge) 
        dec_min=float((coords.dec/u.degree)-deg_edge)
        
        x_cen,y_cen= wcs.world_to_pixel(coords)
        
        xmin=int(x_cen-npix_edge)
        xmax=int(x_cen+npix_edge)
        ymin=int(y_cen-npix_edge)
        ymax=int(y_cen+npix_edge)
        
        #print(ra_max,ra_min,dec_max,dec_min)

        # find which DES and WISE tiles to use
        indices=np.where( (DEStiles[:,2].astype(float)>=ra_min) & (DEStiles[:,1].astype(float)<=ra_max) & (DEStiles[:,3].astype(float)>=dec_min) & (DEStiles[:,4].astype(float)<=dec_max))
        DES_tiles_to_use=DEStiles[indices,0]
        DES_tiles_to_use_coords=DEStiles[indices,:]

        wiseindices=np.where( (WISEtiles[:,6].astype(float)>=ra_min) & (WISEtiles[:,8].astype(float)<=ra_max) & (WISEtiles[:,13].astype(float)>=dec_min) & (WISEtiles[:,9].astype(float)<=dec_max))
        wise_tiles_to_use=WISEtiles[wiseindices,16]

        
        # get radio cutout
        radio_cutout = Cutout2D(image, position=(x_cen,y_cen), size=(2*npix_edge), wcs=wcs, mode='trim')
        # calculate contours
        # TAKE NOTE THIS IS WHERE TO TWEAK RADIO CONTOUR LEVELS
        contourexps=np.arange(start=0,stop=32,step=0.5) 
        #use step=1 for contours doubling each time, 0.5 for a factor of root 2 etc...
        contourmults=np.power(2,contourexps)
        #basecont=3.*background_noise/1000.
        #OR
        basecont=0.00012
        radio_contours = [basecont * i for i in contourmults]
        
        radio_max=np.nanmax(radio_cutout.data)
        nconts=np.nanmax(np.flatnonzero(radio_contours < radio_max))
        
        contcolors=[]
        for c in range(0,nconts+1):
            contcolors.append(greycmap(0.5+(0.5*c/(nconts+1))))  


        R_list=[]
        G_list=[]
        B_list=[]
        
        for j in range(0,len(DES_tiles_to_use[0])):
            #for j in ind:
            Rhdu=fits.open(glob.glob('/mnt/shared/des/DEStiles/'+DES_tiles_to_use[0][j]+'*_i.fits*')[0])
            R=Rhdu[1].data
            des_wcs=WCS(Rhdu[1].header)
            Rhdu.close()

            Ghdu=fits.open(glob.glob('/mnt/shared/des/DEStiles/'+DES_tiles_to_use[0][j]+'*_r.fits*')[0])
            G=Ghdu[1].data
            Ghdu.close()

            Bhdu=fits.open(glob.glob('/mnt/shared/des/DEStiles/'+DES_tiles_to_use[0][j]+'*_g.fits*')[0])
            B=Bhdu[1].data
            Bhdu.close()

            R_cutout=Cutout2D(R,position=coords,size=1.05*arcmins*u.arcmin,wcs=des_wcs,mode='trim')
            G_cutout=Cutout2D(G,position=coords,size=1.05*arcmins*u.arcmin,wcs=des_wcs,mode='trim')
            B_cutout=Cutout2D(B,position=coords,size=1.05*arcmins*u.arcmin,wcs=des_wcs,mode='trim')

            R_hdu=fits.PrimaryHDU(data=R_cutout.data, header=R_cutout.wcs.to_header())
            G_hdu=fits.PrimaryHDU(data=G_cutout.data, header=G_cutout.wcs.to_header())
            B_hdu=fits.PrimaryHDU(data=B_cutout.data, header=B_cutout.wcs.to_header())

            R_list.append(R_hdu)
            G_list.append(G_hdu)
            B_list.append(B_hdu)  

        if len(R_list)==0:
            print("something wrong")
        elif len(R_list)==0:
            #only one image so no need to mosaic
            R=R_list[0].data
            G=G_list[0].data
            B=B_list[0].data
            des_wcs=WCS(R_list[0].header)
        else:
            #need to combine them
            des_wcs, shape_out = find_optimal_celestial_wcs(R_list)
            R, footprint_R = reproject_and_coadd(R_list,des_wcs,shape_out=shape_out,reproject_function=reproject_interp)
            G, footprint_G = reproject_and_coadd(G_list,des_wcs,shape_out=shape_out,reproject_function=reproject_interp)
            B, footprint_B = reproject_and_coadd(B_list,des_wcs,shape_out=shape_out,reproject_function=reproject_interp)
            
        img=lupton_rgb.make_lupton_rgb(R,G,B,Q=10,stretch=50,minimum=1)

        my_dpi = 100
        sub_font_size = 10
        fig = plt.figure(constrained_layout=False,figsize=(1024/my_dpi, 1024/my_dpi),dpi=my_dpi*0.55)
        # Set figure background as white
        fig.patch.set_facecolor('w')        

        ax1=plt.subplot(331,projection=radio_cutout.wcs,fc='grey')            
        ax2=plt.subplot(332,projection=radio_cutout.wcs,fc='grey') 
        ax3=plt.subplot(333,projection=radio_cutout.wcs,fc='grey')    
        ax4=plt.subplot(334,projection=radio_cutout.wcs,fc='grey')    
        ax5=plt.subplot(335,projection=radio_cutout.wcs,fc='grey')    
        ax6=plt.subplot(336,projection=radio_cutout.wcs,fc='grey')    
        ax7=plt.subplot(337,projection=radio_cutout.wcs,fc='grey')    
        ax8=plt.subplot(338,projection=radio_cutout.wcs,fc='grey')    
        ax9=plt.subplot(339,projection=radio_cutout.wcs,fc='grey')    

        ax1.imshow(radio_cutout.data,origin='lower',cmap=magmacmap,norm=colors.LogNorm(vmin=basecont/5, vmax=radio_max))
        ax1.contour(radio_cutout.data,levels=radio_contours,colors='grey')  
        ax2.imshow(radio_cutout.data,origin='lower',cmap=magmacmap,norm=colors.LogNorm(vmin=basecont/5, vmax=radio_max))
        ax2.contour(radio_cutout.data,levels=radio_contours,colors='grey')
        ax3.imshow(radio_cutout.data,origin='lower',cmap=magmacmap,norm=colors.LogNorm(vmin=basecont/5, vmax=radio_max))
        ax3.contour(radio_cutout.data,levels=radio_contours,colors='grey')

        ax4.imshow(img,transform=ax4.get_transform(des_wcs),origin='lower') 
        ax5.imshow(img,transform=ax5.get_transform(des_wcs),origin='lower') 
        ax6.imshow(img,transform=ax6.get_transform(des_wcs),origin='lower')
        ax4.contour(radio_cutout.data,levels=radio_contours,colors=contcolors)
        ax5.contour(radio_cutout.data,levels=radio_contours,colors=contcolors)
        ax6.contour(radio_cutout.data,levels=radio_contours,colors=contcolors)
        
        try:
            wise_list=[]
            for k in range(0,len(wise_tiles_to_use[0])):
                #get wise cutout 
                wise_im='/mnt/shared/des/WISEtiles/'+str(wise_tiles_to_use[0][k])+'-w1-int-3.fits'
                wise_hdu=fits.open(wise_im)
                wise_data= wise_hdu[0].data
                wise_wcs= WCS(wise_hdu[0].header)
                wise_hdu.close()
                wise_cutout=Cutout2D(wise_data,position=coords,size=1.05*arcmins*u.arcmin,wcs=wise_wcs,mode='trim')
                wise_hdu=fits.PrimaryHDU(data=wise_cutout.data, header=wise_cutout.wcs.to_header())               
                wise_list.append(wise_hdu)

            if len(wise_list)==0:
                print("something wrong")
            elif len(wise_list)==1:
                print("only one wise")
                #only one image so no need to mosaic
                wise_data=wise_list[0].data
                wise_wcs=WCS(wise_list[0].header)

            else:
                print("combining wise")
                #need to combine them
                wise_wcs, wise_shape_out = find_optimal_celestial_wcs(wise_list)
                wise_data, footprint_wise = reproject_and_coadd(wise_list,wise_wcs,shape_out=wise_shape_out,reproject_function=reproject_interp)

            #wise_cutout_hdu=fits.PrimaryHDU(data=wise_data, header=wise_wcs.to_header()) 
            #wise_cutout_hdu.writeto(filename)
            
            ax7.imshow(ashinh_scale(wise_data,zeropoint=2,scale=100),transform=ax7.get_transform(wise_wcs),origin='lower',cmap=gist_heat)
            ax8.imshow(ashinh_scale(wise_data,zeropoint=2,scale=100),transform=ax8.get_transform(wise_wcs),origin='lower',cmap=gist_heat)
            ax9.imshow(ashinh_scale(wise_data,zeropoint=2,scale=100),transform=ax9.get_transform(wise_wcs),origin='lower',cmap=gist_heat)
            
        except:
            print("wise failed")
            for k in range(0,len(wise_tiles_to_use[0])):
                #get wise cutout 
                wise_im='/mnt/shared/des/WISEtiles/'+str(wise_tiles_to_use[0][k])+'-w1-int-3.fits'
                wise_hdu=fits.open(wise_im)
                wise_data= wise_hdu[0].data
                wise_wcs= WCS(wise_hdu[0].header)
                wise_hdu.close()
                ax7.imshow(ashinh_scale(wise_data,zeropoint=2,scale=100),transform=ax7.get_transform(wise_wcs),origin='lower',cmap=gist_heat)
                ax8.imshow(ashinh_scale(wise_data,zeropoint=2,scale=100),transform=ax8.get_transform(wise_wcs),origin='lower',cmap=gist_heat)
                ax9.imshow(ashinh_scale(wise_data,zeropoint=2,scale=100),transform=ax9.get_transform(wise_wcs),origin='lower',cmap=gist_heat)

        ax7.contour(radio_cutout.data,levels=radio_contours,colors=contcolors)
        ax8.contour(radio_cutout.data,levels=radio_contours,colors=contcolors)
        ax9.contour(radio_cutout.data,levels=radio_contours,colors=contcolors)


        ax1.set_xlim(0.75*npix_edge,(1.25*npix_edge)-1)
        ax1.set_ylim(0.75*npix_edge,(1.25*npix_edge)-1)
        ax2.set_xlim(0.5*npix_edge,(1.5*npix_edge)-1)
        ax2.set_ylim(0.5*npix_edge,(1.5*npix_edge)-1)
        ax3.set_xlim(0,(2*npix_edge)-1)
        ax3.set_ylim(0,(2*npix_edge)-1)
        
        ax4.set_xlim(0.75*npix_edge,(1.25*npix_edge)-1)
        ax4.set_ylim(0.75*npix_edge,(1.25*npix_edge)-1)
        ax5.set_xlim(0.5*npix_edge,(1.5*npix_edge)-1)
        ax5.set_ylim(0.5*npix_edge,(1.5*npix_edge)-1)
        ax6.set_xlim(0,(2*npix_edge)-1)
        ax6.set_ylim(0,(2*npix_edge)-1)
        
        ax7.set_xlim(0.75*npix_edge,(1.25*npix_edge)-1)
        ax7.set_ylim(0.75*npix_edge,(1.25*npix_edge)-1)
        ax8.set_xlim(0.5*npix_edge,(1.5*npix_edge)-1)
        ax8.set_ylim(0.5*npix_edge,(1.5*npix_edge)-1)
        ax9.set_xlim(0,(2*npix_edge)-1)
        ax9.set_ylim(0,(2*npix_edge)-1)


        ax1.axis('off')
        ax2.axis('off')
        ax3.axis('off')
        ax4.axis('off')
        ax5.axis('off')
        ax6.axis('off')
        ax7.axis('off')
        ax8.axis('off')
        ax9.axis('off')

        plt.rc('font', size=10)    
        plt.annotate('Radio',xy=(0.08,0.75),xycoords='figure fraction',ha='center',va='center',rotation=90)
        plt.annotate('Optical',xy=(0.08,0.45),xycoords='figure fraction',ha='center',va='center',rotation=90)
        plt.annotate('Infrared',xy=(0.08,0.15),xycoords='figure fraction',ha='center',va='center',rotation=90)
        plt.annotate('Zoomed in',xy=(0.25,0.91),xycoords='figure fraction',ha='center')
        plt.annotate('Default',xy=(0.55,0.91),xycoords='figure fraction',ha='center')
        plt.annotate('Zoomed out',xy=(0.85,0.91),xycoords='figure fraction',ha='center')
        plt.subplots_adjust(left=0.1, bottom=0.01, right=0.99, top=0.9, hspace=0,wspace=0.02) # control space between figure and whitespace

        plt.savefig(filename)
        
        xp,yp=utils.skycoord_to_pixel(coords,radio_cutout.wcs)

        #plt.rcParams.update({'lines.linewidth':1.2})

        ax1.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
        ax1.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
        ax1.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
        ax1.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')

        ax2.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
        ax2.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
        ax2.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
        ax2.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')

        ax3.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
        ax3.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
        ax3.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
        ax3.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')

        ax4.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
        ax4.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
        ax4.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
        ax4.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')

        ax5.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
        ax5.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
        ax5.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
        ax5.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')

        ax6.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
        ax6.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
        ax6.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
        ax6.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')
        
        ax7.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
        ax7.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
        ax7.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
        ax7.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')
        
        ax8.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
        ax8.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
        ax8.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
        ax8.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')
        
        ax9.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
        ax9.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
        ax9.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
        ax9.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')

        plt.savefig(filename_cross)
        plt.show()


