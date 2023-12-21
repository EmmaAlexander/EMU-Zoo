#run this in CASA

imdirectory='/Volumes/TARDIS/Work/askap/fullfields/'
outdirectory='/Volumes/TARDIS/Work/EMUzoo/alpha_test/radio/'
EMUSB='9351'
emu_filename= imdirectory +'image.i.SB'+EMUSB+'.cont.taylor.0.restored.fits'
casaimname=EMUSB+'_CASAim'

print("Opening full image")
#full_im=importfits(fitsimage=emu_filename,imagename=casaimname,overwrite=True)
ia.open(casaimname)

import numpy as np
import math
#sources=np.loadtxt('SB_9351_fullcat.txt',dtype='str')

catfile='/Volumes/TARDIS/Work/EMUzoo/alpha_test/cat_below5000.txt'
#cat = np.genfromtxt(catfile, delimiter=',',dtype='str')
cat = np.genfromtxt(catfile, delimiter=' ',dtype='str')


#headers=cat[0,:]
#sources=cat[1:,:]
sources=cat[:,:]

component_names=sources[:,7]
ra_deg=sources[:,10]
dec_deg=sources[:,11]

nsrc=sources.shape[0]

for i in range(0,nsrc):
	objectname=component_names[i]
	print('Working on source {}/{}: {}'.format(i+1,nsrc,objectname))

	ra=ra_deg[i]+'deg'
	dec=dec_deg[i]+'deg'

	position=[ra,dec]
	pixel = ia.topixel( position )['numeric']
	xPixel = str(int( round( pixel[ 0 ] ) ) )
	yPixel = str(int( round( pixel[ 1 ] ) ) )
	x1=int( round( pixel[ 0 ] ) - 180 )
	x2=int( round( pixel[ 0 ] ) + 180 )
	y1=int( round( pixel[ 1 ] ) - 180 )
	y2=int( round( pixel[ 1 ] ) + 180 )
	region='box[['+str(x1)+'pix,'+str(y1)+'pix],['+str(x2)+'pix,'+str(y2)+'pix]]'
	#print(region)
	#region='centerbox[['+xPixel+'pix,'+yPixel+'pix], [12arcmin, 12arcmin]]'
	outfile=outdirectory+objectname+'_SB'+EMUSB+'.fits'
	imsubimage(imagename=emu_filename,outfile='EMU_im_temp',region=region,overwrite=True,dropdeg=True)
	exportfits(imagename='EMU_im_temp',fitsimage=outfile,overwrite=True)
	os.system("rm -r EMU_im_temp")

ia.close()


