#run this in CASA

imdirectory='/Volumes/TARDIS/Work/askap/fullfields/'
outdirectory='/Users/emma/Dropbox/Public/EMUZoo/'

import numpy as np
sources=np.loadtxt('emu-zoo_cutouts_sourcelist.txt',dtype='str')
for i in range(0,sources.shape[0]):
	objectname=sources[i,0]
	sourcecentre=sources[i,1]	
	EMUSB=sources[i,2]

	print(objectname)
	region='centerbox[['+sourcecentre+'], [4arcmin, 4arcmin]]'
	outfile=outdirectory+objectname+'_EMU.fits'
	emu_filename= imdirectory +'image.i.SB'+EMUSB+'.cont.taylor.0.restored.fits'
	print(emu_filename)
	imsubimage(imagename=emu_filename,outfile='EMU_im_temp',region=region,overwrite=True,dropdeg=True)
	exportfits(imagename='EMU_im_temp',fitsimage=outfile,overwrite=True)
	os.system("rm -r EMU_im_temp")



