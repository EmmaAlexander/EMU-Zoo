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

#read in the file with the sources
sourcefile = 'cat_below5000.txt'

SB=9351
directory ='/Volumes/TARDIS/Work/EMUzoo/alpha_test/'
filedirectory='/Volumes/TARDIS/Work/EMUzoo/alpha_test/'
sources=np.loadtxt(filedirectory+sourcefile,dtype='str')
dpi=684

def main():

	sources=np.loadtxt(filedirectory+sourcefile,dtype='str')
	baddes=np.loadtxt(filedirectory+'baddes.txt' ,dtype='str')
	badwise=np.loadtxt(filedirectory+'badwise.txt' ,delimiter=' ',dtype='str')
	nsrc=sources.shape[0]
	#bleh=['J203455-565256', 'J203018-534949', 'J204553-565854','J205720-565940', 'J205319-552733', 'J204124-535646']


	for i in range (0,nsrc):
		print(i)
		gc.collect()
		#check file
		src=sources[i,6]
		#if os.path.isdir(filedirectory+'png/{}'.format(src)) == False:
			#os.system("mkdir "+filedirectory+"png/{}".format(src))
		#os.system("mv "+filedirectory+"png/{}* ".format(src)+filedirectory+"png/{}/ ".format(src))

		#if src in bleh:
			#os.system("mv "+filedirectory+"png/{} ".format(src)+filedirectory+"png/redo/ ".format(src))
			#os.system("mkdir "+filedirectory+"png/{} ".format(src))
			#os.system("mv "+filedirectory+"png/{}* ".format(src)+filedirectory+"png/{}/ ".format(src))
		#else:
			#os.system("mv "+filedirectory+"png/{} ".format(src)+filedirectory+"png/final/ ".format(src))

		if os.path.isdir('/Users/emma/OneDrive/share/EMUzoo/alpha_test/png/{}'.format(src)) == False:
			print(src)

if __name__ == "__main__":
	main()