import numpy as np 
import argparse

def main(args):
	# '/Volumes/TARDIS/Work/EMUzoo/DES_tiles.txt'
	DESfile=str(args.infile[0])
	#print(DESfile)
	DEStiles=np.loadtxt(DESfile,dtype='str')

	#print(DEStiles)

	#ra_min=2000 #format HHMM
	#ra_max=2200 #format HHMM
	#dec_max=-4800 #format DDMM
	#dec_min=-6000 #format DDMM

	ra_min=int(args.ra_min[0])
	ra_max=int(args.ra_max[0])
	dec_min=int(args.dec_min[0])
	dec_max=int(args.dec_max[0])

	outfile=str(args.outfile[0])

	DESkeep=[]
	
	for i in range(0,len(DEStiles)):
		tile=DEStiles[i]
		rastr=tile[3:7]
		decstr=tile[7:]

		ra=int(rastr)
		dec=int(decstr)

		if ra>=ra_min and ra<=ra_max and dec>=dec_min and dec<=dec_max:
			#print(tile)
			DESkeep.append(tile)

	output=np.asarray(DESkeep,dtype='str')

	np.savetxt(outfile,output,fmt='%s')
	
	#DES2159-5748_r4575p01_z_nobkg.fits.fz
	#then in file: 
	'''
	DECC1   =    -51.7470927970672 / DEC corner 1                                   
	RAC2    =     312.157362574794 / RA corner 2                                    
	DECC2   =    -51.7470927970672 / DEC corner 2                                   
	RAC3    =     312.166703090363 / RA corner 3                                    
	DECC3   =    -51.0166582939251 / DEC corner 3                                   
	RAC4    =     313.327840909637 / RA corner 4                                    
	DECC4   =    -51.0166582939251 / DEC corner 4 
	'''

if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument('infile', type=str, nargs='+')
	ap.add_argument('ra_min', type=int, nargs='+')
	ap.add_argument('ra_max', type=int, nargs='+')
	ap.add_argument('dec_min', type=int, nargs='+')
	ap.add_argument('dec_max', type=int, nargs='+')
	ap.add_argument('outfile', type=str, nargs='+')
	args = ap.parse_args()
	main(args)