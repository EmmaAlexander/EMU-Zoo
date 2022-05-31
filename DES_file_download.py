import numpy as np 
import argparse
import os

def main(args):
	# '/Volumes/TARDIS/Work/EMUzoo/DES_tiles.txt'
	DESfile=str(args.infile[0])
	dataloc=str(args.dataloc[0])
	downloader=str(args.dataloc[0])

	DEStiles=np.loadtxt(DESfile,dtype='str')

	rooturl='https://desdr-server.ncsa.illinois.edu/despublic/dr2_tiles/'
	bands=['g'] #choose subset of ['Y','g','i','r','z']

	for i in range(0,len(DEStiles)):
		tile=DEStiles[i]
		for b in bands:
			url=rooturl+'{}/{}_r4575p01_{}.fits.fz'.format(tile,tile,b)
			#imagename='{}_r4575p01_{}.fits.fz'.format(tile,b)
			curlstr="curl --header 'Host: desdr-server.ncsa.illinois.edu' --header 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/101.0.4951.64 Safari/537.36' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9' --header 'Accept-Language: en-US,en;q=0.9,en-GB;q=0.8' --header 'Referer: https://desdr-server.ncsa.illinois.edu/despublic/dr2_tiles/{}' 'https://desdr-server.ncsa.illinois.edu/despublic/dr2_tiles/{}/{}_r4575p01_{}.fits.fz' -L -o '{}_r4575p01_{}.fits.fz'".format(tile,tile,tile,b,tile,b)
			os.system(curlstr)


if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument('infile', type=str, nargs='+')
	ap.add_argument('dataloc', type=str, nargs='+')
	args = ap.parse_args()
	main(args)

