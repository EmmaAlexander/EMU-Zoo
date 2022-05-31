#import os
#from glob import iglob
import glob, os

rootdir_glob = '/Volumes/TARDIS/Work/EMUzoo/alpha_test/DES/**/*'

#file_list = [f for f in iglob(rootdir_glob, recursive=True) if os.path.isfile(f)]

#for f in file_list:
	#print(f)

for filename in glob.iglob(rootdir_glob, recursive=True):
	if os.path.isfile(filename) and '.fits' in filename: # filter dirs
		#print(filename)
		os.system("mv {} /Volumes/TARDIS/Work/EMUzoo/alpha_test/DES/".format(filename))
		#os.system("mv {} ../../".format(filename))