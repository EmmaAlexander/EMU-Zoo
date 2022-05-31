import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
import astropy.io.fits as fits

selavyfile='/Volumes/NARNIA/hydra/AS101_Continuum_Component_Catalogue_9501_92.csv'
selavydata = np.genfromtxt(selavyfile, delimiter=',',dtype='str')

hydra_file='/Volumes/NARNIA/hydra/emu_pilot_sample_2x2deg.hydra.clump_catalogue.fits'
hydrahdu=fits.open()