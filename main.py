# -- Import Libraries --
import astropy
from astropy.io import fits
from astropy.visualization import ZScaleInterval, ImageNormalize

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# -- Path to the File --
main_image_cutout = r"FITS\cutout_192.6009_-28.2675.fits" #future use input()?
main_image_psf = r"FITS\copsf_192.6009_-28.2675.fits"

# -- Open File & Plot (cutout) --
with fits.open(main_image_cutout) as hdu_list: #use with to safely close the file
    hdu_list.info() 
    header = hdu_list[0].header
    image_data = hdu_list[0].data #shape: (4, 256, 256)
    r_band = image_data[1] #red is index 1
    #print(r_band)

norm = ImageNormalize(r_band, interval=ZScaleInterval())
plt.figure(figsize=(8,8))
plt.imshow(r_band, origin="lower", norm=norm, cmap='gray')
plt.colorbar()
plt.title('r-band cutout for FITS File coord: RA=192.6009, Dec=-28.2675')
plt.show()

# -- Inspecting the PSF --
with fits.open(main_image_psf) as psf_hdu_list:
    psf_hdu_list.info()
    psf_data = psf_hdu_list[0].data

    print('PSF shape:', psf_data.shape) #the dimentions of the array -> NAXIS1 is 1D Arrays, NAXIS2 is 2D Arrays, NAXIS3 is 3D Arrays
    print('PSF sum:', psf_data.sum())

plt.figure(figsize=(5,5))
plt.imshow(psf_data, origin='lower', cmap='hot')
plt.colorbar
plt.title('PSF Model Graph')
plt.show()

