# -- Import Libraries --
import astropy
from astropy.io import fits
from astropy.visualization import ZScaleInterval, ImageNormalize
from astropy.table import Table
from astropy.wcs import WCS

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# -- Path to the different Files --
main_image_cutout = r"FITS\cutout_192.6009_-28.2675.fits" #future use input()?
main_image_psf = r"FITS\copsf_192.6009_-28.2675.fits"
neighbors_image = r"FITS\cat-ls-dr10.fits"

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
plt.colorbar()
plt.title('PSF Model Graph')
plt.show()

# -- Opening and Inspecting the Neighbors Image --
neighbors = Table.read(neighbors_image) #Table.read already handles the closing of the file so no with block needed

print(np.unique(neighbors['type'])) #1. check the different exitsting types (e.g., PSF, EXP, REX, etc.)

psf_neighbors_sources = neighbors[neighbors['type'] == b'PSF'] #2. Filter for PSF types ONLY
print(f"Total Neighbors: {len(neighbors)}")
print(f"Total PSF Sources: {len(psf_neighbors_sources)}")

# -- WCS Conversion --
#https://docs.astropy.org/en/latest/wcs/index.html
    #this conversion is needed because the neighbors catalog is in RA and Dec in degrees, but the image for the r_band is a 2D numpy array. 
    #The World Coordinate System transformation keywords are stored in the hearder
wcs = WCS(header).celestial #connect to main image cutout header to access those keywords
#the .celestial property automatically removes non-spacial dimentions, dropping the WCS down to a 2D system 

img_ny, img_nx = r_band.shape #Set the variables to the size of our image 

#setting corner pixels to wolrd coordinate system
bottom_left = wcs.pixel_to_world(0,0)
bottom_right = wcs.pixel_to_world(img_nx, 0)
top_left = wcs.pixel_to_world(0, img_ny)
top_right = wcs.pixel_to_world(img_nx, img_ny)

ra_min  = min(bottom_left.ra.deg,  bottom_right.ra.deg,  top_left.ra.deg,  top_right.ra.deg)
ra_max  = max(bottom_left.ra.deg,  bottom_right.ra.deg,  top_left.ra.deg,  top_right.ra.deg)
dec_min = min(bottom_left.dec.deg, bottom_right.dec.deg, top_left.dec.deg, top_right.dec.deg)
dec_max = max(bottom_left.dec.deg, bottom_right.dec.deg, top_left.dec.deg, top_right.dec.deg)
#this outputs sky coordinates: RA, Dec (Right Accrension, Declination)

print(f"RA covering: {ra_min:.4f} -> {ra_max:.4f} degrees") 
print(f"DEC covering: {dec_min:.4f} -> {dec_max:.4f} degrees") 


# Filtering the neighbor catalog to only sources in the specified RA and Dec coverage
filtered_neighbors_list = (
    (psf_neighbors_sources['ra'] >= ra_min) & (psf_neighbors_sources['ra'] <= ra_max) & (psf_neighbors_sources['dec'] >= dec_min) & (psf_neighbors_sources['dec'] <= dec_max)    
)
filtered_psf_neighbors = psf_neighbors_sources[filtered_neighbors_list]
print(f"PSF sources in this image's field: {len(filtered_psf_neighbors)}")

xs, ys = wcs.world_to_pixel_values(filtered_psf_neighbors['ra'], filtered_psf_neighbors['dec'])
#xs and ys are numpy arrays, one pixel coordinate per source
print(f"xs= {xs}")
print(f"ys = {ys}") #testing to check whether the arrays actually have something stored inside


# WCS Verification by OVERPLOTTING the different PSF sources on the r-band image
norm = ImageNormalize(r_band, interval=ZScaleInterval())
fig, ax = plt.subplots(figsize=(10, 10))
ax.imshow(r_band, origin='lower', norm=norm, cmap='gray')
ax.scatter(xs, ys, s=120, facecolors='none', edgecolors='lime', linewidths=1.5)
ax.set_title(f'{len(filtered_psf_neighbors)} PSF neighbors overlaid on r-band image')
plt.show()