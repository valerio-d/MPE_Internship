# -- Import Libraries --
import astropy
from astropy.io import fits
from astropy.visualization import ZScaleInterval, ImageNormalize
from astropy.table import Table
from astropy.wcs import WCS

import numpy as np

from scipy.ndimage import shift, zoom

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# -- Path to the different Files --
main_image_cutout = r"FITS\cutout_192.6009_-28.2675.fits" #future use input()?
main_image_psf = r"FITS\copsf_192.6009_-28.2675.fits"
neighbors_image = r"FITS\cat-ls-dr10.fits"

#define the object we are interested in
main_object_ra = 192.6009
main_object_dec = -28.2675

# -- Open File & Plot (cutout) --
with fits.open(main_image_cutout) as hdu_list: #use with to safely close the file
    hdu_list.info() 
    header = hdu_list[0].header
    #print(header)
    image_data = hdu_list[0].data #shape: (4, 256, 256)
    r_band = image_data[1] #red is index 1 #this is fluxes
    

norm = ImageNormalize(r_band, interval=ZScaleInterval())
plt.figure(figsize=(8,8))
plt.imshow(r_band, origin="lower", norm=norm, cmap='gray')
plt.colorbar()
plt.title('r-band cutout for FITS File coord: RA=192.6009, Dec=-28.2675')
plt.show()

# -- Inspecting the PSF --
with fits.open(main_image_psf) as psf_hdu_list:
    psf_hdu_list.info()
    psf_data = psf_hdu_list[0].data * 644 #found this multiplication factor by dividing brightest point on og image by the brightest point on psf image

    print('PSF shape:', psf_data.shape) #the dimentions of the array -> NAXIS1 is 1D Arrays, NAXIS2 is 2D Arrays, NAXIS3 is 3D Arrays
    print('PSF sum:', psf_data.sum()) # A PSF that sums to 1.0 means: "if I stamp this at a location and scale it by flux F, the total flux contributed to the image equals F."

plt.figure(figsize=(5,5))
plt.imshow(psf_data, origin='lower', cmap='hot')
plt.colorbar()
plt.title('PSF Model Graph')
plt.show()

# -- Opening and Inspecting the Neighbors Image --
neighbors = Table.read(neighbors_image) #Table.read already handles the closing of the file so no with block needed

print(np.unique(neighbors['type'])) #1. check the different exitsting types (e.g., PSF, EXP, REX, etc.)

psf_neighbors_sources = neighbors #2. Filter for PSF types ONLY
#psf_neighbors_sources = neighbors[neighbors['type'] == b'PSF'] #2. Filter for PSF types ONLY
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
#these two values are the x and y coordinates for each of the filtered psf neighbors

# WCS Verification by OVERPLOTTING the different PSF sources on the r-band image
norm = ImageNormalize(r_band, interval=ZScaleInterval())
fig, ax = plt.subplots(figsize=(10, 10))
ax.imshow(r_band, origin='lower', norm=norm, cmap='bone')
ax.scatter(xs, ys, s=120, facecolors='none', edgecolors='lime', linewidths=1.5)
ax.set_title(f'{len(filtered_psf_neighbors)} PSF neighbors overlaid on r-band image')
#plt.show()

#print(neighbors.colnames)  # see all columns
fluxes = filtered_psf_neighbors['flux_r'] 
print(fluxes) #print all the aperture fluxes in the r band


# -- PSF SYNTHETIC IMAGE CREATION AND PLACEMENT FUNCTION --
#find and exclude the main object 
main_object_distances = np.sqrt((filtered_psf_neighbors['ra'] - main_object_ra)**2 + (filtered_psf_neighbors['dec'] - main_object_dec)**2) #the distance from every source to the object we're interested in
main_obj_index = np.argmin(main_object_distances) #find the index of the main object by finding the smallest distance

psf_object_KEEP = np.arange(len(filtered_psf_neighbors)) != main_obj_index

xs_without_target = xs[psf_object_KEEP]
ys_without_target = ys[psf_object_KEEP]
fluxes_without_target = filtered_psf_neighbors['flux_r'][psf_object_KEEP]


def create_synthetic_image(image_shape, psf, sources_x, sources_y, fluxes): #takes the real size of the image, the psf stamp, the xs and xy arrays from the WCS conversion, and the fluxes (how bright each star is)
    synthetic = np.zeros(image_shape, dtype=np.float64) # creating a BLANK CANVAS, a 256 x 256 grid
    psf_h, psf_w = psf.shape
    psf_cy, psf_cx = psf_h // 2, psf_w // 2 #identifying the center pixel of the PSF Image (height and width / 2)

    img_h, img_w = image_shape

    for x, y, flux in zip(sources_x, sources_y, fluxes): #zip allows to combine the 3 lists in triplets (x1, y1. flux1), (x2, y2, flux2), etc. #this means that we are going through each one of the identified objects

        ix, iy = int(round(x)), int(round(y)) #round the pixel position to the nearest pixel
        print(ix, iy, flux)

        # calculate where the PSF stamp lands in the image, centering it on the star position
        #Say the star is at pixel x=135, and the sticker is 25 pixels wide with center at index 12. If you want the sticker's center pixel to land exactly on x=135, the sticker's left edge (x0) must start at 135 - 12 = 123. Then the sticker's right edge (x1) is 123 + 25 = 148. So the sticker occupies columns 123 through 148 in the big image, and its middle column (135) lines up exactly with the star.
        x0 = ix - psf_cx
        y0 = iy - psf_cy
        x1 = x0 + psf_w
        y1 = y0 + psf_h
        
        # Clip to image boundaries
        img_x0 = max(x0, 0);  img_x1 = min(x1, img_w) #if the column < 0, set it = 0
        img_y0 = max(y0, 0);  img_y1 = min(y1, img_h) #the right edge of the psf cannot go past the image's width
        psf_x0 = img_x0 - x0; psf_x1 = psf_x0 + (img_x1 - img_x0)
        psf_y0 = img_y0 - y0; psf_y1 = psf_y0 + (img_y1 - img_y0)
        
        if img_x1 <= img_x0 or img_y1 <= img_y0:
            continue  # source is fully outside image, skip
        
        synthetic[img_y0:img_y1, img_x0:img_x1] += flux * psf[psf_y0:psf_y1, psf_x0:psf_x1] #equation to create the psf
    
    return synthetic

synthetic_image = create_synthetic_image(r_band.shape, psf_data, xs_without_target, ys_without_target, fluxes_without_target) #create image WITHOUT the main object
synthetic_image_data = synthetic_image.data
#assign the arguments: image shape (the r-band shape), the psf data, the xs and ys coordinates for each psf, and the fluxes for each psf

# -- PLOT THE SYNTHETIC IMAGE --
fig, axes = plt.subplots(1, 2, figsize=(16, 8))

norm_og = ImageNormalize(r_band, interval=ZScaleInterval())
norm_syn  = ImageNormalize(synthetic_image, interval=ZScaleInterval())

axes[0].imshow(r_band, origin='lower', norm=norm_og, cmap='gray')
axes[0].set_title('original r-band image')

axes[1].imshow(synthetic_image, origin='lower', norm=norm_syn, cmap='gray')
axes[1].set_title('synthetic PSF image')

plt.tight_layout()
#plt.show()

# plotting in LOGARTITHMIC scale
fig, axes = plt.subplots(1, 2, figsize=(16, 8))
LogNorm_minimum_value = 0.1 #minumum value to handle fully black
LogNorm_maximum_value = r_band.max() #minumum value to handle fully black

logarithmic_image1 = axes[0].imshow(r_band, origin="lower", norm=LogNorm(vmin=LogNorm_minimum_value, vmax=LogNorm_maximum_value), cmap="gray")
axes[0].set_title('original r-band image')
logarithmic_image2 = axes[1].imshow(synthetic_image, origin="lower", norm=LogNorm(vmin=LogNorm_minimum_value, vmax=LogNorm_maximum_value), cmap="gray")
axes[1].set_title('synthetic PSF image')


low_tick = np.percentile(r_band[r_band > 0], 20) #ticks for image 1
medium_tick = np.percentile(r_band[r_band > 0], 50)#calculate percentiles using only positive pixels to avoid math issues
high_tick = np.percentile(r_band[r_band > 0], 90)
low_tick2 = np.percentile(synthetic_image[synthetic_image > 0], 20) #ticks for image 2
medium_tick2 = np.percentile(synthetic_image[synthetic_image > 0], 50)
high_tick2 = np.percentile(synthetic_image[synthetic_image > 0], 90)


cbar_1 = fig.colorbar(logarithmic_image1, ax=axes[0], ticks=[low_tick, medium_tick, high_tick]) # Create specific colorbars for each of the images
cbar_2 = fig.colorbar(logarithmic_image2, ax=axes[1], ticks=[low_tick, medium_tick, high_tick]) 
# Format the labels
cbar_1.ax.set_yticklabels([f"{int(low_tick):,}", f"{int(medium_tick):,}", f"{int(high_tick):,}"])
cbar_2.ax.set_yticklabels([f"{int(low_tick):,}", f"{int(medium_tick):,}", f"{int(high_tick):,}"])

plt.tight_layout()

# print(header) print(psf_hdu_list[0].header) print(wcs)

# -- IDENTIFICATION OF THE SIZE/BRIGHTNESS FACTOR PROBLEM --
#identify the lowest y-value in the array to find the bottom object (what we will use as reference to find the mismatch)
bottom_object_coord = np.argmin(ys) #get the lowest value in the array

btm_obj_x = int(round(xs[bottom_object_coord])) #round coords to nearest 
btm_obj_y = int(round(ys[bottom_object_coord]))

#slicing the object in the center
slice_center_row = btm_obj_y
slice_start_col = btm_obj_x - 25 #creating a 50x50 around the object (the area that we will graph)
slice_end_col = btm_obj_x + 25

og_line_graph = r_band[slice_center_row, slice_start_col:slice_end_col] #the logic is r_band[50,100] outputs the brightness (the flux value) at row 50, column 100
synthetic_line_graph = synthetic_image[slice_center_row, slice_start_col:slice_end_col]
line_graph_px_positions = np.arange(slice_start_col, slice_end_col)

plt.figure(figsize=(10, 6))
plt.plot(line_graph_px_positions, og_line_graph, label='OG Object Line Graph', color='blue', linewidth=2)
plt.plot(line_graph_px_positions, synthetic_line_graph, label='Synthetic Image Object Line Graph', color='red', linewidth=2)
plt.xlabel('Pixel position (x)')
plt.ylabel('Pixel brightness')
plt.title('inspectig the size/brightness factor problem')
plt.legend()
plt.grid(alpha=0.3)
plt.show()

print(psf_data.shape)

# -- RESIZE THE PSF IMAGE TO THE CORRECT RESOLUTION --
psf_data_enlarged = zoom(psf_data, 4, order=3) #the image to enlarge, the scale factor, the order (cubic interpolation allows for a smooth result)
psf_data_enlarged = psf_data_enlarged / psf_data_enlarged.sum() * psf_data.max() #renormalize everything to ensure that sum = 1

print(psf_data_enlarged.sum()) #verify sum = 1

#replot
synthetic_image = create_synthetic_image(r_band.shape, psf_data_enlarged, xs_without_target, ys_without_target, fluxes)
fig, axes = plt.subplots(1, 2, figsize=(16, 8))
LogNorm_minimum_value = 0.1 #minumum value to handle fully black
LogNorm_maximum_value = r_band.max() #minumum value to handle fully black

logarithmic_image1 = axes[0].imshow(r_band, origin="lower", norm=LogNorm(vmin=LogNorm_minimum_value, vmax=LogNorm_maximum_value), cmap="gray")
axes[0].set_title('original r-band image')
logarithmic_image2 = axes[1].imshow(synthetic_image, origin="lower", norm=LogNorm(vmin=LogNorm_minimum_value, vmax=LogNorm_maximum_value), cmap="gray")
axes[1].set_title('synthetic PSF image')


low_tick = np.percentile(r_band[r_band > 0], 20) #ticks for image 1
medium_tick = np.percentile(r_band[r_band > 0], 50)#calculate percentiles using only positive pixels to avoid math issues
high_tick = np.percentile(r_band[r_band > 0], 90)
low_tick2 = np.percentile(synthetic_image[synthetic_image > 0], 20) #ticks for image 2
medium_tick2 = np.percentile(synthetic_image[synthetic_image > 0], 50)
high_tick2 = np.percentile(synthetic_image[synthetic_image > 0], 90)


cbar_1 = fig.colorbar(logarithmic_image1, ax=axes[0], ticks=[low_tick, medium_tick, high_tick]) # Create specific colorbars for each of the images
cbar_2 = fig.colorbar(logarithmic_image2, ax=axes[1], ticks=[low_tick, medium_tick, high_tick]) 
# Format the labels
cbar_1.ax.set_yticklabels([f"{int(low_tick):,}", f"{int(medium_tick):,}", f"{int(high_tick):,}"])
cbar_2.ax.set_yticklabels([f"{int(low_tick):,}", f"{int(medium_tick):,}", f"{int(high_tick):,}"])

plt.tight_layout()
plt.show()

print(neighbors['shape_r'])