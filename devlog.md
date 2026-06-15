# Devlog

## The Project:

## Day 1

The first phase of this project was importing phase of the different FITS files, which I downloaded and organized succesffully in my project folder. The main issue with this step was that one of the FITS files I downloaded was incorrect or corrupted, and this caused me trouble later on.
After that, I moved on to ensuring that the image cutout opens successfully, and that I can extract the information I need and identify the r-band extension and which idex it corresponds to.
Next, I ensured that the PSF model was correct and that all the information, such as the sum which is ~1.0 or the shape, are right.
I then moved on to filtering the neighbor catalog by finding the PSF sources, and then filtering it even further to identify the PSF sources that are avaible in my specified field. To create this field, I had to use Astropy's World Coordinate System conversion module, so that I could find the Right Ascension and Declination for each pixel that defined my field.
Lastly, I plotted the identified neighboring PSF sources over the original r-band image, to ensure that eveyrthing is correct and corresponds to their actual position in the image.

### What's Next?

Next, I will have to work on creating the actual synthetic image which will then be subtracted from the main image, removing the nearby interfering objects.
In the future, I plan on "generalizing" this project, to ensure that it works with any FITS Files.
