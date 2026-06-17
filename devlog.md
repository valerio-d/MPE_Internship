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

## Day 2

Today, I wrote the function that created the actual Synthetic Image.
Although I had to write code to ensure that the main object (the one which we are interested in) isn't included in the array of neighbors, which took a lot of my time, I managed to achieve a successfilly SYnthetic PSF Image that I will then be able to subtract from my Main Image.
Another problem that I encountered during this process was that the PSF images were in a different resolution than my main image. After creating a simple line graph which showed the size/resolution of both objects at once, I was able to find out that the change in resolution was a factor of 4. Therefore, I used scipy's zoom() function whcih allowed me to resize the image and achieve a precise and effective PSF Synthetic Image.

### What's Next?

Next, I will have to subtract the PSF Synthetic Image I created from the original r-band Image. This will allow me to have a clear image of the quasar only and no interferences from nearby objects. This should be fairly easy, meaning that I will not have to spend much of my time, and that I will be able to immediately start with the Aperture Photometry with the new, cleaned-up image.

## Day 3

Today, I started off by importing the line graph that I performed the other day, to my new Jupyter Notebook. I decided to switch to a Jupyter Notebook for this Porject because of how much more organized the environment is, and how quicker it is to make changes and visualize their results. During the implemtation of the line graph, I noticed that, although I had zoomed the image, the graph still didn't look right. I therefore found the brighstest point on both the original and synthetic image and realized that there was a difference of 1.3. After multiplying the r_band data by this factor to achieve the highest accuracy possible, I moved on to learning how to perform Aperture Photometry, by using the original image so that I could compare my results with the online data. However, I found out that my results were not similar to them.
