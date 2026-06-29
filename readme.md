# Clean Aperture Photometry

## Overview

This project performs aperture photometry on quasar host galaxies by removing contaminating neighbouring sources using PSF modelling.

## Features

- FITS image loading
- World Coordinate System conversion
- Source identification
- Neighboring Sources Removal through Synthetic PSF image creation and subtraction
- Aperture photometry
- Comparison with Legacy Survey fluxes

## Requirements

pip install -r requirements.txt

## Running

Open Project.ipynb.

## Data

The FITS files are not included.
They are automatically downloaded from the Legacy Survey Catalog once Coordinates (RA, Dec) and Band are inserted.
