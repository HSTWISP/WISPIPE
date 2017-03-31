
wispcats
==============

Produce a catalog of UVIS and IR psf-matched photometry using the outputs of the WISP pipeline.


Additional Requirements
-----------------------
The code must be run in astroconda. The following additional packages 
are required:

* photutils: v0.3
* reproject: v0.3


Creating the catalog
--------------------
Options include:
1. Run the code from `PYTHON/wispcats`
2. Add `wispcats` to your `PYTHONPATH` and copy `create_catalog.py` to 
   your desired directory.

To run it on a single field:
```
python create_catalog.py 123
```
To run it on several fields in a row:
```
python create_catalog.py 123 124 125
```
(You can enter the fields however you want, 'Par123', '123', 'field123'.)

Inputs are taken from `ParXXX/DATA/UVIS/IRtoUVIS`

All outputs are written to `ParXXX/DATA/UVIScats`. 
Outputs include several images, a segmentation map, some region files, 
and the final catalog, which is a FITS table called `ParXXX_cat.fits`.
The FITS table also includes some meta data with information about the 
apertures used.

Parameters necessary for photometry and for calculating the convolution
matching kernels are stored in `parameters.cfg`

Steps
-----
* **Image cleaning**
    All edges, chip gaps, and bad pixels are replaced with randomly generated 
    Gaussian noise.

* **Image convolution**
    Each image is convolved with a kernel calculated to match the Hband PSF. 
    All PSFs are approximated as Gaussians and `IRAF`'s `PSFMATCH` 
    was used to calculated the matching kernels. If you want to recalculate 
    the matching kernels using a different threshold, change the thresholds in 
    `parameters.cfg` and set computekernel to `True` in 
    `create_catalog.py`.

* **Segmentation map**
    The segmentation maps in the `SEX/` directory created by the pipeline
    (`JH_combined_seg.fits`, `F110W_*_seg.fits` and `F160W_*_seg.fits`) are 
    combined and regridded onto the UVIS pixel scale. 

    It is possible that this process will produce an incorrect segmentation 
    map. This is very rare. If the process failed, you will see that the 
    background of the segmentation map will be a number > 0. Please check the 
    segmentation map for each field. 

* **Photometry**
The photometry is performed using photutils. 

    * ISO photometry uses the segmentation map. 

    * APER photometry is performed with circular apertures with radii set in 
        `parameters.cfg` (in arcsec). We can add as many as we want. 
        It is currently set to use r=0.2, 0.5, 1.0, 1.5 arcsec.

    * Elliptical apertures are used, and the scaling factor for a and b of the 
        ellipse, defined in `parameters.cfg`, is currently set to 4.

    * A very large circular aperture is used to approximate the "total" 
        magnitude. The circle's radius is `f * a`, with `f` defined in 
        `parameters.cfg` and currently set to 5. For this aperture, all 
        flux from other sources is masked out. These "total" apertures 
        may include too much noise, and one might prefer the elliptical 
        apertures.

    * Local sky subtraction is performed in rectangular apertures around each 
        object. The rectangle size, in `parameters.cfg`, is 10". All source 
        flux is masked out of the sky apertures.

Some Caveats
------------
* The code should be robust for any combination of filters up until the 
    segmentation creation stage. This is really only coded to handle the case 
    where we have F110 and F160. I think this is fine for UVIS fields, but if 
    you have a field that has a different subset of filters, let me know.

* The pipeline's segmentation maps sometimes have two isolated islands of flux 
    identified as the same object. The ISO photometry will include flux from 
    both "islands".

* If there is a very large, bright object in the sky that covers more than 
    10" on a side, the sky apertures will be too small. I have tried to 
    account for this by allowing the sky apertures to increase for a given 
    object to be 15 times the object's semi-major axis. However, I haven't had 
    a chance to rigorously test this feature as I didn't have any fields with 
    bright sources that passed the pipeline.

* The region files that are created show the apertures that were used for 
    photometry. The angles of the ellipses are sometimes incorrect in the 
    region file. This has something to do with how ds9 defines the angles. 
    I have confirmed that the angles are correct internal to photutils.
