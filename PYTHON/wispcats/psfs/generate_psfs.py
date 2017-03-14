#! /usr/bin/env python
import numpy as np
from astropy.io import fits
from astropy.convolution import convolve_fft


class PSF():
    """Class to create PSFS for WFC3 direct images.
    """

    def __init__(self, filt, pixscale=0.04, size=10, savepsf=True, blend_psfs=False):
        """Initialize the class.

        :param filt: WFC3 filter for which to create a PSF
        :type filt: string
        :param pixscale: pixel scale of images (default is UVIS 0.04"/pix)
        :type pixscale: float
        :param size: size of psf in arcsecs (default is 10 arcsec)
        :type size: float
        """
        self.filt = filt
        self.pixscale = pixscale

        # get FWHM of PSF in arcsec
        self.info = self.lookup_fwhms()
        self.fwhm_arcsec = self.info['arcsec']
        # convert to pixels for correct pixscale
        self.fwhm_pix = self.fwhm_arcsec / self.pixscale
        # generate gaussian PSF
        self.gaussian_psf = self.create_gaussian_psf(size)

        # blend Gaussian PSF with TinyTim PSF?
        if blend_psfs:
            self.psf = self.blend_psfs()
        else:
            self.psf = self.gaussian_psf

        if savepsf:
            # filename of output 
            self.output = '%s_psf.fits' % self.filt
            self.save_psf()
        


    def lookup_fwhms(self):
        """Gets info on the filter's PSF.

        Info for the filters can be found in the `UVIS`_ and `IR`_
        handbooks.

        .. _UVIS: http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c06\_uvis07.html
        .. _IR: http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c07\_ir07.html

        :returns: dict -- the wavelength and PSF FWHM in pixels and arcsecs

        """
        info = { \
            'F475X': \
                {'wave':500, 'pix_on_native_scale':1.675, 'arcsec':0.067}, \
            'F606W': \
                {'wave':600, 'pix_on_native_scale':1.681, 'arcsec':0.067}, \
            'F600LP': \
                {'wave':700, 'pix_on_native_scale':1.746, 'arcsec':0.070}, \
            'F814W': \
                {'wave':800, 'pix_on_native_scale':1.844, 'arcsec':0.074}, \
            'F110W': \
                {'wave':1100, 'pix_on_native_scale':1.019, 'arcsec':0.130}, \
            'F140W': \
                {'wave':1400, 'pix_on_native_scale':1.100, 'arcsec':0.141}, \
            'F160W': \
                {'wave':1500, 'pix_on_native_scale':1.136, 'arcsec':0.154}}

        try:
            return info[self.filt]
        except:
            print 'Unknown filter: %s' % self.filt
            exit()

    
    def create_gaussian_psf(self, size):
        """Creates a 2D Gaussian PSF

        :param size: size of psf in arcsecs
        :type size: float

        """
        # convert FWHM to sigma for Gaussian distribution
        sig = self.fwhm_pix / (2. * np.sqrt(2 * np.log(2.)))

        # convert desired PSF size to pixels
        size_pix = size / self.pixscale

        # generate Gaussian of desired size
        x = np.arange(size_pix)
        X,Y = np.meshgrid(x,x)
        center = np.mean(x)
        gauss = np.exp(-((X-center)**2 + (Y-center)**2) / (2.*(sig**2)))
        # normalize by total to match TinyTim PSFs 
        ### (How are TinyTim PSFs normalized?) ###
        gauss = gauss / np.sum(gauss)
        return gauss


    def blend_psfs(self):
        """
        Method to
        http://homepages.inf.ed.ac.uk/rbf/HIPR2/gsmooth.htm
        http://www.site.uottawa.ca/~edubois/courses/CEG4311/plotting_2D.pdf
        """
        pass

    
    def save_psf(self):
        """Saves the PSF as a FITS file."""
        # create HDU with a default header
        hdu = fits.PrimaryHDU(self.psf)
        # add filter 
        hdu.header['FILTER'] = self.filt
        hdu.writeto(self.output, clobber=True)



