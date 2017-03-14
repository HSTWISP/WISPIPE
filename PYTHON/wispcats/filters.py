import os
import re
import shutil
import numpy as np
from astropy.io import fits
from astropy.convolution import convolve_fft

from cleanedges import clean_image


def pyraf_import():
    """Import PSFMATCH to recalculate convolution kernels."""
    from pyraf.iraf import psfmatch


class Filt(object):
    """Class to prepare science image for a given WISP filter.

    Sets all file names necessary for the filter, cleans the images, 
    and convolves the images with higher resolution PSFs to match the 
    PSF of the reddest filter.
        
    Args:
        wisppar (wispcats.wisppar.WISPPar): Instance of 
            :class:`~.wisppar.WISPPar`
        image (str): Filename of image
    """
    def __init__(self, wisppar, image):
        """
        """
        self.wisppar = wisppar
            
        self.image = image
        self.filt = fits.getheader(self.image)['FILTER']

        self.zp = self.wisppar.zps[self.filt]
        self.exptime = fits.getheader(self.image)['EXPTIME']

        # WHT and RMS maps
        base = os.path.basename(self.image)
        self.wht = os.path.join(self.wisppar.outdir, base.replace('sci','wht'))
        self.rms = os.path.join(self.wisppar.outdir, base.replace('sci','rms'))
        shutil.copy(os.path.join(self.wisppar.imdir, base.replace('sci','wht')),
                    self.wht)
        shutil.copy(os.path.join(self.wisppar.imdir, base.replace('sci','rms')),
                    self.rms)

        filenames = self.prepare_images()
        self.cln,self.convim,self.psf,self.ker = filenames
        
    
#    def fix_rms_map(self, rmsfile, output):
#        """Fixes an RMS map for use with SExtractor v2.8+.
#
#        Replaces bad pixels with NaNs. NaNs are interpreted correctly as
#        bad pixels in Sextractor v2.8. If using an earlier version of 
#        Sextractor, photometry errors will probably be NaN if aperture 
#        includes bad pixels. 
#
#        Args:
#            rmsfile (str): Filename of RMS map
#            output (str): Filename of output RMS map. Output RMS map will 
#                be overwritten if it already exists.
#        """
#        im,hdr = fits.getdata(rmsfile, header=True)
#        rms = np.where(im != 0, im, np.nan)
#        fits.writeto(output, rms, header=hdr, clobber=True)


    def prepare_images(self):
        """Sets file names necessary for cleaning, SExtractor, etc.

        Sets the filenames of: 
    
        - the cleaned image: [image]_cln.fits
        - the convolved image and convolution kernel (if necessary): 
            [filter]_convto[red_filter].fits, ker_[filter]to[red_filter].fits
        - the psf corresponding to this filter: [filter]_psf.fits

        Returns: 
            (tuple): Tuple containing:
        
                cln (str): Filename of cleaned image

                convim (str): Filename of convolved image
    
                psf (str): Filename of PSF corresponding to this filter

                kernel (str): Filename of convolution kernel
        """
        # cleaned image
        tmp = os.path.splitext(self.image)[0] + '_cln.fits'
        cln = os.path.join(self.wisppar.outdir, os.path.basename(tmp))

        # convolved image
        if self.filt != self.wisppar.reddest_filt:
            rf = re.search('\d+', self.wisppar.reddest_filt).group(0)
            f = re.search('\d+', self.filt).group(0)
            convim = os.path.join(self.wisppar.outdir, 
                                    '%s_convto%s.fits' % (self.filt, rf))
            kernel = os.path.join(self.wisppar.psfdir, 'ker_%sto%s.fits'%(f,rf))

        else:
            convim = cln
            kernel = None

        # psfs 
        psf = os.path.join(self.wisppar.psfdir, '%s_psf.fits' % self.filt)

        return cln, convim, psf, kernel


    def clean_image(self, conv=True):
        """Cleans images, replacing edges and chip gaps.

        If image is to be convolved, also cleans bad pixels so they 
        are not convolve in the
        These pixels will have weights as zero anyway.

        Args:
            conv (bool): If True, bad pixels are cleaned as well as the edges
                of the image. If False, only the edges are cleaned. 
        """
        print 'Cleaning %s' % os.path.basename(self.image)

        if conv:
            # clean bad pixels as well as edges
            clean_image(self.image, self.cln, cleanbywht=True, whtim=self.wht)
        else:
            # clean only edges
            clean_image(self.image, self.cln, cleanbywht=False, whtim=self.wht)


    def convolve_image(self, computekernel, threshold):
        """Convolve the cleaned image to match the PSF of the reddest filter.

        It is assumed by default that the matching kernel has already been 
        calculated with the appropriate threshold. If :attr:`computekernel`
        is True, the kernel will be recalculated using IRAF's PSFMATCH. 
        In this case, the PSF matching function is computed directly from 
        pre-computed reference and input image PSFs. A filter is used to 
        remove high frequency noise from the matching function. Specifically,
        the high frequency, low signal-to-noise components of the matching 
        function are replaced with a Gaussian model computed from the low 
        frequency, high signal-to-noise components. The low frequency cutoff 
        is set by :attr:`threshold`. This threshold value is read in from 
        the configuration file ``parameters.cfg``.

        Args: 
            computekernel (bool): If True, PSF matching kernel is 
                recalculated. It is therefore possible to change the 
                threshold used in computing the matching function. 
            threshold (float): The low frequency cutoutt in units of the 
                fraction of the total input image spectrum power.
        """
        rfilt = self.wisppar.reddest_filt

        if computekernel:
            # lower resolution filter for convolution is always reddest_filt
            lorespsf = os.path.join(self.wisppar.psfdir, '%s_psf.fits' % rfilt)

            # save newly computed kernel to outdir
            self.ker = os.path.join(self.wisppar.outdir,
                                       os.path.basename(self.ker))

            # compute the matching function and save it to kernel
            pyraf_import()
            psfmatch(self.psf, lorespsf, self.psf, self.ker, filter='replace',
                     convolution='psf', background='none', threshold=threshold)

        print 'Convolving %s - kernel: %s' % (
                   os.path.basename(self.cln), os.path.basename(self.ker))

        # matching function
        k = fits.getdata(self.ker)
        # image to be convolved
        im,hdr = fits.getdata(self.cln, header=True)
        convimg = convolve_fft(im, k, allow_huge=True)
        hdrstr = 'Convolved to the PSF of %s' % rfilt
        hdr['CONVOLV'] = ('',hdrstr)
        fits.writeto(self.convim, convimg, header=hdr, clobber=True)


    def __str__(self):
        if self.ker:
            retstr = '%s\n   Images:   %s, %s\n   Rms:      %s\n   Exptime:  %.1f\n   Zp:       %.4f\n   Convolution using %s\n' % \
                (self.filt, os.path.basename(self.cln), 
                 os.path.basename(self.convim),
                 os.path.basename(self.rms), self.exptime, self.zp, 
                 os.path.basename(self.ker))
        else:
            retstr = '%s\n   Images:   %s, %s\n   Rms:      %s\n   Exptime:  %.1f\n   Zp:       %.4f\n' % (self.filt, os.path.basename(self.cln), 
                 os.path.basename(self.convim),
                 os.path.basename(self.rms), self.exptime, self.zp)
        return retstr


    def __repr__(self):
        return str(self)

