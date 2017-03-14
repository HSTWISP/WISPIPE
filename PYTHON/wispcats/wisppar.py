import os
import re
from glob import glob
import ConfigParser
from datetime import datetime
from astropy.io import fits

import wispcats
from wisperrors import *


class WISPPar(object):
    """Class to set up and compile all necessary information about a WISP field.

    Collects filters and image lists, appropriate zeropoints, input and 
    output directories, etc.

    Args:
        par (str,int): ID of WISP field. Will be reformatted as 'ParXXX'.
        datadir (str): Path to top-level directory containing WISP reductions.
            Input data will be taken from 
            ``[datadir]/ParXXX/DATA/UVIS/IRtoUVIS`` 
        topdir (str): Path to top-level directory for output. Default is 
            same as ``datadir``. Output will be written to 
            ``[topdir]/ParXXX/DATA/UVIScats``
        overwrite (Optional[bool]): If True, files that already exist in 
            the output directory will be overwritten. If False and output 
            directory already exists, an exception will be raised.
            Default is False.
    """

    def __init__(self, par, datadir, topdir=None, overwrite=False):
        # format field ID as Par123, regardless of user input
        self.par = 'Par%s' % re.search('\d+', str(par)).group(0)
        self.datadir = datadir
        if topdir is None:
            self.topdir = datadir
        else:
            self.topdir = topdir

        self.imdir = ''
        self.scriptdir = ''
        self.psfdir = ''
        self.outdir = ''
        
        self.overwrite = overwrite

        self.get_directories()
        self.config, self.thresh = self.get_config()


        self.images = []
        self.irimages = []
        self.uvisimages = []
        self.filts = []
        self.reddest_filt = ''
        self.ircats = []
        self.segmap = ''

        self.get_image_lists()
        self.zps = self.get_zp()


    def create_output_directory(self):
        """Creates output directory for images and catalogs.

        The name of the output directory is ``[topdir]/ParXXX/DATA/UVIScats``
        and is set by :func:`get_directories`. If the output directory 
        already exists (from a previous reduction) and the ``overwrite``
        option if False, an exception will be raised.
        """
        try:
            os.makedirs(self.outdir)
        except OSError:
            if self.overwrite is False:
                msg = '%s already exists. Exiting.' % self.outdir
                raise WISPError(msg)


    def get_config(self):
        """Loads the config file for the reduction.

        The config file includes thresholds for image convolution and 
        parameters for photometry apertures.
        """
        config = ConfigParser.ConfigParser()
        config.read(os.path.join(self.scriptdir, 'parameters.cfg'))

        # get thresholds for convolution
        options = config.options('thresholds')
        thresh = {}
        for option in options:
            thresh[option.upper()] = config.get('thresholds', option)
        
        return config, thresh


    def get_directories(self):
        """Sets relevant input and output directories for the WISP field.

        - The input directory is ``[datadir]/ParXXX/DATA/UVIS/IRtoUVIS`` 
        - The output directory is ``[topdir]/ParXXX/DATA/UVIScats``
        
        Additionally sets the directories containing the PSFs and matching 
        kernels.
        """
        self.imdir = os.path.join(self.datadir, self.par, 'DATA/UVIS/IRtoUVIS')

        self.scriptdir = wispcats.__path__[0]
        self.psfdir = os.path.abspath(os.path.join(self.scriptdir, 'psfs'))
        self.outdir = os.path.abspath(os.path.join(self.topdir, self.par, 
                                                   'DATA/UVIScats'))


    def get_filter(self, image):
        """Returns filter of WISP image.
        
        Args: 
            image (str): Filename of image
        
        Returns:
            filt (str): Filter of image
        """
        return fits.getheader(image)['FILTER']


    def get_image_lists(self):
        """Compiles lists of images and catalogs for the WISP field.

        - Finds all science images and separates them into lists of images
            in the IR and UVIS filters. 
        - Determines the reddest available filter, which will be used for 
            determining the convolution kernels for PSF matching. 
        - Also finds the fin_F*.cat catalogs in DIRECT_GRISM
        - Sets the filename of the segmentation map to be created on the 
            UVIS pixel scale
        """
        self.images = glob(os.path.join(self.imdir, '*_sci.fits'))
        self.filts = [self.get_filter(x) for x in self.images]
        ir = ['F110W', 'F140W', 'F160W']
        uvis = ['F475X', 'F6000LP', 'F606W', 'F814W']

        self.irimages = [x for x in self.images if self.get_filter(x) in ir]
        self.uvisimages = [x for x in self.images if self.get_filter(x) in uvis]

        # determine the reddest filter to use for image convolution
        self.irimages.sort()
        self.reddest_filt = self.get_filter(self.irimages[-1])

        self.ircats = glob(os.path.join(self.datadir, self.par, 
                                        'DATA/DIRECT_GRISM/fin_F*.cat'))
        self.ircats.sort()

        # filename of segmentation map
        self.segmap = os.path.join(self.outdir, 'UVISscale_seg.fits')


    def get_zp(self):
        """Determines zeropoints for WFC3 filters based on observation date.

        Returns:
            zp (dict): A dict of appropriate zeropoints for each WFC3 filter.
        """
        image = self.images[0]
        # obs date for new photometry
        photdate = '2012-03-06'
        HSTdate = datetime.strptime(photdate, '%Y-%m-%d')
        # get DATE-OBS from header
        obsdate = fits.getheader(image)['DATE-OBS']
        date = datetime.strptime(obsdate, '%Y-%m-%d')

        zp = {}
        if date.date() >= HSTdate.date():
            # new zero points
            zp['F110W'] = 26.8223
            zp['F140W'] = 26.4524
            zp['F160W'] = 25.9463
            zp['F475X'] = 26.1579
            zp['F600LP'] = 25.8746
            zp['F606W'] = 26.0691
            zp['F814W'] = 25.0985
        if date.date() < HSTdate.date():
            # old zero points
            zp['F110W'] = 26.83
            zp['F140W'] = 26.46
            zp['F160W'] = 25.96
            zp['F475X'] = 26.15
            zp['F600LP'] = 25.85
            zp['F606W'] = 26.08
            zp['F814W'] = 25.09

        return zp


    def __str__(self):
        print '\n%s:\n   (%s)\n' % (self.par, self.outdir)
        retstr = 'Convolution filter: %s\n' % self.reddest_filt
        return retstr


    def __repr__(self):
        return str(self)

