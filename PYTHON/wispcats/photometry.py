import os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import photutils as pu


class Aperture(object):
    """Class to define apertures for photometry.

    Args:
        wisppar (wispcats.wisppar.WISPPar): Instance of 
            :class:`~.wisppar.WISPPar`

        filters (wispcats.filters.Filt): Instance of 
            :class:`~.filters.Filt`
    """
   
    def __init__(self, wisppar, filters):
        self.wisppar = wisppar
        self.filters = filters
        self.seg = fits.getdata(self.wisppar.segmap)

        self.photparams = self.get_config()
        self.ps = float(self.photparams['ps'])
    

    def get_config(self):
        """ """
        # get values for photometry from the config file
        config = self.wisppar.config
        options = config.options('photometry')
        photparams = {}
        for option in options:
            photparams[option] = config.get('photometry', option)
        return photparams
        

    def get_circular_apertures(self, pix):
        """Defines circular apertures at a variety of radii.

        Args:
            pix (tuple): x,y positions of objects

        Returns:
            (tuple): tuple containing:
                
                aps (list): list of photutils.CircularAperture objects
                radii (float): radii of circular apertures in arcsec
                rdict (dict): dict of aps corresponding to each radii
        """
        radii = np.array([float(x) for x in self.photparams['radii'].split(',')])
        # radii in pixels
        aps = [pu.CircularAperture(pix,r/self.ps) for r in radii]
        # some hacks required to get accurate aperture areas later on
        rdict = {}
        for r in radii:
            rdict['%.1f'%r] = [pu.CircularAperture(p,r/self.ps) for p in pix]

        return aps, radii, rdict


    def get_elliptical_apertures(self, pix, apix, bpix, tpix):
        """Defines elliptical apertures for each object.
        
        Args:
            pix (tuple): x,y positions of objects
            apix (float): semi-major axis of objects in pixels
            bpix (float): semi-minor axis of objects in pixels
            tpix (float): THETA_IMAGE from IR catalog

        Returns
            (tuple): tuple containing:

                aps (list): list of photutils.EllipticalAperture objects
                k (float): approximate isophotal extent, scaling factor for elliptical apertures
        """
        # theta in radians for elliptical apertures
        theta = tpix * np.pi/180.
        k = float(self.photparams['isophotal_factor'])
        aps = [pu.EllipticalAperture(p,a*k,b*k,t) for p,a,b,t in zip(pix,apix,bpix,theta)]
        return aps,k


    def get_total_apertures(self, pix, apix):
        """Defines a circular aperture as an approximate "total" aperture.
        
        x times the 
    
        Args:
            pix (tuple): x,y positions of objects
            apix: semi-major axis of objects in pixels

        Returns:
            (tuple): tuple containing:

                aps (list): list of photutils.CircularAperture objects
                k (float): total aperture scaling factor
        """
        k = float(self.photparams['totalap_factor'])
        aps = [pu.CircularAperture(p,a*k) for p,a in zip(pix,apix)]
        return aps,k


    def get_sky_apertures(self, pix):
        """Defines sky apertures to be used for background subtraction.

        Args:
            pix (tuple): x,y positions of objects

        Returns:
            (tuple): tuple containing:

                aps ():  photutils.RectangularAperture objects
                aplist (list): list of RectangularAperture objects 
        """
        sz = float(self.photparams['skysize'])
        size = sz / self.ps
        aps = pu.RectangularAperture(pix,size,size,0)
        # hack alert: define sky apertures as a list to get the 
        #   correct, masked areas later
        aplist = [pu.RectangularAperture(p,size,size,0.) for p in pix]
        return aps, aplist
        


