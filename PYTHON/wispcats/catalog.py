import os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, Column
import photutils as pu

from segmap import SegmentationMap
from photometry import Aperture
from utils import elliptical_region_pix, circular_region_pix, total_region_pix
from wisperrors import *


class Catalog(object):
    """Class to create the catalog of PSF-matched UVIS and IR photometry.

    Photometry is performed using 
    `photutils <https://photutils.readthedocs.io/en/stable/>`_.

    The IR segmentation map is regridded to the UVIS pixel scale and used 
    for ISO photometry. Circular and elliptical aperture positions are 
    taken from the IR catalog fin_FXXX.cat. The radii of the circular 
    apertures are set in the configuration file ``parameters.cfg``.
    The semi-major and semi-minor axes for each object are taken from the 
    IR catalog, scaled to the new pixel scale, and scaled again by an
    "approximate isophotal extent" that is also set in ``parameters.cfg``.
    A very large circular aperture is used to approximate the 
    "total magnitude" from each object.

    Args:
        wisppar (wispcats.wisppar.WISPPar): Instance of 
            :class:`~.wisppar.WISPPar`

        filters (wispcats.filters.Filt): Instance of 
            :class:`~.filters.Filt`
    """

    def __init__(self, wisppar, filters):
        self.wisppar = wisppar
        self.filters = filters

        # create segmentation map regridded to UVIS pixel scale
        print 'regrid segmentation map'
        segmap = SegmentationMap(wisppar)
        segmap.create_segmap()
        segmap.regrid_segmap()
        self.seg = fits.getdata(self.wisppar.segmap)
        # SegmentationImage object for ISO photometry
        self.segm = pu.SegmentationImage(self.seg)

        # initialize catalog
        self.cat,self.nobj,pix,a,b,t = self.initialize_catalog()

        # define apertures for photometry
        aper = Aperture(self.wisppar, self.filters)
        self.ps = aper.ps
        self.circaps,self.radii,self.rdict = aper.get_circular_apertures(pix)
        self.rad = self.radii / self.ps
        self.ellipaps,self.isophotal_factor = aper.get_elliptical_apertures(pix, a, b, t)
        self.totaps,self.totalap_factor = aper.get_total_apertures(pix,a)
        self.skyaps,self.skylist = aper.get_sky_apertures(pix)

        for f in self.filters:
            self.run_photometry(f, f.cln, conv=False)
            if f.filt != self.wisppar.reddest_filt:
                self.run_photometry(f, f.convim, conv=True)


    def write_cat(self):
        """Writes out the final catalog as a FITS table to ParXXX_cat.fits"""
        outcat = os.path.join(self.wisppar.outdir, 
                              '%s_cat.fits'%self.wisppar.par)
        print 'Writing catalog to %s' % os.path.basename(outcat)
        self.cat.write(outcat, format='fits', overwrite=True)


    def create_region_files(self):
        """Creates region files for the catalog
        
        One region file each for circular apertures, elliptical apertures, 
        and "total" apertures is created.
        """
        circregs = os.path.join(self.wisppar.outdir, 'circular_apertures.reg')
        circular_region_pix(circregs, self.cat['X_IMAGE'], self.cat['Y_IMAGE'],
                            self.rad, self.cat['NUMBER'])

        totregs = os.path.join(self.wisppar.outdir, 'total_apertures.reg')
        k = self.totalap_factor
        total_region_pix(totregs, self.cat['X_IMAGE'], self.cat['Y_IMAGE'],
                         self.cat['A_IMAGE']*k, self.cat['NUMBER'])

        ellregs = os.path.join(self.wisppar.outdir, 'elliptical_apertures.reg')
        k = self.isophotal_factor
        elliptical_region_pix(ellregs, self.cat['X_IMAGE'], self.cat['Y_IMAGE'],
                            self.cat['A_IMAGE']*k, self.cat['B_IMAGE']*k,
                            self.cat['THETA_IMAGE'], self.cat['NUMBER'])

    
    def initialize_catalog(self):
        """Initializes the catalog table of positions and photometry.
        
        Reads in the object positions in RA and Dec from the IR catalog. 
        The WCS coordinates are transformed to pixel coordinates in the 
        new images. The a, b and theta values in WCS units are added 
        directly to the catalog. The a and b values in pixel units are 
        scaled to the new pixel scale.

        Returns:
            (tuple): a tuple containing:
            
                cat (astropy.table.Table): the catalog
                nobj (int): number of objects
                pix (tuple): x,y positions of objects
                a (float): semi-major axis in pixels on new pixel scale
                b (float): semi-minor axis in pixels on new pixel scale
                t (float): THETA_IMAGE from IR catalog
                """
        hdr = fits.getheader(self.filters[0].cln)
        wcs = WCS(hdr)
        wispcat = np.genfromtxt(self.wisppar.ircats[0])
        obj = wispcat[:,1]
        nobj = obj.shape[0]
        ra = wispcat[:,7]
        dec = wispcat[:,8]
        pix = wcs.wcs_world2pix(zip(ra,dec), 1)
        x,y = pix[:,0], pix[:,1]
        awcs = wispcat[:,9]
        bwcs = wispcat[:,10]
        twcs = wispcat[:,11]
        # convert a,b in image coords to the new pixel scale
        scale = 0.08 / 0.04
        a = wispcat[:,4] * scale
        b = wispcat[:,5] * scale
        t = wispcat[:,6]
        cat = Table([obj, ra, dec, x, y, awcs, bwcs, twcs, a, b, t],
                    names=('NUMBER','X_WORLD','Y_WORLD','X_IMAGE','Y_IMAGE',
                           'A_WORLD','B_WORLD','THETA_WORLD',
                           'A_IMAGE','B_IMAGE','THETA_IMAGE'))

        return cat,nobj,pix,a,b,t
            

    def get_aperture_area(self, ap, immask):
        """Find the area of the aperture.
        
        The aperture are is necessary for scaling the mean background flux
        to match the aperture size. The returned area takes
        the subpixelation of the aperture into account.
        
        Returns:
            aparea (float): area of aperture including edge effects 
        """
#        m = ap.to_mask(method='subpixel', subpixels=5)[0]
#        m1 = m.cutout(immask) * m.data
#        mm = np.ma.masked_where(m1 == 0, m1)
#        aparea = mm.count()
        aparea = ap.mask_area(method='subpixel',subpixels=5)[0]
        return aparea


    def append_to_catalog(self, flux, fluxerr, mag, magerr, filt, aptype, conv, extrameta=None):
        """Appends photometry to the catalog.
        
        Columns of fluxes and magnitudes are added to the catalog and 
        labeled by the aperture type and filter. Columns are named as:
            - FLUX_[aptype]_[filt]
            - FLUXERR_[aptype]_[filt]
            - MAG_[aptype]_[filt]
            - MAGERR_[aptype]_[filt]
        
        If the image on which the photometry is performed has been 
        convolved, 'c' is appended to the filter ID. For example, F606Wc. 
        Extra information can be added to the meta data of the catalog 
        using `extrameta`.

        Args:
            flux (float): Array of flux values
            fluxerr (float): Array of flux errors
            mag (float): Array of magnitudes
            magerr (float): Array of magnitude errors
            filt (str): Image filter
            aptype (str): Aperture type, 'SKY','ISO','APER','ELLIP','TOTAL'
            conv (bool): Set to True to indicate that photometry was 
                performed on a convolved image. Set to False to indicate 
                the photometry was performed on an unconvolved image.
            extrameta (Optional[str]): Pass extra information to be written
                to the meta data of the catalog. 
        """
        if conv is True:
            apstr = '%s_%sc'%(aptype,filt)
        else:
            apstr = '%s_%s'%(aptype,filt)

        if aptype == 'SKY':
            self.cat.add_columns([Column(data=flux, name='FLUX_%s'%apstr),
                                 Column(data=fluxerr, name='FLUXERR_%s'%apstr)])
        elif mag is not None:
            self.cat.add_columns([Column(data=flux, name='FLUX_%s'%apstr),
                                  Column(data=fluxerr, name='FLUXERR_%s'%apstr),
                                  Column(data=mag, name='MAG_%s'%apstr),
                                  Column(data=magerr, name='MAGERR_%s'%apstr)])
        else:
            msg = "aptype != 'SKY' and mag is None."
            raise WISPError(msg)

        if extrameta:
            self.cat.meta['%s'%aptype] = extrameta


    def run_photometry(self, f, image, conv=False):
        """Run the photometry

        Args:
            f (wispcats.filters.Filt): Instance of :class:`~.filters.Filt`
                for the given filter
            image (str): Image name 
            conv (Optional[bool]): If True, SExtractor is run on the 
                convolved image. If False, SExtractor is run on the 
                original image.
        """
        print os.path.basename(image)
        im = fits.getdata(image)
        rms = fits.getdata(f.rms)
        # mask out bad pixels
        # wht <= 1/np.sqrt(100.) : rms >= 0.1
        immask = np.ma.masked_array(im, mask=(rms >= 0.1))
        # mask out object flux from sky apertures
        skymask = np.ma.masked_array(im, mask=((rms >= 0.1) & 
                                               (self.seg != 0)))
        # errors
        exptime = f.exptime
        err = pu.utils.calc_total_error(immask, rms, exptime)

        print '    Sky apertures'
        # background photometry for sky subtraction
        skyphot = pu.aperture_photometry(skymask, self.skyaps, error=err,
                                         method='subpixel', subpixels=5)
        meanbkg = np.zeros(self.nobj, dtype=float)
        meanbkgerr = np.zeros(self.nobj, dtype=float)
        # account for the skyap area that is masked out
        for i,ap in enumerate(self.skylist):
            aparea = self.get_aperture_area(ap, immask)
            # for the brightest sources, there might not be enough 
            # source-free sky within the given sky size
            if aparea < 100:
                size = 15. * self.cat['A_IMAGE'][i]
                print i, size
                newap = pu.RectangularAperture(ap.positions, size, size, 0.)
                newphot = pu.aperture_photometry(skymask, newap, error=err,
                                                 method='subpixel', subpixels=5)
                aparea = self.get_aperture_area(newap, immask)
                meanbkg[i] = newphot['aperture_sum'] / aparea
                meanbkgerr[i] = newphot['aperture_sum_err'] / aparea
            else:
                meanbkg[i] = skyphot['aperture_sum'][i] / aparea
                meanbkgerr[i] = skyphot['aperture_sum_err'][i] / aparea
        self.append_to_catalog(meanbkg, meanbkgerr, None, None, f.filt,
                               'SKY', conv=conv)

        # ISO photometry
        print '    ISO apertures'
        props = pu.source_properties(immask, self.segm, error=err)
        isoflux = np.zeros(self.nobj, dtype=float)
        isofluxerr = np.zeros(self.nobj, dtype=float)
        for i,prop in enumerate(props):
            bkg = meanbkg[i] * prop.area.value
            isoflux[i] = prop.source_sum - bkg
            isofluxerr[i] = np.sqrt((prop.source_sum_err)**2 + 
                                    (meanbkgerr[i]*prop.area.value)**2)
        isomag = np.ones_like(isoflux) * 99.
        isomagerr = np.ones_like(isoflux) * -99.
        cond = isoflux > 0.
        isomag[cond] = -2.5*np.log10(isoflux[cond]) + f.zp
        isomagerr[cond] = isofluxerr[cond]/isoflux[cond] * 2.5/np.log(10.)
        self.append_to_catalog(isoflux, isofluxerr, isomag, isomagerr, 
                               f.filt, 'ISO', conv=conv)

        # APER photometry
        print '    APER apertures'
        circphot = pu.aperture_photometry(immask, self.circaps, error=err,
                                          method='subpixel', subpixels=5)
        aperflux = np.zeros((self.nobj,len(self.radii)), dtype=float)
        aperfluxerr = np.zeros((self.nobj,len(self.radii)), dtype=float)
        for j,rad in enumerate(self.radii):
            circaps = self.rdict['%.1f'%self.radii[j]]
            for i,ap in enumerate(circaps):
                aparea = self.get_aperture_area(ap, immask)
                bkg = meanbkg[i] * aparea
                aperflux[i,j] = circphot['aperture_sum_%i'%j][i] - bkg
                aperfluxerr[i,j] = np.sqrt(
                                     (circphot['aperture_sum_err_%i'%j][i])**2 +
                                     (meanbkgerr[i]*aparea)**2)
        apermag = np.ones_like(aperflux) * 99.
        apermagerr = np.ones_like(aperflux) * -99.
        cond = aperflux > 0.
        apermag[cond] = -2.5*np.log10(aperflux[cond]) + f.zp
        apermagerr[cond] = aperfluxerr[cond]/aperflux[cond] * 2.5/np.log(10.)
        radstr = ','.join(['%.1s'%r for r in self.radii])
        apinfo = 'Circular apertures with r=[%s] arcsec' % radstr
        self.append_to_catalog(aperflux, aperfluxerr, apermag, apermagerr, 
                               f.filt, 'APER', conv=conv, extrameta=apinfo)

        # Elliptical Aperture Photometry
        # There is apparently a memory leak somewhere in photutils
        # so use a recarray to pre-allocate memory
#        ellipphot = np.recarray((self.nobj,), dtype=[('id',int),
#                                ('xcenter',float),('ycenter',float),
#                                ('aperture_sum',float),
#                                ('aperture_sum_err',float)])
        print '    Elliptical apertures'
        ellipflux = np.zeros(self.nobj, dtype=float)
        ellipfluxerr = np.zeros(self.nobj, dtype=float)
        for i,ap in enumerate(self.ellipaps):
            p = pu.aperture_photometry(immask, ap, error=err,
                                       method='subpixel', subpixels=5)
#            ellipphot[i] = Table(p).as_array()[0]
            aparea = self.get_aperture_area(ap, immask)
            bkg = meanbkg[i] * aparea
            ellipflux[i] = p['aperture_sum'] - bkg
            ellipfluxerr[i] = np.sqrt(p['aperture_sum_err']**2 + 
                                      (meanbkgerr[i]*aparea)**2)
        ellipmag = np.ones_like(ellipflux) * 99.
        ellipmagerr = np.ones_like(ellipflux) * -99.
        cond = ellipflux > 0.
        ellipmag[cond] = -2.5*np.log10(ellipflux[cond]) + f.zp
        ellipmagerr[cond] = ellipfluxerr[cond]/ellipflux[cond] * 2.5/np.log(10.)
        apinfo = 'Elliptical apertures with isophotal extent r=%.1f' % \
                    self.isophotal_factor
        self.append_to_catalog(ellipflux, ellipfluxerr, ellipmag, ellipmagerr,
                               f.filt, 'ELLIP', conv=conv, extrameta=apinfo)

        # "Total" photometry
        print '    Total apertures\n'
        totflux = np.zeros(self.nobj, dtype=float)
        totfluxerr = np.zeros(self.nobj, dtype=float)
        for i,ap in enumerate(self.totaps):
            m = ap.to_mask(method='subpixel', subpixels=5)[0]
            imcutout = m.cutout(immask) * m.data
            segcutout = m.cutout(self.seg) * m.data
            mask = (imcutout == 0) | ((segcutout != 0) & (segcutout != (i+1)))
            mm = np.ma.masked_where(mask, m.cutout(immask))
            me = np.ma.masked_where(mask, m.cutout(err))
            aparea = mm.count()
            bkg = meanbkg[i] * aparea
            flux = np.sum(mm)
            eflux = np.sqrt(np.sum(me**2))
            totflux[i] = flux - bkg
            totfluxerr[i] = np.sqrt(eflux**2 + (meanbkgerr[i]*aparea)**2)
        totmag = np.ones_like(totflux) * 99.
        totmagerr = np.ones_like(totflux) * -99.
        cond = totflux > 0.
        totmag[cond] = -2.5*np.log10(totflux[cond]) + f.zp
        totmagerr[cond] =totfluxerr[cond]/totflux[cond] * 2.5/np.log(10.)
        apinfo = 'Total apertures are circles with r=%.1f times A_IMAGE' % \
                    self.totalap_factor
        self.append_to_catalog(totflux, totfluxerr, totmag, totmagerr,
                               f.filt, 'TOTAL', conv=conv, extrameta=apinfo)
        


