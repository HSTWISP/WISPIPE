import os
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u


def match_cats(ra, dec, refra, refdec):
    """Matches a catalog of RA's and Dec's to a reference catalog.

    Returns the indices of the reference catalog that mach each
    source in the input catalog, and the on-sky separation
    between each source's closest match in arcsec.
    """
    # create SkyCoord objects to use with the matching
    sccat = SkyCoord(ra=ra, dec=dec, frame='icrs', unit=(u.deg,u.deg))
    screfcat = SkyCoord(ra=refra, dec=refdec, frame='icrs', unit=(u.deg,u.deg))

    # idx    - indices of matched sources in reference cat
    # sep2d  - on-sky angular separation between closest match
    # dist3d - 3D distance between closest matches
    idx,sep2d,dist3d = match_coordinates_sky(sccat, screfcat)

    return (idx, sep2d)


def elliptical_region_pix(filename, xx, yy, a, b, theta, objid, color='green'):
    """Create a ds9 region file with elliptical apertures in image coords. """
    f = open(filename, 'w')
    f.write('global color=green dashlist=8 3 width=3 font="helvetica 10 normal  roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1      source=1\n')
    f.write('image\n')
    for i,obj in enumerate(objid):
        f.write('ellipse(%f,%f,%f,%f,%f) # color=%s text={%i}\n' %
                (xx[i], yy[i], a[i], b[i], 360.-theta[i], color, objid[i]))
    f.close()


def circular_region_pix(filename, xx, yy, radii, objid, color='green'):
    """Create a ds9 region file with elliptical apertures in image coords. """
    f = open(filename, 'w')
    f.write('global color=green dashlist=8 3 width=3 font="helvetica 10 normal  roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1      source=1\n')
    f.write('image\n')
    for i,obj in enumerate(objid):
        for j,r in enumerate(radii):
            if j == 0:
                f.write('circle(%f,%f,%f) # color=%s text={%i}\n' %
                        (xx[i], yy[i], r, color, objid[i]))
            else:
                f.write('circle(%f,%f,%f) # color=%s\n' %
                        (xx[i], yy[i], r, color))
    f.close()


def total_region_pix(filename, xx, yy, rad, objid, color='green'):
    """Create a ds9 region file with elliptical apertures in image coords. """
    f = open(filename, 'w')
    f.write('global color=green dashlist=8 3 width=3 font="helvetica 10 normal  roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1      source=1\n')
    f.write('image\n')
    for i,obj in enumerate(objid):
        f.write('circle(%f,%f,%f) # color=%s text={%i}\n' %
                (xx[i], yy[i], rad[i], color, objid[i]))
    f.close()

