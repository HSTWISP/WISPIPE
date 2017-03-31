import os
import sys
import re
import inspect
import ConfigParser

import wispcats


def createcatalog(field):
    """
    """
    # format par ID as Par123, regardless of user input
    par = 'Par%s' % re.search('\d+', field).group(0)

    datadir = os.path.join(os.environ['WISPDATA'],'aXe')
    topdir = os.path.join(os.environ['WISPDATA'],'aXe')
    overwrite = True
    computekernel = False

    obs = wispcats.observation.Observation(par, datadir, topdir, 
                                           overwrite=overwrite,
                                           computekernel=computekernel)
    obs.process_images()

    try:
        obs.build_catalog()
    except ValueError:
        sys.stdout = obs.orig_stdout
        print 'Segmap generation failed for %s. Skipping %s.'%(par,par)


def create_multiple_catalogs(fields):
    """ """
    for field in fields:
        createcatalog(field)


def main():
    if len(sys.argv) > 2:
        pars = sys.argv[1:]
        print pars
        create_multiple_catalogs(pars)
    elif len(sys.argv) == 2:
        par = sys.argv[1]
        createcatalog(par)
    else:
        print 'Missing argument. Call as:\n   > python create_catalog.py 123\nOr:\n   > python create_catalog.py 123 124 125'
        exit()


if __name__ == '__main__':
    main()
