import sys

from wisppar import *
from filters import *
from catalog import *


class Observation(object):
    """Class to run all necessary steps for creating a WISP PSF-matched photometry catalog.

    Args:
        par (str): WISP field
        datadir (str): Path to WISP field reduction. Input data will be taken 
            from ``[datadir]/ParXXX/DATA/UVIS/IRtoUVIS``
        topdir (Optional[str]): Path to top-level directory for output. 
            Default is same as ``datadir``. Output will be written to
            ``[topdir]/ParXXX/DATA/UVIScats``
        overwrite (Optional[bool]): If True, output directory contents 
            are overwritten. If False and the output directory already exists,
            an exception will be raised.
        computekernel (Optional[bool]): If True, the PSF matching kernels are
            recalculated. In this case, the threshold used in computing 
            the matching function can be changed by editing the
            ``parameters.cfg`` file. If False, the previously calculated
            kernels are used.
    """

    def __init__(self, par, datadir, topdir=None, overwrite=False, computekernel=False):
        print par

        self.wisppar = WISPPar(par, datadir, topdir=topdir, overwrite=overwrite)

        self.wisppar.create_output_directory()

        # write output to file
        orig_stdout = sys.stdout
        myout = file(os.path.join(self.wisppar.outdir, '%s_log.txt'%par), 'w')
        sys.stdout = myout
        print self.wisppar

        self.computekernel = computekernel

        self.filters = []


    def process_images(self):
        """Cleans and convolves (when necessary) images from all filters.

        By default, image convolution assumes that the matching kernel has 
        already been calculated with the appropriate threshold. The 
        kernel can be recalculated by setting :ivar computekernel:=True.
        """
        for image in self.wisppar.images:
            self.filters.append(Filt(self.wisppar, image))

        for f in self.filters:
            # convolve the image to match the psf of the reddest filter
            if f.filt != self.wisppar.reddest_filt:
                # clean the image
                f.clean_image(conv=True)
                threshold = self.wisppar.thresh[f.filt]
                print f.filt, threshold
                f.convolve_image(self.computekernel, threshold)
            else:
                f.clean_image(conv=False)
            print f
            
    
    def build_catalog(self):
        """Create the catalog with PSF-matched UVIS and IR photometry."""
        catalog = Catalog(self.wisppar, self.filters)
        catalog.write_cat()
        catalog.create_region_files()

