import numpy as np
import astropy.io.fits as fits
from astropy.convolution import Gaussian2DKernel,convolve
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

from wisperrors import *

def clean_image(image, output, cleanbywht=False, whtim=None, checkconv=False, showgauss=False):
    """Cleans edges of an image. Also works for images with chip gaps.

    Replaces edges and chip gaps with random noise generated 
    from the distribution of noise in the image.

    Based on cleanedges_general.pro written by Marc Rafelski
    
    Args:
        image (str): Filename of image to be cleaned.
        output (str): Filename of output, cleaned image. Output image will
            be overwritten if it already exists.
        cleanbywht (Optional[bool]): If True, pixels with a weight of 0 
            are also replaced with random noise. Requires a weight image 
            provided as ``whtim``. This option is good for removing bad pixels 
            before image convolution. Default is False.
        whtim (Optional[str]): Filename of weight map associated with ``image``.
            Required if ``cleanbywht`` is True.
        showgauss (Optional[bool]): If True, displays a plot showing the 
            noise distribution and the Gaussian fit. Default is False.
    """
    drz,hdr = fits.getdata(image, header=True)

    # replace all data with 1
    cln = np.copy(drz)
    cln[drz != 0] = 1

    # create Gaussian kernel
    gauss = Gaussian2DKernel(stddev=3.8)#, x_size=30, y_size=30)
    filt = convolve(cln, gauss, boundary=None, normalize_kernel=True)

    # bin data
    dat = drz[drz > 0]
    hist,bins = np.histogram(dat, bins=np.arange(0,0.2,0.001))
    bc = 0.5 * (bins[1:] + bins[:-1])
    binsize = bins[1] - bins[0]
    nbins = hist.shape[0]

    # mirror distribution
    peak = np.argmax(hist)
    right = hist[peak:]
    rbins = bc[peak:]
    left = right[::-1]
    tot = np.append([left], [right])
    ntot = tot.shape[0]
    lbins = np.arange(rbins[0]-right.shape[0]*binsize, rbins[0], binsize)
    bincens = np.append([lbins], [rbins])

    # fit a gaussian
    gaussfunc = lambda x,a,mu,sig: a * np.exp(-(x-mu)**2 / (2.*sig**2))
    p0 = [np.max(tot), bincens[ntot/2.],(np.max(bincens)-np.min(bincens))/4.]
    popt,pcov = curve_fit(gaussfunc, bincens, tot, p0=p0)
    mu = popt[1] if popt[1] >= 0. else 0.
    sigma = np.abs(popt[2])
    
    if showgauss:
        plt.bar(bc, hist, align='center', width=binsize, alpha=0.3, 
                color='b', linewidth=0)
        plt.bar(bincens, tot, align='center', width=binsize, alpha=0.3, 
                color='k', linewidth=0)
        plt.plot(bincens, gaussfunc(bincens, *popt), 'k', linewidth=2)
        plt.show()

    if cleanbywht:
        try:
            wht = fits.getdata(whtim)
        except IOError:
            raise WISPError('cleanedges.py: No weight image provided')

        # replace edges and bad pixels with random noise in the image
        bad = np.where(((filt != 0) & (filt != 1)) | (filt != 0) & (wht == 0))
    else:
        # replace edges with random noise in the image
        bad = np.where((filt != 0) & (filt != 1))

    # random numbers with mean of 0 and sigma of 1*sig
    noise = np.random.normal(0, 1, drz[bad].shape[0])
    drz[bad] = mu + noise * sigma
    fits.writeto(output, drz, hdr, clobber=True)

