import os
from glob import glob
import numpy as np
from astropy.io import fits
from scipy import ndimage
from reproject import reproject_interp


class SegmentationMap(object):
    """ 

    NOTE: this part of the code is still pretty dependent on hardcoded filenames
    All UVIS fields have at least F110
    4 fields have UVIS and neither F140 nor F160 (60, 407, 447, 459)
    In these cases, there would be neither JHseg nor Hseg

    No fields have UVIS and F140 - only need to care about F160

"""

    def __init__(self, wisppar):
        self.wisppar = wisppar
        self.catalog = self.wisppar.ircats[0]
        
        sedir = os.path.join(self.wisppar.datadir, self.wisppar.par, 'SEX')
        self.jhseg = os.path.join(sedir, 'JH_combined_seg.fits')
        # J and H segmap file names reflect the detection threshold used
        try:
            self.jseg = glob(os.path.join(sedir, 'F110W_*seg.fits'))[0]
        except IndexError:
            self.jseg = None
        
        # H band may always be 2.0, but better to be safe
        try:
            self.hseg = glob(os.path.join(sedir, 'F160W_*seg.fits'))[0]
        except IndexError:
            self.hseg = None

        if not os.path.exists(self.jhseg):
           self.jhseg = self.jseg
        print '    Using %s as J+H segmentation map'%os.path.basename(self.jhseg)
        print '    Using %s as J segmentation map'%os.path.basename(self.jseg)
        print '    Using %s as H segmentation map\n'%os.path.basename(self.hseg)


    def object_list(self, cx, cy, seg, filename):
        """ """
        objlist = []
        for xx,yy in zip(cx,cy):
            val = seg[int(yy),int(xx)]
            if val == 0:
                s = seg[int(yy)-2:int(yy)+2,int(xx)-2:int(xx)+2]
                val = np.unique(s[s != 0])
                if val.shape[0] > 1:
                    print 'Wrong ID may be matched to object at '+\
                        '%.1f,%.1f in %s' % (xx,yy,os.path.basename(filename))
            objlist.append(val)
        return objlist


    def create_segmap(self):
        """ """
        cat = np.genfromtxt(self.catalog)
        x = cat[:,2]
        y = cat[:,3]
        obj = cat[:,1]
        # we have at least jhseg and jseg defined
        seg,hdr = fits.getdata(self.jhseg, header=True)
        jseg = fits.getdata(self.jseg)
        if os.path.exists(self.hseg):
            hseg = fits.getdata(self.hseg)

        # start with objects detected in both images
        wjh = np.where(cat[:,1] < 10000)
#        o = [seg[int(yy),int(xx)] for yy,xx in zip(y[wjh],x[wjh])]
        o = self.object_list(x[wjh], y[wjh], seg, self.jhseg)

        # remove objects that are not in the catalog
        b = [i for i in np.unique(seg) if (i not in o) & (i != 0)]
        wb = [np.where(seg == i) for i in b]
        for i in range(len(b)):
            seg[wb[i]] = 0

        # now replace obj IDs with the catalog IDs
        w = [np.where(seg == i) for i in o]
        for i in range(obj[wjh].shape[0]):
            seg[w[i]] = obj[wjh][i]

        # now the Jband-only objects
        wj = np.where((cat[:,1] >= 10000) & (cat[:,1] < 20000))
#        jo = [jseg[int(yy),int(xx)] for yy,xx in zip(y[wj],x[wj])]
        jo = self.object_list(x[wj], y[wj], jseg, self.jseg)
        w = [np.where(jseg == i) for i in jo]
        for i in range(obj[wj].shape[0]):
            seg[w[i]] = obj[wj][i]

        # now the Hband-only objects
        wh = np.where((cat[:,1] >= 20000) & (cat[:,1] < 30000))
#        ho = [hseg[int(yy),int(xx)] for yy,xx in zip(y[wh],x[wh])]
        ho = self.object_list(x[wh], y[wh], hseg, self.hseg)
        w = [np.where(hseg == i) for i in ho]
        for i in range(obj[wh].shape[0]):
            seg[w[i]] = obj[wh][i]

        if np.max(cat[:,1] >= 30000):
            # now the J and H (but not J+H) band objects
            w3 = np.where(cat[:,1] >= 30000)
#            o3 = [jseg[int(yy),int(xx)] for yy,xx in zip(y[w3],x[w3])]
            o3 = self.object_list(x[w3], y[w3], jseg, self.jseg)
            w = [np.where(jseg == i) for i in o3]
            for i in range(obj[w3].shape[0]):
                seg[w[i]] = obj[w3][i]

        fits.writeto(self.wisppar.segmap, seg, header=hdr, clobber=True)

        
    def regrid_segmap(self):
        """ """
        seg = fits.open(self.wisppar.segmap)
        img = fits.open(self.wisppar.uvisimages[0])
        arr,fp = reproject_interp(seg, img[0].header, order='nearest-neighbor')
        arr[np.isnan(arr)] = 0
        fits.writeto(self.wisppar.segmap, arr, img[0].header, clobber=True)
        
