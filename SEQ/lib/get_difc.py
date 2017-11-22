# Difc calculation using scipy from I(d) spectrum
# 
# fit the whole diffraction pattern all together

import os, numpy as np, yaml
from mantid import simpleapi as msa, mtd


def get_difc(
        dvalues = [1.10860231, 1.24596143, 1.357755, 1.63751414, 1.92015553, 3.13560085],
        dmin=1., dmax=3.5,
        I_d_dir = None,
        pack = None,
        fitter = None,
        maxchisq = 3.,
        min_counts = 3000 # at least this many counts in the range
        ):
    nominal_difc_all = np.load(os.path.join(I_d_dir, 'difc-nominal.npy'))
    detIDs = np.load(os.path.join(I_d_dir, 'detIDs.npy'))
    detIDs = list(detIDs)
    packinfo = yaml.load(open(os.path.join(I_d_dir, 'pack-%s.yaml' % pack)))
    pixelIDs = packinfo['pixelIDs']
    nominal_difc = nominal_difc_all[
        detIDs.index(pixelIDs['first']): detIDs.index(pixelIDs['last'])+1]
    Xbb = np.load(os.path.join(I_d_dir, 'I_d-xbb.npy'))
    x = (Xbb[:-1]+Xbb[1:])/2
    y_pack = np.load(os.path.join(I_d_dir, 'I_d-y-%s.npy' % pack))
    inrange = np.logical_and(x>dmin, x<dmax)
    # ## Loop over all pixels in D23
    npixels = pixelIDs['last'] - pixelIDs['first'] + 1
    # arrays for results
    newdifc = nominal_difc.copy()
    mask = np.zeros(npixels, dtype=bool)
    # In[66]:
    for ipixel in range(npixels):
        if (ipixel%100)==0 : print "- Working on pixel", ipixel
        # get data
        y = y_pack[ipixel]
        if np.sum(y[inrange]) < min_counts:
            print "* Not enourgh counts:", ipixel
            mask[ipixel] = 1
            continue
        # fit
        fitres = fitter.fit_pixel(x[inrange],y[inrange], dvalues)
        if fitres is None:
            mask[ipixel] = 1
            continue
        popt, fity, chisq = fitres
        if chisq > maxchisq:
            mask[ipixel] = 1
            continue
        ratio = popt[-4]
        newdifc[ipixel] = nominal_difc[ipixel] / ratio
        # print nominal_difc[ipixel]
        # print ratio
        # print newdifc[ipixel]
        # break
        continue
    return newdifc, mask


class Fitter:

    def __init__(
            self,
            peak_fractional_width=0.1,
            d0_range=(-.05, .05),
            curve_fit_options=dict(maxfev=10000)):
        self.peak_fractional_width = peak_fractional_width
        self.d0_range = d0_range
        self.curve_fit_options = curve_fit_options
        return

    def fit_pixel(self, x, y, dvalues):
        # create fitting function that make sure peak centers are at
        # specifi positions: peaks = dvalues*ratio + d_offset
        ff = create_fit_function(dvalues)
        # gather fitting inputs
        p0 = [];     bgs = [];  lower = []; upper = []
        #
        pfw = self.peak_fractional_width
        for pkindex, d0 in enumerate(dvalues):
            dmin = d0*(1-pfw)
            dmax = d0*(1+pfw)
            startindex, stopindex = np.where(x>dmin)[0][0], np.where(x<dmax)[0][-1]
            x1 = x[startindex:stopindex]
            y1 = y[startindex:stopindex]
            bg = (y1[0]+y1[-1])/2; bgs.append(bg)
            A = np.max(y1)-bg; sigma = d0*pfw/2.
            p0 += [A, sigma]
            lower += [A/10., sigma/10.]; upper+=[A*10., sigma*10.]
            # print A, sigma
            continue
        ave_bg =  np.average(bgs)
        p0 += [1., 0., 0., ave_bg]
        max_bg_slop = ave_bg/(x[-1]-x[0])*3
        lower += [0.95, self.d0_range[0], -max_bg_slop, -ave_bg*3]
        upper += [1.05, self.d0_range[-1], max_bg_slop, ave_bg*3]
        fitopts = self.curve_fit_options
        try:
            popt, pcov = sopt.curve_fit(ff, x, y, p0=p0, bounds=(lower, upper), **fitopts)
        except:
            # raise
            import warnings, traceback as tb
            msg = 'Fitting failed: %s' % (tb.format_exc(),)
            warnings.warn(msg)
            return

        # print "popt=", popt
        fity = ff(x, *popt)
        chisq = np.average((y - fity)**2 / (y+1))
        return popt, fity, chisq
    

# fitting tools
import scipy.optimize as sopt
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))
def gauss_with_bg(x, *p):
    A, mu, sigma, bg = p
    return bg + A*np.exp(-(x-mu)**2/(2.*sigma**2))
def multi_gauss_with_bg(x, *p):
    "p: [(A,mu,sigma),...] + [slope, bg] "
    assert len(p)%3 == 2
    N = len(p)//3
    slope, bg = p[-2:]
    rt = bg + slope*x
    for i in range(N):
        A, mu, sigma = p[i*3: (i+1)*3]
        rt += A*np.exp(-(x-mu)**2/(2.*sigma**2))
        continue
    return rt
def create_fit_function(ds):
    def fit_function(x, *p):
        "p: [(A,sigma),...] + [d_ratio, d0, slope, bg]"
        assert len(p)%2 == 0
        N = (len(p)-4)/2
        d_ratio, d0, slope, bg = p[-4:]
        p2 = []
        for i in range(N):
            A, sigma = p[i*2: (i+1)*2]
            p2 += [ A, d_ratio*ds[i] + d0, sigma ]
            continue
        p2 += [slope, bg]
        return multi_gauss_with_bg(x, *p2)
    return fit_function
