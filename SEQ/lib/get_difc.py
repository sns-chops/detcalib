# Difc calculation using scipy from I(d) spectrum
# 
# fit the whole diffraction pattern all together

import os, numpy as np, yaml, sys
from mantid import simpleapi as msa, mtd


class GetPackDifc:

    def __init__(
            self,
            dvalues = [1.10860231, 1.24596143, 1.357755, 1.63751414, 1.92015553, 3.13560085],
            dmin=1., dmax=3.5,
            I_d_dir = None,
            pack = None,
            fitter = None,
            maxchisq = 3.,
            min_counts = 3000, # at least this many counts in the range
            tofmin=None, tofmax = None,
            max_N_dvalues = 6,
    ):
        """
        dvalues: a list of d values to be considered. the final d list will depend on the pixel
        dmin, dmax: range of d
        I_d_dir: directory with I(d) spectra
        pack: pack name
        fitter: fitting utility
        maxchisq: maximum chi squared
        min_counts: pixels without enough counts will be ignored
        tofmin, tofmax: tof range to be considered
        max_N_dvalues: maximumn number of d values to fit for any pixel
        """
        self.dvalues = dvalues
        self.max_N_dvalues = max_N_dvalues
        # TOF range
        if tofmax is None:
            tofmax = 1e6/60/2
        if tofmin is None:
            tofmin = tofmax/4
        if tofmin >= tofmax:
            raise ValueError("tofmin should be smaller than tofmax. tofmin=%s, tofmax=%s" % (tofmin, tofmax))
        print "* tof range:", tofmin, tofmax
        self.tofmin, self.tofmax = tofmin, tofmax

        nominal_difc_all = np.load(os.path.join(I_d_dir, 'difc-nominal.npy'))
        detIDs = np.load(os.path.join(I_d_dir, 'detIDs.npy'))
        detIDs = list(detIDs)
        packinfo = yaml.load(open(os.path.join(I_d_dir, 'pack-%s.yaml' % pack)))
        pixelIDs = packinfo['pixelIDs']
        self.nominal_difc = nominal_difc = nominal_difc_all[
            detIDs.index(pixelIDs['first']): detIDs.index(pixelIDs['last'])+1]
        self.dmin = dmin; self.dmax = dmax
        # averaged nominal difc
        # ave_nominal_difc = np.average(nominal_difc)
        # compute dmin and dmax from tofmin and tofmax using ave_nominal_difc
        # dmin1 = tofmin/ave_nominal_difc; dmax1 = tofmax/ave_nominal_difc
        # print "* d range computed from tof range:", dmin1, dmax1
        # dvalues should be limited by dmin and dmax
        # self.dvalues = [d for d in dvalues if d>dmin and d<dmax]
        # print "* dvalues=", self.dvalues

        Xbb = np.load(os.path.join(I_d_dir, 'I_d-xbb.npy'))
        self.x = x = (Xbb[:-1]+Xbb[1:])/2
        self.y_pack = np.load(os.path.join(I_d_dir, 'I_d-y-%s.npy' % pack))
        # self.inrange = np.logical_and(x>dmin, x<dmax)
        self.npixels = pixelIDs['last'] - pixelIDs['first'] + 1

        self.fitter = fitter
        self.min_counts = min_counts
        self.maxchisq = maxchisq
        return


    def __call__(self):
        # arrays for results
        npixels = self.npixels
        newdifc = self.nominal_difc.copy()
        mask = np.zeros(npixels, dtype=bool)
        maxchisq = self.maxchisq
        # ## Loop over all pixels
        for ipixel in range(npixels):
            if (ipixel%100)==0 : print "- Working on pixel", ipixel
            sys.stdout.write('.'); sys.stdout.flush()
            res = self.fitOnePixel(ipixel)
            if res is None:
                print "* fit failed: pixel %s" % (ipixel, )
                mask[ipixel] = 1
                continue
            chisq = res[-1][-1]
            if chisq > maxchisq:
                print "* chisq too large: pixel %s, chisq %s" % (ipixel, chisq)
                mask[ipixel] = 1
                continue
            difc = res[0]
            newdifc[ipixel] = difc
            continue
        return newdifc, mask


    def fitOnePixel(self, ipixel):
        # import pdb; pdb.set_trace()
        nominal_difc = self.nominal_difc
        tofmin, tofmax = self.tofmin, self.tofmax
        dmin1 = tofmin/nominal_difc[ipixel]; dmax1 = tofmax/nominal_difc[ipixel]
        # print "* d range computed from tof range:", dmin1, dmax1
        # take the smaller range for d
        dmin = max(self.dmin, dmin1)
        dmax = min(self.dmax, dmax1)
        # dvalues should be limited by dmin and dmax
        dvalues = [d for d in self.dvalues if d>dmin and d<dmax]
        if len(dvalues)>self.max_N_dvalues:
            dmin, dvalues = _limit_number_of_dvalues(dvalues, self.max_N_dvalues)
        # print "* dvalues=", dvalues
        fitter = self.fitter
        min_counts = self.min_counts
        # get data
        x = self.x
        inrange = np.logical_and(x>dmin, x<dmax)
        y = self.y_pack[ipixel]
        if np.sum(y[inrange]) < min_counts:
            print "* Not enough counts:", ipixel
            return
        # fit
        fitres = fitter.fit_pixel(x[inrange],y[inrange], dvalues)
        if fitres is None:
            return
        popt, fitx, fity, chisq, fitfunc = fitres
        ratio = popt[-2-n_bg_terms[fitter.bg_type]]
        return self.nominal_difc[ipixel] / ratio, (x,y, fitx, fity), (fitfunc, popt, chisq)


def _limit_number_of_dvalues(dvalues, N):
    # the last peak (higher d) is usually a higher peak, so keep that,
    # then go backward and retain N peaks
    new_ds = dvalues[-N:]
    # print dvalues
    # print new_ds
    while new_ds[0] - dvalues[-len(new_ds)-1] < new_ds[0]*0.05 and len(new_ds)>3:
        # print new_ds
        new_ds = new_ds[1:]
    n = len(new_ds)
    assert dvalues[-n-1] < new_ds[0]
    dmin = (dvalues[-n-1] + new_ds[0])/2
    return dmin, new_ds

    
def get_difc(
        dvalues = [1.10860231, 1.24596143, 1.357755, 1.63751414, 1.92015553, 3.13560085],
        dmin=1., dmax=3.5,
        I_d_dir = None,
        pack = None,
        fitter = None,
        maxchisq = 3.,
        min_counts = 3000, # at least this many counts in the range
        tofmin=None, tofmax = None,
        ):
    """obtain difc from I(d) spectrum by fitting peaks
    """
    gpd = GetPackDifc(dvalues=dvalues, dmin=dmin, dmax=dmax, I_d_dir=I_d_dir, pack=pack,
                fitter=fitter, maxchisq=maxchisq, min_counts=min_counts,
                tofmin=tofmin, tofmax=tofmax)
    return gpd()


class Fitter:

    def __init__(
            self,
            peak_fractional_width=0.1,
            d0_range=(-.05, .05),
            curve_fit_options=dict(maxfev=1000),
            bg_type='quadratic'):
        self.peak_fractional_width = peak_fractional_width
        self.d0_range = d0_range
        self.curve_fit_options = curve_fit_options
        self.bg_type = bg_type
        return

    def fit_pixel(self, x, y, dvalues):
        # create fitting function that make sure peak centers are at
        # specifi positions: peaks = dvalues*ratio + d_offset
        ff = create_fit_function(dvalues, self.bg_type)
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
            A = max(A, 1.)
            p0 += [A, sigma]
            lower += [A/30., sigma/30.]; upper+=[A*3., sigma*2.]
            # print A, sigma
            continue
        # import pdb; pdb.set_trace()
        ave_bg =  np.average(bgs)
        p0 += [1., 0.]
        bg_type = self.bg_type
        N_bg_coeffs = n_bg_terms[bg_type]
        p0 += [0.]*(N_bg_coeffs-1) + [ave_bg]
        # bounds
        lower += [0.95, self.d0_range[0]]
        upper += [1.05, self.d0_range[-1]]
        max_bg_slope = ave_bg/(x[-1]-x[0])*3
        max_bg_a2 = ave_bg/(x[-1]-x[0])**2*3
        if bg_type == 'quadratic':
            lower.append(-max_bg_a2); upper.append(max_bg_a2)
        if bg_type != 'constant':
            lower.append(-max_bg_slope)
            upper.append(max_bg_slope)
        lower.append(-ave_bg*3)
        upper.append(ave_bg*3)
        for i, (t1, t2) in enumerate(zip(lower, upper)):
            if t1>=t2:
                print i, t1, t2
                raise
        fitopts = self.curve_fit_options
        try:
            # popt, pcov = sopt.curve_fit(ff, x, y, p0=p0, bounds=(lower, upper), sigma=(y+1)**.5, **fitopts)
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
        return popt, x, fity, chisq, ff
    

# fitting tools
import scipy.optimize as sopt
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))
def gauss_with_bg(x, *p):
    A, mu, sigma, bg = p
    return bg + A*np.exp(-(x-mu)**2/(2.*sigma**2))
def multi_gauss_with_constant_bg(x, *p):
    "p: [(A,mu,sigma),...] + [bg] "
    assert len(p)%3 == 1
    N = (len(p)-1)//3
    bg = p[-1]
    rt = bg
    for i in range(N):
        A, mu, sigma = p[i*3: (i+1)*3]
        rt += A*np.exp(-(x-mu)**2/(2.*sigma**2))
        continue
    return rt
def multi_gauss_with_linear_bg(x, *p):
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
def multi_gauss_with_quadratic_bg(x, *p):
    "p: [(A,mu,sigma),...] + [quadratic, slope, bg] "
    assert len(p)%3 == 0
    N = (len(p)-3)//3
    a2, slope, bg = p[-3:]
    rt = bg + slope*x + a2*x*x
    for i in range(N):
        A, mu, sigma = p[i*3: (i+1)*3]
        rt += A*np.exp(-(x-mu)**2/(2.*sigma**2))
        continue
    return rt
n_bg_terms = dict(
    constant=1,
    linear=2,
    quadratic=3
)
def create_fit_function(ds, bg_type='quadratic'):
    def fit_function(x, *p):
        "p: [(A,sigma),...] + [d_ratio, d0] + [bg_terms...]"
        N_bg = n_bg_terms[bg_type]
        N = (len(p)-2-N_bg)/2
        assert N*2 + 2 + N_bg == len(p)
        d_ratio, d0 = p[-2-N_bg:-N_bg]
        bg_terms = p[-N_bg:]
        p2 = []
        for i in range(N):
            A, sigma = p[i*2: (i+1)*2]
            p2 += [ A, d_ratio*ds[i] + d0, sigma ]
            continue
        p2 += bg_terms
        fname = 'multi_gauss_with_%s_bg' % bg_type
        return eval(fname)(x, *p2)
    return fit_function
