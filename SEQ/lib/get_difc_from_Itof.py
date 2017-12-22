# Difc calculation using scipy from I(tof) spectrum
# 
# fit the whole diffraction pattern all together

import os, numpy as np, yaml, sys
from mantid import simpleapi as msa, mtd


def get_nominal_difc(nxspath, init_IDF, outdir):
    if not os.path.exists(outdir): os.makedirs(outdir)
    # ## Compute nominal difc 
    ws = msa.LoadEventNexus(nxspath, FilterByTimeStart=0, FilterByTimeStop=1) # load just one second
    #
    msa.LoadInstrument(ws, Filename=init_IDF, RewriteSpectraMap=False)
    import shutil
    shutil.copyfile(init_IDF, os.path.join(outdir, 'init_IDF.xml'))
    #
    difc = msa.CalculateDIFC(InputWorkspace=ws)
    difc = difc.extractY().flatten().copy()
    msa.DeleteWorkspace('difc')
    np.save(os.path.join(outdir, 'difc-nominal.npy'), difc)
    return



class GetPackDifc:

    def __init__(
            self,
            dvalues = [1.10860231, 1.24596143, 1.357755, 1.63751414, 1.92015553, 3.13560085],
            dmin=1., dmax=3.5,
            I_tof_dir = None,
            pack = None,
            fitter = None,
            maxchisq = 3.,
            min_counts = 3000, # at least this many counts in the range
            tofmin=None, tofmax = None,
            max_N_dvalues = 6,
            T0_of_E = None, L1=20.0114, L2=None,
    ):
        """
        dvalues: a list of d values to be considered. the final d list will depend on the pixel
        dmin, dmax: range of d
        I_tof_dir: directory with I(tof) spectra
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

        nominal_difc_all = np.load(os.path.join(I_tof_dir, 'difc-nominal.npy'))
        detIDs = np.load(os.path.join(I_tof_dir, 'detIDs.npy'))
        detIDs = list(detIDs)
        packinfo = yaml.load(open(os.path.join(I_tof_dir, 'pack-%s.yaml' % pack)))
        pixelIDs = packinfo['pixelIDs']
        self.nominal_difc = nominal_difc = nominal_difc_all[
            detIDs.index(pixelIDs['first']): detIDs.index(pixelIDs['last'])+1]
        self.dmin = dmin; self.dmax = dmax
        Xbb = np.load(os.path.join(I_tof_dir, 'I_tof-xbb.npy'))
        self.x = x = (Xbb[:-1]+Xbb[1:])/2
        self.y_pack = np.load(os.path.join(I_tof_dir, 'I_tof-y-%s.npy' % pack))
        self.npixels = pixelIDs['last'] - pixelIDs['first'] + 1
        self.fitter = fitter
        self.min_counts = min_counts
        self.maxchisq = maxchisq
        self.T0_of_E = T0_of_E; self.L1 = L1; self.L2 = L2
        return


    def __call__(self):
        # arrays for results
        npixels = self.npixels
        newdifc = self.nominal_difc.copy()
        mask = np.zeros(npixels, dtype=bool)
        maxchisq = self.maxchisq
        signature_d = np.zeros(npixels)
        # ## Loop over all pixels
        for ipixel in range(npixels):
            if (ipixel%100)==0 : print "- Working on pixel", ipixel
            sys.stdout.write('.'); sys.stdout.flush()
            res = self.fitOnePixel(ipixel)
            if res is None:
                print "* fit failed: pixel %s" % (ipixel, )
                mask[ipixel] = 1
                continue
            chisq = res.chisq
            if chisq > maxchisq:
                print "* chisq too large: pixel %s, chisq %s" % (ipixel, chisq)
                mask[ipixel] = 1
                continue
            newdifc[ipixel] = res.difc
            # dvalues1 = np.array(res.params.dvalues[::-1])
            # peakcounts = res.popt[-2::-2] [:len(dvalues1)]
            # good_d = dvalues1[peakcounts>self.min_counts/50]
            #if len(good_d):
            #    signature_d[ipixel] = good_d[0]
            ds = np.array(res.params.dvalues[::-1])
            peakheight = np.array(res.popt[-2::-2] [:len(ds)])
            peakwidth = np.array(res.popt[-1::-2] [:len(ds)])
            peakcounts = peakheight*peakwidth
            signature_d[ipixel] = np.sum(ds*peakcounts)/np.sum(peakcounts)
            continue
        return newdifc, mask, signature_d


    def fitOnePixel(self, ipixel):
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
        # convert back to tof range
        tofmin = nominal_difc[ipixel]*dmin; tofmax = nominal_difc[ipixel]*dmax
        # print "* dvalues=", dvalues
        fitter = self.fitter
        min_counts = self.min_counts
        # get data
        if self.T0_of_E:
            x = self._correctT0(self.x, L2 = self.L2[ipixel])
        else:
            x = self.x
        inrange = np.logical_and(x>tofmin, x<tofmax)
        y = self.y_pack[ipixel]
        if np.sum(y[inrange]) < min_counts:
            print "* Not enough counts:", ipixel
            return
        # fit
        fitres = fitter.fit_pixel(x[inrange],y[inrange], dvalues, nominal_difc[ipixel])
        if fitres is None:
            return
        popt, (p0, lower, upper), fitx, expy, fity, chisq, fitfunc = fitres
        difc, t0 = popt[:2]
        d = (fitx-t0)/difc
        return PixelFitRes(
            popt, chisq, fitfunc,
            nominal_difc[ipixel], difc, t0,
            PixelFitRes.PlotData(x,y, fitx, expy, fity, d),
            PixelFitRes.Params(p0, lower, upper, dvalues),
            )

    def _correctT0(self, tof, L2=None):
        L = self.L1 + L2
        v = L/tof*1e6
        #print v
        E = v*v*VS2E
        #print E
        T0 = self.T0_of_E(E)
        # print "T0=", T0
        # import pdb; pdb.set_trace()
        return tof-T0

    
class PixelFitRes:
    def __init__(self, popt, chisq, fitfunc, nominal_difc, difc, t0, plotdata, params):
        self.popt = popt
        self.chisq = chisq
        self.fitfunc = fitfunc
        self.nominal_difc = nominal_difc
        self.difc = difc
        self.t0 = t0
        self.plotdata = plotdata
        self.params = params
        return
    class Params:
        def __init__(self, p0, lower, upper, dvalues):
            self.p0 = p0
            self.lower = lower
            self.upper = upper
            self.dvalues = dvalues
    class PlotData:
        def __init__(self, x, y, fitx, expy, fity, d):
            self.x = x
            self.y = y
            self.fitx = fitx
            self.expy = expy
            self.fity = fity
            self.d = d

neutron_mass = 1.6749286e-27;
mN = neutron_mass;
electron_charge = 1.60217733e-19;
e = electron_charge;
VS2E = mN/(2e-3*e)

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

    
class Fitter:

    def __init__(
            self,
            peak_fractional_width=0.1,
            t0_range=(0., 200.),
            curve_fit_options=dict(maxfev=1000),
            bg_type='quadratic'):
        self.peak_fractional_width = peak_fractional_width
        self.t0_range = t0_range
        self.curve_fit_options = curve_fit_options
        self.bg_type = bg_type
        return

    def fit_pixel(self, x, y, dvalues, nominal_difc):
        # create fitting function that make sure peak centers are at
        # specific positions: peak_centers at dvalues
        ff = create_fit_function(dvalues, self.bg_type)
        # gather initial guess of parameters for fitting
        p0 = [];     bgs = [];  lower = []; upper = []
        pfw = self.peak_fractional_width
        # go through all the peaks
        for pkindex, d0 in enumerate(dvalues):
            dmin = d0*(1-pfw);   dmax = d0*(1+pfw)
            tofmin = dmin*nominal_difc; tofmax = dmax*nominal_difc
            startindex, stopindex = np.where(x>tofmin)[0][0], np.where(x<tofmax)[0][-1]
            x1 = x[startindex:stopindex]
            y1 = y[startindex:stopindex]
            bg = (y1[0]+y1[-1])/2; bgs.append(bg)
            A = np.max(y1)-bg; sigma = d0*pfw/2.
            A = max(A, 1.)
            p0 += [A, sigma]
            lower += [A/30., sigma/30.]; upper+=[A*1000., sigma*5.]
            # print A, sigma
            continue
        # import pdb; pdb.set_trace()
        ave_bg =  np.average(bgs)
        bg_type = self.bg_type
        N_bg_coeffs = n_bg_terms[bg_type]
        p0 = [nominal_difc, 0.] + [0.]*(N_bg_coeffs-1) + [ave_bg] + p0
        # bounds
        delta_d = (x[-1]-x[0])/nominal_difc
        max_bg_slope = ave_bg/delta_d*3; max_bg_a2 = ave_bg/delta_d**2*3
        bg_lower = []; bg_upper = []
        if bg_type == 'quadratic':
            bg_lower.append(-max_bg_a2); bg_upper.append(max_bg_a2)
        if bg_type in ['quadratic', 'linear']:
            bg_lower.append(-max_bg_slope); bg_upper.append(max_bg_slope)
        bg_lower.append(-ave_bg*3); bg_upper.append(ave_bg*3)
        lower = [0.97*nominal_difc, self.t0_range[0]] + bg_lower + lower
        upper = [1.03*nominal_difc, self.t0_range[-1]] + bg_upper + upper
        # import pdb; pdb.set_trace()
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
        return popt, (p0, lower, upper), x, y, fity, chisq, ff
    

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
    def fit_function(tof, *p):
        "p: difc, t0, {bg_coeffs}, {(A,sigma),...}"
        N_bg = n_bg_terms[bg_type]
        N = (len(p)-2-N_bg)/2
        assert N*2 + 2 + N_bg == len(p)
        difc, t0 = p[:2]
        x = (tof - t0)/difc
        bg_coeffs = p[2:2+N_bg]
        offset = 2+N_bg
        p2 = []
        for i in range(N):
            A, sigma = p[offset+i*2: offset+(i+1)*2]
            p2 += [ A, ds[i], sigma ]
            continue
        p2 += bg_coeffs
        fname = 'multi_gauss_with_%s_bg' % bg_type
        return eval(fname)(x, *p2)
    return fit_function
