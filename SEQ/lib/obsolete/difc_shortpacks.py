import os, numpy as np
from mantid import simpleapi as msa, mtd

workdir = "/SNS/users/lj7/dv/sns-chops/detcalib/SEQ"
os.chdir(workdir)

# Inputs
dvalues = [2.72509327,  2.89039789, 4.26940063,  5.00631601, 8.17527981]
dvalues = dvalues[4:]
nxspath = '/SNS/SEQ/IPTS-19573/nexus/SEQ_130249.nxs.h5'
initial_idf = './SEQUOIA_Definition_guessshortpacks.xml'
packname = 'C25T'
# packtype = 'eightpack-bottom'
packtype = 'eightpack-top'
x_path = 'C60-C26T-I_d-x.npy'
# y_path = 'C60-I_d-y-C25T.npy'
y_path = 'C60-C25T-I_d-y_all_pixels.npy'

d_spacing_max_mismatch = 0.2 # maximum fractional mismatch of d spacing values allowed.
d_spacing_peak_width = 0.1 # fractional width of d spacing peak.
maxchisq = 3.  # if chisq>maxchisq, mask this pixel
min_counts = 2000  # if total couts of the peak < min_counts, don't count this peak

# Outputs
difc_outpath = "C60-difc-2-C25T.npy"
difc_mask_outpath = 'C60-difc-2-C25T-mask.npy'


# ## Compute nominal difc
ws = msa.LoadEventNexus(nxspath, FilterByTimeStart=0, FilterByTimeStop=1)
msa.LoadInstrument(ws, Filename=initial_idf, RewriteSpectraMap=False)
difc = msa.CalculateDIFC(InputWorkspace=ws)
difc = difc.extractY().flatten().copy()
msa.DeleteWorkspace('difc')


# # Get pack pixel IDs
instrument = ws.getInstrument()
pack = instrument.getComponentByName("%s/%s" % (packname, packtype))
firstpixel = pack[0][0].getID()
lasttube = pack[pack.nelements()-1]
lastpixel = lasttube[lasttube.nelements()-1]
lastpixel = lastpixel.getID()


# Get detID list
detIDs = []
for i in range(ws.getNumberHistograms()):
    sp = ws.getSpectrum(i)
    dets = list(sp.getDetectorIDs())
    assert len(dets)==1
    detIDs.append(dets[0])
    continue
for i in range(len(detIDs)-1):
    assert detIDs[i] < detIDs[i+1]

# Get pack pixel indexes
startindex = detIDs.index(firstpixel)
endindex = detIDs.index(lastpixel)
msa.DeleteWorkspace('ws')


# # load x and y_pack
x = np.load(x_path)
y_pack = np.load(y_path)

# # fitting
import scipy.optimize as sopt

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def gauss_with_bg(x, *p):
    A, mu, sigma, bg = p
    return bg + A*np.exp(-(x-mu)**2/(2.*sigma**2))


def fit_pixel(x, y, dvalues, peaks_thispixel=None, chisq_thispixel=None):
    if peaks_thispixel is None:
        peaks_thispixel = np.zeros(len(dvalues), dtype=float)
    if chisq_thispixel is None:
        chisq_thispixel = np.zeros(len(dvalues), dtype=float)
    for pkindex, d0 in enumerate(dvalues):
        dmin = d0*(1.-d_spacing_max_mismatch)
        dmax = d0*(1.+d_spacing_max_mismatch)
        startindex, stopindex = np.where(x>dmin)[0][0], np.where(x<dmax)[0][-1]
        x1 = x[startindex:stopindex]
        y1 = y[startindex:stopindex]
        # print x1
        # print y1
        guess_center = x1[np.argmax(y1)]
        startindex, stopindex = np.where(x>guess_center-d_spacing_peak_width*d0)[0][0], np.where(x<guess_center+d_spacing_peak_width*d0)[0][-1]
        x1 = x[startindex:stopindex]
        y1 = y[startindex:stopindex]
        if np.sum(y1) < min_counts:
            peaks_thispixel[pkindex] = np.nan
            continue
        bg = (y1[0]+y1[-1])/2
        p0 = np.max(y1)-bg, guess_center, 0.01, bg
        # print "p0=",p0
        # print x1
        # print y1
        try:
            popt, pcov = sopt.curve_fit(gauss_with_bg, x1, y1, p0=p0)
        except:
            peaks_thispixel[pkindex] = np.nan
            continue
        # print "popt=", popt
        chisq1 = np.average((y1 - gauss_with_bg(x1, *popt))**2 / (y1+.1))
        chisq_thispixel[pkindex] = chisq1
        # print "chisq=", chisq1
        if chisq1 > maxchisq:
            peaks_thispixel[pkindex] = np.nan
            continue
        peaks_thispixel[pkindex] = popt[1]
        continue
    return peaks_thispixel


# ## Loop over all pixels in the pack
newdifc = difc.copy()
mask = np.zeros(difc.shape, dtype=bool)
chisq = np.zeros((len(dvalues),) + difc.shape, dtype=float)
peaks = np.zeros((len(dvalues),) + difc.shape, dtype=float)

firstindex, lastindex = startindex, endindex
for ipixel, pixel in enumerate(range(firstindex, lastindex+1)):
    if (pixel%100)==0 : print pixel
    # hack
    # if ipixel!=356: continue
    # print pixel
    y = y_pack[ipixel]
    peaks_thispixel = peaks[:, pixel]
    chisq_thispixel = chisq[:, pixel]
    fit_pixel(x,y, dvalues, peaks_thispixel, chisq_thispixel)
    Ngoodfit = np.sum(np.isfinite(peaks_thispixel))
    # print Ngoodfit
    if not Ngoodfit:
        mask[pixel] = 1
        continue
    ratios = peaks_thispixel / dvalues
    # print "ratios=", ratios
    good_ratios = ratios[np.isfinite(ratios)]
    # print "good_ratios=", good_ratios
    average_ratio = np.average(good_ratios)
    newdifc[pixel] = difc[pixel] * average_ratio
    # print average_ratio
    # print newdifc[pixel]
    continue


np.save(difc_outpath, newdifc)
np.save(difc_mask_outpath, mask)
