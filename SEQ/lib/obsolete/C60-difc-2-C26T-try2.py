
# coding: utf-8

import matplotlib
import matplotlib.pyplot as plt

import os, numpy as np
from mantid import simpleapi as msa, mtd

workdir = "/SNS/users/lj7/dv/sns-chops/detcalib/SEQ"
os.chdir(workdir)

dvalues = [2.72509327,  2.89039789, 4.26940063,  5.00631601, 8.17527981]


# ## Compute nominal difc
nxspath = '/SNS/SEQ/IPTS-19573/nexus/SEQ_130249.nxs.h5'
ws = msa.LoadEventNexus(nxspath, FilterByTimeStart=0, FilterByTimeStop=1)
msa.LoadInstrument(ws, Filename='./SEQUOIA_Definition_guessshortpacks.xml', RewriteSpectraMap=False)
difc = msa.CalculateDIFC(InputWorkspace=ws)
difc = difc.extractY().flatten().copy()
msa.DeleteWorkspace('difc')

# # Get pack index
instrument = ws.getInstrument()
pack = instrument.getComponentByName("C26T/eightpack-top")
firstpixel = pack[0][0].getID()
lasttube = pack[pack.nelements()-1]
lastpixel = lasttube[lasttube.nelements()-1]
lastpixel = lastpixel.getID()
print "first and last pixel ID:", firstpixel, lastpixel

detIDs = []
for i in range(ws.getNumberHistograms()):
    sp = ws.getSpectrum(i)
    dets = list(sp.getDetectorIDs())
    assert len(dets)==1
    detIDs.append(dets[0])
    continue
for i in range(len(detIDs)-1):
    assert detIDs[i] < detIDs[i+1]
msa.DeleteWorkspace('ws')

startindex = detIDs.index(firstpixel)
endindex = detIDs.index(lastpixel)
print "array indexes of first and last pixel", startindex, endindex


# # Get run times of all data files
nxsfiles = [
    '/SNS/SEQ/IPTS-19573/nexus/SEQ_130249.nxs.h5',
    '/SNS/SEQ/IPTS-19573/nexus/SEQ_130250.nxs.h5',
    '/SNS/SEQ/IPTS-19573/nexus/SEQ_130251.nxs.h5',
    '/SNS/SEQ/IPTS-19573/nexus/SEQ_130252.nxs.h5',
]
def getRunTime(p):
    ws = msa.LoadEventNexus(Filename=p, FilterByTimeStart=0, FilterByTimeStop=0)
    run = ws.getRun()
    t = (run.endTime() - run.startTime()).total_seconds()
    msa.DeleteWorkspace('ws')
    return t
runtimes = dict()
for f in nxsfiles:
    runtimes[f] = getRunTime(f)
print "run times:", runtimes


# # all pixels in C25T
dmin, dmax, delta_d=2., 11., 0.02
Nd = int((dmax-dmin)/delta_d)
print "Number of d bins:", Nd
Npixels = 1024
y_pack = np.zeros((Npixels, Nd))
xbb_saved = None
dt = 1000.
for nxsfile in nxsfiles:
    print nxsfile
    t_total = runtimes[nxsfile]
    for tstart in np.arange(0, t_total-dt, dt):
        print tstart
        tend = min(t_total-1, tstart+dt)
        ws = msa.LoadEventNexus(nxsfile, FilterByTimeStart=tstart, FilterByTimeStop=tend)
        msa.LoadInstrument(ws, Filename='./SEQUOIA_Definition_guessshortpacks.xml', RewriteSpectraMap=False)
        I_d = msa.ConvertUnits(InputWorkspace=ws, Target='dSpacing', EMode='Elastic')
        I_d = msa.Rebin(InputWorkspace=I_d, Params='%s,%s,%s' % (dmin, delta_d, dmax))
        for i, pixelindex in enumerate(range(startindex, endindex+1)):
            I_d_pixel = msa.SumSpectra(
                InputWorkspace=I_d, StartWorkspaceIndex=pixelindex, EndWorkspaceIndex=pixelindex)
            y = I_d_pixel.readY(0)
            xbb = I_d_pixel.readX(0)
            if xbb_saved is None: xbb_saved = np.array(xbb, copy=True)
            y_pack[i] += y
            msa.DeleteWorkspace('I_d_pixel')
            continue
        msa.DeleteWorkspaces(['ws', 'I_d'])
        continue

np.save("xbb_saved.npy", xbb_saved)
xbb = np.arange(dmin, dmax+delta_d/2., delta_d)
np.save("xbb.npy", xbb)
np.save("y_pack.npy", y_pack)


"""
# In[150]:


x = (xbb[1:]+xbb[:-1])/2


# In[152]:


plt.figure(figsize=(7,4))
plt.plot( x, y_pack[0])


# In[153]:


plt.figure(figsize=(7,4))
plt.plot(x, y_pack[2])


# ## pixels near the center of tube 1

# In[155]:


plt.figure(figsize=(7,4))
for i in range(55, 65):
    plt.plot(x, y_pack[i]+i*300)


# ## pixels at the edges

# In[156]:


plt.figure(figsize=(7,4))
pixel = 0
for tube in range(8):
    plt.plot(x, y_pack[pixel+tube*128]+tube*10)


# ## pixels #20

# In[159]:


plt.figure(figsize=(7,4))
pixel = 20
for tube in range(8):
    plt.plot(x, y_pack[pixel+tube*128]+tube*200)


# ## pixels #100

# In[160]:


plt.figure(figsize=(7,4))
pixel = 100
for tube in range(8):
    plt.plot(x, y_pack[pixel+tube*128]+tube*300)


# # fitting

# In[161]:


import scipy.optimize as sopt


# In[104]:


def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


# In[105]:


def gauss_with_bg(x, *p):
    A, mu, sigma, bg = p
    return bg + A*np.exp(-(x-mu)**2/(2.*sigma**2))


# In[121]:


def fit_pixel(x, y, dvalues, peaks_thispixel=None, chisq_thispixel=None):
    if peaks_thispixel is None:
        peaks_thispixel = np.zeros(len(dvalues), dtype=float)
    if chisq_thispixel is None:
        chisq_thispixel = np.zeros(len(dvalues), dtype=float)
    for pkindex, d0 in enumerate(dvalues):
        dmin = d0*.96
        dmax = d0*1.04
        startindex, stopindex = np.where(x>dmin)[0][0], np.where(x<dmax)[0][-1]
        x1 = x[startindex:stopindex]
        y1 = y[startindex:stopindex]
        # print x1
        # print y1
        guess_center = x1[np.argmax(y1)]
        startindex, stopindex = np.where(x>guess_center-0.05*d0)[0][0], np.where(x<guess_center+0.05*d0)[0][-1]
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


# In[132]:


maxchisq = 3.
min_counts = 500


# In[170]:


plt.figure(figsize=(7,4))
plt.plot(x, y_pack[100])


# In[165]:


dvalues


# In[191]:


dvalues2 = dvalues[4:]
dvalues2


# In[192]:


peaks = fit_pixel(x, y_pack[100], dvalues=dvalues2)


# In[193]:


peaks


# In[194]:


print fit_pixel(x, y_pack[20], dvalues2)


# In[195]:


plt.figure(figsize=(7,4))
plt.plot(x, y_pack[20])


# In[196]:


print fit_pixel(x, y_pack[0], dvalues2)


# ## Loop over all pixels in C25T

# In[174]:


get_ipython().run_cell_magic(u'time', u'', u'ws = msa.LoadEventNexus(nxspath, FilterByTimeStart=0, FilterByTimeStop=1)')


# In[176]:


N = ws.getNumberHistograms()


# In[197]:


newdifc = difc.copy()
mask = np.zeros(difc.shape, dtype=bool)
chisq = np.zeros((len(dvalues2),) + difc.shape, dtype=float)
peaks = np.zeros((len(dvalues2),) + difc.shape, dtype=float)


# In[198]:


firstindex, lastindex = startindex, endindex


# In[199]:


for ipixel, pixel in enumerate(range(firstindex, lastindex+1)):
    if (pixel%100)==0 : print pixel
    # print pixel
    y = y_pack[ipixel]
    peaks_thispixel = peaks[:, pixel]
    chisq_thispixel = chisq[:, pixel]
    fit_pixel(x,y, dvalues2, peaks_thispixel, chisq_thispixel)
    Ngoodfit = np.sum(np.isfinite(peaks_thispixel))
    # print Ngoodfit
    if not Ngoodfit:
        mask[pixel] = 1
        continue
    ratios = peaks_thispixel / dvalues2
    # print "ratios=", ratios
    good_ratios = ratios[np.isfinite(ratios)]
    # print "good_ratios=", good_ratios
    average_ratio = np.average(good_ratios)
    newdifc[pixel] = difc[pixel] * average_ratio
    # print average_ratio
    # break
    continue


# In[208]:


plt.figure(figsize=(7,4))
plt.plot(x, y_pack[128*2+127])


# In[205]:


print fit_pixel(x, y_pack[128*2+127], dvalues2)


# In[207]:


8.282/8.175


# In[202]:


tube=2
tube_slice = slice(firstindex+128*tube, firstindex+128*(tube+1))
newdifc[tube_slice]/difc[tube_slice]


# In[209]:


mask[tube_slice]


# # Save

# In[210]:


np.save('./C60-difc-2-C25T.npy', newdifc)
np.save('./C60-difc-2-C25T-mask.npy', mask)

"""
