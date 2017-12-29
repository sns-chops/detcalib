# coding: utf-8

import os, numpy as np
from mantid import simpleapi as msa, mtd

workdir = "/SNS/users/lj7/dv/sns-chops/detcalib/SEQ"
os.chdir(workdir)

# ## Compute nominal difc
nxspath = '/SNS/SEQ/IPTS-19573/nexus/SEQ_130249.nxs.h5'
ws = msa.LoadEventNexus(nxspath, FilterByTimeStart=0, FilterByTimeStop=1)
msa.LoadInstrument(ws, Filename='./SEQUOIA_Definition_guessshortpacks.xml', RewriteSpectraMap=False)
difc = msa.CalculateDIFC(InputWorkspace=ws)
difc = difc.extractY().flatten().copy()
msa.DeleteWorkspace('difc')
np.save('difc-guessshortpacks.npy', difc)

# get det ID list
detIDs = []
for i in range(ws.getNumberHistograms()):
    sp = ws.getSpectrum(i)
    dets = list(sp.getDetectorIDs())
    assert len(dets)==1
    detIDs.append(dets[0])
    continue
for i in range(len(detIDs)-1):
    assert detIDs[i] < detIDs[i+1]

# # Get pack index
def getPackIndexes(ws, name="C26T/eightpack-top"):
    instrument = ws.getInstrument()
    pack = instrument.getComponentByName(name)
    firstpixel = pack[0][0].getID()
    lasttube = pack[pack.nelements()-1]
    lastpixel = lasttube[lasttube.nelements()-1]
    lastpixel = lastpixel.getID()
    print "first and last pixel ID:", firstpixel, lastpixel
    return firstpixel, lastpixel
pack_indexes = dict()
packs = [
    #'C25T/eightpack-top',
    #'C26T/eightpack-top',
    'C25B/eightpack-bottom',
    'C26B/eightpack-bottom',
    ]
for name in packs:
    pack_indexes[name] = getPackIndexes(ws, name)
    continue
# clean up
msa.DeleteWorkspace('ws')


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
# time step for loading files
# too large will need too much memory
dt = 1000.


# d axis
dmin, dmax, delta_d=2., 11., 0.02
Nd = int((dmax-dmin)/delta_d)
print "Number of d bins:", Nd

# pixels per pack
Npixels = 1024

#
Npacks = len(packs)

y_matrix = np.zeros((Npacks, Npixels, Nd))
xbb_saved = None
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
        
        # loop over packs
        for ipack, packname in enumerate(packs):
            firstpixel, lastpixel = pack_indexes[packname]
            startindex = detIDs.index(firstpixel)
            endindex = detIDs.index(lastpixel)
            print "array indexes of first and last pixel", startindex, endindex

            y_pack = y_matrix[ipack]
            # loop over pixels in the pack
            for i, pixelindex in enumerate(range(startindex, endindex+1)):
                I_d_pixel = msa.SumSpectra(
                    InputWorkspace=I_d, StartWorkspaceIndex=pixelindex, EndWorkspaceIndex=pixelindex)
                xbb = I_d_pixel.readX(0)
                if xbb_saved is None: xbb_saved = np.array(xbb, copy=True)
                y = I_d_pixel.readY(0)
                y_pack[i] += y
                msa.DeleteWorkspace('I_d_pixel')
                continue
            continue
        
        msa.DeleteWorkspaces(['ws', 'I_d'])
        continue
    continue

xbb = np.arange(dmin, dmax+delta_d/2., delta_d)
np.save("C60-I_d-xbb.npy", xbb)
np.save("C60-I_d-y_matrix.npy", y_matrix)

for ipack, packname in enumerate(packs):
    y_pack = y_matrix[ipack]
    packname1 = packname.split('/')[0] # "C25T"
    np.save("C60-I_d-y-%s.npy" % packname1, y_pack)
    continue
