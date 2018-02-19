#!/usr/bin/env python

import os
import sys

here = os.path.dirname(__file__) or '.'
sys.path.insert(0, os.path.join(here, "../lib"))

workdir = here
os.chdir(workdir)

# setup plotting
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np, os, sys

# mantid
from mantid import simpleapi as msa, mtd

def compare(
        pack="C25B/eightpack-bottom",
        nxspath="/SNS/SEQ/IPTS-19573/nexus/SEQ_130249.nxs.h5", #C60
        detIDs_npy = '../C60-I_d/detIDs.npy',
        newIDF='./SEQUOIA_Definition.xml',
        dmin=2, dmax=11, dd=0.01,
        dvalues = None,
        tmin=0, tmax=2000
        ):
    orig_ws = msa.LoadEventNexus(Filename=nxspath, FilterByTimeStart=tmin, FilterByTimeStop=tmax)
    
    ws = orig_ws
    instrument = ws.getInstrument()
    packnameandtype = pack
    packname, packtype = pack.split('/')
    pack = instrument.getComponentByName(packnameandtype)
    firstpixel = pack[0][0].getID()
    lasttube = pack[pack.nelements()-1]
    lastpixel = lasttube[lasttube.nelements()-1]
    lastpixel = lastpixel.getID()
    print "first and last pixel IDs:", firstpixel, lastpixel
    #
    #
    detIDs = list(np.load(detIDs_npy))
    startindex = detIDs.index(firstpixel)
    endindex = detIDs.index(lastpixel)
    print "first and last pixel indexes:", startindex, endindex
    del ws

    # # Old I(d)
    daxis = "%s,%s,%s" % (dmin, dd, dmax)
    I_d_0 = msa.ConvertUnits(InputWorkspace=orig_ws, Target='dSpacing', EMode='Elastic')
    I_d_0 = msa.Rebin(InputWorkspace=I_d_0, Params=daxis)
    pack_I_d_0 = msa.SumSpectra(InputWorkspace=I_d_0, StartWorkspaceIndex=startindex, EndWorkspaceIndex=endindex)
    xbb0 = pack_I_d_0.readX(0); y0 = pack_I_d_0.readY(0).copy()
    x0 = (xbb0[1:] + xbb0[:-1])/2
    msa.DeleteWorkspace(I_d_0)
    msa.DeleteWorkspace(pack_I_d_0)

    # # New I(d)
    msa.LoadInstrument(orig_ws, Filename=newIDF, RewriteSpectraMap=False)
    I_d_1 = msa.ConvertUnits(InputWorkspace=orig_ws, Target='dSpacing', EMode='Elastic')
    I_d_1 = msa.Rebin(InputWorkspace=I_d_1, Params=daxis)
    pack_I_d_1 = msa.SumSpectra(InputWorkspace=I_d_1, StartWorkspaceIndex=startindex, EndWorkspaceIndex=endindex)
    xbb1 = pack_I_d_1.readX(0); y1 = pack_I_d_1.readY(0).copy()
    x1 = (xbb1[1:] + xbb1[:-1])/2
    msa.DeleteWorkspace(I_d_1)
    msa.DeleteWorkspace(pack_I_d_1)
    msa.DeleteWorkspace(orig_ws)

    data = [x0, y0, x1, y1]
    np.save("%s-I_d.npy" % packname, data)
    plt.figure(figsize=(7,4))
    plt.title("Pack %s" % packname)
    plt.plot(x0, y0, label='original')
    plt.plot(x1, y1, label='after loading new xml')
    for d in dvalues: plt.axvline(x=d, linewidth=1, color='k')
    # plt.xlim(3,3.3)
    plt.legend(loc='upper left')
    outpng = '%s-I_d.png' % packname
    plt.savefig(outpng)
    return


import d_spacing

def test1():
    compare()
    return

def test2():
    compare(
        pack="D23/eightpack",
        nxspath="/SNS/SEQ/IPTS-19573/nexus/SEQ_130273.nxs.h5",  # Si
        detIDs_npy = '../Si-I_d/detIDs.npy',
        newIDF='./SEQUOIA_Definition.xml',
        dmin=0.4, dmax=7, dd=0.005,
        dvalues = d_spacing.values['Si']
        )

def test3():
    compare(
        pack="D27/eightpack",
        nxspath="/SNS/SEQ/IPTS-19573/nexus/SEQ_130273.nxs.h5",  # Si
        detIDs_npy = '../Si-I_d/detIDs.npy',
        newIDF='./SEQUOIA_Definition.xml',
        dmin=0.6, dmax=3.45, dd=0.005,
        dvalues = d_spacing.values['Si']
        )

if __name__ == "__main__": test3()
