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
        pack="D23/eightpack",
        dmin=0.4, dmax=7, dd=0.005,
        dvalues = None,
        I_tof_dir=os.path.join(workdir, '../Si-I_tof'),
        outdir = None
        ):
    if not os.path.exists(outdir): os.makedirs(outdir)
    packnameandtype = pack
    packname, packtype = pack.split('/')

    # old
    L2 = np.load(os.path.join(workdir, "L2_orig.npy"))
    difc = np.load(os.path.join(workdir, "difc_orig.npy"))
    x0,y0 = compute_I_d(pack, dmin, dmax, dd, dvalues, I_tof_dir, L2, difc)

    # new
    L2 = np.load(os.path.join(workdir, "L2.npy"))
    difc = np.load(os.path.join(workdir, "difc.npy"))
    x1,y1 = compute_I_d(pack, dmin, dmax, dd, dvalues, I_tof_dir, L2, difc)

    plt.figure(figsize=(7,4))
    plt.title("Pack %s" % packname)
    plt.plot(x0, y0, label='original')
    plt.plot(x1, y1, label='after loading new xml')
    # plt.plot(x0, y1-y0, label='diff')
    for d in dvalues: plt.axvline(x=d, linewidth=1, color='lightgray')
    # plt.xlim(3,3.3)
    plt.legend(loc='upper left')
    outpng = '%s-I_d.png' % packname
    plt.savefig(os.path.join(outdir, outpng) )
    plt.close()
    return


def compute_I_d(
        pack="D23/eightpack",
        dmin=0.4, dmax=7, dd=0.005,
        dvalues = None,
        I_tof_dir=os.path.join(workdir, '../Si-I_tof'),
        L2=None, difc=None,
        ):
    packnameandtype = pack
    packname, packtype = pack.split('/')
    from get_I_d_accounting_for_T0 import Get_I_d
    from calibrate import SEQ_T0_of_E
    # I(tof) was already calculated and saved in I_tof_dir
    # so we just need to read it, apply T0 correction
    # and then compute d using difc, which is also pre-calcualted.
    I_d_array = Get_I_d(
        I_tof_dir = I_tof_dir,
        pack = packname,
        T0_of_E = SEQ_T0_of_E, L1=20.0114,
        L2=L2, difc = difc,
    )()
    npixels, _two, npoints = I_d_array.shape
    assert _two == 2
    # get the spectrum for the whole pack
    _tmp = np.transpose(I_d_array, (1,0,2)).copy()
    _tmp.shape = 2, -1
    d, I = _tmp
    hist, bin_edges = np.histogram(d, bins=np.arange(dmin, dmax, dd), weights=I)
    x = (bin_edges[1:] + bin_edges[:-1])/2
    y = hist
    return x,y


import d_spacing

def test1():
    compare(
        pack="D19/eightpack",
        dmin=0.4, dmax=7, dd=0.005,
        dvalues = d_spacing.values['Si'],
        I_tof_dir=os.path.join(workdir, '../Si-I_tof')
        )
    return

def test2():
    compare(
        pack="D23/eightpack",
        dmin=2, dmax=11, dd=0.01,
        dvalues = d_spacing.values['C60'],
        I_tof_dir=os.path.join(workdir, '../C60-I_tof')
        )
    return
    # Si not good. why?
    compare(
        pack="D23/eightpack",
        dmin=0.4, dmax=7, dd=0.005,
        dvalues = d_spacing.values['Si'],
        I_tof_dir=os.path.join(workdir, '../Si-I_tof')
        )
    
packs = ['B%s' % n for n in range(1, 38)] \
        + ['D%s' % n for n in range(1, 38)] \
        + ['C%s' % n for n in range(1, 25)] \
        + ['C%s' % n for n in range(27, 38)] \
        + ['C25T', 'C25B', 'C26T', 'C26B']

def main():
    from pack2sample import pack2sample as p2s
    # packs = ['D27']
    # packs = packs[:10]
    for pack in packs:
        sample = p2s(pack)
        packtype = 'eightpack'
        if sample.endswith('B'): packtype += '-bottom'
        elif sample.endswith('T'): packtype += '-top'
        packname = pack
        pack = '%s/%s' % (packname, packtype)
        print pack
        kwds = dict(pack=pack)
        kwds.update(d_spacing.ranges[sample])
        kwds.update(dvalues=d_spacing.values[sample])
        kwds.update(I_tof_dir=os.path.join(workdir, '../%s-I_tof' % sample))
        kwds.update(outdir = os.path.join(workdir, 'out'))
        print kwds
        compare(**kwds)
        continue
    return

# if __name__ == "__main__": test1()
if __name__ == "__main__": main()
