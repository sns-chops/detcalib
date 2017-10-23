#!/usr/bin/env python
# coding: utf-8

from mantid import simpleapi as msa, mtd
import numpy as np, os, math, collections
import scipy.optimize as sopt

import align


# ## Load difc
detID = np.load('./difc-2-detID.npy')

# for short packs: use C60 data for just C26B
difc = np.load("./C60-difc-2-C26B.npy")
mask = np.load("./C60-difc-2-C26B-mask.npy")

# ## Load L2
L2_calib = msa.Load('L2table.nxs')
L2_calib = msa.SortTableWorkspace(L2_calib, Columns='detid')
L2 = np.array(L2_calib.column('L2'))
L2_detID = np.array(L2_calib.column('detid'))
assert np.all(detID==L2_detID)
nodata = np.array(L2_calib.column('nodata'), dtype=bool)
mask = np.logical_or(mask, nodata)

# ## Apply mask
difc = np.ma.masked_array(difc, mask)
L2 = np.ma.masked_array(L2, mask)

# ## Fit source, sample
import logging
logger = logging.getLogger("Align component")
logger.setLevel(logging.DEBUG)
# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# add formatter to ch
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)

#
eulerConvention='YZX'
wks_name = "alignedWorkspace"
# use C60 data for short packs around forward beam
idf_orig = os.path.join('SEQUOIA_Definition_guessshortpacks.xml')


# ## Fit twotheta and L2
# $ DIFC = (L1+L2)/(\pi) \; sin(\theta) \times 0.0015882549421289758 \times 10^6 $
#L1 = 20.0965
L1 = 20.0114
sin_theta = difc/(0.0015882549421289758*1e6) * np.pi / (L1+L2)

msa.LoadEmptyInstrument(idf_orig, OutputWorkspace=wks_name)
instrument_model = align.InstrumentModel(wks_name, detID, mask, eulerConvention="YZX")
options = collections.OrderedDict()
options['Xposition']=(-.3, .3)
options['Yposition']=False
options['Zposition']=(-.3, .3)
options['AlphaRotation']=(-2., 2.)
options['BetaRotation']=False
options['GammaRotation']=False

crow_around_forward_beam = ['C25T', 'C26T', 'C25B', 'C26B']
packs = ['C26B']

# most of the packs are of type "eightpack"
# the short packs are special
pack_types = collections.defaultdict(lambda: 'eightpack')
pack_types['C25T' ] = pack_types['C26T'] = 'eightpack-top'
pack_types['C25B' ] = pack_types['C26B'] = 'eightpack-bottom'

template = """
  <type name="{0}">
    <component type="{1}">
      <location x="{2:.8f}" y="{3:.8f}" z="{4:.8f}">
         <rot axis-z="0" axis-x="0" axis-y="1" val="{5:.8f}"/>
      </location>
    </component>
  </type>
"""
comps = []
ofile = open('new.xml', 'wt')
for pack in packs:
    print "* Working on %s" % (pack,)
    pack_type = pack_types[pack]
    pack_model = instrument_model.component('%s/%s'%(pack, pack_type), type='detpack')
    print pack_model.getParams()
    init_center = pack_model.position()
    estimate = align.estimate_pack_center_position(sin_theta, L2, pack_model, init_center)
    # fit = align.FitPackTwothetaAndL2(pack_model, options, sin_theta, L2, logger)
    fit = align.FitPack_DifcL2(pack_model, options, difc, L2, logger)
    # fit = align.FitPackDifc(pack_model, options, difc, logger)
    try:
        fit.fit()
        comps.append(pack_model.fullname)
    except:
        print "Skipped %s" % pack
    x,y,z, ry,rz,rx = pack_model.getParams()
    print estimate
    s = template.format(pack,pack_type, x,y,z,ry)
    ofile.write(s)
    continue

# msa.ExportGeometry(InputWorkspace=wks_name,Components=comps,EulerConvention=eulerConvention,Filename='new.xml')
