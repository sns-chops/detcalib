#!/usr/bin/env python

"""compute difc and L2 from IDF. It does not matter which nxs file is used.
What matters is the IDF file.
"""

import os
import sys

here = os.path.dirname(__file__) or '.'
sys.path.insert(0, os.path.join(here, "../lib"))

import from_IDF
from_IDF.get_difc(
    nxspath="/SNS/SEQ/IPTS-19573/nexus/SEQ_130273.nxs.h5",  # Si
    outpath='difc_orig.npy')

from_IDF.get_difc(
    nxspath="/SNS/SEQ/IPTS-19573/nexus/SEQ_130273.nxs.h5",  # Si
    idfpath='./SEQUOIA_Definition.xml',
    outpath='difc.npy')

# original, non-calibrated L2
from_IDF.get_L2(
    nxspath="/SNS/SEQ/IPTS-19573/nexus/SEQ_130273.nxs.h5",  # Si
    outpath='L2_orig.npy')

# calibrated L2
from_IDF.get_L2(
    nxspath="/SNS/SEQ/IPTS-19573/nexus/SEQ_130273.nxs.h5",  # Si
    idfpath='./SEQUOIA_Definition.xml',
    outpath='L2.npy')
