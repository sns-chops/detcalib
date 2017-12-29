#!/usr/bin/env python

from mantid import simpleapi as msa
import os
from matplotlib import pyplot as plt

orig_I_d = msa.Load(os.path.expanduser("~/tmp/130273-I_d.nxs"))
I_d = msa.Load(os.path.expanduser("~/tmp/130273-newxml-using-both-difc-and-L2-I_d.nxs"))

start_index = 2048

orig_I_d_s = msa.SumSpectra(InputWorkspace=orig_I_d, StartWorkspaceIndex=start_index, EndWorkspaceIndex=start_index+1023)
I_d_s = msa.SumSpectra(InputWorkspace=I_d, StartWorkspaceIndex=start_index, EndWorkspaceIndex=start_index+1023)

orig_d_bb = orig_I_d_s.readX(0); orig_I = orig_I_d_s.readY(0)
d_bb = I_d_s.readX(0); I = I_d_s.readY(0)
plt.figure(figsize=(7,4))
plt.plot(orig_d_bb[:-1], orig_I, label='original')
plt.plot(d_bb[:-1], I, label='calibrated')
# plt.xlim(3,3.3)
plt.legend(loc='upper left')
plt.show()
