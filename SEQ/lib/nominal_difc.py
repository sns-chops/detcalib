
from mantid import simpleapi as msa

def get_nominal_difc(nxspath, idfpath, outpath=None):
    ws = msa.LoadEventNexus(nxspath, FilterByTimeStart=0, FilterByTimeStop=1)
    msa.LoadInstrument(ws, Filename=idfpath, RewriteSpectraMap=False)
    difc = msa.CalculateDIFC(InputWorkspace=ws)
    difc = difc.extractY().flatten().copy()
    msa.DeleteWorkspace('difc')
    if outpath:
        np.save(outpath, difc)
    return difc
