import os, numpy as np
import get_difc_from_Itof


def T0_E(E, bg, A, power, E_decay):
    return bg + A*E**power*np.exp(-E/E_decay)
SEQ_T0_of_E = lambda E: T0_E(E, bg=0, A=54.881, power=-0.21413, E_decay=194.4)
# SEQ_T0_of_E = lambda E: T0_E(E, bg=-330, A=54.881*26, power=-0.21413, E_decay=194.4)


def load_L2_from_nxs(path):
    "load L2 computed using mantid and saved as nxs"
    from mantid import simpleapi as msa
    L2_calib = msa.Load(path)
    L2_calib = msa.SortTableWorkspace(L2_calib, Columns='detid')
    L2 = np.array(L2_calib.column('L2'))
    L2_detID = np.array(L2_calib.column('detid'))
    nodata = np.array(L2_calib.column('nodata'), dtype=bool)
    return L2

def calibrate(
        pack = 'D21/eightpack', I_tof_dir='Si-I_tof', peak_fractional_width=0.02,
        dvalues = [1.10860231, 1.24596143, 1.357755, 1.63751414, 1.92015553, 3.13560085],
        dmin=2.5, dmax=3.5,
        maxchisq = 5., min_counts = 800,
        T0_of_E = SEQ_T0_of_E,
        l2table_nxs = './L2table.nxs',
        geometrical_constraints = None,
        align = True,
):
    """
    geometrical_constraints: ex. dict(Xposition=(-0.005, 0.005))
    """
    packname, packtype = pack.split('/')
    fitter = get_difc_from_Itof.Fitter(
        peak_fractional_width=peak_fractional_width,
        bg_type='linear', t0_range=(0, 0.01))
    
    L2 = load_L2_from_nxs(l2table_nxs)

    detIDs = np.load(os.path.join(I_tof_dir, 'detIDs.npy'))
    detID_list = list(detIDs)
    
    import yaml
    packinfo = yaml.load(open(os.path.join( I_tof_dir, 'pack-%s.yaml' % packname)))
    firstpixelID = packinfo['pixelIDs']['first']
    firstpixel_index = detID_list.index(firstpixelID)
    L2_pack = L2[firstpixel_index: firstpixel_index+1024]

    gpd = get_difc_from_Itof.GetPackDifc(
        pack=packname,
        dvalues=dvalues,
        dmin=dmin, dmax=dmax,
        I_tof_dir = I_tof_dir,
        fitter=fitter,
        maxchisq = maxchisq,
        min_counts = min_counts,
        T0_of_E=T0_of_E,
        L2 = L2_pack,
    )

    difc, mask, signature_d = gpd()
    np.save(os.path.join(I_tof_dir, 'difc-%s.npy' % packname), difc)
    np.save(os.path.join(I_tof_dir, 'mask-%s.npy' % packname), mask)

    if align:
        import align
        alignment = align.Align(I_tof_dir)
        if geometrical_constraints:
            alignment.options.update(geometrical_constraints)
        alignment.load_L2_from_nxs(l2table_nxs)
        alignment.align(difc, mask, packname, ofile=open('new-%s.xml' % packname, 'wt'))
    return difc, mask
