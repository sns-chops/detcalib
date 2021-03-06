# coding: utf-8

import os, numpy as np
from mantid import simpleapi as msa, mtd
from get_I_tof import getDetIDs, getFirstLastPixelIDs, getRunTime, dumpYaml


def get_I_d(nxs_files, init_IDF, outdir, packs, dt = 1000., d_axis=(2.,11.,0.02), Npixels_per_pack=1024):
    """nxs_files: paths of calibration nxs files
    init_IDF: initial IDF path
    outdir: output directory
    packs: list of pack names, e.g. C26B/eightpack-bottom
    dt: time step for loading files. too large will need too much memory
    d_axis: dmin, dmax, delta_d. e.g. 2., 11., 0.02
    Npixels_per_pack: number of pixels per pack

    Output files:
    * difc-nominal.npy
    * detIDs.npy
    * I_d-xbb.npy
    * I_d-y-PACKNAME.npy
    * pack-PACKNAME.yaml

    NOTE:
    * Assumed that the difc array from CalculateDIFC is ordered according to the "spectrrum list" in
      the mantid workspace. See function getDetIDs
    * Different combinations of nxs_files, init_IDF, d_axis should use different outdirs
    """
    if not os.path.exists(outdir): os.makedirs(outdir)
    # ## Compute nominal difc using first file in the list
    nxspath = nxs_files[0]
    ws = msa.LoadEventNexus(nxspath, FilterByTimeStart=0, FilterByTimeStop=1) # load just one second
    #
    msa.LoadInstrument(ws, Filename=init_IDF, RewriteSpectraMap=False)
    import shutil
    shutil.copyfile(init_IDF, os.path.join(outdir, 'init_IDF.xml'))
    #
    difc = msa.CalculateDIFC(InputWorkspace=ws)
    difc = difc.extractY().flatten().copy()
    msa.DeleteWorkspace('difc')
    np.save(os.path.join(outdir, 'difc-nominal.npy'), difc)
    # IDs of all pixels
    detIDs = getDetIDs(ws)
    np.save(os.path.join(outdir, 'detIDs.npy'), detIDs)
    #
    # map pack name to (start_pixelID, stop_pixelID)
    pack2pixelID_start_stop = dict()
    for name in packs:
        pack2pixelID_start_stop[name] = getFirstLastPixelIDs(ws, name)
        continue
    # clean up
    msa.DeleteWorkspace('ws')

    runtimes = dict()
    for f in nxs_files:
        runtimes[f] = getRunTime(f)
    print "* run times:", runtimes

    dmin, dmax, delta_d = d_axis
    Nd = int((dmax-dmin)/delta_d)
    print "* Number of d bins:", Nd

    #
    Npacks = len(packs)

    y_matrix = np.zeros((Npacks, Npixels_per_pack, Nd))
    xbb_saved = None
    for nxsfile in nxs_files:
        print "* Working on", nxsfile
        t_total = runtimes[nxsfile]
        for tstart in np.arange(0, t_total-dt, dt):
            print "* tstart", tstart
            tend = min(t_total-1, tstart+dt)
            ws = msa.LoadEventNexus(nxsfile, FilterByTimeStart=tstart, FilterByTimeStop=tend)
            msa.LoadInstrument(ws, Filename=init_IDF, RewriteSpectraMap=False)
            I_d = msa.ConvertUnits(InputWorkspace=ws, Target='dSpacing', EMode='Elastic')
            I_d = msa.Rebin(InputWorkspace=I_d, Params='%s,%s,%s' % (dmin, delta_d, dmax))

            # loop over packs
            for ipack, packname in enumerate(packs):
                firstpixel, lastpixel = pack2pixelID_start_stop[packname]
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
    np.save(os.path.join(outdir, "I_d-xbb.npy"), xbb)
    # for debugging
    np.save(os.path.join(outdir, "I_d-y_matrix.npy"), y_matrix)

    for ipack, packname in enumerate(packs):
        y_pack = y_matrix[ipack]
        packname1 = packname.split('/')[0] # "C25T"
        # save y values of I(d) for the pack
        np.save(os.path.join(outdir, "I_d-y-%s.npy" % packname1), y_pack)
        # save pack info
        first, last = pack2pixelID_start_stop[packname]
        pixelIDs = dict(first=first, last=last)
        pack_info = dict(pixelIDs=pixelIDs)
        dumpYaml(pack_info, os.path.join(outdir, 'pack-%s.yaml' % packname1))
        continue
    return



