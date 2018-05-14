# coding: utf-8

import os, numpy as np
from mantid import simpleapi as msa, mtd

def get_I_tof(nxs_files, outdir, packs, dt = 1000., tofaxis=None, Npixels_per_pack=1024):
    """nxs_files: paths of calibration nxs files
    outdir: output directory
    packs: list of pack names, e.g. C26B/eightpack-bottom
    dt: time step for loading files. too large will need too much memory
    Npixels_per_pack: number of pixels per pack
    tofaxis: tofmin, tofmax, dtof

    Output files:
    * detIDs.npy
    * I_tof-xbb.npy
    * I_tof-y-PACKNAME.npy
    * pack-PACKNAME.yaml

    NOTE:
    * Different combinations of nxs_files, init_IDF, d_axis should use different outdirs
    """
    tofmin, tofmax, dtof = tofaxis
    if not os.path.exists(outdir): os.makedirs(outdir)
    # ## Compute nominal difc using first file in the list
    nxspath = nxs_files[0]
    ws = msa.LoadEventNexus(nxspath, FilterByTimeStart=0, FilterByTimeStop=1) # load just one second
    #
    # IDs of all pixels
    detIDs = getDetIDs(ws)
    np.save(os.path.join(outdir, 'detIDs.npy'), detIDs)
    #
    # map pack name to (start_pixelID, stop_pixelID)
    pack2pixelID_start_stop = dict()
    for name in packs:
        pack2pixelID_start_stop[name] = getFirstLastPixelIDs(ws, name)
        continue
    # get tof axis
    I_tof = msa.Rebin(InputWorkspace=ws, Params='%s,%s,%s' % (tofmin, dtof, tofmax))
    I_tof = msa.SumSpectra(InputWorkspace=I_tof)
    xbb = np.array(I_tof.readX(0), copy=True)
    print xbb[0], xbb[-1], len(xbb)
    # clean up
    msa.DeleteWorkspaces(['ws', 'I_tof'])

    runtimes = dict()
    for f in nxs_files:
        runtimes[f] = getRunTime(f)
    print "* run times:", runtimes

    Ntof = len(xbb) - 1
    print "* Number of TOF bins:", Ntof

    #
    Npacks = len(packs)

    y_matrix = np.zeros((Npacks, Npixels_per_pack, Ntof))
    for nxsfile in nxs_files:
        print "* Working on", nxsfile
        t_total = runtimes[nxsfile]
        for tstart in np.arange(0, t_total-dt, dt):
            print "* tstart", tstart
            tend = min(t_total-1, tstart+dt)
            ws = msa.LoadEventNexus(nxsfile, FilterByTimeStart=tstart, FilterByTimeStop=tend)
            I_tof = msa.Rebin(InputWorkspace=ws, Params='%s,%s,%s' % (tofmin, dtof, tofmax))

            # loop over packs
            for ipack, packname in enumerate(packs):
                firstpixel, lastpixel = pack2pixelID_start_stop[packname]
                startindex = detIDs.index(firstpixel)
                endindex = detIDs.index(lastpixel)
                print "array indexes of first and last pixel", startindex, endindex

                y_pack = y_matrix[ipack]
                # loop over pixels in the pack
                for i, pixelindex in enumerate(range(startindex, endindex+1)):
                    I_tof_pixel = msa.SumSpectra(
                        InputWorkspace=I_tof, StartWorkspaceIndex=pixelindex, EndWorkspaceIndex=pixelindex)
                    y = I_tof_pixel.readY(0)
                    y_pack[i] += y
                    msa.DeleteWorkspace('I_tof_pixel')
                    continue
                continue

            msa.DeleteWorkspaces(['ws', 'I_tof'])
            continue
        continue

    #xbb = np.arange(tofmin, tofmax+dtof/2., dtof)
    # print xbb
    np.save(os.path.join(outdir, "I_tof-xbb.npy"), xbb)
    # for debugging
    np.save(os.path.join(outdir, "I_tof-y_matrix.npy"), y_matrix)

    for ipack, packname in enumerate(packs):
        y_pack = y_matrix[ipack]
        packname1 = packname.split('/')[0] # "C25T"
        # save y values of I(d) for the pack
        np.save(os.path.join(outdir, "I_tof-y-%s.npy" % packname1), y_pack)
        # save pack info
        first, last = pack2pixelID_start_stop[packname]
        pixelIDs = dict(first=first, last=last)
        pack_info = dict(pixelIDs=pixelIDs)
        dumpYaml(pack_info, os.path.join(outdir, 'pack-%s.yaml' % packname1))
        continue
    return


def dumpYaml(d, path):
    import yaml
    with open(path, 'wt') as ostream:
        yaml.dump(d, ostream, default_flow_style=False)
    return


def getDetIDs(ws):
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
    return detIDs


## Get detectorID of the first and last pixels of a pack
def getFirstLastPixelIDs(ws, name="C26T/eightpack-top"):
    instrument = ws.getInstrument()
    pack = instrument.getComponentByName(name)
    if not pack:
        raise RuntimeError("Failed to get Component %s" % (name,))
    firstpixel = pack[0][0].getID()
    lasttube = pack[pack.nelements()-1]
    lastpixel = lasttube[lasttube.nelements()-1]
    lastpixel = lastpixel.getID()
    # print "first and last pixel ID:", firstpixel, lastpixel
    return firstpixel, lastpixel


# # Get run times of all data files
def getRunTime(p):
    __ws = msa.LoadEventNexus(Filename=p, FilterByTimeStart=0, FilterByTimeStop=0)
    run = __ws.getRun()
    t = (run.endTime() - run.startTime()).total_seconds()
    msa.DeleteWorkspace('__ws')
    return t
