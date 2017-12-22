#!/usr/bin/env python
# coding: utf-8

from mantid import simpleapi as msa, mtd
import numpy as np, os, math, collections, yaml
import align_utils

class Align:

    def __init__(
            self, 
            I_d_dir = None,
            logger = None,
            L1 = 20.0114, #L1 = 20.0965,
            eulerConvention='YZX',
            ):
        self.I_d_dir = I_d_dir
        if not logger: logger = createLogger()
        self.logger = logger
        self.detIDs = np.load(os.path.join(I_d_dir, 'detIDs.npy'))
        self.detID_list = list(self.detIDs)
        self.difc_all = np.load(os.path.join(I_d_dir, 'difc-nominal.npy'))
        self.L1 = L1
        self.init_IDF = os.path.join(I_d_dir, 'init_IDF.xml')
        self.eulerConvention = eulerConvention
        self.options = self.default_options()
        return

    def default_options(self):
        options = collections.OrderedDict()
        options['Xposition']=(-.3, .3)
        options['Yposition']=False
        options['Zposition']=(-.3, .3)
        options['AlphaRotation']=(-2., 2.)
        options['BetaRotation']=False
        options['GammaRotation']=False
        return options


    def load_L2_from_nxs(self, path):
        "load L2 computed using mantid and saved as nxs"
        L2_calib = msa.Load(path)
        L2_calib = msa.SortTableWorkspace(L2_calib, Columns='detid')
        L2 = np.array(L2_calib.column('L2'))
        L2_detID = np.array(L2_calib.column('detid'))
        assert np.all(self.detIDs==L2_detID)
        nodata = np.array(L2_calib.column('nodata'), dtype=bool)
        self.L2 = L2
        self.L2_mask = nodata
        return L2, nodata
    
    def align(
            self,
            difc, mask, pack = 'C25T',
            ofile = open('new.xml', 'wt'),
            logger = None,
            params0 = None,
            ):
        """
        difc, mask: for a pack
        pack: pack name
        """
        logger = logger or self.logger
        # ## Fit twotheta and L2
        # $ DIFC = (L1+L2)/(\pi) \; sin(\theta) \times 0.0015882549421289758 \times 10^6 $
        difc, mask = self._prepare_difc_and_mask_from_pack_difc(difc, mask, pack)
        # combine with L2 mask
        mask = np.logical_or(mask, self.L2_mask)
        #
        L1 = self.L1
        L2 = self.L2
        # import pdb; pdb.set_trace()
        # Apply mask
        difc = np.ma.masked_array(difc, mask)
        L2 = np.ma.masked_array(L2, mask)
        # calculate sin(theta)
        self.sin_theta = sin_theta = difc/(0.0015882549421289758*1e6) * np.pi / (L1+L2)
        # 
        wks_name = "alignedWorkspace"
        msa.LoadEmptyInstrument(self.init_IDF, OutputWorkspace=wks_name)
        instrument_model = align_utils.InstrumentModel(
            wks_name, self.detIDs, mask, eulerConvention=self.eulerConvention)
        #
        options = self.options
        #
        print "- Working on %s" % (pack,)
        pack_type = pack_types[pack]
        self.pack_model = pack_model = instrument_model.component('%s/%s'%(pack, pack_type), type='detpack')
        print "- pack params:", pack_model.getParams()
        init_center = pack_model.position()
        estimate = align_utils.estimate_pack_center_position(sin_theta, L2, pack_model, init_center)
        # fit = align_utils.FitPackTwothetaAndL2(pack_model, options, sin_theta, L2, logger)
        fit = align_utils.FitPack_DifcL2(pack_model, options, difc, L2, logger=logger, params0=params0)
        # fit = align_utils.FitPackDifc(pack_model, options, difc, logger=logger)
        fit.fit()
        print "- Estimate:", estimate
        x,y,z, ry,rz,rx = new_params = pack_model.getParams()
        print "- New:", new_params
        s = template.format(pack,pack_type, x,y,z,ry)
        print s
        ofile.write(s)
        return new_params, fit

    def _prepare_difc_and_mask_from_pack_difc(self, difc, mask, pack):
        "prepare full difc and mask array for all pixels from a pack difc array"
        packinfo = yaml.load(open(os.path.join(self.I_d_dir, 'pack-%s.yaml' % pack)))
        detID_list = self.detID_list
        pixelIDs = packinfo['pixelIDs']
        pack_slice = slice(detID_list.index(pixelIDs['first']), detID_list.index(pixelIDs['last']) + 1)
        difc1 = self.difc_all.copy(); difc1[pack_slice] = difc
        mask1 = np.zeros(difc1.shape, dtype=bool);
        mask1[pack_slice] = mask
        self.pack_difc = difc; self.pack_mask = mask
        self.pack_L2 = self.L2[pack_slice]
        self.pack_sin_theta = self.pack_difc/(0.0015882549421289758*1e6) * np.pi / (self.L1+self.pack_L2)
        return difc1, mask1


def createLogger():
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
    return logger
    
# map pack name to type
crow_around_forward_beam = ['C25T', 'C26T', 'C25B', 'C26B']
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

