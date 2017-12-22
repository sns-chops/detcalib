from mantid import simpleapi as msa, mtd
import numpy as np, os, math
import scipy.optimize as sopt

def opts2fullparams(x, params0, options):
    """ 
    Creates an array combining the refining and constant variables
    This is required because scipy.optimise.minimise expect a constant
    number of variable, so need to be able to maps any number of      
    inputs to six outputs. 
    """
    x0_index = 0
    out = []
    for opt, p0 in zip(options, params0):
        if options[opt]:
            v = x[x0_index]
            x0_index += 1
        else:
            v = p0
        out.append(v)
    return out


def eulerToQuat(alpha, beta, gamma, convention):
    """                    
    Convert Euler angles to a quaternion
    """
    from mantid.kernel import Quat, V3D
    getV3D = {'X': V3D(1, 0, 0), 'Y': V3D(0, 1, 0), 'Z': V3D(0, 0, 1)}
    return (Quat(alpha, getV3D[convention[0]]) * Quat(beta, getV3D[convention[1]]) *
            Quat(gamma, getV3D[convention[2]]))

def eulerToAngleAxis(alpha, beta, gamma, convention):
    """
    Convert Euler angles to a angle rotation around an axis
    """
    quat = eulerToQuat(alpha, beta, gamma, convention)
    if quat[0] == 1:
        return 0, 0, 0, 1
    deg = math.acos(quat[0])
    scale = math.sin(deg)
    deg *= 360.0 / math.pi
    ax0 = quat[1] / scale
    ax1 = quat[2] / scale
    ax2 = quat[3] / scale
    return deg, ax0, ax1, ax2


class InstrumentModel(object):
    
    def __init__(self, wks_name, detID, mask, eulerConvention="YZX"):
        self.wks_name = wks_name
        self.detID = detID
        self.detID_list = list(detID)
        self.mask = mask
        self.eulerConvention = eulerConvention
        return
    
    def component(self, name, type='component'):
        klsname = "%sModel" % type.capitalize()
        kls = eval(klsname)
        return kls(self, name)
    
    def adjustComponent(self, component, pos, rot):
        wks_name = self.wks_name
        if pos:
            msa.MoveInstrumentComponent(wks_name, component, X=pos[0], Y=pos[1], Z=pos[2], RelativePosition=False)

        if rot:
            (rotw, rotx, roty, rotz) = eulerToAngleAxis(rot[0], rot[1], rot[2], self.eulerConvention)
            msa.RotateInstrumentComponent(
                wks_name, component, X=rotx, Y=roty, Z=rotz, Angle=rotw, RelativeRotation=False)
        return
    
    def difc(self):
        wks_name = self.wks_name
        msa.CalculateDIFC(InputWorkspace=wks_name, OutputWorkspace='difc')
        return mtd['difc'].extractY().flatten()

    def L2(self):
        wks_name = self.wks_name
        pp = msa.PreprocessDetectorsToMD(InputWorkspace=wks_name, OutputWorkspace='pp')
        pp = msa.SortTableWorkspace(pp, Columns='DetectorID')
        L2 = np.array(pp.column('L2'))
        # pp may contain monitors
        nmon = (np.array(pp.column('DetectorID'))<=0).sum()
        assert nmon + self.detID.size == L2.size
        return L2[nmon:]

    def twotheta_and_L2(self):
        wks_name = self.wks_name
        pp = msa.PreprocessDetectorsToMD(InputWorkspace=wks_name, OutputWorkspace='pp')
        pp = msa.SortTableWorkspace(pp, Columns='DetectorID')
        tt, L2 = np.array( pp.column('TwoTheta') ), np.array(pp.column('L2'))
        # pp may contain monitors
        nmon = (np.array(pp.column('DetectorID'))<=0).sum()
        assert nmon + self.detID.size == tt.size
        return tt[nmon:], L2[nmon:]
    
    
class ComponentModel(object):
    
    def __init__(self, instrument_model, name):
        self.name = name
        self.instrument_model = instrument_model
        comp = self.getComponent()
        self.fullname = comp.getFullName()
        self.init_pos = self.position()
        self.init_rot = self.eulerAngles()
        return
    
    def adjust(self, pos, rot):
        return self.instrument_model.adjustComponent(self.fullname, pos, rot)
    
    def position(self):
        pos = self.comp.getPos()
        return [pos.getX(), pos.getY(), pos.getZ()]
    
    def eulerAngles(self):
        comp = self.comp
        eulerConvention = self.instrument_model.eulerConvention
        res = comp.getRotation().getEulerAngles(eulerConvention)
        return list(res)
    
    def getComponent(self):
        wks_name = self.instrument_model.wks_name
        instrument = mtd[wks_name].getInstrument()
        # self.comp = instrument.getComponentByName(self.name)
        self.comp = getattr(instrument, 'get%s' % self.name.capitalize())()
        return self.comp

    def getParams(self):
        # need to grab the component again
        self.getComponent()
        return self.position() + self.eulerAngles()


class FitSourceSample(object):
    
    def __init__(self, component_model, options, difc, logger=None):
        self.comp_model = component_model
        self.params0 = component_model.getParams()
        self.options = options
        self.move = True
        self.rotate = False
        self.difc = difc
        self.logger = logger
        return
    
    def adjust_model(self, params):
        pos = params[:3]
        rot = False
        self.comp_model.adjust(pos, rot)
        return
    
    def fit(self):
        comp_model = self.comp_model
        options = self.options
        logger = self.logger
        logger.info(
            "Working on %s. Starting position is %s, Starting rotation is %s" % (
                comp_model.fullname, comp_model.position(), comp_model.eulerAngles())
            )
        x0List = []
        boundsList = []
        params0 = self.params0
        for iopt,opt in enumerate(options.keys()):
            if not options[opt]: continue
            min, max = options[opt]
            x0List.append(params0[iopt])
            boundsList.append((params0[iopt]+min, params0[iopt]+max))            

        results = sopt.minimize(
        # results = sopt.differential_evolution(
            self.cost, 
            bounds=boundsList,
            x0=x0List,
            method='L-BFGS-B',
            )

        # Apply the results to the output workspace
        params = opts2fullparams(results.x, params0, options)
        self.adjust_model(params)
        # Need to grab the component again, as things have changed   
        comp_model.getComponent()
        logger.info(
            "Finished optimizing %s. Final position is %s, Final rotation is %s" % (
                comp_model.fullname, comp_model.position(), comp_model.eulerAngles())
            )        
        return
    
    
    def cost(self, x):
        params0 = self.params0
        options = self.options
        params = opts2fullparams(x, params0, options)
        self.adjust_model(params)
        wks_name = self.comp_model.instrument_model.wks_name
        msa.CalculateDIFC(InputWorkspace=wks_name, OutputWorkspace='difc')
        difc_new = mtd['difc'].extractY().flatten()
        difc_new = np.ma.masked_array(difc_new, self.comp_model.instrument_model.mask)
        difc = self.difc
        return np.sum((difc_new-difc)**2)
    
class DetpackModel(ComponentModel):
    
    def __init__(self, instrument_model, name):
        super(DetpackModel, self).__init__(instrument_model, name)
        detID_list = instrument_model.detID_list
        comp = self.comp
        firstDetID = getFirstDetID(comp)
        firstIndex = detID_list.index(firstDetID)
        lastDetID = getLastDetID(comp)
        lastIndex = detID_list.index(lastDetID)
        if lastDetID - firstDetID != lastIndex - firstIndex:
            raise RuntimeError("Calibration detid doesn't match instrument")
        self.firstIndex = firstIndex
        self.lastIndex = lastIndex
        self.mask = instrument_model.mask[firstIndex:lastIndex + 1]
        return
    
    def difc(self):
        return np.ma.masked_array(
            self.instrument_model.difc()[self.firstIndex:self.lastIndex+1], self.mask)

    def L2(self):
        return np.ma.masked_array(
            self.instrument_model.L2()[self.firstIndex:self.lastIndex+1], self.mask)

    def twotheta_and_L2(self):
        tt, L2 = self.instrument_model.twotheta_and_L2()
        return (np.ma.masked_array(tt[self.firstIndex:self.lastIndex+1], self.mask),
                np.ma.masked_array(L2[self.firstIndex:self.lastIndex+1], self.mask))
    
    def getComponent(self):
        wks_name = self.instrument_model.wks_name
        instrument = mtd[wks_name].getInstrument()
        self.comp = instrument.getComponentByName(self.name)
        return self.comp


class FitRow(object):

    """fit a row of packs"""
    
    def __init__(self, difc, pack_models, logger=None):
        self.difc = difc
        self.pack_models = pack_models
        self.original_positions = dict()
        for pm in pack_models:
            self.original_positions[pm] = pm.position()
            continue
        self.logger = logger
        return
    
    def adjust_model(self, vertical_offset):
        for pm in self.pack_models:
            x,y,z = self.original_positions[pm]
            pos = map(float, [x,y+vertical_offset,z])
            print vertical_offset, pos
            pm.adjust(pos, rot=None)
            continue
        return

    def cost(self, dy):
        self.adjust_model(dy)
        r = 0.
        for pm in self.pack_models:
            pack_difc = self.difc[pm.firstIndex:pm.lastIndex+1]
            difc_new = pm.difc()
            r += np.sum((difc_new-pack_difc)**2)# /difc_new.size
            continue
        print r/len(self.pack_models)
        return r
    
    def fit(self, bound=(-0.05, 0.05)):
        pack_models = self.pack_models
        min,max = bound
        bounds = (bound,)
        results = sopt.minimize(
            self.cost, 
            bounds=bounds,
            x0=[0.],
            # method='L-BFGS-B',
            )
        return results


class FitBundleXZ(object):

    """fit a bundle of packs"""
    
    def __init__(self, difc, L2, pack_models, min_pack_distance=0.22, y_offset=0, logger=None):
        """pack_models must be in a sequence
        """
        self.difc = difc
        self.L2 = L2
        self.pack_models = pack_models
        self.original_positions = dict()
        self.original_rotations = dict()
        for pm in pack_models:
            self.original_positions[pm] = pm.position()
            self.original_rotations[pm] = pm.eulerAngles()
            continue
        self.min_pack_distance = min_pack_distance
        self.y_offset = y_offset
        self.logger = logger
        return
    
    def adjust_model(self, xzry_arr):
        for i, pm in enumerate(self.pack_models):
            x,y,z = self.original_positions[pm]
            dx1,dz1,dry = xzry_arr[3*i: 3*(i+1)]
            newpos = map(float, [x+dx1,y+self.y_offset,z+dz1])
            r1,r2,r3 = self.original_rotations[pm]
            newrot = map(float, [r1+dry, r2, r3])
            print "%s, %s\n -> %s, %s" % ((x,y,z),(r1,r2,r3), newpos, newrot)
            pm.adjust(newpos, rot=newrot)
            continue
        return

    def cost(self, x):
        self.adjust_model(x)
        r = 0.
        for pm in self.pack_models:
            pack_difc = self.difc[pm.firstIndex:pm.lastIndex+1]
            difc_new = pm.difc()
            pack_L2 = self.L2[pm.firstIndex:pm.lastIndex+1]
            L2_new = pm.L2()
            residual = np.sum((difc_new/pack_difc-1)**2) + np.sum((L2_new/pack_L2-1)**2)
            residual/=pack_L2.size
            r += residual
            continue
        r/=len(self.pack_models)
        # distance between models are larger than pack width
        from scipy.special import erf
        r2 = 0.
        for i in range(len(self.pack_models)-1):
            pm1 = self.pack_models[i]
            pm2 = self.pack_models[i+1]
            x1,y1,z1 = self.original_positions[pm1]
            x2,y2,z2 = self.original_positions[pm2]
            dx1,dz1,dry1 = x[3*i: 3*(i+1)]
            dx2,dz2,dry2 = x[3*(i+1): 3*(i+2)]
            dist = ((x1+dx1-x2-dx2)**2 + (z1+dz1-z2-dz2)**2)**.5
            print "dist=",dist
            r2+= (erf((-dist+self.min_pack_distance)/0.001) +1)/2. * 1e-2
            continue
        r2/=len(self.pack_models)-1
        r+= r2
        print "cost:", r
        print
        return r
    
    def fit(self, bounds, x0=None):
        # results = sopt.minimize(
        results = sopt.differential_evolution(
            self.cost, 
            bounds=bounds,
            # x0=x0,
            # method='L-BFGS-B',
            )
        return results
    

class FitPack(object):
    
    def __init__(self, pack_model, options, logger=None, params0 = None):
        self.pack_model = pack_model
        self.params0 = params0 or pack_model.getParams()
        print "* initial parameters: ", self.params0
        self.options = options
        self.move = False
        for v in options.values()[:3]:
            if v: self.move = True
        self.rotate = False
        for v in options.values()[3:]:
            if v: self.rotate = True
        self.logger = logger
        return
    
    def adjust_model(self, params):
        pos = self.move and params[:3] # boolean or 3-vector
        rot = self.rotate and params[3:]
        self.pack_model.adjust(pos, rot)
        return
    
    def fit(self):
        pack_model = self.pack_model
        options = self.options
        logger = self.logger
        mask_out = pack_model.mask
        if mask_out.sum() == mask_out.size:
            logger.warning("All pixels in '%s' are masked. Skipping calibration." % pack_model.name)
            return
        # print "move: ", move
        # print "rotate: ", rotate
        logger.info(
            "Working on %s. Starting position is %s, Starting rotation is %s" % (
                pack_model.fullname, pack_model.position(), pack_model.eulerAngles())
            )
        x0List = []
        boundsList = []
        params0 = self.params0
        for iopt,opt in enumerate(options.keys()):
            if not options[opt]: continue
            min, max = options[opt]
            x0List.append(params0[iopt])
            boundsList.append((params0[iopt]+min, params0[iopt]+max))            

        # results = sopt.differential_evolution(
        results = sopt.minimize(
            self.cost, 
            bounds=boundsList,
            x0=x0List,
            method='L-BFGS-B',
            )
        # Apply the results to the output workspace
        params = opts2fullparams(results.x, params0, options)
        self.adjust_model(params)
        # Need to grab the component again, as things have changed   
        pack_model.getComponent()
        logger.info(
            "Finished optimizing %s. Final position is %s, Final rotation is %s" % (
                pack_model.fullname, pack_model.position(), pack_model.eulerAngles())
            )        
        return


class FitPackDifc(FitPack):
    
    def __init__(self, pack_model, options, difc, logger=None):
        self.difc = difc[pack_model.firstIndex:pack_model.lastIndex+1]
        super(FitPackDifc, self).__init__(pack_model, options, logger)
        return

    def cost(self, x):
        params0 = self.params0
        options = self.options
        params = opts2fullparams(x, params0, options)
        self.adjust_model(params)
        wks_name = self.pack_model.instrument_model.wks_name
        difc_new = self.pack_model.difc()
        difc1 = self.difc
        return np.sum((difc_new-difc1)**2)

    
class FitPackTwothetaAndL2(FitPack):
    
    def __init__(self, pack_model, options, sin_theta, L2, logger=None):
        self.sin_theta = np.ma.masked_array(sin_theta[pack_model.firstIndex:pack_model.lastIndex+1], pack_model.mask)
        self.L2 = np.ma.masked_array(L2[pack_model.firstIndex:pack_model.lastIndex+1], pack_model.mask)
        super(FitPackTwothetaAndL2, self).__init__(pack_model, options, logger)
        return

    def cost(self, x):
        params0 = self.params0
        options = self.options
        params = opts2fullparams(x, params0, options)
        self.adjust_model(params)
        pack_model = self.pack_model
        wks_name = pack_model.instrument_model.wks_name

        twotheta_new, L2_new = pack_model.twotheta_and_L2()
        sin_theta_new = np.sin(twotheta_new/2.)
        
        sin_theta = self.sin_theta
        L2 = self.L2
        
        residual = np.sum((L2_new/L2-1)**2) + np.sum((sin_theta_new/sin_theta-1)**2) * 10
        residual /= L2.size
        print residual
        return residual
    
    
class FitPack_DifcL2(FitPack):
    
    def __init__(self, pack_model, options, difc, L2, logger=None, params0=None):
        self.difc = np.ma.masked_array(difc[pack_model.firstIndex:pack_model.lastIndex+1], pack_model.mask)
        self.L2 = np.ma.masked_array(L2[pack_model.firstIndex:pack_model.lastIndex+1], pack_model.mask)
        super(FitPack_DifcL2, self).__init__(pack_model, options, logger=logger, params0=params0)
        return

    def cost(self, x):
        params0 = self.params0
        options = self.options
        params = opts2fullparams(x, params0, options)
        self.adjust_model(params)
        pack_model = self.pack_model

        L2_new = pack_model.L2()
        difc_new = pack_model.difc()

        L2 = self.L2
        difc = self.difc
        
        # residual = 0*np.sum((L2_new/L2-1)**2) + np.sum((difc_new/difc-1)**2)
        # residual = np.sum((difc_new-difc)**2) + np.sum((L2_new-L2)**2)*1000000
        residual = np.sum((difc_new/difc-1)**2) + np.sum((L2_new/L2-1)**2)
        residual /= L2.size
        print residual
        return residual


def compute_initial_guess_of_pixel_positions(
        sin_theta, L2, mask, init_center, init_beta,
        pixel_height=0.0046875*2, tube_spacing = 0.013716*2,
):
    center1 = _estimate_pack_center_position(sin_theta, L2, mask, init_center)
    # y of pixels
    pack_pixel_indexes = np.repeat(np.arange(128), 8)
    pack_pixel_indexes.shape = 128,8
    pack_pixel_indexes = pack_pixel_indexes.T
    pack_y1 = center1[1]+(pack_pixel_indexes-63.5)*pixel_height
    # x and z
    pack_tube_indexes = np.repeat(np.arange(8), 128)
    pack_tube_indexes.shape = 8, 128
    # with no rotation the dx, dz of tubes from center of pack is
    dx_tubes = (pack_tube_indexes-3.5)*tube_spacing
    dz_tubes = pack_tube_indexes*0.
    # after rotation of angle beta around y axis, dx1 = dx*cos(beta), dz1 = -dx*sin(beta)
    beta = init_beta/180.*np.pi
    dx1_tubes = dx_tubes*np.cos(beta)
    dz1_tubes = -dx_tubes*np.sin(beta)
    pack_x1 = center1[0]+dx1_tubes
    pack_z1 = center1[2]+dz1_tubes
    return pack_x1, pack_y1, pack_z1


def compute_pixel_positions(sin_theta, L2, x0, y0, plot=False):
    # calculate positions
    if plot:
        from matplotlib import pyplot as plt
        plt.figure(figsize=(7, 5))
    sin_theta = sin_theta.view(); sin_theta.shape = 8, 128
    if plot:
        plt.subplot(2,3,1)
        plt.title(r"$\sin ( \theta )$")
        plt.imshow(sin_theta.T)
        plt.colorbar()
    L2 = L2.view(); L2.shape = 8, 128
    if plot:
        plt.subplot(2,3,2)
        plt.title(r"L2")
        plt.imshow(L2.T)
        plt.colorbar()
    twotheta = np.arcsin(sin_theta)*2
    if plot:
        plt.subplot(2,3,3)
        plt.title(r"$2\theta$")
        plt.imshow(twotheta.T*180./np.pi)
        plt.colorbar()
    z = L2*np.cos(twotheta)
    y = y0
    # sign: what if x is very close to zero and it could flip sign?
    x = (L2 **2 - z**2 - y**2)**.5 * np.sign(x0)
    if plot:
        plt.subplot(2,3,4)
        plt.title(r"$x$")
        plt.imshow(x.T)
        plt.colorbar()
        plt.subplot(2,3,5)
        plt.title(r"$y$")
        plt.imshow(y.T)
        plt.colorbar()
        plt.subplot(2,3,6)
        plt.title(r"$z$")
        plt.imshow(z.T)
        plt.colorbar()
        plt.tight_layout()
    return x,y,z


def estimate_pack_center_position(sin_theta, L2, pack, init_center, mask=None):
    """estimate pack center position. take the full data arrays and pack info
    """
    # problem: what if pixels near center are all masked?
    # get pack data
    pack_sin_theta = sin_theta[pack.firstIndex:pack.lastIndex+1]
    pack_L2 = L2[pack.firstIndex:pack.lastIndex+1]
    if mask is not None:
        pack_mask = mask[pack.firstIndex:pack.lastIndex+1]
    else:
        pack_mask = np.ma.getmask(pack_L2)
    pack_sin_theta.shape = pack_L2.shape = pack_mask.shape = 8, 128
    # and call the actual algorithm
    return _estimate_pack_center_position(pack_sin_theta, pack_L2, pack_mask, init_center)
    
def _estimate_pack_center_position(sin_theta, L2, mask, init_center):
    "the actual algorithm to find the pack center"
    # import pdb; pdb.set_trace()
    sin_theta = sin_theta.view(); L2 = L2.view(); mask = mask.view()
    sin_theta.shape = L2.shape = mask.shape = 8, 128
    mask2 = mask.copy()
    # symmetrize
    mask2 = np.logical_or(mask[::-1, :], mask2); mask2 = np.logical_or(mask2, mask2[:, ::-1])
    masked2_sin_theta = np.ma.masked_array(sin_theta, mask2)
    center_sin_theta = np.ma.mean(masked2_sin_theta[3:5, 60:68])
    center_twotheta = np.arcsin(center_sin_theta)*2
    print mask2[3:5, 60:68]
    print "center twotheta: ", center_twotheta
    masked2_L2 = np.ma.masked_array(L2, mask2)
    center_L2 = np.ma.mean(masked2_L2[3:5, 60::68])
    y_center = init_center[1]
    z_center = center_L2*np.cos(center_twotheta)
    # sign: what if x is very close to zero and it could flip sign?
    # import pdb; pdb.set_trace()
    x_center = (center_L2 **2 - z_center**2 - y_center**2)**.5 * np.sign(init_center[0])
    return x_center, y_center, z_center

    
def getFirstDetID(component):
    """
    recursive search to find first detID of a component 
    """
    if component.type() == 'DetectorComponent' or component.type() == 'RectangularDetectorPixel':
        return component.getID()
    else:
        return getFirstDetID(component[0])

    
def getLastDetID(component):
    """  
    recursive search to find last detID of a component 
    """
    if component.type() == 'DetectorComponent' or component.type() == 'RectangularDetectorPixel':
        return component.getID()
    else:
        return getLastDetID(component[component.nelements() - 1])
