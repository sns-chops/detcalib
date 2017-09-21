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
            (rotw, rotx, roty, rotz) = eulerToAngleAxis(rot[0], rot[1], rot[2], eulerConvention)
            msa.RotateInstrumentComponent(
                wks_name, component, X=rotx, Y=roty, Z=rotz, Angle=rotw, RelativeRotation=False)
        return
    
    def difc(self):
        wks_name = self.wks_name
        msa.CalculateDIFC(InputWorkspace=wks_name, OutputWorkspace='difc')
        return mtd['difc'].extractY().flatten()

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
            self.cost, x0=x0List,
            method='L-BFGS-B',
            bounds=boundsList)

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

    def twotheta_and_L2(self):
        tt, L2 = self.instrument_model.twotheta_and_L2()
        return (np.ma.masked_array(tt[self.firstIndex:self.lastIndex+1], self.mask),
                np.ma.masked_array(L2[self.firstIndex:self.lastIndex+1], self.mask))
    
    def getComponent(self):
        wks_name = self.instrument_model.wks_name
        instrument = mtd[wks_name].getInstrument()
        self.comp = instrument.getComponentByName(self.name)
        return self.comp

class FitPack(object):
    
    def __init__(self, pack_model, options, logger=None):
        self.pack_model = pack_model
        self.params0 = pack_model.getParams()
        self.options = options
        self.move = np.any(options.values()[:3])
        self.rotate = np.any(options.values()[3:])
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

        results = sopt.differential_evolution(
        # results = sopt.minimize(
            self.cost, 
            # x0=x0List, method='L-BFGS-B',
            bounds=boundsList)

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
        firstIndex = pack_model.firstIndex
        lastIndex = pack_model.lastIndex
        mask = pack_model.mask

        twotheta_new, L2_new = pack_model.twotheta_and_L2()
        sin_theta_new = np.sin(twotheta_new/2.)
        
        sin_theta = self.sin_theta
        L2 = self.L2
        
        residual = np.sum((L2_new/L2-1)**2) + np.sum((sin_theta_new/sin_theta-1)**2)
        residual /= L2.size
        print residual
        return residual
    
    
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
