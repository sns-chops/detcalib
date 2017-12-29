import sys, os, numpy as np
from matplotlib import pyplot as plt

workdir = os.path.abspath('.')
sys.path.insert(0, "%s/lib"%workdir)

import get_difc

pack = 'B5/eightpack'
packname, packtype = pack.split('/')

fitter = get_difc.Fitter(peak_fractional_width=0.02, bg_type='linear')
# dvalues = [1.10860231, 1.24596143, 1.357755, 1.63751414, 1.92015553, 3.13560085]
dvalues = [
    0.4003801 ,  0.40480433,  0.40593349,  0.40937853,  0.41411139,
    0.41532047,  0.41901228,  0.42409141,  0.42539031,  0.42935983,
    0.43623002,  0.44051389,  0.44794297,  0.452585  ,  0.45900496,
    0.4606531 ,  0.46570612,  0.47270961,  0.47451042,  0.48003888,
    0.48969858,  0.49578202,  0.50425754,  0.50644522,  0.52260014,
    0.52503652,  0.53255532,  0.54583804,  0.55430115,  0.56932559,
    0.57894867,  0.59257286,  0.59613189,  0.60720649,  0.62298071,
    0.62712017,  0.64005184,  0.65860791,  0.66350469,  0.6788775 ,]
dvalues = [
    0.70705857,  0.72575057,  0.76049491,  0.78390021,  0.81875707,
    0.82822278,  0.85871966,  0.90517   ,  0.91800993,  0.96007776,
    1.04520028,  1.10860231,  1.24596143,  1.357755  ,  1.56780042,
    1.63751414,  1.92015553,  3.13560085]
dmin=.7; dmax=3.5

I_d_dir = 'Si-I_d'

gpd = get_difc.GetPackDifc(
    pack=packname,
    dvalues=dvalues,
    dmin=dmin, dmax=dmax,
    I_d_dir = I_d_dir,
    fitter=fitter,
    maxchisq = 10.,
    min_counts = 3000)

newdifc_pixel, (x,y,fitx, fity), (fitfunc, popt, chisq) = gpd.fitOnePixel(100)

print chisq

plt.figure()
plt.plot(x,y, '+')
# plt.plot(fitx, fity)
plt.plot(fitx, fitfunc(fitx, *popt))
plt.xlim(0.6, 1.6)
plt.show()
