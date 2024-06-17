from astropy.io import fits, ascii
import glob
import numpy as np
import gettemp
import h5py
from scipy.optimize import minimize_scalar
from astropy import cosmology
astrocosmo = cosmology.FlatLambdaCDM(H0=67.7,Om0=0.307)

zs=np.array([20.05,14.99,11.98,10.976,10.0,9.389,9.002,8.449,8.012,7.595,7.236,7.005,6.491,6.011,5.847,5.530,5.228,4.996,4.665,4.428,4.177,4.008,3.709,3.491,3.283,3.008,2.896,2.733,2.578,2.444,2.316,2.208,2.103,2.002,1.904,1.823,1.743,1.667,1.604,1.531,1.496,1.414,1.358,1.302,1.248,1.206,1.155,1.114,1.074,1.036,0.997,0.951,0.923,0.887,0.851,0.817,0.791,0.757,0.733,0.700,0.676,0.644,0.621,0.598,0.576,0.546,0.525,0.503,0.482,0.461,0.440,0.420,0.40,0.380,0.361,0.348,0.329,0.310,0.298,0.273,0.261,0.244,0.226,0.214,0.1973,0.1804,0.1693,0.1527,0.1419,0.1258,0.1099,0.0994,0.0839,0.0737,0.0585,0.0485,0.0337,0.0240,0.0095,0])
age=astrocosmo.age(zs)

f=fits.open('alltnggal_withr_wmetals_wvmax.fits')
obj=np.loadtxt('sampleselection_sfr5_notide_nomass.txt',skiprows=1,delimiter=',')

wg=np.where((f[1].data.mstar<1E12) & (f[1].data.mstar>3E7) & (f[1].data.sfr>0) & (f[1].data.half_mass_rad_star>0))[0]

#scale factor for z=0.05, 0.1, z=0.5 -> ages of .7 Gyr, 1.3 Gyr, 5Gyr (ages not right cosmology, but should be fine)
agecut1=0.95
agecut2=0.91
agecut3=.67

allm1=[]
allm2=[]
allm3=[]
allm4=[]
allr1=[]
allr2=[]
allr3=[]
allr4=[]
alldist1=[]
alldist2=[]
alldist3=[]
alldist4=[]
alllum1=[]
alllum2=[]
alllum3=[]
alllum4=[]
allrall=[]
allmall=[]
alllumall=[]
mformedintime=np.zeros([len(obj),len(zs)])
metalage=[]

#wdo=np.where(obj[:,2]==4)[0]
wdo=np.arange(len(obj))
#wdo=np.array([0])
#wdo=np.arange(5)
h=0.6774
def massr(rv,masses,dists):
    w=np.where(dists<rv)[0]
    return np.sum(masses[w])

#for i in range(len(obj)):
for i in wdo:
    if (i%100)==0:
        print(i)

#for i in range(15):
    g=glob.glob('/Volumes/Elements/illustris/cutouts/'+str(int(obj[i,1]))+'/*hdf5')
    if len(g)==0:
        print('no files')
        mformedintime[i]=np.nan
        continue

    hf=h5py.File('/Volumes/Elements/illustris/cutouts/'+str(int(obj[i,1]))+'/'+str(int(obj[i,1]))+'_99.hdf5','r')

    if 'PartType4' not in hf:
        mformedintime[i]=np.nan
        continue

    indx=np.digitize(hf['PartType4']['GFM_StellarFormationTime'][:],1/(1+zs))
    for j in range(len(indx)):
        mformedintime[i,indx[j]-1]+=np.array(hf['PartType4']['Masses'][j]*1E10/.6774)

mformedintime=np.cumsum(mformedintime,axis=1)

from scipy import stats
metal=f[1].data.gas_metal_sfr[obj[:,0].astype(np.int)]
corr=[]
times=[]

for i in range(len(zs)):
    wf=np.where(np.isfinite(mformedintime[wdo,-1]-mformedintime[wdo,i]))[0]
    corr.append(stats.spearmanr(mformedintime[wdo[wf],-1]-mformedintime[wdo[wf],i],metal[wdo[wf]])[0])

np.savetxt('corr_all.txt',np.array(corr).T)

corrwgas=[]
for i in range(len(zs)):
    wf=np.where(np.isfinite(mformedintime[wdo,-1]-mformedintime[wdo,i]))[0]
    corrwgas.append(stats.spearmanr(f[1].data.mgas[wdo[wf]]/(mformedintime[wdo[wf],-1]-mformedintime[wdo[wf],i]),metal[wdo[wf]])[0])
    if ~np.isfinite(stats.spearmanr(f[1].data.mgas[wdo[wf]]/(mformedintime[wdo[wf],-1]-mformedintime[wdo[wf],i]),metal[wdo[wf]])[0]):
        print(f[1].data.mgas[wdo[wf]]/(mformedintime[wdo[wf],-1]-mformedintime[wdo[wf],i]))
np.savetxt('corrwgas_all.txt',np.array(corrwgas).T)
