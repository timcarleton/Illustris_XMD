from astropy.io import fits, ascii
import glob
import numpy as np
import gettemp
import h5py
from scipy.optimize import minimize_scalar
import downloadcutout_snap

f=fits.open('tng50/alltng50gal_withr_withspin.fits')
obj=np.loadtxt('sampleselection_tng50_sfr5.txt',skiprows=1,delimiter=',')

tempcut=50

wg=np.where((f[1].data.mstar<1E12) & (f[1].data.mstar>3E7) & (f[1].data.sfr>0) & (f[1].data.half_mass_rad_star>0))[0]

allmcold=[]
allrcold=[]

h=0.6774
def massr(rv,masses,dists):
    w=np.where(dists<rv)[0]
    return np.sum(masses[w])

for i in range(len(obj)):
#for i in range(5):
    g=glob.glob('/Volumes/Elements/illustris/cutouts50/'+str(int(obj[i,1]))+'/*hdf5')
    if len(g)==0:
        print('no files')
        allmcold.append(0)
        allrcold.append(0)
        continue

    hf=h5py.File('/Volumes/Elements/illustris/cutouts50/'+str(int(obj[i,1]))+'/'+str(int(obj[i,1]))+'_99.hdf5','r')

    if 'PartType0' not in hf:
        allmcold.append(0)
        allrcold.append(0)
        continue
    
    temps=gettemp.gettemp('/Volumes/Elements/illustris/cutouts50/'+str(int(obj[i,1]))+'/'+str(int(obj[i,1]))+'_99.hdf5')

    w=np.where(temps<tempcut)[0]

    allmcold.append(np.sum(np.array(hf['PartType0']['Masses'])[w])*1E10/h)
    xcold=np.array(hf['PartType0']['Coordinates'][:,0])[w]
    ycold=np.array(hf['PartType0']['Coordinates'][:,1])[w]
    zcold=np.array(hf['PartType0']['Coordinates'][:,2])[w]
    mcold=np.array(hf['PartType0']['Masses'][:])[w]
    
    centerxcold=np.mean(xcold*mcold)/np.mean(mcold)
    centerycold=np.mean(ycold*mcold)/np.mean(mcold)
    centerzcold=np.mean(zcold*mcold)/np.mean(mcold)

    distcold=np.sqrt((xcold-centerxcold)**2+(ycold-centerycold)**2+(zcold-centerzcold)**2)

    m12i=minimize_scalar(lambda r12i: abs(massr(r12i,mcold,distcold)-0.5*np.sum(mcold)),bounds=[0.1,2*np.median(distcold)],method='bounded')
    
    
    allrcold.append(m12i['x'])


from astropy import table
tab=table.Table([obj[:,0],obj[:,1],obj[:,2],allmcold,allrcold],names=['allgalindex','haloid','type','mcold','rcold'])
tab.write('coldgasinfo50_notide.txt',format='csv')
