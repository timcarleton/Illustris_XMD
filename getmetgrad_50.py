from astropy.io import ascii
import h5py
import getbryan98
from astropy import cosmology,units,constants
from scipy.optimize import minimize_scalar
import glob
from astropy import table
from astropy.io import fits
import numpy as np
import bindata

astrocosmo=cosmology.Planck15

hconst=0.6774

zs=np.array([20.05,14.99,11.98,10.976,10.0,9.389,9.002,8.449,8.012,7.595,7.236,7.005,6.491,6.011,5.847,5.530,5.228,4.996,4.665,4.428,4.177,4.008,3.709,3.491,3.283,3.008,2.896,2.733,2.578,2.444,2.316,2.208,2.103,2.002,1.904,1.823,1.743,1.667,1.604,1.531,1.496,1.414,1.358,1.302,1.248,1.206,1.155,1.114,1.074,1.036,0.997,0.951,0.923,0.887,0.851,0.817,0.791,0.757,0.733,0.700,0.676,0.644,0.621,0.598,0.576,0.546,0.525,0.503,0.482,0.461,0.440,0.420,0.40,0.380,0.361,0.348,0.329,0.310,0.298,0.273,0.261,0.244,0.226,0.214,0.1973,0.1804,0.1693,0.1527,0.1419,0.1258,0.1099,0.0994,0.0839,0.0737,0.0585,0.0485,0.0337,0.0240,0.0095,0])


#obj=np.loadtxt('sampleselection_tng100_sfr5.txt',delimiter=',',skiprows=1)
#obj=np.loadtxt('sampleselection_sfr5_notide_nomass.txt',delimiter=',',skiprows=1)
obj=np.loadtxt('sampleselection_tng50_sfr5_notide_nomass.txt',delimiter=',',skiprows=1)

f=fits.open('tng50/alltng50gal_withr_withspin.fits')
#f=fits.open('alltnggal_withr_wmetals_wvmax.fits')

def massr(rv,masses,dists):
    w=np.where(dists<rv)[0]
    return np.sum(masses[w])

allmetr=[]
allmgr=[]
allmgrb=[]
allmgzr=[]
allzratio=[]
allmetrobj=[]
allzobj=[]
allsize=[]

rbins=np.linspace(0,10,num=50)

for typei in [1,2,3,4]:
    wcat=np.where(obj[:,2]==typei)[0]

    alldistg=[]
    allmassg=[]
    allmet=[]

    #for i in wcat[0:5]:
    for i in wcat:

        #h0=h5py.File('/Volumes/Elements/illustris/cutouts50/'+str(int(obj[i,1]))+'/'+str(int(obj[i,1]))+'_'+str(99)+'.hdf5','r')
        h0=h5py.File('/Volumes/Elements 1/illustris/cutouts50/'+str(int(obj[i,1]))+'/'+str(int(obj[i,1]))+'_'+str(99)+'.hdf5','r')

        centerx=np.mean(h0['PartType4']['Coordinates'][:,0])
        centery=np.mean(h0['PartType4']['Coordinates'][:,1])
        centerz=np.mean(h0['PartType4']['Coordinates'][:,2])

        distst=np.sqrt((h0['PartType4']['Coordinates'][:,0]-centerx)**2+
                       (h0['PartType4']['Coordinates'][:,1]-centery)**2+
                       (h0['PartType4']['Coordinates'][:,2]-centerz)**2)


        distgas=np.sqrt((h0['PartType0']['Coordinates'][:,0]-centerx)**2+
                        (h0['PartType0']['Coordinates'][:,1]-centery)**2+
                        (h0['PartType0']['Coordinates'][:,2]-centerz)**2)

        distdm=np.sqrt((h0['PartType1']['Coordinates'][:,0]-centerx)**2+
                       (h0['PartType1']['Coordinates'][:,1]-centery)**2+
                       (h0['PartType1']['Coordinates'][:,2]-centerz)**2)

        alldist=np.hstack([distdm,distgas,distst])/hconst
        allmass=np.hstack([np.array([5.1E6]*len(h0['PartType1']['ParticleIDs'])),h0['PartType0']['Masses'][:]*1E10,h0['PartType4']['Masses'][:]*1E10])


        dvir=getbryan98.getbryan98(astrocosmo,0)*astrocosmo.critical_density(0).to(units.M_sun/units.kpc**3).value
    
        vir=minimize_scalar(lambda rv:(massr(rv,allmass,alldist)/hconst-(4/3*np.pi*dvir*rv**3))**2,bounds=[.5,200],method='bounded')

        rvir=vir['x']
        mvir=massr(vir['x'],allmass,alldist)
        rhalf=np.median(distst)

        allmet+=list(h0['PartType0']['GFM_Metallicity'][:])
        allmassg+=list(h0['PartType0']['Masses'][:])
        alldistg+=list(np.array(distgas)/rhalf)
        w1=np.where(np.array(distgas)/rhalf<1)[0]
        w2=np.where((np.array(distgas)/rhalf>4) & (np.array(distgas)/rhalf<5))[0]
        z1=np.sum(np.array(h0['PartType0']['GFM_Metallicity'][:])[w1]*np.array(h0['PartType0']['Masses'][:])[w1])/np.sum(np.array(h0['PartType0']['Masses'][:])[w1])
        z2=np.sum(np.array(h0['PartType0']['GFM_Metallicity'][:])[w2]*np.array(h0['PartType0']['Masses'][:])[w2])/np.sum(np.array(h0['PartType0']['Masses'][:])[w2])
        allzratio.append(z1/z2)
        #alldistg+=list(np.array(distgas))
        bi=bindata.bindata(np.array(distgas)/rhalf,np.array(h0['PartType0']['GFM_Metallicity'][:])*np.array(h0['PartType0']['Masses'][:]),rbins,median=1)
        allmetrobj.append(bi[0]/np.median(h0['PartType0']['Masses'][:]))
        allsize.append(rhalf)

    b=bindata.bindata(np.array(alldistg),np.array(allmassg)*np.array(allmet),rbins,median=1)
    bg=bindata.bindata(np.array(alldistg),np.array(allmassg),rbins,median=1)
    
    
    allmetr.append(b[0]/np.median(allmassg))
    allmgr.append(bg[0]*1E10/.6774)
    allmgzr.append(b[0]*1E10/.6774)
    

import matplotlib.pyplot as plt
plt.style.use('../python/stylesheet.txt')

plt.clf()
plt.plot(allsize,allzratio,'o',alpha=.2)
plt.xlim(0,10)
plt.xlabel('size')
plt.ylabel('met grad')
plt.savefig('metrvssize_50.png')

plt.clf()
colors=['C0','C1','C2','C3']
lst=['-','--',':','-.']
for typei in [1,2,3,4]:
    wcat=np.where(obj[:,2]==typei)[0]
    #for i in wcat[0:10]:
        #plt.plot((rbins[0:-1]+rbins[1:])/2,allmetrobj[i]/.0127,color=colors[typei-1],ls=lst[typei-1],alpha=.2)
    
plt.plot((rbins[0:-1]+rbins[1:])/2,allmetr[3]/.0127,color='C3',ls='-.',label='LMD Analogs')
plt.plot((rbins[0:-1]+rbins[1:])/2,allmetr[2]/.0127,color='C2',ls='--',label='XMD Analogs')
plt.plot((rbins[0:-1]+rbins[1:])/2,allmetr[1]/.0127,color='C1',ls=':',label='LMDs')
plt.plot((rbins[0:-1]+rbins[1:])/2,allmetr[0]/.0127,color='C0',ls='-',label='XMDs')

plt.legend()
plt.ylim(1E-2,1)
plt.title('TNG50',y=.92)
plt.yscale('log')
plt.xlabel(r'$r/r_h$')
plt.ylabel(r'$Z_p/Z_{\rm \odot}$')
plt.savefig('metgrad_noscale_50.png')

plt.clf()
colors=['C0','C1','C2','C3']
lst=['-','--',':','-.']
allnormgrad=[]
for typei in [1,2,3,4]:
    allnormgradi=[]
    wcat=np.where(obj[:,2]==typei)[0]
    for i in range(len(wcat)):
        #if i>10 and i<20:
            #print(typei,i,allmetrobj[wcat[i]]/allmetrobj[wcat[i]][-1])
            #plt.plot((rbins[0:-1]+rbins[1:])/2,allmetrobj[wcat[i]]/allmetrobj[wcat[i]][-1],color=colors[typei-1],ls=lst[typei-1],alpha=.2)
        allnormgradi.append(allmetrobj[wcat[i]]/allmetrobj[wcat[i]][-1])
    allnormgrad.append(np.nanmedian(allnormgradi,axis=0))

print(allnormgrad)
plt.plot((rbins[0:-1]+rbins[1:])/2,allnormgrad[3],color='C3',ls='-.',label='LMD Analogs')
plt.plot((rbins[0:-1]+rbins[1:])/2,allnormgrad[2],color='C2',ls=':',label='XMD Analogs')
plt.plot((rbins[0:-1]+rbins[1:])/2,allnormgrad[1],color='C1',ls='--',label='LMDs')
plt.plot((rbins[0:-1]+rbins[1:])/2,allnormgrad[0],color='C0',ls='-',label='XMDs')

plt.legend()
plt.ylim(.5,40)
plt.yscale('log')
plt.xlabel(r'$r/rh$')
plt.ylabel(r'$Z/Z_{\rm 5rh}$')
plt.savefig('normmetgrad_50.png')

plt.clf()
plt.plot((rbins[0:-1]+rbins[1:])/2,allmgr[3],color='C3',ls='-.',label='LMD Analogs',alpha=.5)
plt.plot((rbins[0:-1]+rbins[1:])/2,allmgr[2],color='C2',ls=':',label='XMD Analogs',alpha=.5)
plt.plot((rbins[0:-1]+rbins[1:])/2,allmgr[1],color='C1',ls='--',label='LMDs',alpha=.5)
plt.plot((rbins[0:-1]+rbins[1:])/2,allmgr[0],color='C0',ls='-',label='XMDs',alpha=.5)
plt.yscale('log')

ax2=plt.twinx()
ax2.plot((rbins[0:-1]+rbins[1:])/2,allmgzr[3],color='C3',ls='-.',label='LMD Analogs')
ax2.plot((rbins[0:-1]+rbins[1:])/2,allmgzr[2],color='C2',ls=':',label='XMD Analogs')
ax2.plot((rbins[0:-1]+rbins[1:])/2,allmgzr[1],color='C1',ls='--',label='LMDs')
ax2.plot((rbins[0:-1]+rbins[1:])/2,allmgzr[0],color='C0',ls='-',label='XMDs')

ax2.set_yscale('log')
plt.legend()

plt.xlabel(r'$r/r_h$')
plt.ylabel(r'$Z/Z_{\rm \odot}$')
plt.savefig('gasmetgrad_50.png')
