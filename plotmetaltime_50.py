import matplotlib.pyplot as plt
import os
import numpy as np
import h5py
from astropy import cosmology
astrocosmo = cosmology.FlatLambdaCDM(H0=67.7,Om0=0.307)
from astropy.io import ascii
plt.style.use('../python/stylesheet.txt')

plt.rcParams['xtick.top']=False
plt.rcParams['xtick.labeltop']=False
x=ascii.read('sampleselection_tng50_sfr5_notide_nomass.txt')

wxmd=np.where(x['type']==1)[0]
wlmd=np.where(x['type']==2)[0]
wxmda=np.where(x['type']==3)[0]
wlmda=np.where(x['type']==4)[0]

ts=np.linspace(0,13.8)

zs=np.array([20.05,14.99,11.98,10.976,10.0,9.389,9.002,8.449,8.012,7.595,7.236,7.005,6.491,6.011,5.847,5.530,5.228,4.996,4.665,4.428,4.177,4.008,3.709,3.491,3.283,3.008,2.896,2.733,2.578,2.444,2.316,2.208,2.103,2.002,1.904,1.823,1.743,1.667,1.604,1.531,1.496,1.414,1.358,1.302,1.248,1.206,1.155,1.114,1.074,1.036,0.997,0.951,0.923,0.887,0.851,0.817,0.791,0.757,0.733,0.700,0.676,0.644,0.621,0.598,0.576,0.546,0.525,0.503,0.482,0.461,0.440,0.420,0.40,0.380,0.361,0.348,0.329,0.310,0.298,0.273,0.261,0.244,0.226,0.214,0.1973,0.1804,0.1693,0.1527,0.1419,0.1258,0.1099,0.0994,0.0839,0.0737,0.0585,0.0485,0.0337,0.0240,0.0095,0])
age=astrocosmo.age(zs)
snaps=np.arange(99)
ox_mass = 15.999
fe_mass = 55.845
h_mass = 1.008

n_mass=14.0067
c_mass=12.011
ne_mass=20.180
mg_mass=24.305
si_mass=28.085

o_nfrac_sun = 10**(8.69-12)
fe_nfrac_sun = 10**(7.5-12)
n_nfrac_sun = 10**(7.83-12)
c_nfrac_sun = 10**(8.43-12)
ne_nfrac_sun = 10**(7.93-12)
mg_nfrac_sun = 10**(7.6-12)
si_nfrac_sun = 10**(7.51-12)

allmetsx=np.zeros([len(wxmd),100])+np.nan
allmetsmx=np.zeros([len(wxmd),100])+np.nan
allsfrx=np.zeros([len(wxmd),100])+np.nan
allmstx=np.zeros([len(wxmd),100])+np.nan
allmhx=np.zeros([len(wxmd),100])+np.nan
allmgx=np.zeros([len(wxmd),100])+np.nan
for i in range(len(wxmd)):
    h=h5py.File('/Volumes/Elements/illustris/tng50trees/'+str(x['haloid'][wxmd[i]])+'_tree.hdf5','r')
    snapsi=np.append(h['SnapNum'],99)
    wm=np.where(snapsi[1:]-snapsi[0:-1]>0)[0]

    allmetsx[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],0]/h_mass)/o_nfrac_sun
    allmetsmx[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractions'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractions'][0:wm[0],0]/h_mass)/o_nfrac_sun
    
    #allmetsx[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetallicitySfrWeighted'][0:wm[0]]
    #allmetsmx[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetallicity'][0:wm[0]]
    allsfrx[i,h['SnapNum'][0:wm[0]]]=h['SubhaloSFR'][0:wm[0]]
    allmstx[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMassType'][0:wm[0],4]/.6774*1E10
    allmgx[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMassType'][0:wm[0],0]/.6774*1E10
    allmhx[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMass'][0:wm[0]]/.6774*1E10

allmetsl=np.zeros([len(wlmd),100])+np.nan
allmetsml=np.zeros([len(wlmd),100])+np.nan
allsfrl=np.zeros([len(wlmd),100])+np.nan
allmstl=np.zeros([len(wlmd),100])+np.nan
allmhl=np.zeros([len(wlmd),100])+np.nan
allmgl=np.zeros([len(wlmd),100])+np.nan
for i in range(len(wlmd)):
    h=h5py.File('/Volumes/Elements/illustris/tng50trees/'+str(x['haloid'][wlmd[i]])+'_tree.hdf5','r')

    wm=np.where(h['SnapNum'][1:]-h['SnapNum'][0:-1]>0)[0]
    if len(wm)>0:
        allmetsl[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],0]/h_mass)/o_nfrac_sun
        allmetsml[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractions'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractions'][0:wm[0],0]/h_mass)/o_nfrac_sun
        #allmetsl[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetallicitySfrWeighted'][0:wm[0]]
        #allmetsml[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetallicity'][0:wm[0]]
        allsfrl[i,h['SnapNum'][0:wm[0]]]=h['SubhaloSFR'][0:wm[0]]
        allmstl[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMassType'][0:wm[0],4]/.6774*1E10
        allmgl[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMassType'][0:wm[0],0]/.6774*1E10
        allmhl[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMass'][0:wm[0]]/.6774*1E10


allmetsa=np.zeros([len(wxmda),100])+np.nan
allmetsma=np.zeros([len(wxmda),100])+np.nan
allsfra=np.zeros([len(wxmda),100])+np.nan
allmsta=np.zeros([len(wxmda),100])+np.nan
allmha=np.zeros([len(wxmda),100])+np.nan
allmga=np.zeros([len(wxmda),100])+np.nan
for i in range(len(wxmda)):
    if os.path.exists('/Volumes/Elements/illustris/tng50trees/'+str(x['haloid'][wxmda[i]])+'_tree.hdf5'):
        h=h5py.File('/Volumes/Elements/illustris/tng50trees/'+str(x['haloid'][wxmda[i]])+'_tree.hdf5','r')
        snapsi=np.append(h['SnapNum'],99)
        wm=np.where(snapsi[1:]-snapsi[0:-1]>0)[0]

        if len(wm)>0:
            if h['SubhaloSFR'][0]>0:
                allmetsa[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],0]/h_mass)/o_nfrac_sun
                allmetsma[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractions'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractions'][0:wm[0],0]/h_mass)/o_nfrac_sun

                #allmetsa[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetallicitySfrWeighted'][0:wm[0]]
                #allmetsma[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetallicity'][0:wm[0]]
                allsfra[i,h['SnapNum'][0:wm[0]]]=h['SubhaloSFR'][0:wm[0]]
                allmsta[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMassType'][0:wm[0],4]/.6774*1E10
                allmga[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMassType'][0:wm[0],0]/.6774*1E10
                allmha[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMass'][0:wm[0]]/.6774*1E10
            else:
                allmetsa[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmetsma[i,h['SnapNum'][0:wm[0]]]=np.nan
                allsfra[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmsta[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmga[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmha[i,h['SnapNum'][0:wm[0]]]=np.nan
    else:
        print('ml',i)

allmetsla=np.zeros([len(wlmda),100])+np.nan
allmetsmla=np.zeros([len(wlmda),100])+np.nan
allsfrla=np.zeros([len(wlmda),100])+np.nan
allmstla=np.zeros([len(wlmda),100])+np.nan
allmhla=np.zeros([len(wlmda),100])+np.nan
allmgla=np.zeros([len(wlmda),100])+np.nan
for i in range(len(wlmda)):

    if os.path.exists('/Volumes/Elements/illustris/tng50trees/'+str(x['haloid'][wlmda[i]])+'_tree.hdf5'):
        try:
            h=h5py.File('/Volumes/Elements/illustris/tng50trees/'+str(x['haloid'][wlmda[i]])+'_tree.hdf5','r')
        except:
            print('/Volumes/Elements/illustris/tng50trees/'+str(x['haloid'][wlmda[i]])+'_tree.hdf5')
            continue
        snapsi=np.append(h['SnapNum'],99)
        wm=np.where(snapsi[1:]-snapsi[0:-1]>0)[0]
        if len(wm)>0:
            if h['SubhaloSFR'][0]>0:
                allmetsla[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],0]/h_mass)/o_nfrac_sun
                allmetsmla[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractions'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractions'][0:wm[0],0]/h_mass)/o_nfrac_sun

                #allmetsla[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetallicitySfrWeighted'][0:wm[0]]
                #allmetsmla[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetallicity'][0:wm[0]]
                allsfrla[i,h['SnapNum'][0:wm[0]]]=h['SubhaloSFR'][0:wm[0]]
                allmstla[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMassType'][0:wm[0],4]/.6774*1E10
                allmgla[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMassType'][0:wm[0],0]/.6774*1E10
                allmhla[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMass'][0:wm[0]]/.6774*1E10

            else:
                allmetsla[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmetsmla[i,h['SnapNum'][0:wm[0]]]=np.nan
                allsfrla[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmstla[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmgla[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmhla[i,h['SnapNum'][0:wm[0]]]=np.nan
    else:
        print('mla')    

zg=np.array([0,.1,.25,.5,.75,1,1.5,2,3,7])
zglabels=zg.astype(np.str)
agepoints=astrocosmo.age(zg)
toplotx=np.random.choice(len(allmetsx),size=5)
toplota=np.random.choice(len(allmetsa),size=5)
toplotl=np.random.choice(len(allmetsl),size=5)
toplotla=np.random.choice(len(allmetsla),size=5)


plt.clf()
plt.plot(age,np.nanmedian(allmetsla,axis=0),label='SFR-weighted LMD Analogs',lw=5,color='C3',ls='-.')
plt.plot(age,np.nanmedian(allmetsa,axis=0),label='SFR-weighted XMD Analogs',lw=5,color='C2',ls='--')
plt.plot(age,np.nanmedian(allmetsl,axis=0),label='SFR-weighted LMDs',lw=5,color='C1',ls=':')
plt.plot(age,np.nanmedian(allmetsx,axis=0),label='SFR-weighted XMDs',lw=5,color='C0')

plt.plot(age,np.nanmedian(allmetsmla,axis=0),label='Mass-weighted LMD Analogs',lw=5,color='C3',ls=(0, (3, 1, 1, 1)))
plt.plot(age,np.nanmedian(allmetsma,axis=0),label='Mass-weighted XMD Analogs',lw=5,color='C2',ls=(0, (5, 1)))
plt.plot(age,np.nanmedian(allmetsml,axis=0),label='Mass-weighted LMDs',lw=5,color='C1',ls= (0, (1, 1)))
plt.plot(age,np.nanmedian(allmetsmx,axis=0),label='Mass-weighted XMDs',lw=5,color='C0',ls=(5, (10, 3)))

plt.legend(loc='lower left',ncol=2,fontsize=25)
plt.yscale('log')
plt.ylim(1E-2,2)
plt.xlim(.5,13.9)
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$ or $Z_{\rm mass}/Z_\odot$')
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlim(.5,13.9)
ax.set_xlabel(r'$z$')
plt.title('TNG50',y=.92)
plt.savefig('allsfrmetvsz_mvss_50.png')
         
plt.clf()
for i in toplotla:
    plt.plot(age,allmetsla[i],color='C3',alpha=.3,ls='-.')

for i in toplota:
    plt.plot(age,allmetsa[i],color='C2',alpha=.3,ls=':')


for i in toplotl:
    plt.plot(age,allmetsl[i],color='C1',alpha=.3,ls='--')

for i in toplotx:
    plt.plot(age,allmetsx[i],color='C0',alpha=.3,ls='-')




plt.plot(age,np.nanmedian(allmetsla,axis=0),label='LMD Analogs',lw=5,color='C3',ls='-.')
plt.plot(age,np.nanmedian(allmetsa,axis=0),label='XMD Analogs',lw=5,color='C2',ls='--')
plt.plot(age,np.nanmedian(allmetsl,axis=0),label='LMDs',lw=5,color='C1',ls=':')
plt.plot(age,np.nanmedian(allmetsx,axis=0),label='XMDs',lw=5,color='C0',ls='-')
#plt.fill_between(age,np.nanmedian(allmetsx,axis=0)-np.nanstd(allmetsx,axis=0),np.nanmedian(allmetsx,axis=0)+np.nanstd(allmetsx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsa,axis=0)-np.nanstd(allmetsa,axis=0),np.nanmedian(allmetsa,axis=0)+np.nanstd(allmetsa,axis=0),color='C1',alpha=.3)
plt.legend()
plt.yscale('log')
#plt.ylim(5E-4,5E-3)
plt.ylim(1E-2,2)
#plt.xlim(2,13.9)
plt.xlim(.5,13.9)
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.ylabel(r'$Z_{\rm SFR}/{\rm Z_\odot}$')
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlim(.5,13.9)
ax.set_xlabel(r'$z$')
plt.title('TNG50',y=0.92)
plt.savefig('allsfrmetvsz_50.png')

plt.clf()
for i in toplotx:
    plt.plot(age,allmetsmx[i],color='C0',alpha=.3)

for i in toplota:
    plt.plot(age,allmetsma[i],color='C1',alpha=.3)

for i in toplotl:
    plt.plot(age,allmetsml[i],color='C2',alpha=.3)

for i in toplotla:
    plt.plot(age,allmetsmla[i],color='C3',alpha=.3)

plt.plot(age,np.nanmedian(allmetsmla,axis=0),label='LMD Analogs',lw=5,color='C3')    
plt.plot(age,np.nanmedian(allmetsma,axis=0),label='XMD Analogs',lw=5,color='C1')
plt.plot(age,np.nanmedian(allmetsmx,axis=0),label='XMDs',lw=5,color='C0')
plt.plot(age,np.nanmedian(allmetsml,axis=0),label='LMDs',lw=5,color='C2')


#plt.fill_between(age,np.nanmedian(allmetsmx,axis=0)-np.nanstd(allmetsmx,axis=0),np.nanmedian(allmetsmx,axis=0)+np.nanstd(allmetsmx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsma,axis=0)-np.nanstd(allmetsma,axis=0),np.nanmedian(allmetsma,axis=0)+np.nanstd(allmetsma,axis=0),color='C1',alpha=.3)
plt.legend()
plt.yscale('log')
plt.ylim(5E-4,5E-3)
plt.xlim(2,13.9)
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.ylabel('Mass-weighted Metal fraction')
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlim(2,13.9)
ax.set_xlabel(r'$z$')
plt.savefig('allmassmetvsz_50.png')

plt.clf()
for i in toplotx:
    plt.plot(age,allmetsx[i]/allmetsx[i][-1],color='C0',alpha=.3)

for i in toplota:
    plt.plot(age,allmetsa[i]/allmetsa[i][-1],color='C1',alpha=.3)

for i in toplotl:
    plt.plot(age,allmetsl[i]/allmetsl[i][-1],color='C2',alpha=.3)

for i in toplotla:
    plt.plot(age,allmetsla[i]/allmetsla[i][-1],color='C3',alpha=.3)

z0metsx=np.repeat(allmetsx[:,-1],100).reshape(len(allmetsx),100)
z0metsa=np.repeat(allmetsa[:,-1],100).reshape(len(allmetsa),100)
z0metsl=np.repeat(allmetsl[:,-1],100).reshape(len(allmetsl),100)
z0metsla=np.repeat(allmetsla[:,-1],100).reshape(len(allmetsla),100)
print(allmetsx.shape,z0metsx.shape)
plt.plot(age,np.nanmedian(allmetsx/z0metsx,axis=0),label='XMDs',lw=5,color='C0')
plt.plot(age,np.nanmedian(allmetsa/z0metsa,axis=0),label='XMD Analogs',lw=5,color='C1')
plt.plot(age,np.nanmedian(allmetsl/z0metsl,axis=0),label='LMDs',lw=5,color='C2')
plt.plot(age,np.nanmedian(allmetsla/z0metsla,axis=0),label='LMD Analaogs',lw=5,color='C3')
plt.legend()
plt.yscale('log')
plt.ylim(.1,10)
plt.xlim(.5,13.9)
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.ylabel(r'$Z_{\rm SFR}/Z_{\rm SFR, z=0}$')
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlim(.5,13.9)
ax.set_xlabel(r'$z$')
    
plt.savefig('allsfrmetfracgrowth_50.png')

plt.clf()

# for i in toplotla:
#     plt.plot(age,allsfrla[i],color='C3',alpha=.3,ls='-.')

# for i in toplota:
#     plt.plot(age,allsfra[i],color='C2',alpha=.3,ls=':')


# for i in toplotl:
#     plt.plot(age,allsfrl[i],color='C1',alpha=.3,ls='--')

# for i in toplotx:
#     plt.plot(age,allsfrx[i],color='C0',alpha=.3,ls='-')


plt.plot(age,np.nanmedian(allsfrla,axis=0),label='LMD Analogs',lw=5,color='C3',ls='-.')
plt.plot(age,np.nanmedian(allsfra,axis=0),label='XMD Analogs',lw=5,color='C2',ls='--')
plt.plot(age,np.nanmedian(allsfrl,axis=0),label='LMDs',lw=5,color='C1',ls=':')
plt.plot(age,np.nanmedian(allsfrx,axis=0),label='XMDs',lw=5,color='C0',ls='-')
#plt.fill_between(age,np.nanmedian(allmetsx,axis=0)-np.nanstd(allmetsx,axis=0),np.nanmedian(allmetsx,axis=0)+np.nanstd(allmetsx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsa,axis=0)-np.nanstd(allmetsa,axis=0),np.nanmedian(allmetsa,axis=0)+np.nanstd(allmetsa,axis=0),color='C1',alpha=.3)
plt.legend()
print('sfr0 xmda',np.nanmedian(allsfra,axis=0)[-1])
print('sfr0 xmd',np.nanmedian(allsfrx,axis=0)[-1])
plt.yscale('log')
plt.ylim(5E-4,1E-2)
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.xlim(.5,13.9)
plt.ylabel(r'SFR $({\rm M_\odot~yr^{-1}}$)')
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlabel(r'$z$')
ax.set_xlim(.5,13.9)
plt.title('TNG50',y=0.92)
plt.savefig('allsfrvsz_50.png')

plt.clf()
plt.plot(age,np.nanmedian(allmstx,axis=0),label='XMDs',lw=5,color='C0')
plt.plot(age,np.nanmedian(allmsta,axis=0),label='XMD Analogs',lw=5,color='C1')
plt.plot(age,np.nanmedian(allmstl,axis=0),label='LMDs',lw=5,color='C2')
plt.plot(age,np.nanmedian(allmstla,axis=0),label='LMD Analogs',lw=5,color='C3')
#plt.fill_between(age,np.nanmedian(allmetsx,axis=0)-np.nanstd(allmetsx,axis=0),np.nanmedian(allmetsx,axis=0)+np.nanstd(allmetsx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsa,axis=0)-np.nanstd(allmetsa,axis=0),np.nanmedian(allmetsa,axis=0)+np.nanstd(allmetsa,axis=0),color='C1',alpha=.3)
plt.legend()
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.yscale('log')
plt.ylabel(r'$M_*$ $({\rm M_\odot})$')
plt.xlim(2,13.9)
#plt.ylim(5E-4,5E-3)
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlabel(r'$z$')
ax.set_xlim(2,13.9)
plt.savefig('allmstz_50.png')

plt.clf()
plt.plot(age,np.nanmedian((allmstx)/allmhx,axis=0),label='XMDs',lw=5,color='C0')
plt.plot(age,np.nanmedian((allmsta)/allmha,axis=0),label='XMD Analogs',lw=5,color='C1')
plt.plot(age,np.nanmedian((allmstl)/allmhl,axis=0),label='LMDs',lw=5,color='C2')
plt.plot(age,np.nanmedian((allmstla)/allmhla,axis=0),label='LMD Analogs',lw=5,color='C3')
#plt.fill_between(age,np.nanmedian(allmetsx,axis=0)-np.nanstd(allmetsx,axis=0),np.nanmedian(allmetsx,axis=0)+np.nanstd(allmetsx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsa,axis=0)-np.nanstd(allmetsa,axis=0),np.nanmedian(allmetsa,axis=0)+np.nanstd(allmetsa,axis=0),color='C1',alpha=.3)
plt.legend()
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.ylabel(r'SFR $({\rm M_\odot~yr{-1}}$)')
plt.yscale('log')
plt.ylabel(r'$M_*/M_h$')
plt.xlim(2,13.9)
#plt.ylim(5E-4,5E-3)
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlabel(r'$z$')
ax.set_xlim(2,13.9)
plt.savefig('allmhz_50.png')
print(astrocosmo.age(0))


plt.clf()
plt.plot(age,np.nanmedian(allmhx,axis=0),label='XMDs',lw=5,color='C0')
plt.plot(age,np.nanmedian(allmha,axis=0),label='XMD Analogs',lw=5,color='C1')
plt.plot(age,np.nanmedian(allmhl,axis=0),label='LMDs',lw=5,color='C2')
plt.plot(age,np.nanmedian(allmhla,axis=0),label='LMD Analogs',lw=5,color='C3')
#plt.fill_between(age,np.nanmedian(allmetsx,axis=0)-np.nanstd(allmetsx,axis=0),np.nanmedian(allmetsx,axis=0)+np.nanstd(allmetsx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsa,axis=0)-np.nanstd(allmetsa,axis=0),np.nanmedian(allmetsa,axis=0)+np.nanstd(allmetsa,axis=0),color='C1',alpha=.3)
plt.legend()
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.yscale('log')
plt.ylabel(r'$M_h ({\rm M_\odot})$')
plt.xlim(0,13.9)
#plt.ylim(5E-4,5E-3)
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlabel(r'$z$')
plt.savefig('allhmass_50.png')

plt.clf()
z0mhx=np.repeat(allmhx[:,-1],100).reshape(len(allmhx),100)
z0mha=np.repeat(allmha[:,-1],100).reshape(len(allmha),100)
z0mhl=np.repeat(allmhl[:,-1],100).reshape(len(allmhl),100)
z0mhla=np.repeat(allmhla[:,-1],100).reshape(len(allmhla),100)

plt.plot(age,np.nanmedian(allmhx/z0mhx,axis=0),label='XMDs',lw=5,color='C0')
plt.plot(age,np.nanmedian(allmha/z0mha,axis=0),label='XMD Analogs',lw=5,color='C1')
plt.plot(age,np.nanmedian(allmhl/z0mhl,axis=0),label='LMDs',lw=5,color='C2')
plt.plot(age,np.nanmedian(allmhla/z0mhla,axis=0),label='LMD Analogs',lw=5,color='C3')
#plt.fill_between(age,np.nanmedian(allmetsx,axis=0)-np.nanstd(allmetsx,axis=0),np.nanmedian(allmetsx,axis=0)+np.nanstd(allmetsx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsa,axis=0)-np.nanstd(allmetsa,axis=0),np.nanmedian(allmetsa,axis=0)+np.nanstd(allmetsa,axis=0),color='C1',alpha=.3)
plt.legend()
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.yscale('log')
plt.ylabel(r'$M_h/M_{h,z=0}$')
plt.xlim(0,13.9)
#plt.ylim(5E-4,5E-3)
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlabel(r'$z$')
plt.savefig('allhgrowth_50.png')

plt.clf()
plt.plot(age,np.nanmedian(allmgx,axis=0),label='XMDs',lw=5,color='C0')
plt.plot(age,np.nanmedian(allmga,axis=0),label='XMD Analogs',lw=5,color='C1')
plt.plot(age,np.nanmedian(allmgl,axis=0),label='LMDs',lw=5,color='C2')
plt.plot(age,np.nanmedian(allmgla,axis=0),label='LMD Analogs',lw=5,color='C3')
#plt.fill_between(age,np.nanmedian(allmetsx,axis=0)-np.nanstd(allmetsx,axis=0),np.nanmedian(allmetsx,axis=0)+np.nanstd(allmetsx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsa,axis=0)-np.nanstd(allmetsa,axis=0),np.nanmedian(allmetsa,axis=0)+np.nanstd(allmetsa,axis=0),color='C1',alpha=.3)
plt.legend()
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.yscale('log')
plt.ylabel(r'$M_{\rm gas} ({\rm M_\odot})$')
plt.xlim(0,13.9)
#plt.ylim(5E-4,5E-3)
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlabel(r'$z$')
plt.savefig('allgasmass_50.png')

plt.clf()
z0mgx=np.repeat(allmgx[:,-1],100).reshape(len(allmgx),100)
z0mga=np.repeat(allmga[:,-1],100).reshape(len(allmga),100)
z0mgl=np.repeat(allmgl[:,-1],100).reshape(len(allmgl),100)
z0mgla=np.repeat(allmgla[:,-1],100).reshape(len(allmgla),100)

plt.plot(age,np.nanmedian(allmgx/z0mgx,axis=0),label='XMDs',lw=5,color='C0')
plt.plot(age,np.nanmedian(allmga/z0mga,axis=0),label='XMD Analogs',lw=5,color='C1')
plt.plot(age,np.nanmedian(allmgl/z0mgl,axis=0),label='LMDs',lw=5,color='C2')
plt.plot(age,np.nanmedian(allmgla/z0mgla,axis=0),label='LMD Analogs',lw=5,color='C3')
#plt.fill_between(age,np.nanmedian(allmetsx,axis=0)-np.nanstd(allmetsx,axis=0),np.nanmedian(allmetsx,axis=0)+np.nanstd(allmetsx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsa,axis=0)-np.nanstd(allmetsa,axis=0),np.nanmedian(allmetsa,axis=0)+np.nanstd(allmetsa,axis=0),color='C1',alpha=.3)
plt.legend()
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.yscale('log')
plt.ylabel(r'$M_{\rm gas}/M_{\rm gas, z=0}$')
plt.xlim(0,13.9)
#plt.ylim(5E-4,5E-3)
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlabel(r'$z$')
plt.savefig('allgasgrowth_50.png')


plt.clf()
plt.plot(age,np.nanmedian(allmstx,axis=0),label='XMDs',lw=5,color='C0')
plt.plot(age,np.nanmedian(allmsta,axis=0),label='XMD Analogs',lw=5,color='C1')
plt.plot(age,np.nanmedian(allmstl,axis=0),label='LMDs',lw=5,color='C2')
plt.plot(age,np.nanmedian(allmstla,axis=0),label='LMD Analogs',lw=5,color='C3')
#plt.fill_between(age,np.nanmedian(allmetsx,axis=0)-np.nanstd(allmetsx,axis=0),np.nanmedian(allmetsx,axis=0)+np.nanstd(allmetsx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsa,axis=0)-np.nanstd(allmetsa,axis=0),np.nanmedian(allmetsa,axis=0)+np.nanstd(allmetsa,axis=0),color='C1',alpha=.3)
plt.legend()
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.yscale('log')
plt.ylabel(r'$M_{\rm star} ({\rm M_\odot})$')
plt.xlim(0,13.9)
#plt.ylim(5E-4,5E-3)
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlabel(r'$z$')
plt.savefig('allstarmass_50.png')

plt.clf()
z0mstx=np.repeat(allmstx[:,-1],100).reshape(len(allmstx),100)
z0msta=np.repeat(allmsta[:,-1],100).reshape(len(allmsta),100)
z0mstl=np.repeat(allmstl[:,-1],100).reshape(len(allmstl),100)
z0mstla=np.repeat(allmstla[:,-1],100).reshape(len(allmstla),100)

print('formedbyzp5',np.interp(8.6,age.value,np.nanmedian(allmstx/z0mstx,axis=0)))
print('formedbyzp5',np.interp(8.6,age.value,np.nanmedian(np.vstack([allmstx/z0mstx,allmsta/z0msta,allmstl/z0mstl,allmstla/z0mstla]),axis=0)))
plt.plot(age,np.nanmedian(allmstx/z0mstx,axis=0),label='XMDs',lw=5,color='C0')
plt.plot(age,np.nanmedian(allmsta/z0msta,axis=0),label='XMD Analogs',lw=5,color='C1')
plt.plot(age,np.nanmedian(allmstl/z0mstl,axis=0),label='LMDs',lw=5,color='C2')
plt.plot(age,np.nanmedian(allmstla/z0mstla,axis=0),label='LMD Analogs',lw=5,color='C3')
#plt.fill_between(age,np.nanmedian(allmetsx,axis=0)-np.nanstd(allmetsx,axis=0),np.nanmedian(allmetsx,axis=0)+np.nanstd(allmetsx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsa,axis=0)-np.nanstd(allmetsa,axis=0),np.nanmedian(allmetsa,axis=0)+np.nanstd(allmetsa,axis=0),color='C1',alpha=.3)
plt.legend()
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.yscale('log')
plt.ylabel(r'$M_{\rm star}/M_{\rm star, z=0}$')
plt.xlim(0,13.9)
#plt.ylim(5E-4,5E-3)
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlabel(r'$z$')
plt.savefig('allstargrowth_50.png')
