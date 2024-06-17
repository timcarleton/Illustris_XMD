from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('../python/stylesheet.txt')

zs=np.append(np.array([0,.1]),np.arange(.25,6.1,.25))
files=['alltnggal_withr_wmetals_wvmax.fits','alltnggal_z0_1_withr_wmetals_wvmax.fits','tng100_z_25.fits','tng100_z_5.fits','tng100_z_75.fits']
for i in range(5,len(zs)):
    if zs[i] % 1 ==0.5:
        files.append('tng100_z{:d}_5.fits'.format(int(zs[i])))
    elif zs[i] %1 !=0:
        files.append('tng100_z{:d}_{:d}.fits'.format((int(zs[i])),int((zs[i]%1)*100)))
    else:
        files.append('tng100_z{:d}.fits'.format(int(zs[i])))

xmdab=[]
lmdab=[]
sfab=[]

masscut=3E7
sfcut=.00001
for i in range(len(files)):
    f=fits.open(files[i])
    wsf=np.where((f[1].data.mstar>masscut) & (f[1].data.sfr>sfcut) & (f[1].data.mhalo-f[1].data.mstar-f[1].data.mgas>1E8))[0]
    wx=np.where(f[1].data.gas_metal_sfr[wsf]<.1)[0]
    wl=np.where((f[1].data.gas_metal_sfr[wsf]>.1) & (f[1].data.gas_metal_sfr[wsf]<.2))[0]

    xmdab.append(len(wx))
    lmdab.append(len(wl))
    sfab.append(len(wsf))

print(xmdab,lmdab,sfab)

xmdab=np.array(xmdab)
lmdab=np.array(lmdab)
sfab=np.array(sfab)

zs50=np.append(np.array([0,.1,.31,.55,.75]),np.arange(1,5.6,.5))
files50=['tng50/alltng50gal_withr_withspin.fits','tng50_z_1.fits','tng50_z_31.fits','tng50_z_55.fits','tng50_z_75.fits']
for i in range(5,len(zs50)):
    if zs50[i] % 1 ==0.5:
        files50.append('tng50_z{:d}_5.fits'.format(int(zs50[i])))
    elif zs50[i] %1 !=0:
        files50.append('tng50_z{:d}_{:d}.fits'.format((int(zs50[i])),int((zs50[i]%1)*100)))
    else:
        files50.append('tng50_z{:d}.fits'.format(int(zs50[i])))

print(files)
xmdab50=[]
lmdab50=[]
sfab50=[]

masscut=3E7
sfcut=.00001
for i in range(len(files50)):

    f=fits.open(files50[i])
    wsf=np.where((f[1].data.mstar>masscut) & (f[1].data.sfr>sfcut) & (f[1].data.mhalo-f[1].data.mstar-f[1].data.mgas>1E8))[0]
    wx=np.where(f[1].data.gas_metal_sfr[wsf]<.21)[0]
    wl=np.where((f[1].data.gas_metal_sfr[wsf]>.21) & (f[1].data.gas_metal_sfr[wsf]<.42))[0]

    xmdab50.append(len(wx))
    lmdab50.append(len(wl))
    sfab50.append(len(wsf))

xmdab50=np.array(xmdab50)
lmdab50=np.array(lmdab50)
sfab50=np.array(sfab50)



plt.clf()
plt.plot(zs,xmdab/sfab,'C0*',label='XMDs')
plt.plot(zs,lmdab/sfab,'C1x',label='LMDs')
plt.xlabel(r'$z$')
plt.ylabel(r'Relative Abundance (N$_{\rm x}$/N$_{\rm SF}$)')
plt.legend()
plt.yscale('log')
plt.savefig('relabvsz.png')

plt.clf()

ax1=plt.gca()
sym1,=ax1.plot(zs,xmdab/(75/.6774)**3,'C0*')
sym2,=ax1.plot(zs,lmdab/(75/.6774)**3,'C1x')
sym3,=ax1.plot(zs,sfab/(75/.6774)**3,'k+')
ax1.set_ylabel(r'Abundance (Mpc$^{-3}$)')
#ax2=plt.twinx()
sym4,=ax1.plot(zs50,xmdab50/(50/.6774)**3,'C0v')
sym5,=ax1.plot(zs50,lmdab50/(50/.6774)**3,'C1^')
sym6,=ax1.plot(zs50,sfab50/(50/.6774)**3,'kP')


#ax2.set_yscale('log')
ax1.set_yscale('log')

ax1.set_ylim(1E-5,.1)
#ax2.set_ylim(1E-5,.1)
ax1.legend(handles=[sym1,sym4,sym2,sym5,sym3,sym6],labels=['XMDs TNG100','XMDs TNG50','LMDs TNG100','LMDs TNG50','SF Galaxies TNG100','SF Galaxies TNG50'])
ax1.set_xlabel(r'$z$')
#ax2.set_xlabel(r'$z$')
plt.savefig('abvsz_w50.png')
