from astropy.io import fits, ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy import cosmology
from scipy.optimize import minimize
astrocosmo = cosmology.FlatLambdaCDM(H0=67.7,Om0=0.307)

plt.style.use('../python/stylesheet.txt')

f=fits.open('alltnggal_withr_wmetals_wvmax.fits')
f50=fits.open('tng50/alltng50gal_withr_withspin.fits')

x=ascii.read('sampleselection_sfr5_notide_nomass.txt')
x50=ascii.read('sampleselection_tng50_sfr5_notide_nomass.txt')

zs=np.array([20.05,14.99,11.98,10.976,10.0,9.389,9.002,8.449,8.012,7.595,7.236,7.005,6.491,6.011,5.847,5.530,5.228,4.996,4.665,4.428,4.177,4.008,3.709,3.491,3.283,3.008,2.896,2.733,2.578,2.444,2.316,2.208,2.103,2.002,1.904,1.823,1.743,1.667,1.604,1.531,1.496,1.414,1.358,1.302,1.248,1.206,1.155,1.114,1.074,1.036,0.997,0.951,0.923,0.887,0.851,0.817,0.791,0.757,0.733,0.700,0.676,0.644,0.621,0.598,0.576,0.546,0.525,0.503,0.482,0.461,0.440,0.420,0.40,0.380,0.361,0.348,0.329,0.310,0.298,0.273,0.261,0.244,0.226,0.214,0.1973,0.1804,0.1693,0.1527,0.1419,0.1258,0.1099,0.0994,0.0839,0.0737,0.0585,0.0485,0.0337,0.0240,0.0095,0])
age=astrocosmo.age(zs)
dt=age[1:]-age[0:-1]

obs1=ascii.read('hunter2012.dat',delimiter=';')
obs2=ascii.read('realxmddata.csv')
#wobs1=np.where(obs1['VMag']<-13.7)[0]
mtol=1
obsmst=mtol*10**(-.4*(obs1['VMag']-4.8))
wobs1=np.where((obsmst>3E7) & (obsmst<3e8))[0]
obs1=obs1[wobs1]
obsmst1=obsmst[wobs1]

wxmdhunter=np.where(10**(obs1['[O/H]']-8.69)<.1)[0]
wnmdhunter=np.where(10**(obs1['[O/H]']-8.69)>=.1)[0]

wg=np.where((f[1].data.mstar<1E12) & (f[1].data.mstar>3E7) & (f[1].data.sfr>0) & (f[1].data.mhalo-f[1].data.mstar-f[1].data.mgas>1E8))[0]
wg50=np.where((f50[1].data.mstar<1E12) & (f50[1].data.mstar>3E7) & (f50[1].data.sfr>0))[0]


wxmd=x['allgalindex'][np.where(x['type']==1)[0]]
wlmd=x['allgalindex'][np.where(x['type']==2)[0]]
wxmda=x['allgalindex'][np.where(x['type']==3)[0]]
wlmda=x['allgalindex'][np.where(x['type']==4)[0]]

wxmd50=x50['allgalindex'][np.where(x50['type']==1)[0]]
wlmd50=x50['allgalindex'][np.where(x50['type']==2)[0]]
wxmda50=x50['allgalindex'][np.where(x50['type']==3)[0]]
wlmda50=x50['allgalindex'][np.where(x50['type']==4)[0]]

sfrsxt=np.loadtxt('allsfrx.txt')
sfrsxat=np.loadtxt('allsfrxa.txt')
sfrslt=np.loadtxt('allsfrl.txt')
sfrslat=np.loadtxt('allsfrla.txt')

mgxt=np.loadtxt('allmgx.txt')
mgxat=np.loadtxt('allmgxa.txt')
mglt=np.loadtxt('allmgl.txt')
mglat=np.loadtxt('allmgla.txt')

alleta=np.loadtxt('alleta_100.txt')
alleta50=np.loadtxt('alleta_50.txt')
avgsfrx=np.zeros([len(sfrsxt),len(sfrsxt[0])-1])
avgsfrxa=np.zeros([len(sfrsxat),len(sfrsxat[0])-1])
avgsfrl=np.zeros([len(sfrslt),len(sfrslt[0])-1])
avgsfrla=np.zeros([len(sfrslat),len(sfrslat[0])-1])

for i in range(len(sfrsxt[0])-1):
    avgsfrx[:,i]=np.nansum(sfrsxt[:,i:-1]*np.repeat(dt,len(sfrsxt)).reshape(len(sfrsxt),len(dt))[:,i:],axis=1)
    avgsfrxa[:,i]=np.nansum(sfrsxat[:,i:-1]*np.repeat(dt,len(sfrsxat)).reshape(len(sfrsxat),len(dt))[:,i:],axis=1)
    avgsfrl[:,i]=np.nansum(sfrslt[:,i:-1]*np.repeat(dt,len(sfrslt)).reshape(len(sfrslt),len(dt))[:,i:],axis=1)
    avgsfrla[:,i]=np.nansum(sfrslat[:,i:-1]*np.repeat(dt,len(sfrslat)).reshape(len(sfrslat),len(dt))[:,i:],axis=1)
    

#colddat=ascii.read('coldgasinfo_notide.txt')
colddat=ascii.read('densgasinfob_notide_nomass.txt')
colddat50=ascii.read('densgasinfo_50_notide.txt')
#colddat=ascii.read('sfgasinfo_notide.txt')
wxmdcold=np.where(colddat['type']==1)[0]
wlmdcold=np.where(colddat['type']==2)[0]
wxmdacold=np.where(colddat['type']==3)[0]
wlmdacold=np.where(colddat['type']==4)[0]

wxmdcold50=np.where(x50['type']==1)[0]
wlmdcold50=np.where(x50['type']==2)[0]
wxmdacold50=np.where(x50['type']==3)[0]
wlmdacold50=np.where(x50['type']==4)[0]


plt.clf()
plt.plot(colddat50['mcold'][wlmdacold50]/f50[1].data.mstar[wlmda50],f50[1].data.gas_metal_sfr[wlmda50],'C3h',label='LMD Analogs')
plt.plot(colddat50['mcold'][wxmdacold50]/f50[1].data.mstar[wxmda50],f50[1].data.gas_metal_sfr[wxmda50],'C2s',label='XMD Analogs')
plt.plot(colddat50['mcold'][wlmdcold50]/f50[1].data.mstar[wlmd50],f50[1].data.gas_metal_sfr[wlmd50],'C1x',label='LMDs')
plt.plot(colddat50['mcold'][wxmdcold50]/f50[1].data.mstar[wxmd50],f50[1].data.gas_metal_sfr[wxmd50],'C0*',label='XMDs')
plt.loglog()
plt.legend()
plt.xlabel(r'$M_{\rm gas}/M_*$')
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.title('TNG50',y=0.92)
plt.savefig('gasfracmetal50.png')

plt.clf()
plt.plot(alleta50[wlmdacold50],f50[1].data.gas_metal_sfr[wlmda50],'C3h',label='LMD Analogs')
plt.plot(alleta50[wxmdacold50],f50[1].data.gas_metal_sfr[wxmda50],'C2s',label='XMD Analogs')
plt.plot(alleta50[wlmdcold50],f50[1].data.gas_metal_sfr[wlmd50],'C1x',label='LMDs')
plt.plot(alleta50[wxmdcold50],f50[1].data.gas_metal_sfr[wxmd50],'C0*',label='XMDs')
plt.loglog()
plt.legend()
plt.xlabel(r'$\eta$')
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.title('TNG50',y=0.92)
plt.savefig('etametal50.png')

allavgsfr=np.vstack([avgsfrla,avgsfrl,avgsfrxa,avgsfrx])

wall=np.hstack([wlmda,wlmd,wxmda,wxmd])
wallcold=np.hstack([wlmdacold,wlmdcold,wxmdacold,wxmdcold])

corrall=[]

corrlmda=[]
corrlmd=[]
corrxmda=[]
corrxmd=[]
from scipy import stats
for i in range(len(avgsfrla[0])):
    corrall.append(stats.spearmanr(allavgsfr[:,i]/colddat['mcold'][wallcold]*1e9,f[1].data.gas_metal_sfr[wall])[0])
    corrlmda.append(stats.spearmanr(avgsfrla[:,i]/colddat['mcold'][wlmdacold]*1e9,f[1].data.gas_metal_sfr[wlmda])[0])
    corrlmd.append(stats.spearmanr(avgsfrl[:,i]/colddat['mcold'][wlmdcold]*1e9,f[1].data.gas_metal_sfr[wlmd])[0])
    corrxmda.append(stats.spearmanr(avgsfrxa[:,i]/colddat['mcold'][wxmdacold]*1e9,f[1].data.gas_metal_sfr[wxmda])[0])
    corrxmd.append(stats.spearmanr(avgsfrx[:,i]/colddat['mcold'][wxmdcold]*1e9,f[1].data.gas_metal_sfr[wxmd])[0])



ftmass=np.polyfit(np.log10(f[1].data.mstar[wall]),np.log10(f[1].data.gas_metal_sfr[wall]),1)
wf=np.where(np.isfinite(f[1].data.sfr[wall]/colddat['mcold'][wallcold]*1E9))[0]
wfg=np.where(np.isfinite(np.log10(allavgsfr[:,91]/colddat['mcold'][wallcold]*1E9)))[0]
print(wf)
print(allavgsfr.shape,allavgsfr[wf,91].shape)
ftsfe=np.polyfit(np.log10(f[1].data.sfr[wall[wf]]/colddat['mcold'][wallcold[wf]]*1E9),np.log10(f[1].data.gas_metal_sfr[wall[wf]]),1)
ftsfegyr=np.polyfit(np.log10(allavgsfr[wfg,91]/colddat['mcold'][wallcold[wfg]]*1E9),np.log10(f[1].data.gas_metal_sfr[wall[wfg]]),1)

ftsfecfrac=minimize(lambda ftparam:np.sum((np.log10(f[1].data.gas_metal_sfr[wall[wfg]])-(ftparam[0]+ftparam[1]*np.log10(allavgsfr[wfg,91]/colddat['mcold'][wallcold[wfg]]*1E9)+ftparam[2]*np.log10(colddat['mcold'][wallcold[wfg]]/f[1].data.mgas[wall[wfg]])))**2),[-1,1,1])['x']
print(ftmass,ftsfe,ftsfegyr,ftsfecfrac)

plt.clf()
plt.hist(np.log10(f[1].data.gas_metal_sfr[wall])-(ftmass[0]*np.log10(f[1].data.mstar[wall])+ftmass[1]),histtype='step',color='C0',lw=5,density=1,range=[-1,1],bins=40,label=r'$\log Z_{\rm SFR}=%.2f\log M_* %+.2f; \sigma=%.2f$' % (ftmass[0],ftmass[1],np.nanstd(np.log10(f[1].data.gas_metal_sfr[wall])-(ftmass[0]*np.log10(f[1].data.mstar[wall])+ftmass[1]))))
plt.hist(np.log10(f[1].data.gas_metal_sfr[wall])-(ftsfe[0]*np.log10(f[1].data.sfr[wall]/colddat['mcold'][wallcold]*1E9)+ftsfe[1]),histtype='step',color='C1',lw=5,density=1,ls='--',range=[-1,1],bins=40,label=r'$\log Z_{\rm SFR}=%.2f\log ({\rm SFR}/M_{\rm cold}) %+.2f; \sigma=%.2f$' % (ftsfe[0],ftsfe[1],np.nanstd(np.log10(f[1].data.gas_metal_sfr[wall[wfg]])-(ftsfe[0]*np.log10(f[1].data.sfr[wall[wfg]]/colddat['mcold'][wallcold[wfg]]*1E9)+ftsfe[1]))))
plt.hist(np.log10(f[1].data.gas_metal_sfr[wall])-(ftsfegyr[0]*np.log10(allavgsfr[:,91]/colddat['mcold'][wallcold]*1E9)+ftsfegyr[1]),histtype='step',color='C2',lw=5,density=1,ls=':',range=[-1,1],bins=40,label=r'$\log Z_{\rm SFR}=%.2f\log (<{\rm SFR}>_{\rm 1.3~Gyr}/M_{\rm cold}) %+.2f; \sigma=%.2f$' % (ftsfegyr[0],ftsfegyr[1],np.nanstd(np.log10(f[1].data.gas_metal_sfr[wall[wfg]])-(ftsfegyr[0]*np.log10(allavgsfr[wfg,91]/colddat['mcold'][wallcold[wfg]]*1E9)+ftsfegyr[1]))))
plt.hist(np.log10(f[1].data.gas_metal_sfr[wall])-(ftsfecfrac[0]+ftsfecfrac[1]*np.log10(allavgsfr[:,91]/colddat['mcold'][wallcold]*1E9)+ftsfecfrac[2]*np.log10(colddat['mcold'][wallcold]/f[1].data.mgas[wall])),histtype='step',color='C3',lw=5,density=1,ls='-.',range=[-1,1],bins=40,label=r'$\log Z_{\rm SFR}=%.2f\log (<{\rm SFR}>_{\rm 1.3~Gyr}/M_{\rm cold}) %+.2f \log(M_{\rm dense gas}/M_{\rm total gas}) %+.2f; \sigma=%.2f$' % (ftsfecfrac[1],ftsfecfrac[2],ftsfecfrac[0],np.nanstd(np.log10(f[1].data.gas_metal_sfr[wall[wfg]])-(ftsfecfrac[0]+ftsfecfrac[1]*np.log10(allavgsfr[wfg,91]/colddat['mcold'][wallcold[wfg]]*1E9)+ftsfecfrac[2]*np.log10(colddat['mcold'][wallcold[wfg]]/f[1].data.mgas[wall[wfg]])))))
plt.legend(fontsize=25)
plt.ylim(0,2.5)
plt.xlabel(r'$\log Z_{\rm SFR}/Z_\odot$ - Model')
plt.ylabel(r'$dN/d\log Z_{\rm SFR}/Z_\odot$')
plt.savefig('modelhist.png')

    
plt.clf()
plt.plot(age[-1]-age[0:-1],corrlmda,'C3h',label='LMD Analog')
plt.plot(age[-1]-age[0:-1],corrlmd,'C1x',label='LMD')
plt.plot(age[-1]-age[0:-1],corrxmda,'C2s',label='XMD Analog')
plt.plot(age[-1]-age[0:-1],corrxmd,'C0*',label='XMD')
plt.legend()
plt.xlabel('SFE timescale')
plt.ylabel('Z-SFE correlation coefficient')
plt.savefig('sfecorr.png')

plt.clf()
plt.plot(age[-1]-age[0:-1],corrall,'C0o')
plt.xlabel(r'$t$ Gyr')
plt.ylabel(r'Correlation Coeff.: $\frac{M_*{\rm ~formed~in~last}~t~{\rm Gyr}}{M_{\rm cold~gas,~z=0}}-Z_{\rm SFR}$')
plt.savefig('sfecorrall.png')

print('all',np.array(age[-1]-age[0:-1])[np.argmax(np.array(corrall))])
print('ax',np.array(age[-1]-age[0:-1])[np.argmax(np.array(corrxmd))])
print('axa',np.array(age[-1]-age[0:-1])[np.argmax(np.array(corrxmda))])
print('al',np.array(age[-1]-age[0:-1])[np.argmax(np.array(corrlmd))])
print('ala',np.array(age[-1]-age[0:-1])[np.argmax(np.array(corrlmda))])

plt.clf()
plt.plot(avgsfrla[:,91]/colddat['mcold'][wlmdacold]*1e9,f[1].data.gas_metal_sfr[wlmda],'C3h',alpha=.35,label='LMD Analog')
plt.plot(avgsfrl[:,91]/colddat['mcold'][wlmdcold]*1e9,f[1].data.gas_metal_sfr[wlmd],'C1x',alpha=.75,label='LMD')
plt.plot(avgsfrxa[:,91]/colddat['mcold'][wxmdacold]*1e9,f[1].data.gas_metal_sfr[wxmda],'C2s',alpha=.35,label='XMD Analog')
plt.plot(avgsfrx[:,91]/colddat['mcold'][wxmdcold]*1e9,f[1].data.gas_metal_sfr[wxmd],'C0*',alpha=.75,label='XMD')
plt.xlabel(r'$<SFR>_{\rm  1.5~Gyr}/M_{\rm cold~gas}$ (Gyr$^{-1})$')
plt.loglog()
plt.legend()
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.savefig('coldgasplot_gyr.png')

plt.clf()
plt.plot(avgsfrla[:,91]/colddat['mcold'][wlmdacold]*1e9,f[1].data.gas_metal_mass[wlmda],'C3h',alpha=.35,label='LMD Analog')
plt.plot(avgsfrl[:,91]/colddat['mcold'][wlmdcold]*1e9,f[1].data.gas_metal_mass[wlmd],'C1x',alpha=.75,label='LMD')
plt.plot(avgsfrxa[:,91]/colddat['mcold'][wxmdacold]*1e9,f[1].data.gas_metal_mass[wxmda],'C2s',alpha=.35,label='XMD Analog')
plt.plot(avgsfrx[:,91]/colddat['mcold'][wxmdcold]*1e9,f[1].data.gas_metal_mass[wxmd],'C0*',alpha=.75,label='XMD')
plt.xlabel(r'$<SFR>_{\rm  1.3~Gyr}/M_{\rm cold~gas}$ (Gyr$^{-1})$')
plt.loglog()
plt.legend()
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.savefig('coldgasplot_mass_gyr.png')


plt.clf()
plt.scatter(f[1].data.sfr[wlmda]/colddat['mcold'][wlmdacold]*1e9,f[1].data.gas_metal_sfr[wlmda],c=np.log10(mglat[:,99]/colddat['mcold'][wlmdacold]),alpha=.35,vmin=1,vmax=2.5)
plt.scatter(f[1].data.sfr[wlmd]/colddat['mcold'][wlmdcold]*1e9,f[1].data.gas_metal_sfr[wlmd],c=np.log10(mglt[:,99]/colddat['mcold'][wlmdcold]),alpha=.35,vmin=1,vmax=2.5)
plt.scatter(f[1].data.sfr[wxmda]/colddat['mcold'][wxmdacold]*1e9,f[1].data.gas_metal_sfr[wxmda],c=np.log10(mgxat[:,99]/colddat['mcold'][wxmdacold]),alpha=.35,vmin=1,vmax=2.5)
plt.scatter(f[1].data.sfr[wxmd]/colddat['mcold'][wxmdcold]*1e9,f[1].data.gas_metal_sfr[wxmd],c=np.log10(mgxt[:,99]/colddat['mcold'][wxmdcold]),alpha=.35,vmin=1,vmax=2.5)
plt.xlabel(r'$SFR/M_{\rm cold~gas}$ (Gyr$^{-1})$')
plt.loglog()
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.colorbar(label=r'$\log(M_{\rm gas;~total}/M_{\rm gas;~dense})$')
plt.savefig('coldgasplot_ratio.png')

plt.clf()
plt.scatter(f[1].data.sfr[wlmda]/colddat['mcold'][wlmdacold]*1e9,f[1].data.gas_metal_sfr[wlmda],c=np.log10(alleta[wlmdacold]),alpha=.35,vmin=1,vmax=2.5)
plt.scatter(f[1].data.sfr[wlmd]/colddat['mcold'][wlmdcold]*1e9,f[1].data.gas_metal_sfr[wlmd],c=np.log10(alleta[wlmdcold]),alpha=.35,vmin=1,vmax=2.5)
plt.scatter(f[1].data.sfr[wxmda]/colddat['mcold'][wxmdacold]*1e9,f[1].data.gas_metal_sfr[wxmda],c=np.log10(alleta[wxmdacold]),alpha=.35,vmin=1,vmax=2.5)
plt.scatter(f[1].data.sfr[wxmd]/colddat['mcold'][wxmdcold]*1e9,f[1].data.gas_metal_sfr[wxmd],c=np.log10(alleta[wxmdcold]),alpha=.35,vmin=1,vmax=2.5)
plt.xlabel(r'$SFR/M_{\rm cold~gas}$ (Gyr$^{-1})$')
plt.loglog()
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.colorbar(label=r'$\eta$')
plt.savefig('coldgasplot_eta.png')

plt.clf()
plt.scatter(f50[1].data.sfr[wlmda50]/f50[1].data.mgas[wlmda50]*1e9,f50[1].data.gas_metal_sfr[wlmda50],c=np.log10(alleta50[wlmdacold50]),alpha=.35,vmin=1,vmax=2.5)
plt.scatter(f50[1].data.sfr[wlmd50]/f50[1].data.mgas[wlmd50]*1e9,f50[1].data.gas_metal_sfr[wlmd50],c=np.log10(alleta50[wlmdcold50]),alpha=.35,vmin=1,vmax=2.5)
plt.scatter(f50[1].data.sfr[wxmda50]/f50[1].data.mgas[wxmda50]*1e9,f50[1].data.gas_metal_sfr[wxmda50],c=np.log10(alleta50[wxmdacold50]),alpha=.35,vmin=1,vmax=2.5)
plt.scatter(f50[1].data.sfr[wxmd50]/f50[1].data.mgas[wxmd50]*1e9,f50[1].data.gas_metal_sfr[wxmd50],c=np.log10(alleta50[wxmdcold50]),alpha=.35,vmin=1,vmax=2.5)
plt.xlabel(r'$SFR/M_{\rm gas}$ (Gyr$^{-1})$')
plt.loglog()
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.colorbar(label=r'$\eta$')
plt.savefig('gaseffplot_eta50.png')

alletaa=np.hstack([alleta50[wlmdacold50],alleta50[wlmdcold50],alleta50[wxmdacold50],alleta50[wxmdcold50]])
allmet50=np.hstack([f50[1].data.gas_metal_sfr[wlmda50],f50[1].data.gas_metal_sfr[wlmd50],f50[1].data.gas_metal_sfr[wxmda50],f50[1].data.gas_metal_sfr[wxmd50]])
wnn=np.where(np.isfinite(np.log10(alletaa)))[0]
print('cov met eta',np.cov(np.array([np.log10(allmet50[wnn]),np.log10(alletaa[wnn])])))

plt.clf()
plt.scatter(alleta50[wlmdacold50],f50[1].data.gas_metal_sfr[wlmda50],c=np.log10(f50[1].data.sfr[wlmda50]/f50[1].data.mgas[wlmda50]*1e9),alpha=.35,vmin=-4,vmax=-1)
plt.scatter(alleta50[wlmdcold50],f50[1].data.gas_metal_sfr[wlmd50],c=np.log10(f50[1].data.sfr[wlmd50]/f50[1].data.mgas[wlmd50]*1e9),alpha=.35,vmin=-4,vmax=-1)
plt.scatter(alleta50[wxmdacold50],f50[1].data.gas_metal_sfr[wxmda50],c=np.log10(f50[1].data.sfr[wxmda50]/f50[1].data.mgas[wxmda50]*1e9),alpha=.35,vmin=-4,vmax=-1)
plt.scatter(alleta50[wxmdcold50],f50[1].data.gas_metal_sfr[wxmd50],c=np.log10(f50[1].data.sfr[wxmd50]/f50[1].data.mgas[wxmd50]*1e9),alpha=.35,vmin=-4,vmax=-1)
plt.xlabel(r'$\eta$')
plt.loglog()
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.colorbar(label=r'$SFR/M_{\rm gas}$ (Gyr$^{-1})$')
plt.savefig('gasetaplot_eff50.png')

plt.clf()
plt.plot(alleta[wlmdacold],f[1].data.sfr[wlmda]/colddat['mcold'][wlmdacold]*1E9,'C3h',alpha=.5)
plt.plot(alleta[wlmdcold],f[1].data.sfr[wlmd]/colddat['mcold'][wlmdcold]*1E9,'C1x',alpha=.5)
plt.plot(alleta[wxmdacold],f[1].data.sfr[wxmda]/colddat['mcold'][wxmdacold]*1E9,'C2s',alpha=.5)
plt.plot(alleta[wxmdcold],f[1].data.sfr[wxmd]/colddat['mcold'][wxmdcold]*1E9,'C0*',alpha=.5)
plt.loglog()
plt.ylabel(r'$SFR/M_{\rm cold~gas}$ (Gyr$^{-1})$')
plt.xlabel(r'$\eta$')
plt.savefig('etasfe.png')

plt.clf()
plt.scatter(alleta[wlmdacold],f[1].data.gas_metal_sfr[wlmda],c=np.log10(f[1].data.sfr[wlmda]/colddat['mcold'][wlmdacold]*1e9),alpha=.35,vmin=-2,vmax=-.3)
plt.scatter(alleta[wlmdcold],f[1].data.gas_metal_sfr[wlmd],c=np.log10(f[1].data.sfr[wlmd]/colddat['mcold'][wlmdcold]*1e9),alpha=.35,vmin=-2,vmax=-.3)
plt.scatter(alleta[wxmdacold],f[1].data.gas_metal_sfr[wxmda],c=np.log10(f[1].data.sfr[wxmda]/colddat['mcold'][wxmdacold]*1e9),alpha=.35,vmin=-2,vmax=-.3)
plt.scatter(alleta[wxmdcold],f[1].data.gas_metal_sfr[wxmd],c=np.log10(f[1].data.sfr[wxmd]/colddat['mcold'][wxmdcold]*1e9),alpha=.35,vmin=-2,vmax=-.3)
plt.xlabel(r'$\eta$')
plt.loglog()
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.colorbar(label=r'$SFR/M_{\rm cold~gas}$ (Gyr$^{-1})$')
plt.savefig('coldgasplot_eta2.png')

plt.clf()
plt.scatter(alleta[wlmdacold],f[1].data.gas_metal_sfr[wlmda],c=np.log10(f[1].data.sfr[wlmda]/colddat['mcold'][wlmdacold]*1e9),alpha=.35,vmin=-2,vmax=-.3)
plt.scatter(alleta[wxmdacold],f[1].data.gas_metal_sfr[wxmda],c=np.log10(f[1].data.sfr[wxmda]/colddat['mcold'][wxmdacold]*1e9),alpha=.35,vmin=-2,vmax=-.3)
plt.xlabel(r'$\eta$')
plt.loglog()
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')

plt.xlim(5,2E4)
plt.ylim(.01,2)
plt.colorbar(label=r'$SFR/M_{\rm cold~gas}$ (Gyr$^{-1})$')
plt.savefig('coldgasplot_eta_sep1.png')

plt.clf()
plt.scatter(alleta[wlmdcold],f[1].data.gas_metal_sfr[wlmd],c=np.log10(f[1].data.sfr[wlmd]/colddat['mcold'][wlmdcold]*1e9),alpha=.35,vmin=-2,vmax=-.3)
plt.scatter(alleta[wxmdcold],f[1].data.gas_metal_sfr[wxmd],c=np.log10(f[1].data.sfr[wxmd]/colddat['mcold'][wxmdcold]*1e9),alpha=.35,vmin=-2,vmax=-.3)
plt.xlabel(r'$\eta$')
plt.loglog()
plt.xlim(5,2E4)
plt.ylim(.01,2)

plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.colorbar(label=r'$SFR/M_{\rm cold~gas}$ (Gyr$^{-1})$')
plt.savefig('coldgasplot_eta_sep2.png')

plt.clf()
#plt.plot(.02/(1+alleta[wlmdacold]/.8+f[1].data.mgas[wlmda]/f[1].data.mstar[wlmda]),f[1].data.gas_metal_sfr[wlmda],'o')
#plt.plot(.02/(1+alleta[wlmdcold]/.8+f[1].data.mgas[wlmd]/f[1].data.mstar[wlmd]),f[1].data.gas_metal_sfr[wlmd],'o')
#plt.plot(.02/(1+alleta[wxmdacold]/.8+f[1].data.mgas[wxmda]/f[1].data.mstar[wxmda]),f[1].data.gas_metal_sfr[wxmda],'o')
#plt.plot(.02/(1+alleta[wxmdcold]/.8+f[1].data.mgas[wxmd]/f[1].data.mstar[wxmd]),f[1].data.gas_metal_sfr[wxmd],'o')

print(f[1].data.mgas[wlmda]/avgsfrla[:,91])
plt.plot(2/(1+alleta[wlmdacold]/.8+f[1].data.mgas[wlmda]/avgsfrla[:,91]/1E9),f[1].data.gas_metal_sfr[wlmda],'o')
plt.plot(2/(1+alleta[wlmdcold]/.8+f[1].data.mgas[wlmd]/avgsfrl[:,91]/1E9),f[1].data.gas_metal_sfr[wlmd],'o')
plt.plot(2/(1+alleta[wxmdacold]/.8+f[1].data.mgas[wxmda]/avgsfrxa[:,91]/1E9),f[1].data.gas_metal_sfr[wxmda],'o')
plt.plot(2/(1+alleta[wxmdcold]/.8+f[1].data.mgas[wxmd]/avgsfrx[:,91]/1E9),f[1].data.gas_metal_sfr[wxmd],'o')
plt.xlabel(r'pred')
plt.loglog()
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.savefig('coldgasplot_model1.png')

plt.clf()

#plt.scatter(colddat['mcold'][wlmdacold]/mglat[:,99],f[1].data.gas_metal_mass[wlmda]*f[1].data.mgas[wlmda],c=np.log10(avgsfrla[:,91]/colddat['mcold'][wlmdacold]*1e9),vmin=-2,vmax=0,alpha=.35)
#plt.scatter(colddat['mcold'][wlmdcold]/mglt[:,99],f[1].data.gas_metal_mass[wlmd]*f[1].data.mgas[wlmd],c=np.log10(avgsfrl[:,91]/colddat['mcold'][wlmdcold]*1e9),vmin=-2,vmax=0,alpha=.35)
#plt.scatter(colddat['mcold'][wxmdacold]/mgxat[:,99],f[1].data.gas_metal_mass[wxmda]*f[1].data.mgas[wxmda],c=np.log10(avgsfrxa[:,91]/colddat['mcold'][wxmdacold]*1e9),vmin=-2,vmax=0,alpha=.35)
#plt.scatter(colddat['mcold'][wxmdcold]/mgxt[:,99],f[1].data.gas_metal_mass[wxmd]*f[1].data.mgas[wxmd],c=np.log10(avgsfrx[:,91]/colddat['mcold'][wxmdcold]*1e9),vmin=-2,vmax=0,alpha=.35)

plt.scatter(colddat['mcold'][wlmdacold]/mglat[:,99],f[1].data.gas_metal_sfr[wlmda],c=np.log10(f[1].data.sfr[wlmda]),vmin=-4,vmax=-2,alpha=.35)
plt.scatter(colddat['mcold'][wlmdcold]/mglt[:,99],f[1].data.gas_metal_sfr[wlmd],c=np.log10(f[1].data.sfr[wlmd]),vmin=-4,vmax=2,alpha=.35)
plt.scatter(colddat['mcold'][wxmdacold]/mgxat[:,99],f[1].data.gas_metal_sfr[wxmda],c=np.log10(f[1].data.sfr[wxmda]),vmin=-4,vmax=-2,alpha=.35)
plt.scatter(colddat['mcold'][wxmdcold]/mgxt[:,99],f[1].data.gas_metal_sfr[wxmd],c=np.log10(f[1].data.sfr[wxmd]),vmin=-4,vmax=-2,alpha=.35)

plt.xlabel('coldgas fraction')
plt.ylabel('metal mass')
plt.loglog()
plt.savefig('metalmassfrac.png')

plt.clf()
plt.scatter(avgsfrla[:,91]/colddat['mcold'][wlmdacold]*1e9,f[1].data.gas_metal_mass[wlmda]*f[1].data.mgas[wlmda],c=np.log10(colddat['mcold'][wlmdacold]/mglat[:,99]),alpha=.35,vmin=-2.5,vmax=-1)
plt.scatter(avgsfrl[:,91]/colddat['mcold'][wlmdcold]*1e9,f[1].data.gas_metal_mass[wlmd]*f[1].data.mgas[wlmd],c=np.log10(colddat['mcold'][wlmdcold]/mglt[:,99]),alpha=.35,vmin=-2.5,vmax=-1)
plt.scatter(avgsfrxa[:,91]/colddat['mcold'][wxmdacold]*1e9,f[1].data.gas_metal_mass[wxmda]*f[1].data.mgas[wxmda],c=np.log10(colddat['mcold'][wxmdacold]/mgxat[:,99]),alpha=.35,vmin=-2.5,vmax=-1)
plt.scatter(avgsfrx[:,91]/colddat['mcold'][wxmdcold]*1e9,f[1].data.gas_metal_mass[wxmd]*f[1].data.mgas[wxmd],c=np.log10(colddat['mcold'][wxmdcold]/mgxt[:,99]),alpha=.35,vmin=-2.5,vmax=-1)
plt.xlabel(r'$SFR/M_{\rm cold~gas}$ (Gyr$^{-1})$')
plt.loglog()
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.colorbar(label=r'$\log(M_{\rm gas;~dense}/M_{\rm gas;~total})$')
plt.savefig('coldgasplot_ratio_mass.png')


plt.clf()
plt.scatter(avgsfrla[:,91]/colddat['mcold'][wlmdacold]*1e9,f[1].data.gas_metal_sfr[wlmda],c=np.log10(colddat['mcold'][wlmdacold]/mglat[:,99]),alpha=.35,vmin=-2.5,vmax=-1)
plt.scatter(avgsfrl[:,91]/colddat['mcold'][wlmdcold]*1e9,f[1].data.gas_metal_sfr[wlmd],c=np.log10(colddat['mcold'][wlmdcold]/mglt[:,99]),alpha=.35,vmin=-2.5,vmax=-1)
plt.scatter(avgsfrxa[:,91]/colddat['mcold'][wxmdacold]*1e9,f[1].data.gas_metal_sfr[wxmda],c=np.log10(colddat['mcold'][wxmdacold]/mgxat[:,99]),alpha=.35,vmin=-2.5,vmax=-1)
plt.scatter(avgsfrx[:,91]/colddat['mcold'][wxmdcold]*1e9,f[1].data.gas_metal_sfr[wxmd],c=np.log10(colddat['mcold'][wxmdcold]/mgxt[:,99]),alpha=.35,vmin=-2.5,vmax=-1)
plt.xlabel(r'$<SFR>_{\rm  1.3~Gyr}/M_{\rm cold~gas}$ (Gyr$^{-1})$')
plt.loglog()
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.colorbar(label=r'$\log(M_{\rm gas;~total}/M_{\rm gas;~dense})$')
plt.savefig('coldgasplot_gyr_ratio.png')

plt.clf()
#plt.scatter(avgsfrla[:,91]/colddat['mcold'][wlmdacold]*1e9,f[1].data.gas_metal_sfr[wlmda],c=np.log10(avgsfrla[:,98]/colddat['mcold'][wlmdacold]*1e9),alpha=.35,vmin=-2.5,vmax=-1)
#plt.scatter(avgsfrl[:,91]/colddat['mcold'][wlmdcold]*1e9,f[1].data.gas_metal_sfr[wlmd],c=np.log10(avgsfrl[:,98]/colddat['mcold'][wlmdcold]*1e9),alpha=.35,vmin=-2.5,vmax=-1)
#plt.scatter(avgsfrxa[:,91]/colddat['mcold'][wxmdacold]*1e9,f[1].data.gas_metal_sfr[wxmda],c=np.log10(avgsfrxa[:,98]/colddat['mcold'][wxmdacold]*1e9),alpha=.35,vmin=-2.5,vmax=-1)
#plt.scatter(avgsfrx[:,91]/colddat['mcold'][wxmdcold]*1e9,f[1].data.gas_metal_sfr[wxmd],c=np.log10(avgsfrx[:,98]/colddat['mcold'][wxmdcold]*1e9),alpha=.35,vmin=-2.5,vmax=-1)

plt.scatter(avgsfrla[:,91]/colddat['mcold'][wlmdacold]*1e9,f[1].data.gas_metal_sfr[wlmda],c=np.log10(avgsfrla[:,91]/avgsfrla[:,98]),alpha=.35,vmin=-1,vmax=1)
plt.scatter(avgsfrl[:,91]/colddat['mcold'][wlmdcold]*1e9,f[1].data.gas_metal_sfr[wlmd],c=np.log10(avgsfrl[:,91]/avgsfrl[:,98]),alpha=.35,vmin=-1,vmax=1)
plt.scatter(avgsfrxa[:,91]/colddat['mcold'][wxmdacold]*1e9,f[1].data.gas_metal_sfr[wxmda],c=np.log10(avgsfrxa[:,91]/avgsfrxa[:,98]),alpha=.35,vmin=-1,vmax=1)
plt.scatter(avgsfrx[:,91]/colddat['mcold'][wxmdcold]*1e9,f[1].data.gas_metal_sfr[wxmd],c=np.log10(avgsfrx[:,91]/avgsfrx[:,98]),alpha=.35,vmin=-1,vmax=1)
plt.xlabel(r'$<SFR>_{\rm  1.3~Gyr}/M_{\rm cold~gas}$ (Gyr$^{-1})$')
plt.loglog()
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.colorbar(label=r'$\log({\rm <SFR>_{1.3~Gyr}/<SFR>_{0~Gyr}})$')
plt.savefig('coldgasplot_gyr_ratio2.png')


plt.clf()
plt.plot(f[1].data.sfr[wlmda]/colddat['mcold'][wlmdacold]*1e9,f[1].data.gas_metal_sfr[wlmda],'C3h',alpha=.35)
plt.plot(f[1].data.sfr[wlmd]/colddat['mcold'][wlmdcold]*1e9,f[1].data.gas_metal_sfr[wlmd],'C1x',alpha=.75)
plt.plot(f[1].data.sfr[wxmda]/colddat['mcold'][wxmdacold]*1e9,f[1].data.gas_metal_sfr[wxmda],'C2s',alpha=.35)
plt.plot(f[1].data.sfr[wxmd]/colddat['mcold'][wxmdcold]*1e9,f[1].data.gas_metal_sfr[wxmd],'C0*',alpha=.75)
plt.errorbar(obs2['SFR_Ha [Msol]']/obs2['HI mass [Msol]']*1E9,10**(obs2['12 +(O/H)']-8.69),color='violet',ms=10,marker='*',xerr=np.sqrt((obs2['SFR err']/obs2['HI mass [Msol]'])**2 + (obs2['HI err']*obs2['SFR_Ha [Msol]']/obs2['HI mass [Msol]']**2)**2)*1E9,yerr=obs2['Z err']*10**(obs2['12 +(O/H)']-8.69)*np.log(10),ls='None',label='Other XMDs')

#plt.errorbar(10**-obs1['MHI'][wnmdhunter]*10**obs1['logSFR2'][wnmdhunter]*1E9,10**(obs1['[O/H]'][wnmdhunter]-8.69),markeredgecolor='purple',ms=10,marker='d',xerr=obs1['e_logSFR2'][wnmdhunter]*10**-obs1['MHI'][wnmdhunter]*10**obs1['logSFR2'][wnmdhunter]*np.log(10)*1E9,ls='None',markerfacecolor='None',markeredgewidth=2)
#plt.errorbar(10**-obs1['MHI'][wxmdhunter]*10**obs1['logSFR2'][wxmdhunter]*1E9,10**(obs1['[O/H]'][wxmdhunter]-8.69),color='purple',ms=10,marker='d',xerr=obs1['e_logSFR2'][wxmdhunter]*10**-obs1['MHI'][wxmdhunter]*10**obs1['logSFR2'][wxmdhunter]*np.log(10)*1E9,ls='None',markeredgewidth=2)

plt.plot(10**-obs1['MHI'][wnmdhunter]*10**obs1['logSFR2'][wnmdhunter]*1E9,10**(obs1['[O/H]'][wnmdhunter]-8.69),markeredgecolor='purple',ms=10,marker='d',ls='None',markerfacecolor='None',markeredgewidth=2,label='Hunter+12 NMDs')
plt.plot(10**-obs1['MHI'][wxmdhunter]*10**obs1['logSFR2'][wxmdhunter]*1E9,10**(obs1['[O/H]'][wxmdhunter]-8.69),color='purple',ms=10,marker='d',ls='None',markeredgewidth=2,label='Hunter+12 XMDs')
plt.legend()
print(len(np.where(colddat['mcold'][wlmdacold]/1e9/f[1].data.sfr[wlmda]>50)[0]))
print(len(np.where(colddat['mcold'][wxmdacold]/1e9/f[1].data.sfr[wxmda]>50)[0]))
print(len(np.where(colddat['mcold'][wlmdcold]/1e9/f[1].data.sfr[wlmd]>50)[0]))
print(len(np.where(colddat['mcold'][wxmdcold]/1e9/f[1].data.sfr[wxmd]>50)[0]))

plt.loglog()
plt.xlabel(r'SFR/$M_{\rm cold~gas}$ (Gyr$^{-1}$)')
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.savefig('coldgasplot.png')

plt.clf()
plt.scatter(colddat['mcold'][wlmdacold]/1e9/f[1].data.sfr[wlmda],f[1].data.gas_metal_sfr[wlmda],c=np.log10(colddat['mcold'][wlmdacold]/f[1].data.mgas[wlmda]),alpha=.35,vmin=-3,vmax=-1)
plt.scatter(colddat['mcold'][wlmdcold]/1e9/f[1].data.sfr[wlmd],f[1].data.gas_metal_sfr[wlmd],c=np.log10(colddat['mcold'][wlmdcold]/f[1].data.mgas[wlmd]),alpha=.75,vmin=-3,vmax=-1)
plt.scatter(colddat['mcold'][wxmdacold]/1e9/f[1].data.sfr[wxmda],f[1].data.gas_metal_sfr[wxmda],c=np.log10(colddat['mcold'][wxmdacold]/f[1].data.mgas[wxmda]),alpha=.35,vmin=-3,vmax=-1)
plt.scatter(colddat['mcold'][wxmdcold]/1e9/f[1].data.sfr[wxmd],f[1].data.gas_metal_sfr[wxmd],c=np.log10(colddat['mcold'][wxmdcold]/f[1].data.mgas[wxmd]),alpha=.75,vmin=-3,vmax=-1)
plt.errorbar(obs2['HI mass [Msol]']/obs2['SFR_Ha [Msol]']/1E9,10**(obs2['12 +(O/H)']-8.69),color='violet',ms=10,marker='^',xerr=np.sqrt((obs2['HI err']/obs2['SFR_Ha [Msol]'])**2 + (obs2['SFR err']*obs2['HI mass [Msol]']/obs2['SFR_Ha [Msol]']**2)**2)/1E9,yerr=obs2['Z err']*10**(obs2['12 +(O/H)']-8.69)*np.log(10),ls='None')
plt.errorbar(10**obs1['MHI']*10**-obs1['logSFR2']/1E9,10**(obs1['[O/H]']-8.69),color='purple',ms=10,marker='d',xerr=obs1['e_logSFR2']*10**obs1['MHI']*10**-obs1['logSFR2']*np.log(10)/1E9,ls='None')

plt.loglog()
plt.xlabel(r'$M_{\rm cold}$/SFR (Gyr)')
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.savefig('coldgasscatterplot.png')

plt.clf()
plt.scatter(colddat['mcold'][wlmdacold]/1e9/f[1].data.sfr[wlmda],f[1].data.gas_metal_sfr[wlmda],c=np.log10(f[1].data.sfr[wlmda]),alpha=.35,vmin=-3,vmax=-1)
plt.scatter(colddat['mcold'][wlmdcold]/1e9/f[1].data.sfr[wlmd],f[1].data.gas_metal_sfr[wlmd],c=np.log10(f[1].data.sfr[wlmd]),alpha=.75,vmin=-3,vmax=-1)
plt.scatter(colddat['mcold'][wxmdacold]/1e9/f[1].data.sfr[wxmda],f[1].data.gas_metal_sfr[wxmda],c=np.log10(f[1].data.sfr[wxmda]),alpha=.35,vmin=-3,vmax=-1)
plt.scatter(colddat['mcold'][wxmdcold]/1e9/f[1].data.sfr[wxmd],f[1].data.gas_metal_sfr[wxmd],c=np.log10(f[1].data.sfr[wxmd]),alpha=.75,vmin=-3,vmax=-1)
plt.errorbar(obs2['HI mass [Msol]']/obs2['SFR_Ha [Msol]']/1E9,10**(obs2['12 +(O/H)']-8.69),color='violet',ms=10,marker='^',xerr=np.sqrt((obs2['HI err']/obs2['SFR_Ha [Msol]'])**2 + (obs2['SFR err']*obs2['HI mass [Msol]']/obs2['SFR_Ha [Msol]']**2)**2)/1E9,yerr=obs2['Z err']*10**(obs2['12 +(O/H)']-8.69)*np.log(10),ls='None')
plt.errorbar(10**obs1['MHI']*10**-obs1['logSFR2']/1E9,10**(obs1['[O/H]']-8.69),color='purple',ms=10,marker='d',xerr=obs1['e_logSFR2']*10**obs1['MHI']*10**-obs1['logSFR2']*np.log(10)/1E9,ls='None')

plt.loglog()
plt.xlabel(r'$M_{\rm cold}$/SFR (Gyr)')
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.savefig('coldgassfrscatterplot.png')

plt.clf()
plt.plot(colddat['mcold'][wlmdacold]/f[1].data.mgas[wlmda],f[1].data.gas_metal_sfr[wlmda],'C3h',alpha=.35)
plt.plot(colddat['mcold'][wlmdcold]/f[1].data.mgas[wlmd],f[1].data.gas_metal_sfr[wlmd],'C1x',alpha=.75)
plt.plot(colddat['mcold'][wxmdacold]/f[1].data.mgas[wxmda],f[1].data.gas_metal_sfr[wxmda],'C2s',alpha=.35)
plt.plot(colddat['mcold'][wxmdcold]/f[1].data.mgas[wxmd],f[1].data.gas_metal_sfr[wxmd],'C0*',alpha=.75)

plt.loglog()
plt.xlabel(r'$M_{\rm cold}/M_{\rm gas}$')
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.savefig('coldgasfracplot.png')

plt.clf()
plt.plot(colddat['mcold'][wlmdacold]/f[1].data.mstar[wlmda],f[1].data.gas_metal_sfr[wlmda],'C3h',alpha=.35)
plt.plot(colddat['mcold'][wlmdcold]/f[1].data.mstar[wlmd],f[1].data.gas_metal_sfr[wlmd],'C1x',alpha=.75)
plt.plot(colddat['mcold'][wxmdacold]/f[1].data.mstar[wxmda],f[1].data.gas_metal_sfr[wxmda],'C2s',alpha=.35)
plt.plot(colddat['mcold'][wxmdcold]/f[1].data.mstar[wxmd],f[1].data.gas_metal_sfr[wxmd],'C0*',alpha=.75)
plt.errorbar(obs2['HI mass [Msol]']/obs2['stellar mass [Msol]'],10**(obs2['12 +(O/H)']-8.69),color='violet',ms=10,marker='^',xerr=np.sqrt((obs2['HI err']/obs2['stellar mass [Msol]'])**2 + (obs2['M* err']*obs2['HI mass [Msol]']/obs2['stellar mass [Msol]']**2)**2),yerr=obs2['Z err']*10**(obs2['12 +(O/H)']-8.69)*np.log(10),ls='None')
plt.plot(10**obs1['MHI'][wnmdhunter]/obsmst1[wnmdhunter],10**(obs1['[O/H]'][wnmdhunter]-8.69),color='purple',markeredgecolor='purple',ms=10,marker='d',ls='None',markerfacecolor='None',markeredgewidth=2,label='Hunter+12 NMDs')
plt.plot(10**obs1['MHI'][wxmdhunter]/obsmst1[wxmdhunter],10**(obs1['[O/H]'][wxmdhunter]-8.69),color='purple',ms=10,marker='d',ls='None',markeredgewidth=2,label='Hunter+12 XMDs')



#plt.plot(10**-obs1['MHI'][wnmdhunter]*10**obs1['logSFR2'][wnmdhunter]*1E9,10**(obs1['[O/H]'][wnmdhunter]-8.69),markeredgecolor='purple',ms=10,marker='d',ls='None',markerfacecolor='None',markeredgewidth=2,label='Hunter+12 NMDs')
#plt.plot(10**-obs1['MHI'][wxmdhunter]*10**obs1['logSFR2'][wxmdhunter]*1E9,10**(obs1['[O/H]'][wxmdhunter]-8.69),color='purple',ms=10,marker='d',ls='None',markeredgewidth=2,label='Hunter+12 XMDs')

plt.loglog()
plt.xlabel(r'$M_{\rm cold}/M_*$')
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.savefig('gasfracplot.png')

import bindata
b1=bindata.bindata(np.log10(f[1].data.mstar[wlmda]),np.log10(colddat['mcold'][wlmdacold]/f[1].data.mstar[wlmda]),np.arange(7.2,8.5,.15),logerr=1,median=1)
b2=bindata.bindata(np.log10(f[1].data.mstar[wlmd]),np.log10(colddat['mcold'][wlmdcold]/f[1].data.mstar[wlmd]),np.arange(7.2,8.5,.15),logerr=1,median=1)
b3=bindata.bindata(np.log10(f[1].data.mstar[wxmda]),np.log10(colddat['mcold'][wxmdacold]/f[1].data.mstar[wxmda]),np.arange(7.2,8.5,.15),logerr=1,median=1)
b4=bindata.bindata(np.log10(f[1].data.mstar[wxmd]),np.log10(colddat['mcold'][wxmdcold]/f[1].data.mstar[wxmd]),np.arange(7.2,8.5,.15),logerr=1,median=1)

plt.clf()
plt.plot(f[1].data.mstar[wlmda],colddat['mcold'][wlmdacold]/f[1].data.mstar[wlmda],'C3h',alpha=.15)
plt.plot(f[1].data.mstar[wlmd],colddat['mcold'][wlmdcold]/f[1].data.mstar[wlmd],'C1x',alpha=.5)
plt.plot(f[1].data.mstar[wxmda],colddat['mcold'][wxmdacold]/f[1].data.mstar[wxmda],'C2s',alpha=.15)
plt.plot(f[1].data.mstar[wxmd],colddat['mcold'][wxmdcold]/f[1].data.mstar[wxmd],'C0*',alpha=.5)

plt.plot(10**b1[2],10**b1[0],color='C3')
plt.plot(10**b2[2],10**b2[0],color='C1')
plt.plot(10**b3[2],10**b3[0],color='C2')
plt.plot(10**b4[2],10**b4[0],color='C0')

plt.loglog()
plt.xlim(2E7,1E9)
plt.xlabel(r'$M_*$~(M$_{\rm \odot}$)')
plt.ylabel(r'$M_{\rm cold}/M_*$')
plt.savefig('coldgasplot_mstar.png')

b1=bindata.bindata(np.log10(colddat['mcold'][wlmdacold]),np.log10(colddat['rcold'][wlmdacold]),np.arange(6,9,.15),logerr=1,median=1)
b2=bindata.bindata(np.log10(colddat['mcold'][wlmdcold]),np.log10(colddat['rcold'][wlmdcold]),np.arange(6,9,.15),logerr=1,median=1)
b3=bindata.bindata(np.log10(colddat['mcold'][wxmdacold]),np.log10(colddat['rcold'][wxmdacold]),np.arange(6,9,.15),logerr=1,median=1)
b4=bindata.bindata(np.log10(colddat['mcold'][wxmdcold]),np.log10(colddat['rcold'][wxmdcold]),np.arange(6,9,.15),logerr=1,median=1)


plt.clf()
plt.plot(colddat['mcold'][wlmdacold],colddat['rcold'][wlmdacold],'C3h',alpha=.2)
plt.plot(colddat['mcold'][wlmdcold],colddat['rcold'][wlmdcold],'C1x',alpha=.2)
plt.plot(colddat['mcold'][wxmdacold],colddat['rcold'][wxmdacold],'C2s',alpha=.35)
plt.plot(colddat['mcold'][wxmdcold],colddat['rcold'][wxmdcold],'C0*',alpha=.35)
plt.plot(10**b1[2],10**b1[0],'C3-.')
plt.plot(10**b2[2],10**b2[0],'C1--')
plt.plot(10**b3[2],10**b3[0],'C2:')
plt.plot(10**b4[2],10**b4[0],'C0-')
plt.loglog()
plt.ylim(.1,100)
plt.xlabel(r'$M_{\rm cold}~({\rm M_\odot})$')
plt.ylabel(r'$r_{\rm 1/2,cold}~({\rm kpc})$')
plt.savefig('masssizecold.png')

b1=bindata.bindata(np.log10(f[1].data.mstar[wlmda]),np.log10(colddat['rcold'][wlmdacold]),np.arange(6,9,.15),logerr=1,median=1)
b2=bindata.bindata(np.log10(f[1].data.mstar[wlmd]),np.log10(colddat['rcold'][wlmdcold]),np.arange(6,9,.15),logerr=1,median=1)
b3=bindata.bindata(np.log10(f[1].data.mstar[wxmda]),np.log10(colddat['rcold'][wxmdacold]),np.arange(6,9,.15),logerr=1,median=1)
b4=bindata.bindata(np.log10(f[1].data.mstar[wxmd]),np.log10(colddat['rcold'][wxmdcold]),np.arange(6,9,.15),logerr=1,median=1)


plt.clf()
plt.plot(f[1].data.mstar[wlmda],colddat['rcold'][wlmdacold],'C3h',alpha=.2)
plt.plot(f[1].data.mstar[wlmd],colddat['rcold'][wlmdcold],'C1x',alpha=.2)
plt.plot(f[1].data.mstar[wxmda],colddat['rcold'][wxmdacold],'C2s',alpha=.35)
plt.plot(f[1].data.mstar[wxmd],colddat['rcold'][wxmdcold],'C0*',alpha=.35)
plt.plot(10**b1[2],10**b1[0],'C3-.')
plt.plot(10**b2[2],10**b2[0],'C1--')
plt.plot(10**b3[2],10**b3[0],'C2:')
plt.plot(10**b4[2],10**b4[0],'C0-')
plt.loglog()
plt.ylim(.1,100)
plt.xlim(2E7,1e9)
plt.xlabel(r'$M_*~({\rm M_\odot})$')
plt.ylabel(r'$r_{\rm 1/2,cold}~({\rm kpc})$')
plt.savefig('masssizecold_star.png')


#allplotcold=np.unique(np.hstack([wlmdacold,wlmdcold,wxmdacold,wxmdcold]))
#allplotcold=np.unique(np.hstack([colddat['mcold'][wlmdacold],colddat['mcold'][wlmdcold],colddat['mcold'][wxmdacold],colddat['mcold'][wxmdcold]]))
#allplotsfr=np.unique(np.hstack([f[1].data.sfr[wlmda],f[1].data.sfr[wlmd],f[1].data.sfr[wxmda],f[1].data.sfr[wxmd]]))
#allplotz=np.unique(np.hstack([f[1].data.gas_metal_sfr[wlmda],f[1].data.gas_metal_sfr[wlmd],f[1].data.gas_metal_sfr[wxmda],f[1].data.gas_metal_sfr[wxmd]]))
#allplot=np.unique(np.hstack([wlmda,wlmd,wxmda,wxmd]))
plt.clf()
#plt.plot(colddat['mcold'][allplotcold]/1e9/f[1].data.sfr[allplot],f[1].data.gas_metal_sfr[allplot],'C0o',alpha=.35)
plt.plot(np.log10(colddat['mcold'][wlmdacold]/1e9/f[1].data.sfr[wlmda]),np.log10(f[1].data.gas_metal_sfr[wlmda]),'C0o',alpha=.35)
plt.plot(np.log10(colddat['mcold'][wxmdacold]/1e9/f[1].data.sfr[wxmda]),np.log10(f[1].data.gas_metal_sfr[wxmda]),'C0o',alpha=.35)

plt.xticks(fontsize=40)
plt.yticks(fontsize=40)

plt.xlabel(r'$\log$ SFR/$M_{\rm cold}$ (Gyr$^{-1}$)',fontsize=50)
plt.ylabel(r'$\log Z_{\rm SFR}/Z_\odot$',fontsize=50)
plt.savefig('coldgasplot_all.png')

plt.clf()

#plt.plot(colddat['mcold'][allplotcold]/1e9/f[1].data.sfr[allplot],f[1].data.ga\s_metal_sfr[allplot],'C0o',alpha=.35)                                           
plt.plot(np.log10(colddat['mcold'][wlmdacold]),np.log10(f[1].data.gas_metal_sfr[wlmda]),'C0o',alpha=.35)
plt.plot(np.log10(colddat['mcold'][wxmdacold]),np.log10(f[1].data.gas_metal_sfr[wxmda]),'C0o',alpha=.35)

plt.xticks(fontsize=40)
plt.yticks(fontsize=40)

plt.xlabel(r'$\log$ $M_{\rm cold}$ (${\rm M_\odot}$)',fontsize=50)
plt.ylabel(r'$\log Z_{\rm SFR}/Z_\odot$',fontsize=50)
plt.savefig('coldgasplot2_all.png')

plt.clf()

#plt.plot(colddat['mcold'][allplotcold]/1e9/f[1].data.sfr[allplot],f[1].data.ga\s_metal_sfr[allplot],'C0o',alpha=.35)                                           
#plt.plot(np.log10(f[1].data.sfr[wlmda]),np.log10(f[1].data.gas_metal_sfr[wlmda]),'C0o',alpha=.35)
#plt.plot(np.log10(f[1].data.sfr[wxmda]),np.log10(f[1].data.gas_metal_sfr[wxmda]),'C0o',alpha=.35)

plt.scatter(np.log10(f[1].data.sfr[wlmda]),np.log10(f[1].data.gas_metal_sfr[wlmda]),c=np.log10(colddat['mcold'][wlmdacold]),alpha=.35)
plt.scatter(np.log10(f[1].data.sfr[wxmda]),np.log10(f[1].data.gas_metal_sfr[wxmda]),c=np.log10(colddat['mcold'][wxmdacold]),alpha=.35)

plt.xticks(fontsize=40)
plt.yticks(fontsize=40)

plt.xlabel(r'$\log$ ${\rm SFR}$ (${\rm M_\odot/yr}$)',fontsize=50)
plt.ylabel(r'$\log Z_{\rm SFR}/Z_\odot$',fontsize=50)
plt.savefig('coldgasplot3_all.png')

