import getgasinfo2_40 as getgasinfo2
import imp
from astropy import cosmology
import numpy as np
import matplotlib.pyplot as plt
import getdeterminationcoef
from astropy.io import fits, ascii
plt.style.use('../python/stylesheet.txt')

astrocosmo=cosmology.Planck15

hconst=0.6774

zs=np.array([20.05,14.99,11.98,10.976,10.0,9.389,9.002,8.449,8.012,7.595,7.236,7.005,6.491,6.011,5.847,5.530,5.228,4.996,4.665,4.428,4.177,4.008,3.709,3.491,3.283,3.008,2.896,2.733,2.578,2.444,2.316,2.208,2.103,2.002,1.904,1.823,1.743,1.667,1.604,1.531,1.496,1.414,1.358,1.302,1.248,1.206,1.155,1.114,1.074,1.036,0.997,0.951,0.923,0.887,0.851,0.817,0.791,0.757,0.733,0.700,0.676,0.644,0.621,0.598,0.576,0.546,0.525,0.503,0.482,0.461,0.440,0.420,0.40,0.380,0.361,0.348,0.329,0.310,0.298,0.273,0.261,0.244,0.226,0.214,0.1973,0.1804,0.1693,0.1527,0.1419,0.1258,0.1099,0.0994,0.0839,0.0737,0.0585,0.0485,0.0337,0.0240,0.0095,0])

time=astrocosmo.age(zs).value
time40=time[-40:]
res1=getgasinfo2.getgasinfo(1)
res3=getgasinfo2.getgasinfo(3)

res2=getgasinfo2.getgasinfo(2)
res4=getgasinfo2.getgasinfo(4)
f=fits.open('tng50/alltng50gal_withr_withspin.fits')
obj=ascii.read('sampleselection_tng50_sfr5_notide_nomass.txt')

wok1=np.array([i for i in range(len(res1[3])) if len(res1[3][i])>1])
wok2=np.array([i for i in range(len(res2[3])) if len(res2[3][i])>1])
wok3=np.array([i for i in range(len(res3[3])) if len(res3[3][i])>1])
wok4=np.array([i for i in range(len(res4[3])) if len(res4[3][i])>1])

wok1=np.array([i for i in wok1 if len(res1[3][i].shape)>1])
wok2=np.array([i for i in wok2 if len(res2[3][i].shape)>1])
wok3=np.array([i for i in wok3 if len(res3[3][i].shape)>1])
wok4=np.array([i for i in wok4 if len(res4[3][i].shape)>1])

wxmd=np.where(obj['type']==1)[0][wok1]
wlmd=np.where(obj['type']==2)[0][wok2]
wxmda=np.where(obj['type']==3)[0][wok3]
wlmda=np.where(obj['type']==4)[0][wok4]


nts=len(res1[0][0][0])
nobj1=len(wok1)
nobj2=len(wok2)
nobj3=len(wok3)
nobj4=len(wok4)

avgmet1=np.array([np.nanmedian(res1[1][i]*res1[2][i],axis=0)/np.nanmedian(res1[2][i],axis=0) for i in wok1])
avgmet2=np.array([np.nanmedian(res2[1][i]*res2[2][i],axis=0)/np.nanmedian(res2[2][i],axis=0) for i in wok2])
avgmet3=np.array([np.nanmedian(res3[1][i]*res3[2][i],axis=0)/np.nanmedian(res3[2][i],axis=0) for i in wok3])
avgmet4=np.array([np.nanmedian(res4[1][i]*res4[2][i],axis=0)/np.nanmedian(res4[2][i],axis=0) for i in wok4])

avgdist1=np.array([np.nanmedian(res1[0][i]*res1[2][i],axis=0)/np.nanmedian(res1[2][i],axis=0) for i in wok1])
avgdist2=np.array([np.nanmedian(res2[0][i]*res2[2][i],axis=0)/np.nanmedian(res2[2][i],axis=0) for i in wok2])
avgdist3=np.array([np.nanmedian(res3[0][i]*res3[2][i],axis=0)/np.nanmedian(res3[2][i],axis=0) for i in wok3])
avgdist4=np.array([np.nanmedian(res4[0][i]*res4[2][i],axis=0)/np.nanmedian(res4[2][i],axis=0) for i in wok4])

avgsfr1=np.array([np.nanmedian(res1[3][i]*res1[2][i],axis=0)/np.nanmedian(res1[2][i],axis=0) for i in wok1])
avgsfr2=np.array([np.nanmedian(res2[3][i]*res2[2][i],axis=0)/np.nanmedian(res2[2][i],axis=0) for i in wok2])
avgsfr3=np.array([np.nanmedian(res3[3][i]*res3[2][i],axis=0)/np.nanmedian(res3[2][i],axis=0) for i in wok3])
avgsfr4=np.array([np.nanmedian(res4[3][i]*res4[2][i],axis=0)/np.nanmedian(res4[2][i],axis=0) for i in wok4])

timeingal1=[]
for i in range(len(wok1)):
    timeingal1i=[]
    for j in range(len(res1[0][wok1[i]])):
        wingal=np.where(res1[0][wok1[i]][j]<f[1].data.half_mass_rad_star[obj['allgalindex'][wxmd[i]]])[0]
        if len(wingal)>0:
            print(wingal,time40[wingal],time40[wingal[0]])
            timeingal1i.append(time[-1]-time40[wingal[0]])
        else:
            timeingal1i.append(0)
    timeingal1.append(np.median(timeingal1i))

timeingal2=[]
for i in range(len(wok2)):
    timeingal2i=[]
    for j in range(len(res2[0][wok2[i]])):
        wingal=np.where(res2[0][wok2[i]][j]<f[1].data.half_mass_rad_star[obj['allgalindex'][wlmd[i]]])[0]
        if len(wingal)>0:
            timeingal2i.append(time[-1]-time40[wingal[0]])
        else:
            timeingal2i.append(0)
    timeingal2.append(np.median(timeingal2i))
    
timeingal3=[]
for i in range(len(wok3)):
    timeingal3i=[]
    for j in range(len(res3[0][wok3[i]])):
        wingal=np.where(res3[0][wok3[i]][j]<f[1].data.half_mass_rad_star[obj['allgalindex'][wxmda[i]]])[0]
        if len(wingal)>0:
            timeingal3i.append(time[-1]-time40[wingal[0]])
        else:
            timeingal3i.append(0)
    timeingal3.append(np.median(timeingal3i))

timeingal4=[]
for i in range(len(wok4)):
    timeingal4i=[]
    for j in range(len(res4[0][wok4[i]])):
        wingal=np.where(res4[0][wok4[i]][j]<f[1].data.half_mass_rad_star[obj['allgalindex'][wlmda[i]]])[0]
        if len(wingal)>0:
            timeingal4i.append(time[-1]-time40[wingal[0]])
        else:
            timeingal4i.append(0)
    timeingal4.append(np.median(timeingal4i))
    
zsp=zs[-nts:]

plt.clf()
#for i in range(1):
#    plt.plot(astrocosmo.age(zsp),avgmet4[i]/.0127,alpha=.2,color='C3',ls='-.')

#for i in range(1):
#    plt.plot(astrocosmo.age(zsp),avgmet3[i]/.0127,alpha=.2,color='C2',ls=':')
    
#for i in range(1):
#    plt.plot(astrocosmo.age(zsp),avgmet2[i]/.0127,alpha=.2,color='C1',ls='--')

#for i in range(1):
#    plt.plot(astrocosmo.age(zsp),avgmet1[i]/.0127,alpha=.2,color='C0',ls='-')

plt.plot(astrocosmo.age(zsp),np.nanmedian(avgmet4,axis=0)/.0127,color='C3',ls='-.',label='LMD Analogs')
plt.plot(astrocosmo.age(zsp),np.nanmedian(avgmet3,axis=0)/.0127,color='C2',ls='--',label='XMD Analogs')
plt.plot(astrocosmo.age(zsp),np.nanmedian(avgmet2,axis=0)/.0127,color='C1',ls=':',label='LMDs')
plt.plot(astrocosmo.age(zsp),np.nanmedian(avgmet1,axis=0)/.0127,color='C0',ls='-',label='XMDs')
plt.xlabel(r'$t$ (Gyr)')
plt.ylabel(r'$Z_p/Z_\odot$')
plt.legend()
plt.yscale('log')
plt.title('TNG50',y=0.92)
plt.savefig('gasevolution_50.png')

plt.clf()

for i in range(1):
    plt.plot(avgdist4[i]/hconst,avgmet4[i]/.0127,alpha=.2,color='C3',ls='-.')

for i in range(1):
    plt.plot(avgdist3[i]/hconst,avgmet3[i]/.0127,alpha=.2,color='C2',ls=':')
    
for i in range(1):
    plt.plot(avgdist2[i]/hconst,avgmet2[i]/.0127,alpha=.2,color='C1',ls='--')

for i in range(1):
    plt.plot(avgdist1[i]/hconst,avgmet1[i]/.0127,alpha=.2,color='C0',ls='-')

plt.plot(np.nanmedian(avgdist4,axis=0)/hconst,np.nanmedian(avgmet4,axis=0)/.0127,color='C3',ls='-.',label='LMD Analogs')
plt.plot(np.nanmedian(avgdist3,axis=0)/hconst,np.nanmedian(avgmet3,axis=0)/.0127,color='C2',ls='--',label='XMD Analogs')
plt.plot(np.nanmedian(avgdist2,axis=0)/hconst,np.nanmedian(avgmet2,axis=0)/.0127,color='C1',ls=':',label='LMDs')
plt.plot(np.nanmedian(avgdist1,axis=0)/hconst,np.nanmedian(avgmet1,axis=0)/.0127,color='C0',ls='-',label='XMDs')
plt.xlabel(r'D (kpc)')
plt.ylabel(r'$Z_p/Z_\odot$')
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.savefig('gasevolution_dist1_50.png')

plt.clf()
plt.plot(np.nanmedian(avgdist4,axis=0)/hconst,np.nanmedian(avgmet4,axis=0)/.0127,ls='-.',marker='h',color='C3',label='LMD Analogs')
plt.plot(np.nanmedian(avgdist3,axis=0)/hconst,np.nanmedian(avgmet3,axis=0)/.0127,ls='--',marker='s',color='C2',label='XMD Analogs')
plt.plot(np.nanmedian(avgdist2,axis=0)/hconst,np.nanmedian(avgmet2,axis=0)/.0127,ls=':',marker='x',color='C1',label='LMDs')
plt.plot(np.nanmedian(avgdist1,axis=0)/hconst,np.nanmedian(avgmet1,axis=0)/.0127,ls='-',marker='*',color='C0',label='XMDs')
plt.xlabel(r'D (kpc)')
plt.ylabel(r'$Z_p/Z_\odot$')
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.savefig('gasevolution_dist_50.png')

plt.clf()
for i in range(1):
    plt.plot(astrocosmo.age(zsp),avgdist4[i]/hconst,alpha=.2,color='C3',ls='-.')

for i in range(1):
    plt.plot(astrocosmo.age(zsp),avgdist3[i]/hconst,alpha=.2,color='C2',ls=':')
    
for i in range(1):
    plt.plot(astrocosmo.age(zsp),avgdist2[i]/hconst,alpha=.2,color='C1',ls='--')

for i in range(1):
    plt.plot(astrocosmo.age(zsp),avgdist1[i]/hconst,alpha=.2,color='C0',ls='-')

plt.plot(astrocosmo.age(zsp),np.nanmedian(avgdist4,axis=0)/hconst,color='C3',ls='-.',label='LMD Analogs')
plt.plot(astrocosmo.age(zsp),np.nanmedian(avgdist3,axis=0)/hconst,color='C2',ls=':',label='XMD Analogs')
plt.plot(astrocosmo.age(zsp),np.nanmedian(avgdist2,axis=0)/hconst,color='C1',ls='--',label='LMDs')
plt.plot(astrocosmo.age(zsp),np.nanmedian(avgdist1,axis=0)/hconst,color='C0',ls='-',label='XMDs')
plt.xlabel(r'$t$ (Gyr)')
plt.ylabel(r'$Z_p/Z_\odot$')
plt.legend()
plt.yscale('log')
plt.savefig('gasdistevolution_50.png')



plt.clf()
for i in range(1):
    plt.plot(avgdist1[i]/hconst,avgsfr1[i],alpha=.2,color='C0',ls='-')

for i in range(1):
    plt.plot(avgdist3[i]/hconst,avgsfr3[i],alpha=.2,color='C2',ls=':')

plt.plot(np.nanmedian(avgdist1,axis=0)/hconst,np.nanmedian(avgsfr1,axis=0),color='C0',ls='-',label='LMDs')
plt.plot(np.nanmedian(avgdist3,axis=0)/hconst,np.nanmedian(avgsfr3,axis=0),color='C2',ls=':',label='LMD Analogs')
plt.xlabel(r'D (kpc)')
plt.xlim(0,10)
plt.ylim(1E-5,1E-2)
plt.ylabel(r'Regional SFR (Msun/yr)')
plt.legend()
plt.yscale('log')
plt.savefig('gasevolution_distsfr_50.png')

plt.clf()
#for i in range(1):
#    plt.plot(astrocosmo.age(zsp),avgsfr1[i],alpha=.2,color='C0',ls='-')

#for i in range(1):
#    plt.plot(astrocosmo.age(zsp),avgsfr3[i],alpha=.2,color='C2',ls=':')

plt.plot(astrocosmo.age(zsp),np.nanmedian(avgsfr4,axis=0),color='C3',ls='-.',label='LMD Analogs')
plt.plot(astrocosmo.age(zsp),np.nanmedian(avgsfr3,axis=0),color='C2',ls=':',label='XMD Analogs')
plt.plot(astrocosmo.age(zsp),np.nanmedian(avgsfr2,axis=0),color='C1',ls='--',label='LMDs')
plt.plot(astrocosmo.age(zsp),np.nanmedian(avgsfr1,axis=0),color='C0',ls='-',label='XMDs')
plt.xlabel(r'Age (Gyr)')
plt.ylabel(r'Regional SFR (M$_\odot$~yr$^{-1}$)')
plt.legend()
plt.ylim(1E-5,1E-2)
plt.yscale('log')
plt.savefig('gasevolution_timesfr_50.png')

plt.clf()
#for i in range(1):
#    plt.plot(astrocosmo.age(zsp),avgsfr1[i],alpha=.2,color='C0',ls='-')

#for i in range(1):
#    plt.plot(astrocosmo.age(zsp),avgsfr3[i],alpha=.2,color='C2',ls=':')

plt.plot(astrocosmo.age(zsp),np.cumsum(np.nanmedian(avgsfr4,axis=0)),color='C3',ls='-.',label='LMD Analogs')
plt.plot(astrocosmo.age(zsp),np.cumsum(np.nanmedian(avgsfr3,axis=0)),color='C2',ls=':',label='XMD Analogs')
plt.plot(astrocosmo.age(zsp),np.cumsum(np.nanmedian(avgsfr2,axis=0)),color='C1',ls='--',label='LMDs')
plt.plot(astrocosmo.age(zsp),np.cumsum(np.nanmedian(avgsfr1,axis=0)),color='C0',ls='-',label='XMDs')
plt.xlabel(r'Age (Gyr)')
plt.ylabel(r'Regional SFR (M$_\odot$~yr$^{-1}$)')
plt.legend()
plt.ylim(1E-5,2E-2)
plt.yscale('log')
plt.savefig('gasevolution_timesfr_50_cum.png')


dzs=zs[-nts-1:]
deltat=(astrocosmo.age(dzs[1:])-astrocosmo.age(dzs[0:-1])).value
print(len(deltat),len(wlmda))
plt.clf()
plt.title('TNG50',y=.92)
plt.plot(np.nansum(avgsfr4*np.tile(deltat,len(wlmda)).reshape(len(wlmda),40),axis=1)*1E9,f[1].data.gas_metal_sfr[obj['allgalindex'][wlmda]],'C3h')
plt.plot(np.nansum(avgsfr3*np.tile(deltat,len(wxmda)).reshape(len(wxmda),40),axis=1)*1E9,f[1].data.gas_metal_sfr[obj['allgalindex'][wxmda]],'C2s')
plt.plot(np.nansum(avgsfr2*np.tile(deltat,len(wlmd)).reshape(len(wlmd),40),axis=1)*1E9,f[1].data.gas_metal_sfr[obj['allgalindex'][wlmd]],'C1x')
plt.plot(np.nansum(avgsfr1*np.tile(deltat,len(wxmd)).reshape(len(wxmd),40),axis=1)*1E9,f[1].data.gas_metal_sfr[obj['allgalindex'][wxmd]],'C0*')
plt.loglog()
plt.xlabel(r'Cumulative SFR (M$_\odot$)')
plt.ylabel(r'$Z_{\rm SFR}/{\rm Z_\odot}$')
plt.savefig('metvscumsfr1.png')

plt.clf()
plt.title('TNG50',y=.92)
plt.plot(np.nansum(avgsfr4,axis=1),avgmet4[:,-1],'C3h')
plt.plot(np.nansum(avgsfr3,axis=1),avgmet3[:,-1],'C2s')
plt.plot(np.nansum(avgsfr2,axis=1),avgmet2[:,-1],'C1x')
plt.plot(np.nansum(avgsfr1,axis=1),avgmet1[:,-1],'C0*')
plt.loglog()
plt.xlabel(r'Cumulative SFR (M$_\odot$)')
plt.ylabel(r'$Z_{\rm SFR}/{\rm Z_\odot}$')
plt.savefig('metvscumsfr.png')


cumsfrall=np.hstack([np.nansum(avgsfr4*np.tile(deltat,len(wlmda)).reshape(len(wlmda),40),axis=1)*1E9,np.nansum(avgsfr3*np.tile(deltat,len(wxmda)).reshape(len(wxmda),40),axis=1)*1E9,np.nansum(avgsfr2*np.tile(deltat,len(wlmd)).reshape(len(wlmd),40),axis=1)*1E9,np.nansum(avgsfr1*np.tile(deltat,len(wxmd)).reshape(len(wxmd),40),axis=1)*1E9])
metall=np.hstack([f[1].data.gas_metal_sfr[obj['allgalindex'][wlmda]],f[1].data.gas_metal_sfr[obj['allgalindex'][wxmda]],f[1].data.gas_metal_sfr[obj['allgalindex'][wlmd]],f[1].data.gas_metal_sfr[obj['allgalindex'][wxmd]]])
#print(getdeterminationcoef.getdeterminationcoef(np.log10(np.nansum(avgsfr4,axis=1)),np.log10(avgmet4[:,-1]),deg=1))
print(getdeterminationcoef.getdeterminationcoef(np.log10(cumsfrall),np.log10(metall),deg=1))

wa=np.where(np.log10(cumsfrall)>5)[0]
print(wa)
print(getdeterminationcoef.getdeterminationcoef(np.log10(cumsfrall)[wa],np.log10(metall)[wa],deg=1))

plt.clf()
plt.scatter(np.nansum(avgsfr4,axis=1),avgmet4[:,-1],c=np.log10(avgmet4[:,-4]))
plt.scatter(np.nansum(avgsfr3,axis=1),avgmet3[:,-1],c=np.log10(avgmet3[:,-4]))
plt.scatter(np.nansum(avgsfr2,axis=1),avgmet2[:,-1],c=np.log10(avgmet2[:,-4]))
plt.scatter(np.nansum(avgsfr1,axis=1),avgmet1[:,-1],c=np.log10(avgmet1[:,-4]))
plt.loglog()
plt.savefig('metvscumsfr_0.png')
