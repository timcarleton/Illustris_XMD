import numpy as np
from astropy.io import ascii, fits
from scipy import stats
from astropy import cosmology
astrocosmo = cosmology.FlatLambdaCDM(H0=67.74,Om0=0.307)
from scipy.optimize import minimize
import getdeterminationcoef

allznowx=[]
allznowa=[]
allznowl=[]
allznowla=[]

allzfa=[]
allz0a=[]
allsfr2a=[]
allsfea=[]
allsfrwa=[]
allmz0a=[]
allznowa=[]
allmg0a=[]
allrha=[]
alletaa=[]
allsfrnowa=[]
alldza=[]
alldlnmu=[]
allmu=[]
allzsfratiox=[]
allzsfratioa=[]
allsfe2x=[]
allsfe2a=[]
allsfe2l=[]
allsfe2la=[]
allmstl=[]
allmstla=[]

allsfrl=[]
allzfl=[]
allzfla=[]
allsfr2l=[]
allsfr2la=[]
allsfrnowla=[]
allsfrnowl=[]
allzsfratiol=[]
allzsfratiola=[]
allrhla=[]
allrhl=[]
allzfla=[]
allsfr2la=[]
allmg0la=[]
allmg0l=[]
allznowl=[]
allznowla=[]
alletal=[]
alletala=[]
allsfel=[]
allsfela=[]
allzsfratiola=[]
allzsfratiol=[]
allz0la=[]
allz0l=[]
allmula=[]
allmua=[]
allmul=[]
allmux=[]

allz0x=[]
alletax=[]
allsfr2x=[]
allsfr2a=[]
allsfex=[]
allzfx=[]
allmg0x=[]
allrhx=[]
allsfrnowx=[]
allmstx=[]
allmsta=[]

obj=ascii.read('sampleselection_tng50_sfr5_notide_nomass.txt')
f=fits.open('tng50/alltng50gal_withr_withspin.fits')
#obj=ascii.read('sampleselection_tng50_sfr5.txt')

wxmd=np.where(obj['type']==1)[0]
wlmd=np.where(obj['type']==2)[0]
wxmda=np.where(obj['type']==3)[0]
wlmda=np.where(obj['type']==4)[0]


allmgi=np.hstack([f[1].data.mgas[obj['allgalindex'][wlmda]],f[1].data.mgas[obj['allgalindex'][wxmda]],f[1].data.mgas[obj['allgalindex'][wlmd]],f[1].data.mgas[obj['allgalindex'][wxmd]]])
allmsti=np.hstack([f[1].data.mstar[obj['allgalindex'][wlmda]],f[1].data.mstar[obj['allgalindex'][wxmda]],f[1].data.mstar[obj['allgalindex'][wlmd]],f[1].data.mstar[obj['allgalindex'][wxmd]]])

print(np.median(allmsti/allmgi),np.median(allmgi))
deltat=1.32E8


for i in wlmda:
    try:
        x=ascii.read('flowinfo/flowinfo_50_4_'+str(int(obj['haloid'][i]))+'.txt')
    except FileNotFoundError:
        continue
    allz0la.append(x['mzinflow_5rh'][-1]/x['minflow_5rh'][-1]/.0127)
    wi=np.where(x['snapshot']==89)[0]
    if len(wi)==0:
        allsfr2la.append(np.nan)
    else:
        wi=wi[0]
        allsfr2la.append(np.median(x['sfr'][wi-2:wi+2]))
    alletala.append(x['moutflow_5rh'][-1]/deltat/x['sfr'][-1])
    allsfela.append(x['sfr'][-1]/x['m_gas_rh'][-1])
    #allsfex.append(x['moutflow_2rh'][-1])
    allmg0la.append(x['mgas_tot'][-1])
    #allzfx.append(x['sfr_weight_metallicity'][-1])
    allzfla.append(f[1].data.gas_metal_sfr[int(obj['allgalindex'][i])])
    allznowla.append(f[1].data.gas_metal_sfr[int(obj['allgalindex'][i])])
    allrhla.append(x['rhstar'][-1])
    allsfrnowla.append(x['sfr'][-1])
    allmula.append(x['m_gas_rh'][-1]/x['mstar'][-1])
    allmstla.append(x['mstar'][-1])
    allzsfratiola.append(x['sfr_weight_metallicity'][-1]/x['mass_weight_metallicity'][-1])

for i in wlmd:
    try:
        x=ascii.read('flowinfo/flowinfo_50_2_'+str(int(obj['haloid'][i]))+'.txt')
    except FileNotFoundError:
        continue
    allz0l.append(x['mzinflow_5rh'][-1]/x['minflow_5rh'][-1]/.0127)
    wi=np.where(x['snapshot']==89)[0]
    if len(wi)==0:
        allsfr2l.append(np.nan)
    else:
        wi=wi[0]
        allsfr2l.append(np.median(x['sfr'][wi-2:wi+2]))
    alletal.append(x['moutflow_5rh'][-1]/deltat/x['sfr'][-1])
    allsfel.append(x['sfr'][-1]/x['m_gas_rh'][-1])
    #allsfex.append(x['moutflow_2rh'][-1])
    allmg0l.append(x['mgas_tot'][-1])
    #allzfx.append(x['sfr_weight_metallicity'][-1])
    allzfl.append(f[1].data.gas_metal_sfr[int(obj['allgalindex'][i])])
    allznowl.append(f[1].data.gas_metal_sfr[int(obj['allgalindex'][i])])
    allrhl.append(x['rhstar'][-1])
    allsfrnowl.append(x['sfr'][-1])
    allmul.append(x['m_gas_rh'][-1]/x['mstar'][-1])
    allmstl.append(x['mstar'][-1])
    allzsfratiol.append(x['sfr_weight_metallicity'][-1]/x['mass_weight_metallicity'][-1])


for i in wxmda:
    x=ascii.read('flowinfo/flowinfo_50_3_'+str(int(obj['haloid'][i]))+'.txt')
    allz0a.append(x['mzinflow_5rh'][-1]/x['minflow_5rh'][-1]/.0127)
    wi=np.where(x['snapshot']==89)[0]
    if len(wi)==0:
        allsfr2a.append(np.nan)
    else:
        wi=wi[0]
        allsfr2a.append(np.median(x['sfr'][wi-2:wi+2]))
    allsfrnowa.append(x['sfr'][-1])
    allsfea.append(x['sfr'][-1]/x['m_gas_rh'][-1])
    #allsfea.append(x['moutflow_2rh'][-1])
    alletaa.append(x['moutflow_5rh'][-1]/deltat/x['sfr'][-1])
    #allznowa.append(x['sfr_weight_metallicity'][-1])
    allzfa.append(f[1].data.gas_metal_sfr[int(obj['allgalindex'][i])])
    allznowa.append(f[1].data.gas_metal_sfr[int(obj['allgalindex'][i])])
    allzsfratioa.append(x['sfr_weight_metallicity'][-1]/x['mass_weight_metallicity'][-1])
    allsfrwa.append(np.mean(x['sfr']/(13.8-astrocosmo.age(x['redshift']).value))/np.mean(1/(13.8-astrocosmo.age(x['redshift']).value)))
    allrha.append(x['rhstar'][-1])
    allmz0a.append(x['mzinflow_vir'][-1])
    allmg0a.append(x['mgas_tot'][-1])
    alldza.append((x['sfr_weight_metallicity'][-1]-x['sfr_weight_metallicity'][-2])/.13)
    allmua.append(x['m_gas_rh'][-1]/x['mstar'][-1])
    allmsta.append(x['mstar'][-1])
    alldlnmu.append((np.log(x['m_gas_rh'][-1]/x['mstar'][-1])-np.log(x['m_gas_rh'][-2]/x['mstar'][-2]))/.13)


for i in wxmd:
    x=ascii.read('flowinfo/flowinfo_50_1_'+str(int(obj['haloid'][i]))+'.txt')
    allz0x.append(x['mzinflow_5rh'][-1]/x['minflow_5rh'][-1]/.0127)
    wi=np.where(x['snapshot']==89)[0]
    if len(wi)==0:
        allsfr2x.append(np.nan)
    else:
        wi=wi[0]
        allsfr2x.append(np.median(x['sfr'][wi-2:wi+2]))
    alletax.append(x['moutflow_5rh'][-1]/deltat/x['sfr'][-1])
    allsfex.append(x['sfr'][-1]/x['m_gas_rh'][-1])
    #allsfex.append(x['moutflow_2rh'][-1])
    allmg0x.append(x['mgas_tot'][-1])
    #allzfx.append(x['sfr_weight_metallicity'][-1])
    allzfx.append(f[1].data.gas_metal_sfr[int(obj['allgalindex'][i])])
    allznowx.append(f[1].data.gas_metal_sfr[int(obj['allgalindex'][i])])
    allrhx.append(x['rhstar'][-1])
    allsfrnowx.append(x['sfr'][-1])
    allmux.append(x['m_gas_rh'][-1]/x['mstar'][-1])
    allmstx.append(x['mstar'][-1])
    allzsfratiox.append(x['sfr_weight_metallicity'][-1]/x['mass_weight_metallicity'][-1])
    
allmu=np.hstack([allmula,allmua,allmul,allmux])
allmg0=np.hstack([allmg0la,allmg0a,allmg0l,allmg0x])
allsfr2=np.hstack([allsfr2la,allsfr2a,allsfr2l,allsfr2x])
allzsfratio=np.hstack([allzsfratiola,allzsfratiol,allzsfratioa,allzsfratiox])
allsfe=np.hstack([allsfela,allsfea,allsfel,allsfex])*1E9
allzf=np.hstack([allznowla,allznowa,allznowl,allznowx])
allz0=np.hstack([allz0la,allz0a,allz0l,allz0x])
allrh=np.hstack([allrhla,allrha,allrhl,allrhx])
alleta=np.hstack([alletala,alletaa,alletal,alletax])
allsfrnow=np.hstack([allsfrnowla,allsfrnowa,allsfrnowl,allsfrnowx])
w=np.where(alleta>0)[0]

print(np.shape(allsfr2),np.shape(allmg0),np.shape(allsfe),np.shape(allzsfratio))
def minfunc(xf):
    tst=10**xf[0]*allsfr2/allmg0+10**xf[1]*allsfe+xf[2]
    #print(xf,tst)
    if np.any(tst<0):
        return np.inf
    else:
        return np.log10(tst)

def minfunc2(xf):
    tst=xf[0]*allz0+10**xf[1]*allsfe+xf[2]
    #print(xf,tst)
    if len(np.where(tst<0)[0])>5:
        print(xf)
        return np.inf
    else:
        return np.log10(tst)

def minfunc3(xf):
    tst=(10**xf[0]*allsfr2+10**xf[1]*allsfrnow+xf[2])/(allmg0/1E9)
    #print(xf,tst)
    if len(np.where(tst<0)[0])>5:
        print('t',xf)
        return np.inf
    else:
        return np.log10(tst)

def minfunc4(xf):
    tst=xf[0]*allz0[w]+10**xf[1]/alleta[w]+xf[2]
    #print(xf,tst)
    if len(np.where(tst<0)[0])>5:
        print('four',xf)
        return np.inf
    else:
        return np.log10(tst)

def minfunc5(xf):
    tst=xf[0]*allz0[w]+10**xf[1]/(alleta[w]+10**xf[2]/allsfe[w])+xf[3]
    #print(xf,tst)
    if len(np.where(tst<0)[0])>5:
        print('four',xf)
        return np.inf
    else:
        return np.log10(tst)

    
print(allsfr2/allmg0,allsfe)
#print(minfunc([1,.01,0]))
#print(minfunc([1,.1,0]))
ft=minimize(lambda xf:np.nansum((np.log10(allzf)-minfunc(xf))**2),[-3,-2,0])
ft2=minimize(lambda xf:np.nansum((np.log10(allzf)-minfunc2(xf))**2),[1,-2,0],bounds=[[.5,2],[-5,0],[-.00000001,.00000001]])
ft3=minimize(lambda xf:np.nansum((np.log10(allzf)-minfunc3(xf))**2),[1,-2,0],bounds=[[-3,0],[-5,0],[0,.0001]])
ft4=minimize(lambda xf:np.nansum((np.log10(allzf[w])-minfunc4(xf))**2),[1,0,0],bounds=[[.5,2],[-11,2],[-.000001,.000001]])
ft5=minimize(lambda xf:np.nansum((np.log10(allzf[w])-minfunc5(xf))**2),[1,0,-2,0],bounds=[[.5,2],[-11,2],[-5,0],[-.000001,.000001]])
print(ft,ft2,ft3,ft4,ft5)
print('3',ft3)
print(np.nanstd(np.log10(allzf)-minfunc(ft['x'])))
print(np.nanstd(np.log10(allzf)-minfunc2(ft2['x'])))
print('t',np.nanstd(np.log10(allzf)-minfunc3(ft3['x'])))
print('t4',np.nanstd(np.log10(allzf[w])-minfunc4(ft4['x'])))
print('t5',np.nanstd(np.log10(allzf[w])-minfunc5(ft5['x'])))
print(np.nanstd(np.log10(allzf)))
#print(np.log10(allzf),np.log10(allz0+10**-2*allsfe),np.log10(allzf)-np.log10(allz0+10**-2*allsfe))
print(minfunc2([1,-2,0]))
print(minfunc2([1.1,-2,0.001]))
#print(minfunc2([1,-2.1,-.001]))
print(x.keys())
import matplotlib.pyplot as plt
plt.style.use('../python/stylesheet.txt')

plt.clf()
plt.scatter(allsfe,allsfr2,c=np.log10(allzf))
plt.loglog()
plt.xlim(1E-3,1)
plt.ylim(1E-4,.1)
plt.savefig('sfesfr2_50.png')


#FSPS Z of-1, tau=1000 results in z=0 M* of 5.96578072e-01, so R=.6
plt.clf()
plt.plot(allz0la+.008/(1+np.array(alletala)/.4+np.array(allmula)),allzfla,'C3h')
plt.plot(allz0a+.008/(1+np.array(alletaa)/.4+np.array(allmua)),allzfa,'C2s')
plt.plot(allz0l+.008/(1+np.array(alletal)/.4+np.array(allmul)),allzfl,'C1x')
plt.plot(allz0x+.008/(1+np.array(alletax)/.4+np.array(allmux)),allzfx,'C0*')
plt.loglog()
plt.xlabel(r'$Z_{\rm model}$')
plt.ylabel(r'$Z_{\rm SFR}$')
plt.title('TNG50',y=0.92)
plt.savefig('zvszmodel.png')
allmod0=np.hstack([allz0la+.008/(1+np.array(alletala)/.4+np.array(allmula)),allz0a+.008/(1+np.array(alletaa)/.4+np.array(allmua)),allz0l+.008/(1+np.array(alletal)/.4+np.array(allmul)),allz0x+.008/(1+np.array(alletax)/.4+np.array(allmux))])
allmetfit=np.hstack([allzfla,allzfa,allzfl,allzfx])
print('coefmodel',getdeterminationcoef.getdeterminationcoef(np.log10(allmod0),np.log10(allmetfit),deg=1))

plt.clf()
plt.scatter(allsfr2,allz0,c=-np.log10(np.array(allzsfratio)))
plt.xlim(1E-4,.1)
plt.xlabel(r'SFR$_{z=0.13}$ (M$_\odot$~yr$^{-1}$)')
plt.ylabel(r'$Z_{\rm infall}$')
plt.colorbar(label=r'$\log Z_{\rm SFR}/Z_{\rm mass}$')
plt.loglog()
plt.savefig('sfr2z0ratio_50.png')


plt.clf()
plt.plot(allsfr2,allzsfratio,'o')
plt.loglog()
plt.savefig('sfr2ratio_50.png')

plt.clf()
plt.scatter(allmz0a,allznowa,c=np.log10(np.array(allsfea)))
plt.loglog()
plt.savefig('mz0zg_50.png')

plt.clf()
plt.scatter(np.array(allmz0a)*np.array(allmg0a),allznowa,c=np.log10(np.array(allsfea)))
plt.loglog()
plt.savefig('mz0mgz_50.png')

plt.clf()
plt.scatter(allz0la,allznowla,c=np.log10(np.array(allzsfratiola)),vmin=-.2,vmax=1.3,alpha=.5)
plt.scatter(allz0a,allznowa,c=np.log10(np.array(allzsfratioa)),vmin=-.2,vmax=1.3,alpha=.5)
plt.scatter(allz0l,allznowl,c=np.log10(np.array(allzsfratiol)),vmin=-.2,vmax=1.3,alpha=.5)
plt.scatter(allz0x,allzfx,c=np.log10(np.array(allzsfratiox)),vmin=-.2,vmax=1.3,alpha=.5)
plt.loglog()
plt.plot([3E-3,2],[3E-3,2],'k--',alpha=.5)
plt.xlim(1E-2,2)
plt.ylim(1E-1,2)

plt.title('TNG50',y=0.92)
plt.xlabel(r'$Z_{{\rm inf,}~p}/Z_\odot$')
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.colorbar(label=r'$\log Z_{\rm SFR}/Z_{\rm mass}$')
plt.savefig('zsfrsinf_50.png')

allzsfratiola=np.array(allzsfratiola)
allzsfratioa=np.array(allzsfratioa)
allzsfratiol=np.array(allzsfratiol)
allzsfratiox=np.array(allzsfratiox)

plt.clf()
plt.scatter(allz0la,allznowla/allzsfratiola,c=np.log10(np.array(allzsfratiola)),vmin=-.2,vmax=1.3,alpha=.5)
plt.scatter(allz0a,allznowa/allzsfratioa,c=np.log10(np.array(allzsfratioa)),vmin=-.2,vmax=1.3,alpha=.5)
plt.scatter(allz0l,allznowl/allzsfratiol,c=np.log10(np.array(allzsfratiol)),vmin=-.2,vmax=1.3,alpha=.5)
plt.scatter(allz0x,allzfx/allzsfratiox,c=np.log10(np.array(allzsfratiox)),vmin=-.2,vmax=1.3,alpha=.5)
plt.loglog()
plt.plot([3E-3,2],[3E-3,2],'k--',alpha=.5)
plt.xlim(1E-2,2)
plt.ylim(1E-2,2)

plt.xlabel(r'$Z_{\rm inf,~p}/Z_\odot$')
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.colorbar(label=r'$\log Z_{\rm SFR}/Z_{\rm mass}$')
plt.savefig('zmzinf_50.png')

plt.clf()
#plt.scatter(allz0a,allznowa,c=(.008-np.array(alldza)/np.array(allsfea)/1E9)/(np.array(alletaa)+np.array(allmu)),vmin=-1E-3,vmax=1E-2)
plt.scatter(allz0a,allznowa,c=(.008-np.array(alldza)),vmin=-1E-3,vmax=1E-2)
plt.loglog()
plt.xlim(1E-4,2E-2)
plt.ylim(1E-4,2E-2)
plt.xlabel(r'$Z_{\rm inf}$')
plt.ylabel(r'$Z_{\rm gas}$')
plt.colorbar(label='Delta Z')
plt.savefig('z0zgdz_50.png')

plt.clf()
plt.scatter(allz0,allzf,c=np.log10(allmu))
plt.loglog()
plt.xlim(1E-4,2E-2)
plt.ylim(1E-4,2E-2)
plt.xlabel(r'$Z_{\rm inf}$')
plt.ylabel(r'$Z_{\rm gas}$')
plt.colorbar(label=r'$\log$ SFR/$M_{\rm dense~gas}$ (Gyr$^{-1}$)')
plt.savefig('z0zgas_50.png')

plt.clf()
plt.scatter(allz0,10**minfunc2(ft2['x']),c=np.log10(allsfe))
plt.loglog()
plt.xlim(1E-4,2E-2)
plt.ylim(1E-4,2E-2)
plt.xlabel(r'$Z_{\rm inf}$')
plt.ylabel(r'$Z_{\rm gas}$')
plt.colorbar(label=r'$\log$ SFR/$M_{\rm dense~gas}$ (Gyr$^{-1}$)')
plt.savefig('z0zgasmod2_50.png')

plt.clf()
plt.scatter(allz0,10**minfunc3(ft3['x']),c=np.log10(allsfe))
plt.loglog()
plt.xlim(1E-4,2E-2)
plt.ylim(1E-4,2E-2)
plt.xlabel(r'$Z_{\rm inf}$')
plt.ylabel(r'$Z_{\rm gas}$')
plt.colorbar(label=r'$\log$ SFR/$M_{\rm dense~gas}$ (Gyr$^{-1}$)')
plt.savefig('z0zgasmod3_50.png')

plt.clf()
plt.scatter(allsfr2/allmg0*1E9,allzf,c=np.log10(allsfe))
plt.loglog()
plt.xlim(1E-3,3)

plt.savefig('modtst0_50.png')

plt.clf()
plt.scatter(allsfr2/allmg0*1E9,10**minfunc3(ft3['x']),c=np.log10(allsfe))
plt.loglog()
plt.xlim(1E-3,3)
plt.savefig('modtst1_50.png')

from scipy import stats
plt.clf()
print(stats.spearmanr(np.log10(alleta),np.log10(allzf)-minfunc2(ft2['x'])))
plt.plot(np.log10(alleta),np.log10(allzf)-minfunc2(ft2['x']),'o')
plt.savefig('etavsresidual_50.png')

plt.clf()
plt.scatter(allsfe,allzf,c=np.log10(allz0))
plt.loglog()
plt.savefig('sfezz0_50.png')

plt.clf()
plt.scatter(allz0,allzf,c=np.log10(allsfe))
plt.loglog()
plt.savefig('z0sfezf_50.png')

plt.clf()
plt.scatter(allz0,10**minfunc2(ft2['x']),c=np.log10(allsfe))
plt.loglog()
plt.savefig('z0sfezfm_50.png')

plt.clf()
plt.plot(allsfr2,allzf,'o')
wsfr=np.where(allsfr2>0)[0]
print('coef sfr',getdeterminationcoef.getdeterminationcoef(np.log10(allsfr2)[wsfr],np.log10(allzf)[wsfr],deg=1))
plt.xlabel(r'SFR$_{z=0.13}$ (M$_\odot$~yr$^{-1}$)')
plt.ylabel(r'$Z_{\rm infall}/Z_\odot$')
plt.savefig('zvssfr.png')

plt.clf()
#plt.plot(allsfr2a,allz0a,'s',color='C2',label='XMD Analogs')
#plt.scatter(allsfr2a,allz0a,c=np.log10(allsfea),vmin=-1,vmax=0)
plt.scatter(allsfr2a,allz0a,c=np.log10(allsfea),vmin=-3,vmax=0)
plt.plot(allsfr2x,allz0x,'*',color='C0',label='XMDs')
plt.loglog()
plt.xlabel(r'SFR$_{z=0.13}$ (M$_\odot$~yr$^{-1}$)')
plt.ylabel(r'$Z_{\rm infall}/Z_\odot$')
plt.savefig('z0sfr_50.png')

plt.clf()
plt.scatter(allsfr2a,allznowa,c=np.log10(allmg0a))
plt.loglog()
plt.xlim(1E-4,1E-1)
plt.savefig('znowsfr_50.png')

plt.clf()
plt.plot(allsfr2a,allmz0a,'o')
plt.loglog()
plt.savefig('mz0sfr_50.png')

#plt.clf()
#plt.plot(allsfr2a,allsfea,'o')
#plt.loglog()
#plt.savefig('sfesfr2_50.png')

plt.clf()
plt.plot(allsfrwa,allz0a,'o')
plt.loglog()
plt.savefig('z0sfrw_50.png')

