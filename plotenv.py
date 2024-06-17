from scipy import stats
from astropy.io import fits, ascii
import numpy as np
import h5py
from astropy import cosmology
astrocosmo = cosmology.FlatLambdaCDM(H0=67.74,Om0=0.307)
from astropy.io import ascii
import matplotlib.pyplot as plt
import getperiodicdist
plt.style.use('../python/stylesheet.txt')

plt.rcParams['xtick.top']=True
plt.rcParams['xtick.labeltop']=False
f=fits.open('alltnggal_withr_wmetals_wvmax.fits')
x=ascii.read('sampleselection_sfr5_notide_nomass.txt')

wxmd=np.where(x['type']==1)[0]
wlmd=np.where(x['type']==2)[0]
wxmda=np.where(x['type']==3)[0]
wlmda=np.where(x['type']==4)[0]

mstarlim=3E7
sfrlim=.00001

wmst=np.where(f[1].data.mstar>mstarlim)[0]

d5=[]
dens5=[]

for i in range(len(x)):
    #dist=np.sqrt((f[1].data.x[wmst]-f[1].data.x[x['allgalindex'][i]])**2+(f[1].data.y[wmst]-f[1].data.y[x['allgalindex'][i]])**2)
    dist=getperiodicdist.getperiodicdist(np.array([f[1].data.x[x['allgalindex'][i]],f[1].data.y[x['allgalindex'][i]],f[1].data.z[x['allgalindex'][i]]]),np.array([f[1].data.x[wmst],f[1].data.y[wmst],f[1].data.z[wmst]]).T,boxsize=75E3/.6774)
    s=np.argsort(dist)
    d5.append(dist[s[6]])
    w5=np.where(dist<dist[s[6]])[0]
    dens5.append(np.sum(f[1].data.mstar[wmst[w5]])/(np.pi*4/3*d5[-1]**3))
    
d5=np.array(d5)
plt.clf()

fig,ax=plt.subplots(nrows=2,ncols=1,gridspec_kw={'hspace':0},sharey=True,figsize=(16,28))

h40=ax[0].hist([],bins=50,range=[0,8000],ec='C3',hatch='|',label='LMD Analogs',density=1,fc='None',lw=3)
h30=ax[1].hist([],bins=50,range=[0,8000],ec='C2',hatch='-',label='XMD Analogs',density=1,fc='None',lw=3)
h20=ax[0].hist([],bins=50,range=[0,8000],ec='C1',hatch='/',label='LMDs',density=1,fc='None',lw=3)
h10=ax[1].hist([],bins=50,range=[0,8000],ec='C0',hatch='\\',label='XMDs',density=1,fc='None',lw=3)

print('meds',np.median(d5[wlmda]),np.median(d5[wlmd]),np.median(d5[wxmda]),np.median(d5[wxmd]))
h4=ax[0].hist(d5[wlmda],bins=50,range=[0,8000],ec='C3',hatch='|',label='LMD Analogs',density=1,fc='None')
h3=ax[1].hist(d5[wxmda],bins=50,range=[0,8000],ec='C2',hatch='-',label='XMD Analogs',density=1,fc='None')
h2=ax[0].hist(d5[wlmd],bins=50,range=[0,8000],ec='C1',hatch='/',label='LMDs',density=1,fc='None')
h1=ax[1].hist(d5[wxmd],bins=50,range=[0,8000],ec='C0',hatch='\\',label='XMDs',density=1,fc='None')

h4a=ax[0].hist(d5[wlmda],bins=50,range=[0,8000],color='C3',histtype='step',label='LMD Analogs',density=1,lw=3)
h3a=ax[1].hist(d5[wxmda],bins=50,range=[0,8000],color='C2',histtype='step',label='XMD Analogs',density=1,lw=3)
h2a=ax[0].hist(d5[wlmd],bins=50,range=[0,8000],color='C1',histtype='step',label='LMDs',density=1,lw=3)
h1a=ax[1].hist(d5[wxmd],bins=50,range=[0,8000],color='C0',histtype='step',label='XMDs',density=1,lw=3)

ax[0].legend([h40[2],h20[2]],['LMD Analogs','LMDs'])
ax[1].legend([h30[2],h10[2]],['XMD Analogs','XMDs'])
#print(h4[2])

#plt.legend([(h4[2],h4a[2]),(h3[2],h3a[2]),(h2[2],h2a[2]),(h1[2],h1a[2])],['LMD Analogs','XMD Analogs','LMDs','XMDs'])
#plt.legend([h4[2],h3[2],h2[2],h1[2]],['LMD Analogs','XMD Analogs','LMDs','XMDs'])

ax[1].set_xlabel(r'$d_5$ (kpc)')
ax[0].set_ylabel(r'$dN/d d_5$ (kpc$^{-1}$)')
ax[1].set_ylabel(r'$dN/d d_5$ (kpc$^{-1}$)')
ax[0].set_xticklabels(['']*9)
from matplotlib.ticker import MultipleLocator
mly = MultipleLocator(200)
ax[0].xaxis.set_minor_locator(mly)
ax[1].xaxis.set_minor_locator(mly)
ax[0].set_xlim(0,8000)
ax[1].set_xlim(0,8000)
plt.savefig('d5env.png',dpi=300)

import bindata
plt.clf()
plt.plot(d5,f[1].data.gas_metal_sfr[x['allgalindex']],'o',alpha=.2)
plt.yscale('log')
b=bindata.bindata(d5,np.log10(f[1].data.gas_metal_sfr[x['allgalindex']]),np.linspace(0,7000,num=50))
plt.plot(b[2],10**b[0],'k')
print(stats.spearmanr(d5,f[1].data.gas_metal_sfr[x['allgalindex']]))
wv=np.where(d5>2000)[0]
print(len(wv),len(d5))
print(stats.spearmanr(d5[wv],f[1].data.gas_metal_sfr[x['allgalindex'][wv]]))
plt.savefig('metvsd5.png')

wtopenv=np.where(d5>np.percentile(d5,95))[0]
wbotenv=np.where(d5<np.percentile(d5,5))[0]
print('mettop',np.median(f[1].data.gas_metal_sfr[x['allgalindex'][wtopenv]]))
print('metbot',np.median(f[1].data.gas_metal_sfr[x['allgalindex'][wbotenv]]))


plt.clf()
plt.plot(np.log10(dens5),f[1].data.gas_metal_sfr[x['allgalindex']],'o')
plt.yscale('log')
plt.savefig('metvsdens5.png')

np.savetxt('tng100_d5.txt',np.array([x['allgalindex'],x['haloid'],x['type'],d5]).T)
