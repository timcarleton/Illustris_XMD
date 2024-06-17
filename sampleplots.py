import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import massmetallicity
import bindata
import getdeterminationcoef

plt.style.use('../python/stylesheet.txt')


f=fits.open('alltnggal_withr_wmetals_wvmax.fits')
f50=fits.open('tng50/tng50_wvmax.fits')

x=ascii.read('sampleselection_sfr5_notide_nomass.txt')
x50=ascii.read('sampleselection_tng50_sfr5_notide_nomass.txt')

obs1=ascii.read('littlethings_oh.csv')
obs2=ascii.read('realxmddata.csv')

wobs1=np.where((obs1['mstar']>.5) & (obs1['mstar']<30))[0]
obs1=obs1[wobs1]
obs1mst=obs1['mstar']*1E7
obs1mg=obs1['mgas']*1E7
obs1mh=10**obs1['m200']
print(obs1['logOH'])
wxmdhunter=np.where(10**(obs1['logOH']-8.69)<.1)[0]
wnmdhunter=np.where(10**(obs1['logOH']-8.69)>=.1)[0]
print(wxmdhunter,wnmdhunter)

wg=np.where((f[1].data.mstar<1E12) & (f[1].data.mstar>3E7) & (f[1].data.sfr>0) & (f[1].data.mhalo-f[1].data.mstar-f[1].data.mgas>1E8))[0]
wg50=np.where((f50[1].data.mstar<1E12) & (f50[1].data.mstar>3E7) & (f50[1].data.sfr>0))[0]


o_nfrac_sun = 10**(8.69-12)
fe_nfrac_sun = 10**(7.5-12)
n_nfrac_sun = 10**(7.83-12)
c_nfrac_sun = 10**(8.43-12)
ne_nfrac_sun = 10**(7.93-12)
mg_nfrac_sun = 10**(7.6-12)
si_nfrac_sun = 10**(7.51-12)

alpha_nfrac_sun=c_nfrac_sun+n_nfrac_sun+o_nfrac_sun+ne_nfrac_sun+mg_nfrac_sun+si_nfrac_sun

wxmd=x['allgalindex'][np.where(x['type']==1)[0]]
wlmd=x['allgalindex'][np.where(x['type']==2)[0]]
wxmda=x['allgalindex'][np.where(x['type']==3)[0]]
wlmda=x['allgalindex'][np.where(x['type']==4)[0]]

wxmd50=x50['allgalindex'][np.where(x50['type']==1)[0]]
wlmd50=x50['allgalindex'][np.where(x50['type']==2)[0]]
wxmda50=x50['allgalindex'][np.where(x50['type']==3)[0]]
wlmda50=x50['allgalindex'][np.where(x50['type']==4)[0]]

def getbindat(xdat,ydat,xscale,yscale,bins=[],nbins=20):
    if xscale=='log':
        xtmp=np.log10(xdat)
        if bins==[]:
            bins=np.logspace(min(xdat),max(xdat),num=nbins)
    else:
        xtmp=xdat
        if bins==[]:
            bins=np.linspace(min(xdat),max(xdat),num=nbins)
    if yscale=='log':
        ytmp=np.log10(ydat)
    else:
        ytmp=ydat
    b=bindata.bindata(xdat,ytmp,bins,median=1,logerr=yscale=='log')
    return2=b[2]
    if yscale=='log':
        return0=10**b[0]
        return1=10**b[1]
    else:
        return0=b[0]
        return1=b[1]

    return return0,return1,return2
    


def plotsamples(xdat,ydat,xlabel,ylabel,xscale,yscale,lmd=False,xlim=[],ylim=[],tng50=False,xmdbins=False,lmdbins=False,xmdabins=False,lmdabins=False,bins=[],nbins=20):

    if tng50:
        wgplot=wg50
        wxmdplot=wxmd50
        wlmdplot=wlmd50
        wxmdaplot=wxmda50
        wlmdaplot=wlmda50
    else:
        wgplot=wg
        wxmdplot=wxmd
        wlmdplot=wlmd
        wxmdaplot=wxmda
        wlmdaplot=wlmda        
    
    if xlim==[]:
        plt.hexbin(xdat[wgplot],ydat[wgplot],mincnt=3,alpha=.1,xscale=xscale,yscale=yscale)
    else:
        plt.hexbin(xdat[wgplot],ydat[wgplot],mincnt=3,alpha=.1,xscale=xscale,yscale=yscale,extent=[np.log10(xlim[0]),np.log10(xlim[1]),np.log10(ylim[0]),np.log10(ylim[1])])
        
    if lmd:
        plt.plot(xdat[wlmdaplot],ydat[wlmdaplot],'C3h',alpha=.25,label='LMD Analogs')

        if xscale=='log':
            allx=np.hstack([np.log10(xdat[wlmdaplot]),np.log10(xdat[wxmdaplot]),np.log10(xdat[wlmdplot]),np.log10(xdat[wxmdplot])])
        else:
            allx=np.hstack([xdat[wlmdaplot],xdat[wxmdaplot],xdat[wlmdplot],xdat[wxmdplot]])

        if xscale=='log':
            ally=np.hstack([np.log10(ydat[wlmdaplot]),np.log10(ydat[wxmdaplot]),np.log10(ydat[wlmdplot]),np.log10(ydat[wxmdplot])])
        else:
            ally=np.hstack([ydat[wlmdaplot],ydat[wxmdaplot],ydat[wlmdplot],ydat[wxmdplot]])

        print('coef',getdeterminationcoef.getdeterminationcoef(allx,ally,deg=1))

    plt.plot(xdat[wxmdaplot],ydat[wxmdaplot],'C2s',alpha=.25,label='XMD Analogs')
    if lmd:
        plt.plot(xdat[wlmdplot],ydat[wlmdplot],'C1x',alpha=.65,label='LMDs')
    plt.plot(xdat[wxmdplot],ydat[wxmdplot],'C0*',alpha=.65,label='XMDs')

    #plt.plot(xdat[wxmdaplot],ydat[wxmdaplot],'s',alpha=1,label='XMD Analogs',mfc='None',mec='C2')
    #plt.plot(xdat[wxmdplot],ydat[wxmdplot],'*',alpha=1,label='XMDs',mfc='None',mec='C0')


    if xlim==[]:
        if xscale=='log':
            bins=np.logspace(np.log10(min(xdat)),np.log10(max(xdat)),num=nbins)
        else:
            bins=np.linspace(min(xdat),max(xdat),num=nbins)
    else:
        if xscale=='log':
            bins=np.logspace(np.log10(xlim[0]),np.log10(xlim[1]),num=nbins)
        else:
            bins=np.linspace(xlim[0],xlim[1],num=nbins)

    if lmdabins:
        labins=getbindat(xdat[wlmdaplot],ydat[wlmdaplot],xscale,yscale,bins=bins,nbins=nbins)
        plt.plot(labins[2],labins[0],'C3-.')

    if xmdabins:
        xabins=getbindat(xdat[wxmdaplot],ydat[wxmdaplot],xscale,yscale,bins=bins,nbins=nbins)
        plt.plot(xabins[2],xabins[0],'C2--')

    if lmdbins:
        lbins=getbindat(xdat[wlmdplot],ydat[wlmdplot],xscale,yscale,bins=bins,nbins=nbins)
        plt.plot(lbins[2],lbins[0],'C1:')
    
    if xmdbins:
        xbins=getbindat(xdat[wxmdplot],ydat[wxmdplot],xscale,yscale,bins=bins,nbins=nbins)
        print('b',xbins)
        plt.plot(xbins[2],xbins[0],'C0')

    plt.xscale(xscale)
    plt.yscale(yscale)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if len(xlim)!=0:
        plt.xlim(xlim)
        plt.ylim(ylim)

    if tng50:
        plt.title('TNG50',y=0.92)
    else:
        plt.title('TNG100',y=0.92)

plt.clf()
plotsamples(f[1].data.mstar,f[1].data.gas_metal_sfr,r'$M_*~({\rm M_\odot})$',r'$Z_{\rm SFR}/Z_\odot$','log','log',lmd=True,xlim=[1E7,5E10],ylim=[.02,5])
xplot=np.logspace(7,12)
plt.plot(xplot,massmetallicity.getmzr(np.log10(xplot),solarunits=1,d02=0),color='k',label='Zahid+14',ls='--',alpha=.75)
#plt.errorbar(obs1mst[wnmdhunter],10**(obs1['logOH'][wnmdhunter]-8.69),mec='purple',mfc='None',marker='d',color='purple',yerr=obs1['logOHerr'][wnmdhunter]*10**(obs1['logOH'][wnmdhunter]-8.69),ls='None',label='LITTLE THINGS NMDs',mew=2,lw=2)
#plt.errorbar(obs1mst[wxmdhunter],10**(obs1['logOH'][wxmdhunter]-8.69),mec='purple',mfc='purple',marker='d',color='purple',yerr=obs1['logOHerr'][wxmdhunter]*10**(obs1['logOH'][wxmdhunter]-8.69),ls='None',label='LITTLE THINGS XMDs',mew=2,lw=2)
#plt.errorbar(obs2['stellar mass [Msol]'],10**(obs2['12 +(O/H)']-8.69),color='violet',ms=10,marker='*',xerr=obs2['M* err'],yerr=obs2['Z err']*10**(obs2['12 +(O/H)']-8.69)*np.log(10),label='Other XMDs',ls='None',lw=2)
plt.legend()
plt.ylim(.02,5)
#plt.title('TNG100')
plt.savefig('sampleplots/mzr.png')

plt.clf()
plotsamples(f50[1].data.mstar,f50[1].data.gas_metal_sfr,r'$M_*~({\rm M_\odot})$',r'$Z_{\rm SFR}/Z_\odot$','log','log',lmd=True,tng50=True,xlim=[1E7,5E10],ylim=[.02,5])
xplot=np.logspace(7,12)
plt.plot(xplot,massmetallicity.getmzr(np.log10(xplot),solarunits=1,d02=0),color='k',label='Zahid+14',ls='--',alpha=.75)
#plt.errorbar(obs1mst[wnmdhunter],10**(obs1['logOH'][wnmdhunter]-8.69),mec='purple',mfc='None',marker='d',color='purple',yerr=obs1['logOHerr'][wnmdhunter]*10**(obs1['logOH'][wnmdhunter]-8.69),ls='None',label='LITTLE THINGS NMDs',mew=2,lw=2)
#plt.errorbar(obs1mst[wxmdhunter],10**(obs1['logOH'][wxmdhunter]-8.69),mec='purple',mfc='purple',marker='d',color='purple',yerr=obs1['logOHerr'][wxmdhunter]*10**(obs1['logOH'][wxmdhunter]-8.69),ls='None',label='LITTLE THINGS XMDs',mew=2,lw=2)
#plt.errorbar(obs2['stellar mass [Msol]'],10**(obs2['12 +(O/H)']-8.69),color='violet',ms=10,marker='*',xerr=obs2['M* err'],yerr=obs2['Z err']*10**(obs2['12 +(O/H)']-8.69)*np.log(10),label='Other XMDs',ls='None',lw=2)

#plt.legend()
plt.ylim(.02,5)
#plt.title('TNG50')
plt.savefig('sampleplots/mzrtng50.png')

plt.clf()
plotsamples(f[1].data.gas_metal_sfr,f[1].data.gas_metal_mass,r'$Z_{\rm SFR}/Z_\odot$',r'$Z_{\rm mass}/Z_\odot$','log','log',lmd=True)
plt.plot([.01,10],[.01,10],ls='--',alpha=.75,color='k')
plt.legend()
plt.savefig('sampleplots/zsfrvszmass.png')

plt.clf()
plotsamples(f50[1].data.gas_metal_sfr,f50[1].data.gas_metal_mass,r'$Z_{\rm SFR}/Z_\odot$',r'$Z_{\rm mass}/Z_\odot$','log','log',lmd=True,tng50=True)
xplot=np.logspace(7,12)
plt.plot([.01,10],[.01,10],ls='--',alpha=.75,color='k')
plt.legend()
plt.savefig('sampleplots/zsfrvszmass50.png')
        
plt.clf()
plotsamples(f[1].data.mstar,f[1].data.mstar/f[1].data.mhalo,r'$M_{\rm *}~({\rm M_\odot})$',r'$M_*/M_{\rm halo}$','log','log',xmdbins=True,xmdabins=True,xlim=[3E7,3E8],ylim=[3E-4,.01])
#plt.plot(obs1mh[wnmdhunter],obs1mst[wnmdhunter]/obs1mh[wnmdhunter],'d',ms=10,mew=2,mec='purple',mfc='None',ls='None')
#plt.plot(obs1mh[wxmdhunter],obs1mst[wxmdhunter]/obs1mh[wxmdhunter],'d',ms=10,mew=2,mec='purple',mfc='purple',ls='None')
plt.savefig('sampleplots/mstarhalo.png')

plt.clf()
plotsamples(f50[1].data.mstar,f50[1].data.mstar/f50[1].data.mhalo,r'$M_{\rm *}~({\rm M_\odot})$',r'$M_*/M_{\rm halo}$','log','log',xmdbins=True,xmdabins=True,xlim=[3E7,3E8],ylim=[3E-4,.01],tng50=True,lmd=True,lmdbins=True,lmdabins=True)
#plt.plot(obs1mh[wnmdhunter],obs1mst[wnmdhunter]/obs1mh[wnmdhunter],'d',ms=10,mew=2,mec='purple',mfc='None',ls='None')
#plt.plot(obs1mh[wxmdhunter],obs1mst[wxmdhunter]/obs1mh[wxmdhunter],'d',ms=10,mew=2,mec='purple',mfc='purple',ls='None')
plt.savefig('sampleplots/mstarhalo_50.png')

plt.clf()
plotsamples(f[1].data.mhalo,f[1].data.mstar/f[1].data.mhalo,r'$M_{\rm halo}~({\rm M_\odot})$',r'$M_*/M_{\rm halo}$','log','log',xmdbins=True,xmdabins=True)
#plt.plot(obs1mh[wnmdhunter],obs1mst[wnmdhunter]/obs1mh[wnmdhunter],'d',ms=10,mew=2,mec='purple',mfc='None',ls='None')
#plt.plot(obs1mh[wxmdhunter],obs1mst[wxmdhunter]/obs1mh[wxmdhunter],'d',ms=10,mew=2,mec='purple',mfc='purple',ls='None')
plt.savefig('sampleplots/mhalostar.png')

plt.clf()
plotsamples(f[1].data.vmax,f[1].data.mstar,r'$V_{\rm max}~({\rm km~s^{-1}})$',r'$M_*~({\rm M_\odot})$','log','log',xmdbins=True,xmdabins=True)
plt.errorbar(obs2['W50 [km/s]'],obs2['stellar mass [Msol]'],color='violet',ms=10,marker='*',xerr=obs2['W50 err'],yerr=obs2['M* err'],ls='None')
plt.plot(obs1['vmax'][wnmdhunter],obs1mst[wnmdhunter],'d',ms=10,mew=2,mec='purple',mfc='None',ls='None')
plt.plot(obs1['vmax'][wxmdhunter],obs1mst[wxmdhunter],'d',ms=10,mew=2,mec='purple',mfc='purple',ls='None')

plt.savefig('sampleplots/vmaxmstar.png')

plt.clf()
plotsamples(f[1].data.mstar,f[1].data.vmax,r'$M_*~({\rm M_\odot})$',r'$V_{\rm max}~({\rm km~s^{-1}})$','log','log',xmdbins=True,xmdabins=True,xlim=[3E7,3E8],ylim=[25,100])
#plt.errorbar(obs2['W50 [km/s]'],obs2['stellar mass [Msol]'],color='violet',ms=10,marker='*',xerr=obs2['W50 err'],yerr=obs2['M* err'],ls='None')
#plt.plot(obs1['vmax'][wnmdhunter],obs1mst[wnmdhunter],'d',ms=10,mew=2,mec='purple',mfc='None',ls='None')
#plt.plot(obs1['vmax'][wxmdhunter],obs1mst[wxmdhunter],'d',ms=10,mew=2,mec='purple',mfc='purple',ls='None')

plt.savefig('sampleplots/mstarvmax.png')

plt.clf()
plotsamples(f50[1].data.mstar,f50[1].data.vmax,r'$M_*~({\rm M_\odot})$',r'$V_{\rm max}~({\rm km~s^{-1}})$','log','log',xmdbins=True,xmdabins=True,xlim=[3E7,3E8],ylim=[25,100],tng50=True)
#plt.errorbar(obs2['W50 [km/s]'],obs2['stellar mass [Msol]'],color='violet',ms=10,marker='*',xerr=obs2['W50 err'],yerr=obs2['M* err'],ls='None')
#plt.plot(obs1['vmax'][wnmdhunter],obs1mst[wnmdhunter],'d',ms=10,mew=2,mec='purple',mfc='None',ls='None')
#plt.plot(obs1['vmax'][wxmdhunter],obs1mst[wxmdhunter],'d',ms=10,mew=2,mec='purple',mfc='purple',ls='None')

plt.savefig('sampleplots/mstarvmax_50.png')


plt.clf()
plotsamples(f[1].data.mstar,f[1].data.vmax,r'$M_*~({\rm M_\odot})$',r'$V_{\rm max}~({\rm km~s^{-1}})$','log','log',xmdbins=True,xmdabins=True,xlim=[3E7,3E8],ylim=[25,100])
#plt.errorbar(obs2['W50 [km/s]'],obs2['stellar mass [Msol]'],color='violet',ms=10,marker='*',xerr=obs2['W50 err'],yerr=obs2['M* err'],ls='None')
#plt.plot(obs1['vmax'][wnmdhunter],obs1mst[wnmdhunter],'d',ms=10,mew=2,mec='purple',mfc='None',ls='None')
#plt.plot(obs1['vmax'][wxmdhunter],obs1mst[wxmdhunter],'d',ms=10,mew=2,mec='purple',mfc='purple',ls='None')

plt.savefig('sampleplots/mstarvmax.png')

#plt.clf()
#plotsamples(f50[1].data.mstar,f50[1].data.vmax,r'$M_*~({\rm M_\odot})$',r'$V_{\rm max}~({\rm km~s^{-1}})$','log','log',xmdbins=True,xmdabins=True,xlim=[3E7,3E8],ylim=[25,100],tng50=True)
#plt.errorbar(obs2['W50 [km/s]'],obs2['stellar mass [Msol]'],color='violet',ms=10,marker='*',xerr=obs2['W50 err'],yerr=obs2['M* err'],ls='None')
#plt.plot(obs1['vmax'][wnmdhunter],obs1mst[wnmdhunter],'d',ms=10,mew=2,mec='purple',mfc='None',ls='None')
#plt.plot(obs1['vmax'][wxmdhunter],obs1mst[wxmdhunter],'d',ms=10,mew=2,mec='purple',mfc='purple',ls='None')

#plt.savefig('sampleplots/mstarvmax_50.png')

plt.clf()
plotsamples(f50[1].data.mstar,f50[1].data.half_mass_rad_star,r'$M_{\rm *}~({\rm M_\odot})$',r'$r_{1/2,*}$ (kpc)','log','log',xlim=[3E7,3E8],ylim=[.5,10],xmdbins=True,xmdabins=True,tng50=True)
#plt.plot(obsmst1[wnmdhunter],obs1['Rd'][wnmdhunter],markeredgecolor='purple',ms=10,marker='d',color='purple',markerfacecolor='None',markeredgewidth=2,ls='None')
#plt.plot(obsmst1[wxmdhunter],obs1['Rd'][wxmdhunter],markeredgecolor='purple',ms=10,marker='d',color='purple',markerfacecolor='purple',markeredgewidth=2,ls='None')
plt.savefig('sampleplots/masssize_50.png')

plt.clf()
plotsamples(f[1].data.mstar,f[1].data.half_mass_rad_star,r'$M_{\rm *}~({\rm M_\odot})$',r'$r_{1/2,*}$ (kpc)','log','log',xlim=[3E7,3E8],ylim=[.5,10],xmdbins=True,xmdabins=True,tng50=False)
#plt.plot(obsmst1[wnmdhunter],obs1['Rd'][wnmdhunter],markeredgecolor='purple',ms=10,marker='d',color='purple',markerfacecolor='None',markeredgewidth=2,ls='None')
#plt.plot(obsmst1[wxmdhunter],obs1['Rd'][wxmdhunter],markeredgecolor='purple',ms=10,marker='d',color='purple',markerfacecolor='purple',markeredgewidth=2,ls='None')
plt.savefig('sampleplots/masssize.png')

plt.clf()
plotsamples(f[1].data.mhalo,(f[1].data.mstar+f[1].data.mgas)/f[1].data.mhalo,r'$M_{\rm halo}~({\rm M_\odot})$',r'$M_{\rm bar}/M_{\rm halo}$','log','log',xmdbins=True,xmdabins=True,xlim=[1E10,2E11],ylim=[.01,.2])
plt.plot(obs1mh[wnmdhunter],(obs1mst[wnmdhunter]+obs1mg[wnmdhunter])/obs1mh[wnmdhunter],'d',ms=10,mew=2,mec='purple',mfc='None')
plt.plot(obs1mh[wxmdhunter],(obs1mst[wxmdhunter]+obs1mg[wxmdhunter])/obs1mh[wxmdhunter],'d',ms=10,mew=2,mec='purple',mfc='purple')

plt.savefig('sampleplots/mbarmhalo.png')

plt.clf()
plotsamples(f50[1].data.mhalo,(f50[1].data.mstar+f50[1].data.mgas)/f50[1].data.mhalo,r'$M_{\rm halo}~({\rm M_\odot})$',r'$M_{\rm bar}/M_{\rm halo}$','log','log',xmdbins=True,xmdabins=True,xlim=[1E10,2E11],ylim=[.01,.2],tng50=True)
plt.plot(obs1mh[wnmdhunter],(obs1mst[wnmdhunter]+obs1mg[wnmdhunter])/obs1mh[wnmdhunter],'d',ms=10,mew=2,mec='purple',mfc='None')
plt.plot(obs1mh[wxmdhunter],(obs1mst[wxmdhunter]+obs1mg[wxmdhunter])/obs1mh[wxmdhunter],'d',ms=10,mew=2,mec='purple',mfc='purple')

plt.savefig('sampleplots/mbarmhalo_50.png')

plt.clf()
plotsamples(f[1].data.mstar,f[1].data.gas_metal_mass,r'$M_{\rm *}~({\rm M_\odot})$',r'$Z_{\rm mass}/Z_\odot$','log','log')
plt.savefig('sampleplots/mstargasmetalmass.png')

plt.clf()
plotsamples(f[1].data.mstar,f[1].data.mgas/f[1].data.sfr,r'$M_{\rm *}~({\rm M_\odot})$',r'$M_{\rm gas}/{\rm SFR}~({\rm yr^{-1}})$','log','log')
#plt.errorbar(obs2['stellar mass [Msol]'],obs2['HI mass [Msol]']/obs2['SFR_Ha [Msol]'],color='violet',ms=10,marker='^',xerr=obs2['M* err'],yerr=np.sqrt((obs2['HI err']/obs2['SFR_Ha [Msol]'])**2 + (obs2['SFR err']*obs2['HI mass [Msol]']/obs2['SFR_Ha [Msol]']**2)**2),ls='None')

plt.savefig('sampleplots/mstarsfe.png')

plt.clf()
print('sfe')
plotsamples(f[1].data.mgas/f[1].data.sfr,f[1].data.gas_metal_mass,r'$M_{\rm gas}/{\rm SFR}~({\rm yr^{-1}})$',r'$Z_{\rm mass}/Z_\odot$','log','log')
plt.errorbar(obs2['HI mass [Msol]']/obs2['SFR_Ha [Msol]'],10**(obs2['12 +(O/H)']-8.69),color='violet',ms=10,marker='*',xerr=np.sqrt((obs2['HI err']/obs2['SFR_Ha [Msol]'])**2 + (obs2['SFR err']*obs2['HI mass [Msol]']/obs2['SFR_Ha [Msol]']**2)**2),yerr=obs2['Z err']*10**(obs2['12 +(O/H)']-8.69)*np.log(10),ls='None')
plt.errorbar(obs1mg*10**-obs1['logSFRHa'],10**(obs1['logOH']-8.69),color='purple',ms=10,marker='d',xerr=obs1['logSFRHaerr']*obs1mg*10**-obs1['logSFRHa']*np.log(10),ls='None')
plt.savefig('sampleplots/massmetalsfe.png')

plt.clf()
plotsamples(f50[1].data.mgas/f50[1].data.sfr,f50[1].data.gas_metal_mass,r'$Z_{\rm SFR}/Z_\odot$',r'$Z_{\rm mass}/Z_\odot$','log','log',lmd=True,tng50=True)
plt.errorbar(obs1mg*10**-obs1['logSFRHa'],10**(obs1['logOH']-8.69),ms=10,marker='d',color='purple',xerr=obs1['logSFRHaerr']*obs1mg*10**-obs1['logSFRHa']*np.log(10),ls='None')
plt.errorbar(obs2['HI mass [Msol]']/obs2['SFR_Ha [Msol]'],10**(obs2['12 +(O/H)']-8.69),color='violet',ms=10,marker='*',xerr=np.sqrt((obs2['HI err']/obs2['SFR_Ha [Msol]'])**2 + (obs2['SFR err']*obs2['HI mass [Msol]']/obs2['SFR_Ha [Msol]']**2)**2),yerr=obs2['Z err']*10**(obs2['12 +(O/H)']-8.69)*np.log(10),ls='None')
plt.savefig('sampleplots/massmetalsfe50.png')

plt.clf()
plotsamples(f[1].data.sfr/f[1].data.mgas*1E9,f[1].data.gas_metal_sfr,r'${\rm SFR}/M_{\rm gas}~({\rm Gyr^{-1}})$',r'$Z_{\rm SFR}/Z_\odot$','log','log',lmd=True)
#plt.errorbar(obs2['SFR_Ha [Msol]']/obs2['HI mass [Msol]'],10**(obs2['12 +(O/H)']-8.69),color='violet',ms=10,marker='*',xerr=np.sqrt((obs2['SFR err']/obs2['HI mass [Msol]'])**2 + (obs2['HI err']*obs2['SFR_Ha [Msol]']/obs2['HI Mass [Msol]']**2)**2),yerr=obs2['Z err']*10**(obs2['12 +(O/H)']-8.69)*np.log(10),ls='None')
#plt.errorbar(obs1mg*10**-obs1['logSFRHa'],10**(obs1['logOH']-8.69),color='purple',ms=10,marker='d',xerr=obs1['logSFRHaerr']*obs1mg*10**-obs1['logSFRHa']*np.log(10),ls='None')

plt.savefig('sampleplots/metalsfe.png')

plt.clf()
plotsamples(f[1].data.mgas/f[1].data.mstar,f[1].data.gas_metal_sfr,r'$M_{\rm gas}/{\rm M_*}$',r'$Z_{\rm SFR}/Z_\odot$','log','log')
plt.errorbar(obs2['HI mass [Msol]']/obs2['stellar mass [Msol]'],10**(obs2['12 +(O/H)']-8.69),color='violet',ms=10,marker='*',xerr=np.sqrt((obs2['HI err']/obs2['stellar mass [Msol]'])**2 + (obs2['M* err']*obs2['HI mass [Msol]']/obs2['stellar mass [Msol]']**2)**2),yerr=obs2['Z err']*10**(obs2['12 +(O/H)']-8.69)*np.log(10),ls='None')
plt.plot(obs1mg/obs1mst,10**(obs1['logOH']-8.69),color='purple',ms=10,marker='d',ls='None')
plt.savefig('sampleplots/metalgasfrac.png')

plt.clf()
plotsamples(f[1].data.half_mass_rad_star,f[1].data.mgas/f[1].data.sfr,r'$r_{\rm 1/2,*}~({\rm kpc})$',r'$M_{\rm gas}/{\rm SFR}~({\rm yr^{-1}})$','log','log')
plt.savefig('sampleplots/rstarsfe.png')

plt.clf()
plotsamples(f[1].data.sfr/(2*np.pi*f[1].data.half_mass_rad_star**2),f[1].data.gas_metal_sfr,r'$\Sigma_{\rm SFR}$',r'$logzson$','log','log')
plt.savefig('sampleplots/metalsfrdens.png')

plt.clf()
plotsamples(f[1].data.gas_metal_sfr,f[1].data.star_metal,r'$Z_{\rm SFR}/Z_\odot$',r'$Z_{*}/Z_\odot$','log','log',xlim=[.01,2],ylim=[.01,2])
plt.savefig('sampleplots/gasmetalstarmetal.png')

plt.clf()
plotsamples(f50[1].data.gas_metal_sfr,f50[1].data.star_metal,r'$Z_{\rm SFR}/Z_\odot$',r'$Z_{*}/Z_\odot$','log','log',xlim=[.01,10],ylim=[1E-3,5],tng50=True)
plt.savefig('sampleplots/gasmetalstarmetal_50.png')

plt.clf()
rvir=(f[1].data.mhalo/(4/3*np.pi*135.99294735*200))**(1/3)
plotsamples(rvir/f[1].data.rmax,f[1].data.vmax,r'$R_{\rm 200}/R_{\rm max}$',r'$V_{\rm max}$ (km/s)','log','log')
plt.savefig('sampleplots/vmaxrmaxplot.png')

plt.clf()
plotsamples(f[1].data.mhalo,np.sqrt(f[1].data.xspin**2+f[1].data.yspin**2+f[1].data.zspin**2),r'$M_{\rm halo}$ (M$_\odot$)',r'$\lambda$ (kpc km/s)','log','log')
plt.savefig('sampleplots/mhalospin.png')


plt.clf()
#b1=bindata.bindata(np.log10(f[1].data.fe_h_star/fe_nfrac_sun)[wxmd],np.log10((f[1].data.c_h_star+f[1].data.o_h_star+f[1].data.ne_h_star+f[1].data.mg_h_star+f[1].data.si_h_star)/(f[1].data.fe_h_star)*(fe_nfrac_sun/alpha_nfrac_sun))[wxmd],np.linspace(-2,0,num=15))
#b2=bindata.bindata(np.log10(f[1].data.fe_h_star/fe_nfrac_sun)[wxmda],np.log10((f[1].data.c_h_star+f[1].data.o_h_star+f[1].data.ne_h_star+f[1].data.mg_h_star+f[1].data.si_h_star)/(f[1].data.fe_h_star)*(fe_nfrac_sun/alpha_nfrac_sun))[wxmda],np.linspace(-2,0,num=15))
plotsamples(f[1].data.fe_h_star/fe_nfrac_sun,(f[1].data.c_h_star+f[1].data.o_h_star+f[1].data.ne_h_star+f[1].data.mg_h_star+f[1].data.si_h_star)/(f[1].data.fe_h_star)*(fe_nfrac_sun/alpha_nfrac_sun),r'[Fe/H]$_*$',r'[$\alpha$/Fe]$_*$','log','log',xlim=[1E-2,.5],ylim=[1.8,4],xmdbins=True,xmdabins=True)
#plt.plot(10**b1[2],10**b1[0],color='C0')
#plt.plot(10**b2[2],10**b2[0],color='C2')
plt.savefig('sampleplots/alphafe.png')

plt.clf()
#b1=bindata.bindata(np.log10(f[1].data.fe_h_star/fe_nfrac_sun)[wxmd],np.log10((f[1].data.c_h_star+f[1].data.o_h_star+f[1].data.ne_h_star+f[1].data.mg_h_star+f[1].data.si_h_star)/(f[1].data.fe_h_star)*(fe_nfrac_sun/alpha_nfrac_sun))[wxmd],np.linspace(-2,0,num=15))
#b2=bindata.bindata(np.log10(f[1].data.fe_h_star/fe_nfrac_sun)[wxmda],np.log10((f[1].data.c_h_star+f[1].data.o_h_star+f[1].data.ne_h_star+f[1].data.mg_h_star+f[1].data.si_h_star)/(f[1].data.fe_h_star)*(fe_nfrac_sun/alpha_nfrac_sun))[wxmda],np.linspace(-2,0,num=15))
plotsamples(f50[1].data.fe_h_star/fe_nfrac_sun,(f50[1].data.c_h_star+f50[1].data.o_h_star+f50[1].data.ne_h_star+f50[1].data.mg_h_star+f50[1].data.si_h_star)/(f50[1].data.fe_h_star)*(fe_nfrac_sun/alpha_nfrac_sun),r'[Fe/H]$_*$',r'[$\alpha$/Fe]$_*$','log','log',xlim=[1E-2,.5],ylim=[1.8,4],xmdbins=True,xmdabins=True,tng50=True)
#plt.plot(10**b1[2],10**b1[0],color='C0')
#plt.plot(10**b2[2],10**b2[0],color='C2')
plt.savefig('sampleplots/alphafe_50.png')


plt.clf()
b1=bindata.bindata(np.log10(f[1].data.sfr)[wxmd],np.log10((f[1].data.c_h_star+f[1].data.o_h_star+f[1].data.ne_h_star+f[1].data.mg_h_star+f[1].data.si_h_star)/(f[1].data.fe_h_star)*(fe_nfrac_sun/alpha_nfrac_sun))[wxmd],np.linspace(-3,-1,num=15))
b2=bindata.bindata(np.log10(f[1].data.sfr)[wxmda],np.log10((f[1].data.c_h_star+f[1].data.o_h_star+f[1].data.ne_h_star+f[1].data.mg_h_star+f[1].data.si_h_star)/(f[1].data.fe_h_star)*(fe_nfrac_sun/alpha_nfrac_sun))[wxmda],np.linspace(-3,-1,num=15))
plotsamples(f[1].data.sfr,(f[1].data.c_h_star+f[1].data.o_h_star+f[1].data.ne_h_star+f[1].data.mg_h_star+f[1].data.si_h_star)/(f[1].data.fe_h_star)*(fe_nfrac_sun/alpha_nfrac_sun),r'SFR (M$_\odot$~yr$^{-1}$)',r'[$\alpha$/Fe]$_*$','log','log',xlim=[1E-3,.2],ylim=[1.3,4])
plt.plot(10**b1[2],10**b1[0],color='C0')
plt.plot(10**b2[2],10**b2[0],color='C2')
plt.savefig('sampleplots/sfralphafe.png')


plt.clf()
#b1=bindata.bindata(np.log10(f[1].data.sfr)[wxmd],np.log10((f[1].data.c_h_star+f[1].data.o_h_star+f[1].data.ne_h_star+f[1].data.mg_h_star+f[1].data.si_h_star)/(f[1].data.fe_h_star)*(fe_nfrac_sun/alpha_nfrac_sun))[wxmd],np.linspace(-3,-1,num=15))
#b2=bindata.bindata(np.log10(f[1].data.sfr)[wxmda],np.log10((f[1].data.c_h_star+f[1].data.o_h_star+f[1].data.ne_h_star+f[1].data.mg_h_star+f[1].data.si_h_star)/(f[1].data.fe_h_star)*(fe_nfrac_sun/alpha_nfrac_sun))[wxmda],np.linspace(-3,-1,num=15))
plotsamples(f[1].data.sfr,f[1].data.fe_h_star/fe_nfrac_sun,r'SFR (M$_\odot$~yr$^{-1}$)',r'[Fe/H]$_*$','log','log',xlim=[1E-3,.2],ylim=[1E-2,3])
#plt.plot(10**b1[2],10**b1[0],color='C0')
#plt.plot(10**b2[2],10**b2[0],color='C2')
plt.savefig('sampleplots/sfrfeh.png')



plt.clf()
b1=bindata.bindata(np.log10(f[1].data.fe_h_star/fe_nfrac_sun)[wxmd],np.log10((f[1].data.mg_h_star)/(f[1].data.fe_h_star)*(fe_nfrac_sun/alpha_nfrac_sun))[wxmd],np.linspace(-2,0,num=15),logerr=1)
b2=bindata.bindata(np.log10(f[1].data.fe_h_star/fe_nfrac_sun)[wxmda],np.log10((f[1].data.mg_h_star)/(f[1].data.fe_h_star)*(fe_nfrac_sun/alpha_nfrac_sun))[wxmda],np.linspace(-2,0,num=15),logerr=1)
plotsamples(f[1].data.fe_h_star/fe_nfrac_sun,(f[1].data.mg_h_star)/(f[1].data.fe_h_star)*(fe_nfrac_sun/alpha_nfrac_sun),r'[Fe/H]$_*$',r'[Mg/Fe]$_*$','log','log',xlim=[1E-2,3],ylim=[.07,.13])
plt.errorbar(10**b1[2],10**b1[0],color='C0',yerr=b1[1])
plt.errorbar(10**b2[2],10**b2[0],color='C2',yerr=b2[1])
plt.savefig('sampleplots/mgfe.png')

plt.clf()
plotsamples(f[1].data.gas_metal_sfr,f[1].data.sdss_g-f[1].data.sdss_r,r'$Z_{\rm SFR}$',r'$g-r$','log','linear')
plt.savefig('sampleplots/colormet.png')

plt.clf()
plotsamples(f[1].data.bessel_u-f[1].data.bessel_v,f[1].data.sdss_z-f[1].data.bessel_k,r'$U-V$',r'$z-K$','linear','linear')
plt.savefig('sampleplots/colorcolor.png')
