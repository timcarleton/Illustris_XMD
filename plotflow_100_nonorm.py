import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from astropy import cosmology
import os
import getdeterminationcoef
from mpl_toolkits.mplot3d import Axes3D
astrocosmo=cosmology.FlatLambdaCDM(H0=67.7,Om0=0.307)

plt.style.use('../python/stylesheet.txt')
plt.rcParams['xtick.top']=False
plt.rcParams['xtick.labeltop']=False

obj=np.loadtxt('sampleselection_sfr5_notide_nomass.txt',delimiter=',',skiprows=1)
zg=np.array([0,.1,.25,.5,.75,1,1.5,2,3,7])
zglabels=zg.astype(str)
agepoints=astrocosmo.age(zg)

doiplots=False

#envinfo=ascii.read('xmd_env.txt')

wcat_x=np.where(obj[:,2]==1)[0]
wcat_xa=np.where(obj[:,2]==3)[0]
wcat_l=np.where(obj[:,2]==2)[0]
wcat_la=np.where(obj[:,2]==4)[0]

allzsxa=np.zeros([len(wcat_xa),100])+np.nan
allzsx=np.zeros([len(wcat_x),100])+np.nan
allzsl=np.zeros([len(wcat_l),100])+np.nan
allzsla=np.zeros([len(wcat_la),100])+np.nan

allzinfoxa=np.zeros([len(wcat_xa),100])+np.nan
allzinfox=np.zeros([len(wcat_x),100])+np.nan
allzinfol=np.zeros([len(wcat_l),100])+np.nan
allzinfola=np.zeros([len(wcat_la),100])+np.nan

allzinfovxa=np.zeros([len(wcat_xa),100])+np.nan
allzinfovx=np.zeros([len(wcat_x),100])+np.nan
allzinfovl=np.zeros([len(wcat_l),100])+np.nan
allzinfovla=np.zeros([len(wcat_la),100])+np.nan

allflowinfoxa=np.zeros([len(wcat_xa),100])+np.nan
allflowinfox=np.zeros([len(wcat_x),100])+np.nan
allflowinfol=np.zeros([len(wcat_l),100])+np.nan
allflowinfola=np.zeros([len(wcat_la),100])+np.nan

alloflowinfoxa=np.zeros([len(wcat_xa),100])+np.nan
alloflowinfox=np.zeros([len(wcat_x),100])+np.nan
alloflowinfol=np.zeros([len(wcat_l),100])+np.nan
alloflowinfola=np.zeros([len(wcat_la),100])+np.nan

allzoflowinfoxa=np.zeros([len(wcat_xa),100])+np.nan
allzoflowinfox=np.zeros([len(wcat_x),100])+np.nan
allzoflowinfol=np.zeros([len(wcat_l),100])+np.nan
allzoflowinfola=np.zeros([len(wcat_la),100])+np.nan

allsfrxa=np.zeros([len(wcat_xa),100])+np.nan
allsfrx=np.zeros([len(wcat_x),100])+np.nan
allsfrl=np.zeros([len(wcat_l),100])+np.nan
allsfrla=np.zeros([len(wcat_la),100])+np.nan

alletaxa=np.zeros([len(wcat_xa),100])+np.nan
alletax=np.zeros([len(wcat_x),100])+np.nan
alletal=np.zeros([len(wcat_l),100])+np.nan
alletala=np.zeros([len(wcat_la),100])+np.nan

alltdxa=np.zeros([len(wcat_xa),100])+np.nan
alltdx=np.zeros([len(wcat_x),100])+np.nan
alltdl=np.zeros([len(wcat_l),100])+np.nan
alltdla=np.zeros([len(wcat_la),100])+np.nan

allvflowinfoxa=np.zeros([len(wcat_xa),100])+np.nan
allvflowinfox=np.zeros([len(wcat_x),100])+np.nan
allvflowinfol=np.zeros([len(wcat_l),100])+np.nan
allvflowinfola=np.zeros([len(wcat_la),100])+np.nan

alldflowinfoxa=np.zeros([len(wcat_xa),100])+np.nan
alldflowinfox=np.zeros([len(wcat_x),100])+np.nan
alldflowinfol=np.zeros([len(wcat_l),100])+np.nan
alldflowinfola=np.zeros([len(wcat_la),100])+np.nan

allmgxa=np.zeros([len(wcat_xa),100])+np.nan
allmgx=np.zeros([len(wcat_x),100])+np.nan
allmgl=np.zeros([len(wcat_l),100])+np.nan
allmgla=np.zeros([len(wcat_la),100])+np.nan

allrhxa=np.zeros([len(wcat_xa),100])+np.nan
allrhx=np.zeros([len(wcat_x),100])+np.nan
allrhl=np.zeros([len(wcat_l),100])+np.nan
allrhla=np.zeros([len(wcat_la),100])+np.nan

allmetxa=np.zeros([len(wcat_xa),100])+np.nan
allmetx=np.zeros([len(wcat_x),100])+np.nan
allmetl=np.zeros([len(wcat_l),100])+np.nan
allmetla=np.zeros([len(wcat_la),100])+np.nan

env=np.loadtxt('tng100_d5.txt')

wenv_x=np.where(env[:,2]==1)[0]
wenv_xa=np.where(env[:,2]==3)[0]
wenv_l=np.where(env[:,2]==2)[0]
wenv_la=np.where(env[:,2]==4)[0]

alltdep=[]

def flowplot1(x,savename,title):
    ax0=plt.axes()
    ax1=plt.twinx()

    deltat=np.append(astrocosmo.age(x['redshift'])[1],astrocosmo.age(x['redshift'])[1:]-astrocosmo.age(x['redshift'])[0:-1])*1E9
    
    infline =ax0.plot(astrocosmo.age(x['redshift']),x['minflow_5rh']/deltat,'--',color='C5')
    outline=ax0.plot(astrocosmo.age(x['redshift']),x['moutflow_5rh']/deltat,'-',color='C6')

    zline=ax1.plot(astrocosmo.age(x['redshift']),x['sfr_weight_metallicity'],":",color='C7')
    ax0.set_ylabel(r'Flow Rate within $5 r_h$ (M$_\odot$ yr$^{-1}$)')
    ax1.set_yscale('log')
    #w=np.where(x['redshift']<2)[0]
    #ax0.set_ylim(0,1.25*max(x['moutflow_5rh'][w]/x['mstar'][w]))
    ax0.set_ylim(0,5)
    ax1.set_ylabel('SFR-Weighted Metallicity')
    ax1.legend([infline[0],outline[0],zline[0]],[r'$M_{\rm inflow,5rh}$',r'$M_{\rm outflow,5rh}$',r'$Z_{\rm SFR}/Z_\odot$'],loc='lower left')
    ax1.set_ylim(1E-4,2E-2)

    plt.title(title,y=1.1)
    ax1.set_xlim(0.5,13.9)
    ax0.set_xlabel(r'$t_{\rm age}$ (Gyr)')
    ax=plt.twiny()
    ax.set_xticks(astrocosmo.age(zg[::-1]).value)
    ax.set_xticklabels(zglabels[::-1])
    ax.set_xlim(0.5,13.9)
    ax.set_xlabel(r'$z$')

    plt.savefig(savename)

def flowplot2(x,savename,title):
    ax0=plt.axes()
    ax1=plt.twinx()

    deltat=np.append(astrocosmo.age(x['redshift'])[1],astrocosmo.age(x['redshift'])[1:]-astrocosmo.age(x['redshift'])[0:-1])*1E9
    
    infline =ax0.plot(astrocosmo.age(x['redshift']),x['mzinflow_5rh']/x['minflow_5rh'],'--',color='C5')
    outline=ax0.plot(astrocosmo.age(x['redshift']),x['mzoutflow_5rh']/x['moutflow_5rh'],'-',color='C6')

    zline=ax1.plot(astrocosmo.age(x['redshift']),x['sfr_weight_metallicity'],":",color='C7')
    ax0.set_ylabel(r'Flow Metallicity within $5 r_h$')
    ax1.set_yscale('log')
    ax0.set_yscale('log')
    ax0.set_ylim(1E-4,1E-2)
    ax1.set_ylabel('SFR-Weighted Metallicity')
    ax1.legend([infline[0],outline[0],zline[0]],[r'$Z_{\rm inflow,5rh}/Z_\odot$',r'$Z_{\rm outflow,5rh}/Z_\odot$',r'$Z_{\rm SFR}/Z_\odot$'],loc='upper right')
    ax1.set_ylim(1E-4,1E-2)

    plt.title(title,y=1.1)
    ax1.set_xlim(0.5,13.9)
    ax0.set_xlabel(r'$t_{\rm age}$ (Gyr)')
    ax=plt.twiny()
    ax.set_xticks(astrocosmo.age(zg[::-1]).value)
    ax.set_xticklabels(zglabels[::-1])
    ax.set_xlim(0.5,13.9)
    ax.set_xlabel(r'$z$')

    plt.savefig(savename)


def flowplot1stack(x,y,savename,labels,colors,lstyles,ylabel,ylog=False,ylim=[]):
    ax0=plt.axes()

    for i in range(len(x)):
        plt.plot(x[i],y[i],label=labels[i],color=colors[i],ls=lstyles[i])
    
    ax0.set_ylabel(ylabel)
    if ylog:
        ax0.set_yscale('log')

    if ylim!=[]:
        ax0.set_ylim(ylim)
    ax0.legend()
    ax0.set_xlim(0.5,13.9)
    ax0.set_xlabel(r'$t_{\rm age}$ (Gyr)')
    ax=plt.twiny()
    ax.set_xticks(astrocosmo.age(zg[::-1]).value)
    ax.set_xticklabels(zglabels[::-1])
    ax.set_xlim(0.5,13.9)
    ax.set_xlabel(r'$z$')
    ax.set_title("TNG100",y=.92)
    
    plt.savefig(savename)

    
for i in range(len(wcat_x)):
    plt.clf()
    x=ascii.read('flowinfo/flowinfo_100_1_'+str(int(obj[wcat_x[i],1]))+'.txt')

    deltat=np.append(astrocosmo.age(x['redshift'])[1],astrocosmo.age(x['redshift'])[1:]-astrocosmo.age(x['redshift'])[0:-1])*1E9
    deltat=deltat.value

    allzsx[i,x['snapshot']]=x['redshift']
    allflowinfox[i,x['snapshot']]=x['minflow_5rh']/deltat
    allsfrx[i,x['snapshot']]=x['sfr']
    alletax[i,x['snapshot']]=x['moutflow_5rh']/deltat/x['sfr']
    alltdx[i,x['snapshot']]=x['m_dens_gas']/(x['minflow_5rh']/deltat+x['sfr'])

    #alldinfox[i,x['snapshot']]=x['minflow_5rh']/deltat
    alloflowinfox[i,x['snapshot']]=x['moutflow_5rh']/deltat
    allzoflowinfox[i,x['snapshot']]=x['mzoutflow_5rh']/deltat
    allvflowinfox[i,x['snapshot']]=x['minflow_vir']/deltat
    alldflowinfox[i,x['snapshot']]=x['minflow_vir']/deltat-x['minflow_5rh']/deltat
    allzinfox[i,x['snapshot']]=x['mzinflow_5rh']/x['minflow_5rh']
    allzinfovx[i,x['snapshot']]=x['mzinflow_vir']/x['minflow_vir']
    allmgx[i,x['snapshot']]=x['m_gas_rh']
    allrhx[i,x['snapshot']]=x['rhstar']
    allmetx[i,x['snapshot']]=x['sfr_weight_metallicity']

    if doiplots:
        flowplot1(x,'flowinfo/flowplot_100_xmd_unnorm_'+str(int(obj[wcat_x[i],1]))+'.png','TNG100 XMD: %d' % int(obj[i,1]))
        flowplot2(x,'flowinfo/zflowplot_100_xmd_unnorm_'+str(int(obj[wcat_x[i],1]))+'.png','TNG100 XMD: %d' % int(obj[i,1]))

    
for i in range(len(wcat_l[0:150])):
    plt.clf()
    if os.path.exists('flowinfo/flowinfo_100_2_'+str(int(obj[wcat_l[i],1]))+'.txt'):
        x=ascii.read('flowinfo/flowinfo_100_2_'+str(int(obj[wcat_l[i],1]))+'.txt')
    else:
        continue
    deltat=np.append(astrocosmo.age(x['redshift'])[1],astrocosmo.age(x['redshift'])[1:]-astrocosmo.age(x['redshift'])[0:-1])*1E9
    deltat=deltat.value

    allzsl[i,x['snapshot']]=x['redshift']
    allflowinfol[i,x['snapshot']]=x['minflow_5rh']/deltat
    allsfrl[i,x['snapshot']]=x['sfr']
    alletal[i,x['snapshot']]=x['moutflow_5rh']/deltat/x['sfr']
    alltdl[i,x['snapshot']]=x['m_dens_gas']/(x['minflow_5rh']/deltat+x['sfr'])
    alloflowinfol[i,x['snapshot']]=x['moutflow_5rh']/deltat
    allzoflowinfol[i,x['snapshot']]=x['mzoutflow_5rh']/deltat
    allvflowinfol[i,x['snapshot']]=x['minflow_vir']/deltat
    alldflowinfol[i,x['snapshot']]=x['minflow_vir']/deltat-x['minflow_5rh']/deltat
    allzinfol[i,x['snapshot']]=x['mzinflow_5rh']/x['minflow_5rh']
    allzinfovl[i,x['snapshot']]=x['mzinflow_vir']/x['minflow_vir']
    allmgl[i,x['snapshot']]=x['m_gas_rh']
    allrhl[i,x['snapshot']]=x['rhstar']
    allmetl[i,x['snapshot']]=x['sfr_weight_metallicity']


    if doiplots:
        flowplot1(x,'flowinfo/flowplot_100_lmd_unnorm_'+str(int(obj[wcat_l[i],1]))+'.png','TNG100 LMD: %d' % int(obj[i,1]))
        flowplot2(x,'flowinfo/zflowplot_100_lmd_unnorm_'+str(int(obj[wcat_l[i],1]))+'.png','TNG100 LMD: %d' % int(obj[i,1]))

for i in range(len(wcat_xa)):
    plt.clf()
    if os.path.exists('flowinfo/flowinfo_100_3_'+str(int(obj[wcat_xa[i],1]))+'.txt'):
        x=ascii.read('flowinfo/flowinfo_100_3_'+str(int(obj[wcat_xa[i],1]))+'.txt')
    else:
        continue
    deltat=np.append(astrocosmo.age(x['redshift'])[1],astrocosmo.age(x['redshift'])[1:]-astrocosmo.age(x['redshift'])[0:-1])*1E9
    deltat=deltat.value

    allzsxa[i,x['snapshot']]=x['redshift']
    allflowinfoxa[i,x['snapshot']]=x['minflow_5rh']/deltat
    alloflowinfoxa[i,x['snapshot']]=x['moutflow_5rh']/deltat
    allzoflowinfoxa[i,x['snapshot']]=x['mzoutflow_5rh']/deltat
    allvflowinfoxa[i,x['snapshot']]=x['minflow_vir']/deltat
    alldflowinfoxa[i,x['snapshot']]=x['minflow_vir']/deltat-x['minflow_5rh']/deltat
    
    allzinfovxa[i,x['snapshot']]=x['mzinflow_vir']/x['minflow_vir']
    allzinfoxa[i,x['snapshot']]=x['mzinflow_5rh']/x['minflow_5rh']
    allsfrxa[i,x['snapshot']]=x['sfr']
    alletaxa[i,x['snapshot']]=x['moutflow_5rh']/deltat/x['sfr']
    alltdxa[i,x['snapshot']]=x['m_dens_gas']/(x['minflow_5rh']/deltat+x['sfr'])
    allmgxa[i,x['snapshot']]=x['m_gas_rh']
    allrhxa[i,x['snapshot']]=x['rhstar']
    allmetxa[i,x['snapshot']]=x['sfr_weight_metallicity']

    if doiplots:
        flowplot1(x,'flowinfo/flowplot_100_nmd_unnorm_'+str(int(obj[wcat_xa[i],1]))+'.png','TNG100 XMD Analog: %d' % int(obj[i,1]))
        flowplot2(x,'flowinfo/zflowplot_100_nmd_unnorm_'+str(int(obj[wcat_xa[i],1]))+'.png','TNG100 XMD Analog: %d' % int(obj[i,1]))

for i in range(len(wcat_la[0:200])):
    plt.clf()
    if os.path.exists('flowinfo/flowinfo_100_4_'+str(int(obj[wcat_la[i],1]))+'.txt'):
        x=ascii.read('flowinfo/flowinfo_100_4_'+str(int(obj[wcat_la[i],1]))+'.txt')
    else:
        continue
    deltat=np.append(astrocosmo.age(x['redshift'])[1],astrocosmo.age(x['redshift'])[1:]-astrocosmo.age(x['redshift'])[0:-1])*1E9
    deltat=deltat.value
    
    allzsla[i,x['snapshot']]=x['redshift']
    allsfrla[i,x['snapshot']]=x['sfr']
    alletala[i,x['snapshot']]=x['moutflow_5rh']/deltat/x['sfr']
    alltdla[i,x['snapshot']]=x['m_dens_gas']/(x['minflow_5rh']/deltat+x['sfr'])

    allflowinfola[i,x['snapshot']]=x['minflow_5rh']/deltat
    alloflowinfola[i,x['snapshot']]=x['moutflow_5rh']/deltat
    allzoflowinfola[i,x['snapshot']]=x['mzoutflow_5rh']/deltat
    allvflowinfola[i,x['snapshot']]=x['minflow_vir']/deltat
    alldflowinfola[i,x['snapshot']]=x['minflow_vir']/deltat-x['minflow_5rh']/deltat
    allzinfovla[i,x['snapshot']]=x['mzinflow_vir']/x['minflow_vir']
    allzinfola[i,x['snapshot']]=x['mzinflow_5rh']/x['minflow_5rh']
    allmgla[i,x['snapshot']]=x['m_gas_rh']
    allrhla[i,x['snapshot']]=x['rhstar']
    allmetla[i,x['snapshot']]=x['sfr_weight_metallicity']


    if doiplots:
        flowplot1(x,'flowinfo/flowplot_100_lmda_unnorm_'+str(int(obj[wcat_la[i],1]))+'.png','TNG100 LMD Analog: %d' % int(obj[i,1]))
        flowplot2(x,'flowinfo/zflowplot_100_lmda_unnorm_'+str(int(obj[wcat_la[i],1]))+'.png','TNG100 LMD Analog: %d' % int(obj[i,1]))

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allzoflowinfola,axis=0),np.nanmedian(allzoflowinfoxa,axis=0),np.nanmedian(allzoflowinfol,axis=0),np.nanmedian(allzoflowinfox,axis=0)],'zoflowstack_unnorm_5rh_100.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.',':','--','-'],r'Mass of outflowing metals (${\rm M_\odot~yr^{-1}}$)',ylim=[0,2E-3])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allzoflowinfola/(allmetla*allmgla),axis=0),np.nanmedian(allzoflowinfoxa/(allmetxa*allmgxa),axis=0),np.nanmedian(allzoflowinfol/(allmetl*allmgl),axis=0),np.nanmedian(allzoflowinfox/(allmetx*allmgx),axis=0)],'zfracoflowstack_unnorm_5rh_100.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.',':','--','-'],r'frac of outflowing metals (${\rm M_\odot~yr^{-1}}$)')


plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allflowinfola,axis=0),np.nanmedian(allflowinfoxa,axis=0),np.nanmedian(allflowinfol,axis=0),np.nanmedian(allflowinfox,axis=0)],'flowstack_unnorm_5rh_100.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.',':','--','-'],r'Inflow Rate within $5 r_h$ (${\rm M_\odot~yr^{-1}}$)',ylim=[0,2])

allmet=np.hstack([allmetla[:,-1],allmetxa[:,-1],allmetl[:,-1],allmetx[:,-1]])
allinflow=np.hstack([allflowinfola[:,-1],allflowinfoxa[:,-1],allflowinfol[:,-1],allflowinfox[:,-1]])
allavginflow=np.hstack([np.nanmean(allflowinfola,axis=1),np.nanmean(allflowinfoxa,axis=1),np.nanmean(allflowinfol,axis=1),np.nanmean(allflowinfox,axis=1)])
allmaxinflow=np.hstack([np.nanmax(allflowinfola,axis=1),np.nanmax(allflowinfoxa,axis=1),np.nanmax(allflowinfol,axis=1),np.nanmax(allflowinfox,axis=1)])

alloutflow=np.hstack([alloflowinfola[:,-1],alloflowinfoxa[:,-1],alloflowinfol[:,-1],alloflowinfox[:,-1]])
allavgoutflow=np.hstack([np.nanmean(alloflowinfola,axis=1),np.nanmean(alloflowinfoxa,axis=1),np.nanmean(alloflowinfol,axis=1),np.nanmean(alloflowinfox,axis=1)])
allmaxoutflow=np.hstack([np.nanmax(alloflowinfola,axis=1),np.nanmax(alloflowinfoxa,axis=1),np.nanmax(alloflowinfol,axis=1),np.nanmax(alloflowinfox,axis=1)])

alleta=np.hstack([alletala[:,-1],alletaxa[:,-1],alletal[:,-1],alletax[:,-1]])
allavgeta=np.hstack([np.nanmean(alletala,axis=1),np.nanmean(alletaxa,axis=1),np.nanmean(alletal,axis=1),np.nanmean(alletax,axis=1)])
allmaxeta=np.hstack([np.nanmax(alletala,axis=1),np.nanmax(alletaxa,axis=1),np.nanmax(alletal,axis=1),np.nanmax(alletax,axis=1)])

wok=np.where(np.isfinite(allmet))[0]
print('cov z0 inflow',np.cov(np.array([np.log10(allmet[wok]),np.log10(allinflow[wok])])))
print('cov avg inflow',np.cov(np.array([np.log10(allmet[wok]),np.log10(allavginflow[wok])])))
print('cov max inflow',np.cov(np.array([np.log10(allmet[wok]),np.log10(allmaxinflow[wok])])))

print('cov z0 outflow',np.cov(np.array([np.log10(allmet[wok]),np.log10(alloutflow[wok])])))
print('cov avg outflow',np.cov(np.array([np.log10(allmet[wok]),np.log10(allavgoutflow[wok])])))
print('cov max outflow',np.cov(np.array([np.log10(allmet[wok]),np.log10(allmaxoutflow[wok])])))

#print(allmet,allinflow)
#print('varmet',np.var(allmet))
#print('covinf',np.cov(np.array([allmet,allinflow])))

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(alletala,axis=0),np.nanmedian(alletaxa,axis=0),np.nanmedian(alletal,axis=0),np.nanmedian(alletax,axis=0)],'etastack_unnorm_5rh_100.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.',':','--','-'],r'$\eta$',ylim=[0,500])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allflowinfola/allsfrla,axis=0),np.nanmedian(allflowinfoxa/allsfrxa,axis=0),np.nanmedian(allflowinfol/allsfrl,axis=0),np.nanmedian(allflowinfox/allsfrx,axis=0)],'infsfrstack_unnorm_5rh_100.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.',':','--','-'],r'Inflow over SFR',ylim=[0,2000])


allinfsfr=np.hstack([np.sum(allflowinfola[:,-10:],axis=1)/np.sum(allsfrla[:,-10:],axis=1),np.sum(allflowinfoxa[:,-10:],axis=1)/np.sum(allsfrxa[:,-10:],axis=1),np.sum(allflowinfol[:,-10:],axis=1)/np.sum(allsfrl[:,-10:],axis=1),np.sum(allflowinfox[:,-10:],axis=1)/np.sum(allsfrx[:,-10:],axis=1)])
allinfsfr2=np.hstack([np.sum(allflowinfola[:,-5:],axis=1)/np.sum(allsfrla[:,-5:],axis=1),np.sum(allflowinfoxa[:,-5:],axis=1)/np.sum(allsfrxa[:,-5:],axis=1),np.sum(allflowinfol[:,-5:],axis=1)/np.sum(allsfrl[:,-5:],axis=1),np.sum(allflowinfox[:,-5:],axis=1)/np.sum(allsfrx[:,-5:],axis=1)])
allavginfsfr=np.hstack([np.nanmedian(allflowinfola/allsfrla,axis=1),np.nanmedian(allflowinfoxa/allsfrxa,axis=1),np.nanmedian(allflowinfol/allsfrl,axis=1),np.nanmedian(allflowinfox/allsfrx,axis=1)])
allmmaxinfsfr=np.hstack([np.nanmax(allflowinfola/allsfrla,axis=1),np.nanmax(allflowinfoxa/allsfrxa,axis=1),np.nanmax(allflowinfol/allsfrl,axis=1),np.nanmax(allflowinfox/allsfrx,axis=1)])

allrh=np.hstack([allrhla[:,-1],allrhxa[:,-1],allrhl[:,-1],allrhx[:,-1]])
allmg=np.hstack([allmgla[:,-1],allmgxa[:,-1],allmgl[:,-1],allmgx[:,-1]])
allz0=np.hstack([allzinfola[:,-1],allzinfoxa[:,-1],allzinfol[:,-1],allzinfox[:,-1]])
allz0v=np.hstack([allzinfovla[:,-1],allzinfovxa[:,-1],allzinfovl[:,-1],allzinfovx[:,-1]])

plt.clf()
plt.plot(np.hstack([env[wenv_la,3],env[wenv_xa,3],env[wenv_l,3],env[wenv_x,3]]),allz0,'o')
plt.yscale('log')
plt.savefig('envz0.png')
plt.clf()

#plt.scatter(allinfsfr,allmet,c=np.log10(alloutflow),alpha=.5,vmin=-.5,vmax=1.5)
#plt.scatter(allinfsfr,allmet,c=np.log10(alleta),alpha=.5,vmin=2,vmax=3.5)
#plt.scatter(allinfsfr,allmet,c=np.log10(allmg),alpha=.5,vmin=-3,vmax=-.5)
#plt.scatter(allinfsfr,allmet,c=allrh,alpha=.5,vmin=1,vmax=5)
#plt.scatter(allinfsfr,allmet,c=np.log10(allinfsfr2),alpha=.5,vmin=1,vmax=4)
plt.scatter(allinfsfr,allmet,c=np.log10(allz0),alpha=.5,vmin=-4,vmax=-2)
plt.loglog()
plt.savefig('infsfrmet.png')

plt.clf()
plt.scatter(allz0,allmet,c=np.log10(allinfsfr),alpha=.5,vmin=1,vmax=4)
plt.loglog()
plt.ylim(1E-3,2E-2)
plt.savefig('infmetinfsfr.png')

plt.clf()
plt.scatter(allz0v,allmet,c=np.log10(allinfsfr),alpha=.5,vmin=1,vmax=4)
plt.loglog()
plt.ylim(1E-3,2E-2)
plt.savefig('infmetvinfsfr.png')

#print('e')
#plt.show()
#plt.clf()
#fig = plt.figure()
#ax=fig.add_subplot(111,projection='3d')
#ax.scatter(np.log10(allz0),np.log10(allinfsfr),np.log10(allmet),c=np.log10(allinfsfr))
#ax.set_xlabel('z0')
#ax.set_ylabel('infsfr')
#ax.set_zlabel('met')
#plt.show()

wok=np.where((np.isfinite(allz0v)) & (np.isfinite(np.log10(allmet))) & (np.isfinite(np.log10(allinfsfr))))[0]
print(np.log10(allz0v[wok]),np.log10(allmet[wok]))
print('dallout',getdeterminationcoef.getdeterminationcoef(np.log10(alloutflow[wok]),np.log10(allmet[wok]),deg=1))
print('davgout',getdeterminationcoef.getdeterminationcoef(np.log10(allavgoutflow[wok]),np.log10(allmet[wok]),deg=1))
print('dmaxout',getdeterminationcoef.getdeterminationcoef(np.log10(allmaxoutflow[wok]),np.log10(allmet[wok]),deg=1))

print('diallout',getdeterminationcoef.getdeterminationcoef(np.log10(allinflow[wok]),np.log10(allmet[wok]),deg=1))
print('diavgout',getdeterminationcoef.getdeterminationcoef(np.log10(allavginflow[wok]),np.log10(allmet[wok]),deg=1))
print('dimaxout',getdeterminationcoef.getdeterminationcoef(np.log10(allmaxinflow[wok]),np.log10(allmet[wok]),deg=1))
print('d',getdeterminationcoef.getdeterminationcoef(np.log10(allz0v[wok]),np.log10(allmet[wok]),deg=1))
print('d',getdeterminationcoef.getdeterminationcoef(np.log10(allinfsfr[wok]),np.log10(allmet[wok]),deg=1))
plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(alltdla,axis=0)/1E9,np.nanmedian(alltdxa,axis=0)/1E9,np.nanmedian(alltdl,axis=0)/1E9,np.nanmedian(alltdx,axis=0)/1E9],'tdstack_unnorm_5rh_100.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.',':','--','-'],r'$t_d$',ylim=[0,10])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(alloflowinfola,axis=0),np.nanmedian(alloflowinfoxa,axis=0),np.nanmedian(alloflowinfol,axis=0),np.nanmedian(alloflowinfox,axis=0)],'oflowstack_unnorm_5rh_100.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.',':','--','-'],r'Outflow Rate within $5 r_h$ (${\rm M_\odot~yr^{-1}}$)',ylim=[0,2])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(alloflowinfola/allflowinfola,axis=0),np.nanmedian(alloflowinfoxa/allflowinfoxa,axis=0),np.nanmedian(alloflowinfol/allflowinfol,axis=0),np.nanmedian(alloflowinfox/allflowinfox,axis=0)],'ratiostack_unnorm_5rh_100.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.',':','--','-'],r'Outflow Rate/Inflow Rate',ylim=[0,2])
 
plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allzinfola,axis=0),np.nanmedian(allzinfoxa,axis=0),np.nanmedian(allzinfol,axis=0),np.nanmedian(allzinfox,axis=0)],'zflowstack_5rh_100.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.',':','--','-'],r'$Z_{\rm inflow}/Z_\odot$',ylog=True,ylim=[1E-4,1E-2])


