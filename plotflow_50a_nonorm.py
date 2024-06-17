import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from astropy import cosmology
import os
import getdeterminationcoef
astrocosmo=cosmology.FlatLambdaCDM(H0=70,Om0=0.3)

plt.style.use('../python/stylesheet.txt')
plt.rcParams['xtick.top']=False
plt.rcParams['xtick.labeltop']=False

obj=np.loadtxt('sampleselection_tng50_sfr5_notide_nomass.txt',delimiter=',',skiprows=1)
#zg=np.array([0,.1,.25,.5,.75,1,1.5,2,3,7])
zg=np.array([0,.1,.2,.3,.4,.5,.6])
zglabels=zg.astype(np.str)
agepoints=astrocosmo.age(zg)

#envinfo=ascii.read('xmd_env.txt')

wcat_x=np.where(obj[:,2]==1)[0]
wcat_xa=np.where(obj[:,2]==3)[0]
wcat_l=np.where(obj[:,2]==2)[0]
wcat_la=np.where(obj[:,2]==4)[0]


allzsxa=np.zeros([len(wcat_xa),100])+np.nan
allzsx=np.zeros([len(wcat_x),100])+np.nan
allzsl=np.zeros([len(wcat_l),100])+np.nan
allzsla=np.zeros([len(wcat_la),100])+np.nan

alldxa=np.zeros([len(wcat_xa),100])+np.nan
alldx=np.zeros([len(wcat_x),100])+np.nan
alldl=np.zeros([len(wcat_l),100])+np.nan
alldla=np.zeros([len(wcat_la),100])+np.nan

all2zinfoxa=np.zeros([len(wcat_xa),100])+np.nan
all2zinfox=np.zeros([len(wcat_x),100])+np.nan
all2zinfol=np.zeros([len(wcat_l),100])+np.nan
all2zinfola=np.zeros([len(wcat_la),100])+np.nan

allzinfoxa=np.zeros([len(wcat_xa),100])+np.nan
allzinfox=np.zeros([len(wcat_x),100])+np.nan
allzinfol=np.zeros([len(wcat_l),100])+np.nan
allzinfola=np.zeros([len(wcat_la),100])+np.nan

allvzinfoxa=np.zeros([len(wcat_xa),100])+np.nan
allvzinfox=np.zeros([len(wcat_x),100])+np.nan
allvzinfol=np.zeros([len(wcat_l),100])+np.nan
allvzinfola=np.zeros([len(wcat_la),100])+np.nan

allzoutxa=np.zeros([len(wcat_xa),100])+np.nan
allzoutx=np.zeros([len(wcat_x),100])+np.nan
allzoutl=np.zeros([len(wcat_l),100])+np.nan
allzoutla=np.zeros([len(wcat_la),100])+np.nan

allmetxa=np.zeros([len(wcat_xa),100])+np.nan
allmetx=np.zeros([len(wcat_x),100])+np.nan
allmetl=np.zeros([len(wcat_l),100])+np.nan
allmetla=np.zeros([len(wcat_la),100])+np.nan

allsfrxa=np.zeros([len(wcat_xa),100])+np.nan
allsfrx=np.zeros([len(wcat_x),100])+np.nan
allsfrl=np.zeros([len(wcat_l),100])+np.nan
allsfrla=np.zeros([len(wcat_la),100])+np.nan

allzoutfoxa=np.zeros([len(wcat_xa),100])+np.nan
allzoutfox=np.zeros([len(wcat_x),100])+np.nan
allzoutfol=np.zeros([len(wcat_l),100])+np.nan
allzoutfola=np.zeros([len(wcat_la),100])+np.nan

allzstartxa=np.zeros([len(wcat_xa),100])+np.nan
allzstartx=np.zeros([len(wcat_x),100])+np.nan
allzstartl=np.zeros([len(wcat_l),100])+np.nan
allzstartla=np.zeros([len(wcat_la),100])+np.nan

allmetoutfoxa=np.zeros([len(wcat_xa),100])+np.nan
allmetoutfox=np.zeros([len(wcat_x),100])+np.nan
allmetoutfol=np.zeros([len(wcat_l),100])+np.nan
allmetoutfola=np.zeros([len(wcat_la),100])+np.nan

allgasxa=np.zeros([len(wcat_xa),100])+np.nan
allgasx=np.zeros([len(wcat_x),100])+np.nan
allgasl=np.zeros([len(wcat_l),100])+np.nan
allgasla=np.zeros([len(wcat_la),100])+np.nan

allflowinfoxa=np.zeros([len(wcat_xa),100])+np.nan
allflowinfox=np.zeros([len(wcat_x),100])+np.nan
allflowinfol=np.zeros([len(wcat_l),100])+np.nan
allflowinfola=np.zeros([len(wcat_la),100])+np.nan

alloflowinfoxa=np.zeros([len(wcat_xa),100])+np.nan
alloflowinfox=np.zeros([len(wcat_x),100])+np.nan
alloflowinfol=np.zeros([len(wcat_l),100])+np.nan
alloflowinfola=np.zeros([len(wcat_la),100])+np.nan

allvflowinfoxa=np.zeros([len(wcat_xa),100])+np.nan
allvflowinfox=np.zeros([len(wcat_x),100])+np.nan
allvflowinfol=np.zeros([len(wcat_l),100])+np.nan
allvflowinfola=np.zeros([len(wcat_la),100])+np.nan

alldflowinfoxa=np.zeros([len(wcat_xa),100])+np.nan
alldflowinfox=np.zeros([len(wcat_x),100])+np.nan
alldflowinfol=np.zeros([len(wcat_l),100])+np.nan
alldflowinfola=np.zeros([len(wcat_la),100])+np.nan

alltdxa=np.zeros([len(wcat_xa),100])+np.nan
alltdx=np.zeros([len(wcat_x),100])+np.nan
alltdl=np.zeros([len(wcat_l),100])+np.nan
alltdla=np.zeros([len(wcat_la),100])+np.nan

allmgxa=np.zeros([len(wcat_xa),100])+np.nan
allmgx=np.zeros([len(wcat_x),100])+np.nan
allmgl=np.zeros([len(wcat_l),100])+np.nan
allmgla=np.zeros([len(wcat_la),100])+np.nan

allrhxa=np.zeros([len(wcat_xa),100])+np.nan
allrhx=np.zeros([len(wcat_x),100])+np.nan
allrhl=np.zeros([len(wcat_l),100])+np.nan
allrhla=np.zeros([len(wcat_la),100])+np.nan


indplot=False

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
    #ax0.set_xlim(0.5,13.9)
    ax0.set_xlim(7.5,13.9)
    ax0.set_xlabel(r'$t_{\rm age}$ (Gyr)')
    ax=plt.twiny()
    ax.set_xticks(astrocosmo.age(zg[::-1]).value)
    ax.set_xticklabels(zglabels[::-1])
    #ax.set_xlim(0.5,13.9)
    ax.set_xlim(7.5,13.9)
    ax.set_xlabel(r'$z$')
    ax.set_title('TNG50',y=.92)

    plt.savefig(savename)

    
for i in range(len(wcat_x)):
    plt.clf()
    x=ascii.read('flowinfo/flowinfo_50_1_'+str(int(obj[wcat_x[i],1]))+'.txt')

    deltat=np.append(astrocosmo.age(x['redshift'])[1],astrocosmo.age(x['redshift'])[1:]-astrocosmo.age(x['redshift'])[0:-1])*1E9

    allzsx[i,x['snapshot']]=x['redshift']
    allflowinfox[i,x['snapshot']]=x['minflow_5rh']/deltat
    alldx[i,x['snapshot']]=x['mendsf']/deltat
    allmetx[i,x['snapshot']]=x['sfr_weight_metallicity']
    allsfrx[i,x['snapshot']]=x['sfr']
    allgasx[i,x['snapshot']]=x['m_gas_rh']
    alloflowinfox[i,x['snapshot']]=x['moutflow_5rh']/deltat
    allmetoutfox[i,x['snapshot']]=x['mzoutflow_5rh']/deltat
    allvflowinfox[i,x['snapshot']]=x['minflow_vir']/deltat
    alldflowinfox[i,x['snapshot']]=x['minflow_vir']/deltat-x['minflow_5rh']/deltat
    etai=x['moutflow_5rh']/deltat/x['sfr']
    alltdep.append(x['moutflow_5rh']/deltat)
    allzstartx[i,x['snapshot']]=x['zstartsf']/x['mstartsf']
    allzinfox[i,x['snapshot']]=x['mzinflow_5rh']/x['minflow_5rh']/.0127
    all2zinfox[i,x['snapshot']]=x['mzinflow_5rh']/x['minflow_5rh']
    allvzinfox[i,x['snapshot']]=x['mzinflow_vir']/x['minflow_vir']/.0127
    allzoutfox[i,x['snapshot']]=x['mzoutflow_5rh']/x['moutflow_5rh']

    allmgx[i,x['snapshot']]=x['mgas_tot']
    allrhx[i,x['snapshot']]=x['rhstar']
    
    alltdx[i,x['snapshot']]=x['m_gas_rh']*1E10/.6774/(x['moutflow_5rh']/deltat.value+x['sfr'])
    
    if indplot:
        flowplot1(x,'flowinfo/flowplot_50_xmd_unnorm_'+str(int(obj[wcat_x[i],1]))+'.png','TNG50 XMD: %d' % int(obj[i,1]))
        flowplot2(x,'flowinfo/zflowplot_50_xmd_unnorm_'+str(int(obj[wcat_x[i],1]))+'.png','TNG50 XMD: %d' % int(obj[i,1]))
    
print(allzoutfox,allmetx)

for i in range(len(wcat_l)):
    plt.clf()
    if os.path.exists('flowinfo/flowinfo_50_2_'+str(int(obj[wcat_l[i],1]))+'.txt'):
        x=ascii.read('flowinfo/flowinfo_50_2_'+str(int(obj[wcat_l[i],1]))+'.txt')
        print(x['mstar'][-1])
        if x['mstar'][-1]>1E9:
            continue
    else:
        continue
    deltat=np.append(astrocosmo.age(x['redshift'])[1],astrocosmo.age(x['redshift'])[1:]-astrocosmo.age(x['redshift'])[0:-1])*1E9

    allzsl[i,x['snapshot']]=x['redshift']
    allflowinfol[i,x['snapshot']]=x['minflow_5rh']/deltat
    alldl[i,x['snapshot']]=x['mendsf']/deltat
    allsfrl[i,x['snapshot']]=x['sfr']
    allmetl[i,x['snapshot']]=x['sfr_weight_metallicity']
    allgasl[i,x['snapshot']]=x['m_gas_rh']
    allmetoutfol[i,x['snapshot']]=x['mzoutflow_5rh']/deltat
    alloflowinfol[i,x['snapshot']]=x['moutflow_5rh']/deltat
    allvflowinfol[i,x['snapshot']]=x['minflow_vir']/deltat
    alldflowinfol[i,x['snapshot']]=x['minflow_vir']/deltat-x['minflow_5rh']/deltat
    allzstartl[i,x['snapshot']]=x['zstartsf']/x['mstartsf']
    allzinfol[i,x['snapshot']]=x['mzinflow_5rh']/x['minflow_5rh']/.0127
    all2zinfol[i,x['snapshot']]=x['mzinflow_5rh']/x['minflow_5rh']
    allvzinfol[i,x['snapshot']]=x['mzinflow_vir']/x['minflow_vir']/.0127
    allzoutfol[i,x['snapshot']]=x['mzoutflow_5rh']/x['moutflow_5rh']
    alltdl[i,x['snapshot']]=x['m_gas_rh']*1E10/.6774/(x['moutflow_5rh']/deltat.value+x['sfr'])

    allmgl[i,x['snapshot']]=x['mgas_tot']
    allrhl[i,x['snapshot']]=x['rhstar']

    if indplot:
        flowplot1(x,'flowinfo/flowplot_50_lmd_unnorm_'+str(int(obj[wcat_l[i],1]))+'.png','TNG50 LMD: %d' % int(obj[i,1]))
        flowplot2(x,'flowinfo/zflowplot_50_lmd_unnorm_'+str(int(obj[wcat_l[i],1]))+'.png','TNG50 LMD: %d' % int(obj[i,1]))
    
for i in range(len(wcat_xa)):
    plt.clf()
    if os.path.exists('flowinfo/flowinfo_50_3_'+str(int(obj[wcat_xa[i],1]))+'.txt'):
        x=ascii.read('flowinfo/flowinfo_50_3_'+str(int(obj[wcat_xa[i],1]))+'.txt')
    else:
        continue
    deltat=np.append(astrocosmo.age(x['redshift'])[1],astrocosmo.age(x['redshift'])[1:]-astrocosmo.age(x['redshift'])[0:-1])*1E9

    allzsxa[i,x['snapshot']]=x['redshift']
    alldxa[i,x['snapshot']]=x['mendsf']/deltat
    allmetxa[i,x['snapshot']]=x['sfr_weight_metallicity']
    allsfrxa[i,x['snapshot']]=x['sfr']
    allgasxa[i,x['snapshot']]=x['m_gas_rh']
    allflowinfoxa[i,x['snapshot']]=x['minflow_5rh']/deltat
    allmetoutfoxa[i,x['snapshot']]=x['mzoutflow_5rh']/deltat
    alloflowinfoxa[i,x['snapshot']]=x['moutflow_5rh']/deltat
    allvflowinfoxa[i,x['snapshot']]=x['minflow_vir']/deltat
    alldflowinfoxa[i,x['snapshot']]=x['minflow_vir']/deltat-x['minflow_5rh']/deltat
    allzstartxa[i,x['snapshot']]=x['zstartsf']/x['mstartsf']
    allzinfoxa[i,x['snapshot']]=x['mzinflow_5rh']/x['minflow_5rh']/.0127
    all2zinfoxa[i,x['snapshot']]=x['mzinflow_5rh']/x['minflow_5rh']
    allvzinfoxa[i,x['snapshot']]=x['mzinflow_vir']/x['minflow_vir']/.0127
    alltdxa[i,x['snapshot']]=x['m_gas_rh']*1E10/.6774/(x['moutflow_5rh']/deltat.value+x['sfr'])
    allzoutfoxa[i,x['snapshot']]=x['mzoutflow_5rh']/x['moutflow_5rh']

    allmgxa[i,x['snapshot']]=x['mgas_tot']
    allrhxa[i,x['snapshot']]=x['rhstar']

    if indplot:
        flowplot1(x,'flowinfo/flowplot_50_nmd_unnorm_'+str(int(obj[wcat_xa[i],1]))+'.png','TNG50 XMD Analog: %d' % int(obj[i,1]))
        flowplot2(x,'flowinfo/zflowplot_50_nmd_unnorm_'+str(int(obj[wcat_xa[i],1]))+'.png','TNG50 XMD Analog: %d' % int(obj[i,1]))

for i in range(len(wcat_la)):
    plt.clf()
    if os.path.exists('flowinfo/flowinfo_50_4_'+str(int(obj[wcat_la[i],1]))+'.txt'):
        x=ascii.read('flowinfo/flowinfo_50_4_'+str(int(obj[wcat_la[i],1]))+'.txt')
        if x['mstar'][-1]>1E9:
            continue

    else:
        continue
    deltat=np.append(astrocosmo.age(x['redshift'])[1],astrocosmo.age(x['redshift'])[1:]-astrocosmo.age(x['redshift'])[0:-1])*1E9

    allzsla[i,x['snapshot']]=x['redshift']
    alldla[i,x['snapshot']]=x['mendsf']/deltat
    allsfrla[i,x['snapshot']]=x['sfr']
    allmetla[i,x['snapshot']]=x['sfr_weight_metallicity']
    allgasla[i,x['snapshot']]=x['m_gas_rh']
    allflowinfola[i,x['snapshot']]=x['minflow_5rh']/deltat
    allmetoutfola[i,x['snapshot']]=x['mzoutflow_5rh']/deltat
    alloflowinfola[i,x['snapshot']]=x['moutflow_5rh']/deltat
    allvflowinfola[i,x['snapshot']]=x['minflow_vir']/deltat
    alldflowinfola[i,x['snapshot']]=x['minflow_vir']/deltat-x['minflow_5rh']/deltat
    allzstartla[i,x['snapshot']]=x['zstartsf']/x['mstartsf']
    allzinfola[i,x['snapshot']]=x['mzinflow_5rh']/x['minflow_5rh']/.0127
    all2zinfola[i,x['snapshot']]=x['mzinflow_5rh']/x['minflow_5rh']
    allvzinfola[i,x['snapshot']]=x['mzinflow_5rh']/x['minflow_vir']/.0127
    allzoutfola[i,x['snapshot']]=x['mzoutflow_5rh']/x['moutflow_5rh']

    alltdla[i,x['snapshot']]=x['m_gas_rh']*1E10/.6774/(x['moutflow_5rh']/deltat.value+x['sfr'])

    allmgla[i,x['snapshot']]=x['mgas_tot']
    allrhla[i,x['snapshot']]=x['rhstar']

    
    if indplot:
        flowplot1(x,'flowinfo/flowplot_50_lmda_unnorm_'+str(int(obj[wcat_la[i],1]))+'.png','TNG50 LMD Analog: %d' % int(obj[i,1]))
        flowplot2(x,'flowinfo/zflowplot_50_lmda_unnorm_'+str(int(obj[wcat_la[i],1]))+'.png','TNG50 LMD Analog: %d' % int(obj[i,1]))

    
plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allflowinfola,axis=0),np.nanmedian(allflowinfoxa,axis=0),np.nanmedian(allflowinfol,axis=0),np.nanmedian(allflowinfox,axis=0)],'flowstack_unnorm_5rh.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'Inflow Rate within $5 r_h$ (${\rm M_\odot~yr^{-1}}$)',ylim=[0,1])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allvflowinfola,axis=0),np.nanmedian(allvflowinfoxa,axis=0),np.nanmedian(allvflowinfol,axis=0),np.nanmedian(allvflowinfox,axis=0)],'flowstack_unnorm_vir.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'Inflow Rate within $rv$ (${\rm M_\odot~yr^{-1}}$)',ylim=[0,1])

print(np.nanmedian(allgasxa/allmetoutfoxa,axis=0))
plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allgasla*allmetla/allmetoutfola,axis=0),np.nanmedian(allgasxa*allmetxa/allmetoutfoxa,axis=0),np.nanmedian(allgasl*allmetl/allmetoutfol,axis=0),np.nanmedian(allgasx*allmetx/allmetoutfox,axis=0)],'zej_unnorm_5rh.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'Fraction of metals ejected (Gyr$^{-1}$)',ylim=[0,.2])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allmetoutfola,axis=0),np.nanmedian(allmetoutfoxa,axis=0),np.nanmedian(allmetoutfol,axis=0),np.nanmedian(allmetoutfox,axis=0)],'mzej_unnorm_5rh.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'Rate of metal ejection (${\rm M_\odot~yr^{-1}}$)',ylim=[0,2E-3])


print(np.nanmedian(allmetoutfola,axis=0))
plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allflowinfola,axis=0),np.nanmedian(allmetoutfoxa,axis=0),np.nanmedian(allmetoutfol,axis=0),np.nanmedian(allmetoutfox,axis=0)],'metflowstack_unnorm_5rh.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'Metal Outflow Rate within $5 r_h$ (${\rm M_\odot~yr^{-1}}$)',ylim=[0,1E-3])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allmetoutfola,axis=0),np.nanmedian(allmetoutfoxa,axis=0),np.nanmedian(allmetoutfol,axis=0),np.nanmedian(allmetoutfox,axis=0)],'metflowstack_unnorm_5rh.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'Metal Outflow Rate within $5 r_h$ (${\rm M_\odot~yr^{-1}}$)',ylim=[0,1E-3])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allmetoutfola,axis=0),np.nanmedian(allmetoutfoxa,axis=0),np.nanmedian(allmetoutfol,axis=0),np.nanmedian(allmetoutfox,axis=0)],'metflowstack_unnorm_5rh.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'Metal Outflow Rate within $5 r_h$ (${\rm M_\odot~yr^{-1}}$)',ylim=[0,1E-3])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allflowinfola/allsfrla,axis=0),np.nanmedian(allflowinfoxa/allsfrxa,axis=0),np.nanmedian(allflowinfol/allsfrl,axis=0),np.nanmedian(allflowinfox/allsfrx,axis=0)],'infsfrstack_unnorm_5rh.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'Inflow over SFR',ylim=[0,1000])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(alloflowinfola/allsfrla,axis=0),np.nanmedian(alloflowinfoxa/allsfrxa,axis=0),np.nanmedian(alloflowinfol/allsfrl,axis=0),np.nanmedian(alloflowinfox/allsfrx,axis=0)],'etastack_unnorm_5rh.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'Outflow over SFR',ylim=[0,1000])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(alldla,axis=0),np.nanmedian(alldxa,axis=0),np.nanmedian(alldl,axis=0),np.nanmedian(alldx,axis=0)],'diffusestack_unnorm_5rh.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'Star formation stoppiing rate (${\rm M_\odot~yr^{-1}}$)',ylim=[0,1])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(alloflowinfola,axis=0),np.nanmedian(alloflowinfoxa,axis=0),np.nanmedian(alloflowinfol,axis=0),np.nanmedian(alloflowinfox,axis=0)],'oflowstack_unnorm_5rh.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'Outflow Rate within $5 r_h$ (${\rm M_\odot~yr^{-1}}$)',ylim=[0,1])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(alloflowinfola/allflowinfola,axis=0),np.nanmedian(alloflowinfoxa/allflowinfoxa,axis=0),np.nanmedian(alloflowinfol/allflowinfol,axis=0),np.nanmedian(alloflowinfox/allflowinfox,axis=0)],'ratiostack_unnorm_5rh.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'Outflow Rate/Inflow Rate',ylim=[0,2])
 
plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allzinfola,axis=0),np.nanmedian(allzinfoxa,axis=0),np.nanmedian(allzinfol,axis=0),np.nanmedian(allzinfox,axis=0)],'zflowstack_5rh.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'$Z_{\rm inflow,~p}/Z_\odot$',ylog=True,ylim=[1E-2,1])
print([np.nanmedian(allzinfola,axis=0),np.nanmedian(allzinfoxa,axis=0),np.nanmedian(allzinfol,axis=0),np.nanmedian(allzinfox,axis=0)])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allzstartla,axis=0),np.nanmedian(allzstartxa,axis=0),np.nanmedian(allzstartl,axis=0),np.nanmedian(allzstartx,axis=0)],'zflowstartstack.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'$Z_{\rm start sf}/Z_\odot$',ylog=True,ylim=[1E-4,1E-2])

#plt.clf()
#flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(all2zinfola,axis=0),np.nanmedian(all2zinfoxa,axis=0),np.nanmedian(all2zinfol,axis=0),np.nanmedian(all2zinfox,axis=0)],'zflowstack_5rh.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'$Z_{\rm inflow}/Z_\odot$',ylog=True,ylim=[1E-2,2])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allvzinfola,axis=0),np.nanmedian(allvzinfoxa,axis=0),np.nanmedian(allvzinfol,axis=0),np.nanmedian(allvzinfox,axis=0)],'zflowstack_vir.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':' ,'-'],r'$Z_{\rm inflow,~p}/Z_\odot$',ylog=True,ylim=[1E-2,1])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allvzinfola,axis=0),np.nanmedian(allvzinfoxa,axis=0),np.nanmedian(allvzinfol,axis=0),np.nanmedian(allvzinfox,axis=0)],'zflowstack_out.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'$Z_{\rm outer}/Z_\odot$',ylog=True,ylim=[1E-4,1E-2])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allzoutfola,axis=0),np.nanmedian(allzoutfoxa,axis=0),np.nanmedian(allzoutfol,axis=0),np.nanmedian(allzoutfox,axis=0)],'zoutflowstack_5rh.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'$Z_{\rm outflow}/Z_\odot$',ylog=True,ylim=[1E-4,1E-2])


plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allzoutfola/allmetla,axis=0),np.nanmedian(allzoutfoxa/allmetxa,axis=0),np.nanmedian(allzoutfol/allmetl,axis=0),np.nanmedian(allzoutfox/allmetx,axis=0)],'zcompstack_5rh.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.','--',':','-'],r'$Z_{\rm outflow}/Z_{\rm SFR}$',ylog=True,ylim=[1E-1,10])

allmet=np.hstack([allmetla[:,-1],allmetxa[:,-1],allmetl[:,-1],allmetx[:,-1]])
allinflow=np.hstack([allflowinfola[:,-1],allflowinfoxa[:,-1],allflowinfol[:,-1],allflowinfox[:,-1]])
allavginflow=np.hstack([np.nanmean(allflowinfola,axis=1),np.nanmean(allflowinfoxa,axis=1),np.nanmean(allflowinfol,axis=1),np.nanmean(allflowinfox,axis=1)])
allmaxinflow=np.hstack([np.nanmax(allflowinfola,axis=1),np.nanmax(allflowinfoxa,axis=1),np.nanmax(allflowinfol,axis=1),np.nanmax(allflowinfox,axis=1)])

alloutflow=np.hstack([alloflowinfola[:,-1],alloflowinfoxa[:,-1],alloflowinfol[:,-1],alloflowinfox[:,-1]])
allavgoutflow=np.hstack([np.nanmean(alloflowinfola,axis=1),np.nanmean(alloflowinfoxa,axis=1),np.nanmean(alloflowinfol,axis=1),np.nanmean(alloflowinfox,axis=1)])
allmaxoutflow=np.hstack([np.nanmax(alloflowinfola,axis=1),np.nanmax(alloflowinfoxa,axis=1),np.nanmax(alloflowinfol,axis=1),np.nanmax(alloflowinfox,axis=1)])

alletala=alloflowinfola/allsfrla
alletaxa=alloflowinfoxa/allsfrxa
alletal=alloflowinfol/allsfrl
alletax=alloflowinfox/allsfrx

alleta=np.hstack([alletala[:,-1],alletaxa[:,-1],alletal[:,-1],alletax[:,-1]])
allavgeta=np.hstack([np.nanmean(alletala,axis=1),np.nanmean(alletaxa,axis=1),np.nanmean(alletal,axis=1),np.nanmean(alletax,axis=1)])
allmaxeta=np.hstack([np.nanmax(alletala,axis=1),np.nanmax(alletaxa,axis=1),np.nanmax(alletal,axis=1),np.nanmax(alletax,axis=1)])

allmet=np.hstack([allmetla[:,-1],allmetxa[:,-1],allmetl[:,-1],allmetx[:,-1]])
allz0=np.hstack([allzinfola[:,-1],allzinfoxa[:,-1],allzinfol[:,-1],allzinfox[:,-1]])
wok=np.where(np.isfinite(allmet))[0]
print('cov z0 inflow',np.cov(np.array([np.log10(allmet[wok]),np.log10(allinflow[wok])])))
print('cov avg inflow',np.cov(np.array([np.log10(allmet[wok]),np.log10(allavginflow[wok])])))
print('cov max inflow',np.cov(np.array([np.log10(allmet[wok]),np.log10(allmaxinflow[wok])])))

print('cov z0 outflow',np.cov(np.array([np.log10(allmet[wok]),np.log10(alloutflow[wok])])))
print('cov avg outflow',np.cov(np.array([np.log10(allmet[wok]),np.log10(allavgoutflow[wok])])))
print('cov max outflow',np.cov(np.array([np.log10(allmet[wok]),np.log10(allmaxoutflow[wok])])))

allinfsfr=np.hstack([np.nansum(allflowinfola[:,-5:],axis=1)/np.nansum(allsfrla[:,-5:],axis=1),np.nansum(allflowinfoxa[:,-5:],axis=1)/np.nansum(allsfrxa[:,-5:],axis=1),np.nansum(allflowinfol[:,-5:],axis=1)/np.nansum(allsfrl[:,-5:],axis=1),np.nansum(allflowinfox[:,-5:],axis=1)/np.nansum(allsfrx[:,-5:],axis=1)])
#print(np.nansum(allflowinfola[:,-5:],axis=1),np.nansum(allsfrla[:,-5:],axis=1))
allavginfsfr=np.hstack([np.nanmedian(alloflowinfola/allsfrla,axis=1),np.nanmedian(alloflowinfoxa/allsfrxa,axis=1),np.nanmedian(alloflowinfol/allsfrl,axis=1),np.nanmedian(alloflowinfox/allsfrx,axis=1)])
allmmaxinfsfr=np.hstack([np.nanmax(alloflowinfola/allsfrla,axis=1),np.nanmax(alloflowinfoxa/allsfrxa,axis=1),np.nanmax(alloflowinfol/allsfrl,axis=1),np.nanmax(alloflowinfox/allsfrx,axis=1)])

plt.clf()
plt.plot(np.nansum(allflowinfola[:,-5:],axis=1)/np.nansum(allsfrla[:,-5:],axis=1),allmetla[:,-1],'C3h',label='LMD Analogs')
plt.plot(np.nansum(allflowinfoxa[:,-5:],axis=1)/np.nansum(allsfrxa[:,-5:],axis=1),allmetxa[:,-1],'C2s',label='XMD Analogs')
plt.plot(np.nansum(allflowinfol[:,-5:],axis=1)/np.nansum(allsfrl[:,-5:],axis=1),allmetl[:,-1],'C1x',label='LMDs')
plt.plot(np.nansum(allflowinfox[:,-5:],axis=1)/np.nansum(allsfrx[:,-5:],axis=1),allmetx[:,-1],'C0*',label='XMDs')
plt.loglog()
plt.legend()
plt.xlabel(r'$M_{\rm infall}/$SFR')
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.title('TNG50',y=0.92)
plt.savefig('infsfrmet50.png')


allmg=np.hstack([allmgla[:,-1],allmgxa[:,-1],allmgl[:,-1],allmgx[:,-1]])
allmgr=np.hstack([allgasla[:,-1],allgasxa[:,-1],allgasl[:,-1],allgasx[:,-1]])
allinfsfr=np.hstack([np.nansum(allflowinfola[:,-5:],axis=1)/np.nansum(allsfrla[:,-5:],axis=1),np.nansum(allflowinfoxa[:,-5:],axis=1)/np.nansum(allsfrxa[:,-5:],axis=1),np.nansum(allflowinfol[:,-5:],axis=1)/np.nansum(allsfrl[:,-5:],axis=1),np.nansum(allflowinfox[:,-5:],axis=1)/np.nansum(allsfrx[:,-5:],axis=1)])
allinffrac=np.hstack([np.nansum(allflowinfola[:,-10:],axis=1)/allmgla[:,-10],np.nansum(allflowinfoxa[:,-10:],axis=1)/allmgxa[:,-10],np.nansum(allflowinfol[:,-10:],axis=1)/allmgl[:,-10],np.nansum(allflowinfox[:,-10:],axis=1)/allmgx[:,-10]])
print(getdeterminationcoef.getdeterminationcoef(np.log10(allinfsfr),np.log10(allmet),deg=1))
#plt.clf()
#plt.scatter(allinfsfr,allmet,c=np.log10(allmg),alpha=.5)
#plt.scatter(np.sum(allflowinfox[:,-10:],axis=1)/np.sum(allsfrx[:,-10:],axis=1),allmetx[:,-1],c=np.log10(alletax[:,-1]),alpha=.5)
#plt.scatter(allinfsfr,allmet,c=np.log10(allz0),alpha=.5,vmin=-2,vmax=0)
#plt.loglog()
#plt.ylim(1E-3,2E-2)
#plt.savefig('infsfrmet50.png')

#plt.clf()
#plt.scatter(allz0,allmet,c=np.log10(allinfsfr),alpha=.5,vmin=1,vmax=4)
#plt.loglog()
#plt.ylim(1E-3,2E-2)
#plt.savefig('infmetinfsfr50.png')

plt.clf()
plt.scatter(allinffrac,allmet,c=np.log10(allmg),alpha=.5)
plt.loglog()
plt.savefig('inffracmet50.png')
