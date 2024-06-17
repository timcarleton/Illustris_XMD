import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from astropy import cosmology
import os

astrocosmo = cosmology.FlatLambdaCDM(H0=67.74,Om0=0.307)

plt.style.use('../python/stylesheet.txt')
plt.rcParams['xtick.top']=False
plt.rcParams['xtick.labeltop']=False

obj=np.loadtxt('sampleselection_tng50_sfr5.txt',delimiter=',',skiprows=1)
zg=np.array([0,.1,.25,.5,.75,1,1.5,2,3,7])
zglabels=zg.astype(np.str)
agepoints=astrocosmo.age(zg)

#envinfo=ascii.read('xmd_env.txt')

allzsxa=np.zeros([10,100])+np.nan
allzsx=np.zeros([6,100])+np.nan
allzsl=np.zeros([20,100])+np.nan
allzsla=np.zeros([20,100])+np.nan

alldxa=np.zeros([10,100])+np.nan
alldx=np.zeros([6,100])+np.nan
alldl=np.zeros([20,100])+np.nan
alldla=np.zeros([20,100])+np.nan

allzinfoxa=np.zeros([10,100])+np.nan
allzinfox=np.zeros([6,100])+np.nan
allzinfol=np.zeros([20,100])+np.nan
allzinfola=np.zeros([20,100])+np.nan

allz2infoxa=np.zeros([10,100])+np.nan
allz2infox=np.zeros([6,100])+np.nan
allz2infol=np.zeros([20,100])+np.nan
allz2infola=np.zeros([20,100])+np.nan

allvzinfoxa=np.zeros([10,100])+np.nan
allvzinfox=np.zeros([6,100])+np.nan
allvzinfol=np.zeros([20,100])+np.nan
allvzinfola=np.zeros([20,100])+np.nan

allmetxa=np.zeros([10,100])+np.nan
allmetx=np.zeros([6,100])+np.nan
allmetl=np.zeros([20,100])+np.nan
allmetla=np.zeros([20,100])+np.nan

allsfrxa=np.zeros([10,100])+np.nan
allsfrx=np.zeros([6,100])+np.nan
allsfrl=np.zeros([20,100])+np.nan
allsfrla=np.zeros([20,100])+np.nan

allzoutfoxa=np.zeros([10,100])+np.nan
allzoutfox=np.zeros([6,100])+np.nan
allzoutfol=np.zeros([20,100])+np.nan
allzoutfola=np.zeros([20,100])+np.nan

allmetoutfoxa=np.zeros([10,100])+np.nan
allmetoutfox=np.zeros([6,100])+np.nan
allmetoutfol=np.zeros([20,100])+np.nan
allmetoutfola=np.zeros([20,100])+np.nan

allgasxa=np.zeros([10,100])+np.nan
allgasx=np.zeros([6,100])+np.nan
allgasl=np.zeros([20,100])+np.nan
allgasla=np.zeros([20,100])+np.nan

allflowinfoxa=np.zeros([10,100])+np.nan
allflowinfox=np.zeros([6,100])+np.nan
allflowinfol=np.zeros([20,100])+np.nan
allflowinfola=np.zeros([20,100])+np.nan

alloflowinfoxa=np.zeros([10,100])+np.nan
alloflowinfox=np.zeros([6,100])+np.nan
alloflowinfol=np.zeros([20,100])+np.nan
alloflowinfola=np.zeros([20,100])+np.nan

allvflowinfoxa=np.zeros([10,100])+np.nan
allvflowinfox=np.zeros([6,100])+np.nan
allvflowinfol=np.zeros([20,100])+np.nan
allvflowinfola=np.zeros([20,100])+np.nan

alldflowinfoxa=np.zeros([10,100])+np.nan
alldflowinfox=np.zeros([6,100])+np.nan
alldflowinfol=np.zeros([20,100])+np.nan
alldflowinfola=np.zeros([20,100])+np.nan

alltdxa=np.zeros([10,100])+np.nan
alltdx=np.zeros([6,100])+np.nan
alltdl=np.zeros([20,100])+np.nan
alltdla=np.zeros([20,100])+np.nan

wcat_x=np.where(obj[:,2]==1)[0]
wcat_xa=np.where(obj[:,2]==3)[0]
wcat_l=np.where(obj[:,2]==2)[0]
wcat_la=np.where(obj[:,2]==4)[0]

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
    ax0.set_xlim(0.5,13.9)
    ax0.set_xlabel(r'$t_{\rm age}$ (Gyr)')
    ax=plt.twiny()
    ax.set_xticks(astrocosmo.age(zg[::-1]).value)
    ax.set_xticklabels(zglabels[::-1])
    ax.set_xlim(0.5,13.9)
    ax.set_xlabel(r'$z$')

    plt.savefig(savename)

    
for i in range(len(wcat_x)):
    plt.clf()
    x=ascii.read('flowinfo/flowinfo_cos_'+str(int(obj[wcat_x[i],1]))+'.txt')

    deltat=np.append(astrocosmo.age(x['redshift'])[1],astrocosmo.age(x['redshift'])[1:]-astrocosmo.age(x['redshift'])[0:-1])*1E9

    allzsx[i,x['snapshot']]=x['redshift']
    allflowinfox[i,x['snapshot']]=x['mcosinflow_5rh']/deltat
    allsfrx[i,x['snapshot']]=x['sfr']

    allvflowinfox[i,x['snapshot']]=x['mcosinflow_vir']/deltat
    alldflowinfox[i,x['snapshot']]=x['mcosinflow_vir']/deltat-x['mcosinflow_5rh']/deltat

    allzinfox[i,x['snapshot']]=x['mzcosinflow_5rh']/x['mcosinflow_5rh']
    allz2infox[i,x['snapshot']]=x['mzcosinflow_2rh']/x['mcosinflow_2rh']
    allvzinfox[i,x['snapshot']]=x['mzcosinflow_vir']/x['mcosinflow_5rh']
    
    if indplot:
        flowplot1(x,'flowinfo/flowplot_cos_unnorm_'+str(int(obj[wcat_x[i],1]))+'.png','TNG50 XMD: %d' % int(obj[i,1]))
        flowplot2(x,'flowinfo/zflowplot_cos_unnorm_'+str(int(obj[wcat_x[i],1]))+'.png','TNG50 XMD: %d' % int(obj[i,1]))
    
print(allzoutfox,allmetx)

for i in range(len(wcat_l[0:20])):
    plt.clf()
    if os.path.exists('flowinfo/flowinfo_cos_'+str(int(obj[wcat_l[i],1]))+'.txt'):
        x=ascii.read('flowinfo/flowinfo_cos_'+str(int(obj[wcat_l[i],1]))+'.txt')
    else:
        continue
    deltat=np.append(astrocosmo.age(x['redshift'])[1],astrocosmo.age(x['redshift'])[1:]-astrocosmo.age(x['redshift'])[0:-1])*1E9

    allzsl[i,x['snapshot']]=x['redshift']
    allflowinfol[i,x['snapshot']]=x['mcosinflow_5rh']/deltat
    allsfrl[i,x['snapshot']]=x['sfr']

    allvflowinfol[i,x['snapshot']]=x['mcosinflow_vir']/deltat
    alldflowinfol[i,x['snapshot']]=x['mcosinflow_vir']/deltat-x['mcosinflow_5rh']/deltat
    allz2infol[i,x['snapshot']]=x['mzcosinflow_2rh']/x['mcosinflow_2rh']
    allzinfol[i,x['snapshot']]=x['mzcosinflow_5rh']/x['mcosinflow_5rh']
    allvzinfol[i,x['snapshot']]=x['mzcosinflow_vir']/x['mcosinflow_vir']


    if indplot:
        flowplot1(x,'flowinfo/flowplot_cos_unnorm_'+str(int(obj[wcat_l[i],1]))+'.png','TNG50 LMD: %d' % int(obj[i,1]))
        flowplot2(x,'flowinfo/zflowplot_cos_unnorm_'+str(int(obj[wcat_l[i],1]))+'.png','TNG50 LMD: %d' % int(obj[i,1]))
    
for i in range(len(wcat_xa[0:10])):
    plt.clf()
    if os.path.exists('flowinfo/flowinfo_cos_'+str(int(obj[wcat_xa[i],1]))+'.txt'):
        x=ascii.read('flowinfo/flowinfo_cos_'+str(int(obj[wcat_xa[i],1]))+'.txt')
    else:
        continue
    deltat=np.append(astrocosmo.age(x['redshift'])[1],astrocosmo.age(x['redshift'])[1:]-astrocosmo.age(x['redshift'])[0:-1])*1E9

    allzsxa[i,x['snapshot']]=x['redshift']
    allsfrxa[i,x['snapshot']]=x['sfr']
    allflowinfoxa[i,x['snapshot']]=x['mcosinflow_5rh']/deltat


    allvflowinfoxa[i,x['snapshot']]=x['mcosinflow_vir']/deltat
    alldflowinfoxa[i,x['snapshot']]=x['mcosinflow_vir']/deltat-x['mcosinflow_5rh']/deltat
    allzinfoxa[i,x['snapshot']]=x['mzcosinflow_5rh']/x['mcosinflow_5rh']
    allz2infoxa[i,x['snapshot']]=x['mzcosinflow_2rh']/x['mcosinflow_2rh']
    allvzinfoxa[i,x['snapshot']]=x['mzcosinflow_vir']/x['mcosinflow_vir']

    
    if indplot:
        flowplot1(x,'flowinfo/flowplot_cos_unnorm_'+str(int(obj[wcat_xa[i],1]))+'.png','TNG50 XMD Analog: %d' % int(obj[i,1]))
        flowplot2(x,'flowinfo/zflowplot_cos_unnorm_'+str(int(obj[wcat_xa[i],1]))+'.png','TNG50 XMD Analog: %d' % int(obj[i,1]))

for i in range(len(wcat_la[0:10])):
    plt.clf()
    if os.path.exists('flowinfo/flowinfo_cos_'+str(int(obj[wcat_la[i],1]))+'.txt'):
        x=ascii.read('flowinfo/flowinfo_cos_'+str(int(obj[wcat_la[i],1]))+'.txt')
    else:
        continue
    deltat=np.append(astrocosmo.age(x['redshift'])[1],astrocosmo.age(x['redshift'])[1:]-astrocosmo.age(x['redshift'])[0:-1])*1E9

    allzsla[i,x['snapshot']]=x['redshift']

    allsfrla[i,x['snapshot']]=x['sfr']
    allflowinfola[i,x['snapshot']]=x['mcosinflow_5rh']/deltat
    allvflowinfola[i,x['snapshot']]=x['mcosinflow_vir']/deltat
    alldflowinfola[i,x['snapshot']]=x['mcosinflow_vir']/deltat-x['mcosinflow_5rh']/deltat
    allzinfola[i,x['snapshot']]=x['mzcosinflow_5rh']/x['mcosinflow_5rh']
    allz2infola[i,x['snapshot']]=x['mzcosinflow_2rh']/x['mcosinflow_2rh']
    allvzinfola[i,x['snapshot']]=x['mzcosinflow_vir']/x['mcosinflow_vir']

    
    if indplot:
        flowplot1(x,'flowinfo/flowplot_cos_unnorm_'+str(int(obj[wcat_la[i],1]))+'.png','TNG50 LMD Analog: %d' % int(obj[i,1]))
        flowplot2(x,'flowinfo/zflowplot_cos_unnorm_'+str(int(obj[wcat_la[i],1]))+'.png','TNG50 LMD Analog: %d' % int(obj[i,1]))

    
plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allflowinfola,axis=0),np.nanmedian(allflowinfoxa,axis=0),np.nanmedian(allflowinfol,axis=0),np.nanmedian(allflowinfox,axis=0)],'flowstack_unnorm_5rhcos.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.',':','--','-'],r'Inflow Rate within $5 r_h$ (${\rm M_\odot~yr^{-1}}$)',ylim=[0,1])
 
plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allzinfola,axis=0),np.nanmedian(allzinfoxa,axis=0),np.nanmedian(allzinfol,axis=0),np.nanmedian(allzinfox,axis=0)],'zflowstack_5rhcos.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.',':','--','-'],r'$Z_{\rm inflow}/Z_\odot$',ylog=True,ylim=[1E-4,1E-2])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allz2infola,axis=0),np.nanmedian(allz2infoxa,axis=0),np.nanmedian(allz2infol,axis=0),np.nanmedian(allz2infox,axis=0)],'zflowstack_2rhcos.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.',':','--','-'],r'$Z_{\rm inflow}/Z_\odot$',ylog=True,ylim=[1E-4,1E-2])

plt.clf()
flowplot1stack([astrocosmo.age(np.nanmedian(allzsla,axis=0)),astrocosmo.age(np.nanmedian(allzsxa,axis=0)),astrocosmo.age(np.nanmedian(allzsl,axis=0)),astrocosmo.age(np.nanmedian(allzsx,axis=0))],[np.nanmedian(allvzinfola,axis=0),np.nanmedian(allvzinfoxa,axis=0),np.nanmedian(allvzinfol,axis=0),np.nanmedian(allvzinfox,axis=0)],'zflowstack_vircos.png',['LMD Analogs','XMD Analogs','LMDs','XMDs'],['C3','C2','C1','C0'],['-.',':','--','-'],r'$Z_{\rm inflow}/Z_\odot$',ylog=True,ylim=[1E-4,1E-2])


