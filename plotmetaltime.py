import matplotlib.pyplot as plt
import os
import numpy as np
import h5py
from astropy import cosmology
astrocosmo = cosmology.FlatLambdaCDM(H0=67.74,Om0=0.307)
import matplotlib.ticker as ticker
from astropy.io import ascii
plt.style.use('../python/stylesheet.txt')

plt.rcParams['xtick.top']=False
plt.rcParams['xtick.labeltop']=False
x=ascii.read('sampleselection_sfr5_notide_nomass.txt')

wxmd=np.where(x['type']==1)[0]
wlmd=np.where(x['type']==2)[0]
wxmda=np.where(x['type']==3)[0]
wlmda=np.where(x['type']==4)[0]

ts=np.linspace(0,13.8)

zs=np.array([20.05,14.99,11.98,10.976,10.0,9.389,9.002,8.449,8.012,7.595,7.236,7.005,6.491,6.011,5.847,5.530,5.228,4.996,4.665,4.428,4.177,4.008,3.709,3.491,3.283,3.008,2.896,2.733,2.578,2.444,2.316,2.208,2.103,2.002,1.904,1.823,1.743,1.667,1.604,1.531,1.496,1.414,1.358,1.302,1.248,1.206,1.155,1.114,1.074,1.036,0.997,0.951,0.923,0.887,0.851,0.817,0.791,0.757,0.733,0.700,0.676,0.644,0.621,0.598,0.576,0.546,0.525,0.503,0.482,0.461,0.440,0.420,0.40,0.380,0.361,0.348,0.329,0.310,0.298,0.273,0.261,0.244,0.226,0.214,0.1973,0.1804,0.1693,0.1527,0.1419,0.1258,0.1099,0.0994,0.0839,0.0737,0.0585,0.0485,0.0337,0.0240,0.0095,0])
age=astrocosmo.age(zs)
dt=age[1:]-age[0:-1]
snaps=np.arange(99)
wgyr=np.where(age.value<age.value[-1]-2)[0][-1]
wmyr=np.where(age.value<age.value[-1]-.4)[0][-1]
ylduse=0.2

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
allmetsrx=np.zeros([len(wxmd),100])+np.nan
allstmetsx=np.zeros([len(wxmd),100])+np.nan
allmetsmx=np.zeros([len(wxmd),100])+np.nan
allmassmetsx=np.zeros([len(wxmd),100])+np.nan
allsizex=np.zeros([len(wxmd),100])+np.nan
allsfrx=np.zeros([len(wxmd),100])+np.nan
allmstx=np.zeros([len(wxmd),100])+np.nan
allmhx=np.zeros([len(wxmd),100])+np.nan
allmgx=np.zeros([len(wxmd),100])+np.nan
allmetalholdx=np.zeros([len(wxmd),100])+np.nan
allmetcgmx=np.zeros([len(wxmd),100])+np.nan

for i in range(len(wxmd)):
    h=h5py.File('/Volumes/Elements/illustris/tng100trees/'+str(x['haloid'][wxmd[i]])+'_tree.hdf5','r')
    wm=np.where(h['SnapNum'][1:]-h['SnapNum'][0:-1]>0)[0]
    allmetsx[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],0]/h_mass)/o_nfrac_sun
    allmetsrx[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractionsHalfRad'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractionsHalfRad'][0:wm[0],0]/h_mass)/o_nfrac_sun
    allmetsmx[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractions'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractions'][0:wm[0],0]/h_mass)/o_nfrac_sun
    allmassmetsx[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractions'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],0]/h_mass)/o_nfrac_sun
    allstmetsx[i,h['SnapNum'][0:wm[0]]]=h['SubhaloStarMetalFractions'][0:wm[0],8]/fe_mass/(h['SubhaloStarMetalFractions'][0:wm[0],0]/h_mass)/fe_nfrac_sun
    allsfrx[i,h['SnapNum'][0:wm[0]]]=h['SubhaloSFR'][0:wm[0]]
    allmstx[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMassType'][0:wm[0],4]/.6774*1E10
    allmgx[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMassType'][0:wm[0],0]/.6774*1E10
    allmhx[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMass'][0:wm[0]]/.6774*1E10
    allsizex[i,h['SnapNum'][0:wm[0]]]=h['SubhaloHalfmassRadType'][0:wm[0],4]/.6774
    allmetalholdx[i,h['SnapNum'][0:wm[0]]]=(h['SubhaloStarMetalFractions'][0:wm[0],4]*h['SubhaloMassType'][0:wm[0],4]*1E10/.6774+h['SubhaloGasMetalFractions'][0:wm[0],4]*h['SubhaloMassType'][0:wm[0],0]*1E10/.6774)/(ylduse*h['SubhaloMassType'][0:wm[0],4])
    cgmomass=(h['SubhaloMassType'][0:wm[0],0]*h['SubhaloGasMetalFractions'][0:wm[0],4]-h['SubhaloMassInHalfRadType'][0:wm[0],0]*h['SubhaloGasMetalFractionsHalfRad'][0:wm[0],4])
    cgmhmass=(h['SubhaloMassType'][0:wm[0],0]*h['SubhaloGasMetalFractions'][0:wm[0],0]-h['SubhaloMassInHalfRadType'][0:wm[0],0]*h['SubhaloGasMetalFractionsHalfRad'][0:wm[0],0])
    allmetcgmx[i,h['SnapNum'][0:wm[0]]]=(cgmomass/cgmhmass)/(ox_mass/h_mass)/o_nfrac_sun

allmetsl=np.zeros([len(wlmd),100])+np.nan
allmetsrl=np.zeros([len(wlmd),100])+np.nan
allmetsml=np.zeros([len(wlmd),100])+np.nan
allstmetsl=np.zeros([len(wlmd),100])+np.nan
allmassmetsl=np.zeros([len(wlmd),100])+np.nan
allmassmetsml=np.zeros([len(wlmd),100])+np.nan
allsfrl=np.zeros([len(wlmd),100])+np.nan
allsizel=np.zeros([len(wlmd),100])+np.nan
allmstl=np.zeros([len(wlmd),100])+np.nan
allmhl=np.zeros([len(wlmd),100])+np.nan
allmgl=np.zeros([len(wlmd),100])+np.nan
allmetalholdl=np.zeros([len(wlmd),100])+np.nan
allmetcgml=np.zeros([len(wlmd),100])+np.nan

for i in range(len(wlmd)):
    h=h5py.File('/Volumes/Elements/illustris/tng100trees/'+str(x['haloid'][wlmd[i]])+'_tree.hdf5','r')
    wm=np.where(h['SnapNum'][1:]-h['SnapNum'][0:-1]>0)[0]
    if len(wm)>0:
        allmetsl[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],0]/h_mass)/o_nfrac_sun
        allmetsrl[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractionsHalfRad'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractionsHalfRad'][0:wm[0],0]/h_mass)/o_nfrac_sun
        allmetsml[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractions'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractions'][0:wm[0],0]/h_mass)/o_nfrac_sun
        allmassmetsl[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractions'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],0]/h_mass)/o_nfrac_sun
        allstmetsl[i,h['SnapNum'][0:wm[0]]]=h['SubhaloStarMetalFractions'][0:wm[0],8]/fe_mass/(h['SubhaloStarMetalFractions'][0:wm[0],0]/h_mass)/fe_nfrac_sun

        allsfrl[i,h['SnapNum'][0:wm[0]]]=h['SubhaloSFR'][0:wm[0]]
        allmstl[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMassType'][0:wm[0],4]/.6774*1E10
        allmgl[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMassType'][0:wm[0],0]/.6774*1E10
        allmhl[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMass'][0:wm[0]]/.6774*1E10
        allsizel[i,h['SnapNum'][0:wm[0]]]=h['SubhaloHalfmassRadType'][0:wm[0],4]/.6774
        allmetalholdl[i,h['SnapNum'][0:wm[0]]]=(h['SubhaloStarMetalFractions'][0:wm[0],4]*h['SubhaloMassType'][0:wm[0],4]*1E10/.6774+h['SubhaloGasMetalFractions'][0:wm[0],4]*h['SubhaloMassType'][0:wm[0],0]*1E10/.6774)/(ylduse*h['SubhaloMassType'][0:wm[0],4])
        cgmomass=(h['SubhaloMassType'][0:wm[0],0]*h['SubhaloGasMetalFractions'][0:wm[0],4]-h['SubhaloMassInHalfRadType'][0:wm[0],0]*h['SubhaloGasMetalFractionsHalfRad'][0:wm[0],4])
        cgmhmass=(h['SubhaloMassType'][0:wm[0],0]*h['SubhaloGasMetalFractions'][0:wm[0],0]-h['SubhaloMassInHalfRadType'][0:wm[0],0]*h['SubhaloGasMetalFractionsHalfRad'][0:wm[0],0])
        allmetcgml[i,h['SnapNum'][0:wm[0]]]=(cgmomass/cgmhmass)/(ox_mass/h_mass)/o_nfrac_sun


allstmetsa=np.zeros([len(wxmda),100])+np.nan
allmetsa=np.zeros([len(wxmda),100])+np.nan
allmetsra=np.zeros([len(wxmda),100])+np.nan
allmetsma=np.zeros([len(wxmda),100])+np.nan
allmassmetsa=np.zeros([len(wxmda),100])+np.nan
allsfra=np.zeros([len(wxmda),100])+np.nan
allsizea=np.zeros([len(wxmda),100])+np.nan
allmsta=np.zeros([len(wxmda),100])+np.nan
allmha=np.zeros([len(wxmda),100])+np.nan
allmga=np.zeros([len(wxmda),100])+np.nan
allmetalholda=np.zeros([len(wxmda),100])+np.nan
allmetcgma=np.zeros([len(wxmda),100])+np.nan

for i in range(len(wxmda)):
    if os.path.exists('/Volumes/Elements/illustris/tng100trees/'+str(x['haloid'][wxmda[i]])+'_tree.hdf5'):
        h=h5py.File('/Volumes/Elements/illustris/tng100trees/'+str(x['haloid'][wxmda[i]])+'_tree.hdf5','r')
        wm=np.where(h['SnapNum'][1:]-h['SnapNum'][0:-1]>0)[0]
        if len(wm)>0:
            if h['SubhaloSFR'][0]>.001:
                allmetsa[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],0]/h_mass)/o_nfrac_sun
                allmetsra[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractionsHalfRad'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractionsHalfRad'][0:wm[0],0]/h_mass)/o_nfrac_sun
                allmetsma[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractions'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractions'][0:wm[0],0]/h_mass)/o_nfrac_sun
                allmassmetsa[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractions'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],0]/h_mass)/o_nfrac_sun

                allstmetsa[i,h['SnapNum'][0:wm[0]]]=h['SubhaloStarMetalFractions'][0:wm[0],8]/fe_mass/(h['SubhaloStarMetalFractions'][0:wm[0],0]/h_mass)/fe_nfrac_sun
                
                allsfra[i,h['SnapNum'][0:wm[0]]]=h['SubhaloSFR'][0:wm[0]]
                allmsta[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMassType'][0:wm[0],4]/.6774*1E10
                allmga[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMassType'][0:wm[0],0]/.6774*1E10
                allmha[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMass'][0:wm[0]]/.6774*1E10
                allsizea[i,h['SnapNum'][0:wm[0]]]=h['SubhaloHalfmassRadType'][0:wm[0],4]/.6774
                allmetalholda[i,h['SnapNum'][0:wm[0]]]=(h['SubhaloStarMetalFractions'][0:wm[0],4]*h['SubhaloMassType'][0:wm[0],4]*1E10/.6774+h['SubhaloGasMetalFractions'][0:wm[0],4]*h['SubhaloMassType'][0:wm[0],0]*1E10/.6774)/(ylduse*h['SubhaloMassType'][0:wm[0],4])
                cgmomass=(h['SubhaloMassType'][0:wm[0],0]*h['SubhaloGasMetalFractions'][0:wm[0],4]-h['SubhaloMassInHalfRadType'][0:wm[0],0]*h['SubhaloGasMetalFractionsHalfRad'][0:wm[0],4])
                cgmhmass=(h['SubhaloMassType'][0:wm[0],0]*h['SubhaloGasMetalFractions'][0:wm[0],0]-h['SubhaloMassInHalfRadType'][0:wm[0],0]*h['SubhaloGasMetalFractionsHalfRad'][0:wm[0],0])
                allmetcgma[i,h['SnapNum'][0:wm[0]]]=(cgmomass/cgmhmass)/(ox_mass/h_mass)/o_nfrac_sun


            else:
                allmetsa[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmetsma[i,h['SnapNum'][0:wm[0]]]=np.nan
                allstmetsa[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmassmetsa[i,h['SnapNum'][0:wm[0]]]=np.nan
                allsfra[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmsta[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmga[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmha[i,h['SnapNum'][0:wm[0]]]=np.nan
                allsizea[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmetalholda[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmetcgma[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmetsra[i,h['SnapNum'][0:wm[0]]]=np.nan
                                              
allstmetsla=np.zeros([len(wlmda),100])+np.nan
allmetsla=np.zeros([len(wlmda),100])+np.nan
allmetsrla=np.zeros([len(wlmda),100])+np.nan
allmetsmla=np.zeros([len(wlmda),100])+np.nan
allmassmetsla=np.zeros([len(wlmda),100])+np.nan
allmassmetsmla=np.zeros([len(wlmda),100])+np.nan
allsfrla=np.zeros([len(wlmda),100])+np.nan
allsizela=np.zeros([len(wlmda),100])+np.nan
allmstla=np.zeros([len(wlmda),100])+np.nan
allmhla=np.zeros([len(wlmda),100])+np.nan
allmgla=np.zeros([len(wlmda),100])+np.nan
allmetalholdla=np.zeros([len(wlmda),100])+np.nan
allmetcgmla=np.zeros([len(wlmda),100])+np.nan

for i in range(len(wlmda)):
    if os.path.exists('/Volumes/Elements/illustris/tng100trees/'+str(x['haloid'][wlmda[i]])+'_tree.hdf5'):
        h=h5py.File('/Volumes/Elements/illustris/tng100trees/'+str(x['haloid'][wlmda[i]])+'_tree.hdf5','r')
        wm=np.where(h['SnapNum'][1:]-h['SnapNum'][0:-1]>0)[0]
        if len(wm)>0:
            if h['SubhaloSFR'][0]>.001:
                allmetsla[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],0]/h_mass)/o_nfrac_sun
                allmetsrla[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractionsHalfRad'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractionsHalfRad'][0:wm[0],0]/h_mass)/o_nfrac_sun
                allmetsmla[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractions'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractions'][0:wm[0],0]/h_mass)/o_nfrac_sun
                allmassmetsla[i,h['SnapNum'][0:wm[0]]]=h['SubhaloGasMetalFractions'][0:wm[0],4]/ox_mass/(h['SubhaloGasMetalFractionsSfrWeighted'][0:wm[0],0]/h_mass)/o_nfrac_sun
                allstmetsla[i,h['SnapNum'][0:wm[0]]]=h['SubhaloStarMetalFractions'][0:wm[0],8]/fe_mass/(h['SubhaloStarMetalFractions'][0:wm[0],0]/h_mass)/fe_nfrac_sun

                allsfrla[i,h['SnapNum'][0:wm[0]]]=h['SubhaloSFR'][0:wm[0]]
                allmstla[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMassType'][0:wm[0],4]/.6774*1E10
                allmgla[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMassType'][0:wm[0],0]/.6774*1E10
                allmhla[i,h['SnapNum'][0:wm[0]]]=h['SubhaloMass'][0:wm[0]]/.6774*1E10
                allsizela[i,h['SnapNum'][0:wm[0]]]=h['SubhaloHalfmassRadType'][0:wm[0],4]/.6774
                allmetalholdla[i,h['SnapNum'][0:wm[0]]]=(h['SubhaloStarMetalFractions'][0:wm[0],4]*h['SubhaloMassType'][0:wm[0],4]*1E10/.6774+h['SubhaloGasMetalFractions'][0:wm[0],4]*h['SubhaloMassType'][0:wm[0],0]*1E10/.6774)/(ylduse*h['SubhaloMassType'][0:wm[0],4])
                cgmomass=(h['SubhaloMassType'][0:wm[0],0]*h['SubhaloGasMetalFractions'][0:wm[0],4]-h['SubhaloMassInHalfRadType'][0:wm[0],0]*h['SubhaloGasMetalFractionsHalfRad'][0:wm[0],4])
                cgmhmass=(h['SubhaloMassType'][0:wm[0],0]*h['SubhaloGasMetalFractions'][0:wm[0],0]-h['SubhaloMassInHalfRadType'][0:wm[0],0]*h['SubhaloGasMetalFractionsHalfRad'][0:wm[0],0])
                allmetcgmla[i,h['SnapNum'][0:wm[0]]]=(cgmomass/cgmhmass)/(ox_mass/h_mass)/o_nfrac_sun


            else:
                allmetsla[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmetsmla[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmassmetsla[i,h['SnapNum'][0:wm[0]]]=np.nan
                allstmetsla[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmassmetsmla[i,h['SnapNum'][0:wm[0]]]=np.nan
                allsfrla[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmstla[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmgla[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmhla[i,h['SnapNum'][0:wm[0]]]=np.nan
                allsizela[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmetalholdla[i,h['SnapNum'][0:wm[0]]]=np.nan
                allmetcgmla[i,h['SnapNum'][0:wm[0]]]
                allmetsrla[i,h['SnapNum'][0:wm[0]]]=np.nan

np.savetxt('allsfrx.txt',allsfrx)
np.savetxt('allsfrxa.txt',allsfra)
np.savetxt('allsfrl.txt',allsfrl)
np.savetxt('allsfrla.txt',allsfrla)

np.savetxt('allmetx.txt',allmetsx)
np.savetxt('allmetxa.txt',allmetsa)
np.savetxt('allmetl.txt',allmetsl)
np.savetxt('allmetla.txt',allmetsla)

np.savetxt('allmassmetx.txt',allmetsmx)
np.savetxt('allmassmetxa.txt',allmetsma)
np.savetxt('allmassmetl.txt',allmetsml)
np.savetxt('allmassmetla.txt',allmetsmla)

np.savetxt('allmstx.txt',allmstx)
np.savetxt('allmstxa.txt',allmsta)
np.savetxt('allmstl.txt',allmstl)
np.savetxt('allmstla.txt',allmstla)

np.savetxt('allmgx.txt',allmgx)
np.savetxt('allmgxa.txt',allmga)
np.savetxt('allmgl.txt',allmgl)
np.savetxt('allmgla.txt',allmgla)

np.savetxt('allmhx.txt',allmhx)
np.savetxt('allmhxa.txt',allmha)
np.savetxt('allmhl.txt',allmhl)
np.savetxt('allmhla.txt',allmhla)

avgsfrx=np.zeros([len(allsfrx),len(allsfrx[0])-1])
avgsfrxa=np.zeros([len(allsfra),len(allsfra[0])-1])
avgsfrl=np.zeros([len(allsfrl),len(allsfrl[0])-1])
avgsfrla=np.zeros([len(allsfrla),len(allsfrla[0])-1])

for i in range(len(allsfrx[0])-1):
    avgsfrx[:,i]=np.nansum(allsfrx[:,i:-1]*np.repeat(dt,len(allsfrx)).reshape(len(allsfrx),len(dt))[:,i:],axis=1)
    avgsfrxa[:,i]=np.nansum(allsfra[:,i:-1]*np.repeat(dt,len(allsfra)).reshape(len(allsfra),len(dt))[:,i:],axis=1)
    avgsfrl[:,i]=np.nansum(allsfrl[:,i:-1]*np.repeat(dt,len(allsfrl)).reshape(len(allsfrl),len(dt))[:,i:],axis=1)
    avgsfrla[:,i]=np.nansum(allsfrla[:,i:-1]*np.repeat(dt,len(allsfrla)).reshape(len(allsfrla),len(dt))[:,i:],axis=1)

colddat=ascii.read('densgasinfob_notide_nomass.txt')
wlmdacold=np.where(colddat['type']==4)[0]
wxmdacold=np.where(colddat['type']==3)[0]
wlmdcold=np.where(colddat['type']==2)[0]
wxmdcold=np.where(colddat['type']==1)[0]
whsfe=np.where(avgsfrx[:,90]/colddat['mcold'][wxmdcold]*1e9>.1)[0]
wmsfe=np.where((avgsfrx[:,90]/colddat['mcold'][wxmdcold]*1e9<.1) & (avgsfrx[:,90]/colddat['mcold'][wxmdcold]*1e9>.01))[0]
wlsfe=np.where(avgsfrx[:,90]/colddat['mcold'][wxmdcold]*1e9<.01)[0]

plt.clf()
#plt.plot(avgsfrla[:,90]/colddat['mcold'][wlmdacold],allmetsla[:,99]/allmassmetsla[:,99],'C3D',alpha=.5)
#plt.plot(avgsfrxa[:,90]/colddat['mcold'][wxmdacold],allmetsa[:,99]/allmassmetsa[:,99],'C2s',alpha=.5)
#plt.plot(avgsfrl[:,90]/colddat['mcold'][wlmdcold],allmetsl[:,99]/allmassmetsl[:,99],'C1x',alpha=.5)
#plt.plot(avgsfrx[:,90]/colddat['mcold'][wxmdcold],allmetsx[:,99]/allmassmetsx[:,99],'C0*',alpha=.5)

#plt.plot(avgsfrla[:,77]/allmgla[:,99],allmetsla[:,99]/allmassmetsla[:,99],'C3D',alpha=.5)
#plt.plot(avgsfrxa[:,77]/allmga[:,99],allmetsa[:,99]/allmassmetsa[:,99],'C2s',alpha=.5)
#plt.plot(avgsfrl[:,77]/allmgl[:,99],allmetsl[:,99]/allmassmetsl[:,99],'C1x',alpha=.5)
#plt.plot(avgsfrx[:,77]/allmgx[:,99],allmetsx[:,99]/allmassmetsx[:,99],'C0*',alpha=.5)

#plt.plot(avgsfrla[:,77],allmetsla[:,99]/allmassmetsla[:,99],'C3D',alpha=.5)
#plt.plot(avgsfrxa[:,77],allmetsa[:,99]/allmassmetsa[:,99],'C2s',alpha=.5)
#plt.plot(avgsfrl[:,77],allmetsl[:,99]/allmassmetsl[:,99],'C1x',alpha=.5)
#plt.plot(avgsfrx[:,77],allmetsx[:,99]/allmassmetsx[:,99],'C0*',alpha=.5)

plt.scatter(avgsfrla[:,77]/allmgla[:,99],allmetsrla[:,99]/allmetcgmla[:,99],c=np.log10(avgsfrla[:,98]/allmgla[:,99]),alpha=.5,vmin=-14,vmax=-9)
plt.scatter(avgsfrxa[:,77]/allmga[:,99],allmetsra[:,99]/allmetcgma[:,99],c=np.log10(avgsfrxa[:,98]/allmga[:,99]),alpha=.5,vmin=-14,vmax=-9)
plt.scatter(avgsfrl[:,77]/allmgl[:,99],allmetsrl[:,99]/allmetcgml[:,99],c=np.log10(avgsfrl[:,98]/allmgl[:,99]),alpha=.5,vmin=-14,vmax=-9)
plt.scatter(avgsfrx[:,77]/allmgx[:,99],allmetsrx[:,99]/allmetcgmx[:,99],c=np.log10(avgsfrx[:,98]/allmgx[:,99]),alpha=.5,vmin=-14,vmax=-9)
plt.loglog()
plt.ylim(.1,10)
plt.xlabel('SFE90')
plt.ylabel('met ratio')
plt.savefig('zsfrmasssfe.png')

plt.clf()
#plt.plot(age,np.nanmedian((allmetsla[whsfe]/colddat['mcold'][wxmdcold[whsfe]])/(allmetsla[whsfe]/colddat['mcold'][wxmdcold[whsfe]]+allmassmetsx[whsfe]*allmgx[whsfe]),axis=0),label='High SFE XMDs',lw=5,color='C0',ls='--')
#plt.plot(age,np.nanmedian((allmetsla[wmsfe]/colddat['mcold'][wxmdcold[wmsfe]])/(allmetsla[wmsfe]/colddat['mcold'][wxmdcold[wmsfe]]+allmassmetsx[wmsfe]*allmgx[wmsfe]),axis=0),label='Mid SFE XMDs',lw=5,color='C1',ls='-')
#plt.plot(age,np.nanmedian((allmetsla[wlsfe]/colddat['mcold'][wxmdcold[wlsfe]])/(allmetsla[wlsfe]/colddat['mcold'][wxmdcold[wlsfe]]+allmassmetsx[wlsfe]*allmgx[wlsfe]),axis=0),label='Low SFE XMDs',lw=5,color='C2',ls=':')
plt.plot(age,np.nanmedian(allmetsx[whsfe]/allmassmetsx[whsfe],axis=0),label='High SFE XMDs',lw=5,color='C0',ls='--')
plt.plot(age,np.nanmedian(allmetsx[wmsfe]/allmassmetsx[wmsfe],axis=0),label='Mid SFE XMDs',lw=5,color='C1',ls='-')
plt.plot(age,np.nanmedian(allmetsx[wlsfe]/allmassmetsx[wlsfe],axis=0),label='Low SFE XMDs',lw=5,color='C2',ls=':')
plt.legend()
plt.yscale('log')
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.ylabel(r'$Z_{\rm SFR}/Z_{\rm mass}$')
plt.savefig('metalsffrac.png')
print('a')
allsfrratiox=(allmstx[:,99]-allmstx[:,wgyr])/(allmstx[:,99]-allmstx[:,wmyr])/10
allsfrratioa=(allmsta[:,99]-allmsta[:,wgyr])/(allmsta[:,99]-allmsta[:,wmyr])/10
allsfrratiol=(allmstl[:,99]-allmstl[:,wgyr])/(allmstl[:,99]-allmstl[:,wmyr])/10
allsfrratiola=(allmstla[:,99]-allmstla[:,wgyr])/(allmstla[:,99]-allmstla[:,wmyr])/10

plt.clf()
plt.plot((allmstla[:,99]-allmstla[:,wgyr])/allmgla[:,-1],allmetsla[:,99],'C3D',alpha=.5)
plt.plot((allmstl[:,99]-allmstl[:,wgyr])/allmgl[:,-1],allmetsl[:,99],'C3D',alpha=.5)
plt.plot((allmsta[:,99]-allmsta[:,wgyr])/allmga[:,-1],allmetsa[:,99],'C3D',alpha=.5)
plt.plot((allmstx[:,99]-allmstx[:,wgyr])/allmgx[:,-1],allmetsx[:,99],'C3D',alpha=.5)
plt.loglog()
plt.xlabel(r'$M_{\rm *,~last~Gyr}~({\rm M_\odot})$')
plt.ylabel(r'$Z_{\rm SFR}/Z_{\rm \odot}$')
plt.savefig('massgyrmet.png')

import bindata

w1r=np.where(np.isfinite(np.log10(allsfrratiola)))[0]
w2r=np.where(np.isfinite(np.log10(allsfrratiol)))[0]
w3r=np.where(np.isfinite(np.log10(allsfrratioa)))[0]
w4r=np.where(np.isfinite(np.log10(allsfrratiox)))[0]

b1=bindata.bindata(np.log10(allmetsla[w1r,99]),np.log10(allsfrratiola[w1r]),np.arange(-2,0,.15),logerr=1,median=1)
b2=bindata.bindata(np.log10(allmetsl[w2r,99]),np.log10(allsfrratiol[w2r]),np.arange(-2,0,.15),logerr=1,median=1)
b3=bindata.bindata(np.log10(allmetsa[w3r,99]),np.log10(allsfrratioa[w3r]),np.arange(-2,0,.15),logerr=1,median=1)
b4=bindata.bindata(np.log10(allmetsx[w4r,99]),np.log10(allsfrratiox[w4r]),np.arange(-2,0,.15),logerr=1,median=1)
plt.clf()
plt.plot(allmetsla[:,99],allsfrratiola,'C3D',alpha=.5)
plt.plot(allmetsl[:,99],allsfrratiol,'C1x',alpha=.5)
plt.plot(allmetsa[:,99],allsfrratioa,'C2s',alpha=.5)
plt.plot(allmetsx[:,99],allsfrratiox,'C0*',alpha=.5)
plt.plot(10**b1[2],10**b1[0],'C3-.')
plt.plot(10**b2[2],10**b2[0],'C1--')
plt.plot(10**b3[2],10**b3[0],'C2:')
plt.plot(10**b4[2],10**b4[0],'C0-')
plt.ylabel(r'${\rm SFR_{1~Gyr}/SFR_{100~Myr}}$')
plt.xlabel(r'$Z_{\rm SFR}/Z_\odot$')
plt.loglog()
plt.savefig('sfrratiovsmet.png')


zg=np.array([0,.1,.25,.5,.75,1,1.5,2,3,7])
zglabels=zg.astype(np.str)
agepoints=astrocosmo.age(zg)
toplotx=np.random.choice(len(allmetsx),size=10)
toplota=np.random.choice(len(allmetsa),size=10)
toplotl=np.random.choice(len(allmetsl),size=10)
toplotla=np.random.choice(len(allmetsla),size=10)


plt.clf()
plt.plot(age,np.nanmedian(allmetcgmla,axis=0),label='SFR-weighted LMD Analogs',lw=5,color='C3',ls='-.')
plt.plot(age,np.nanmedian(allmetcgml,axis=0),label='SFR-weighted LMDs',lw=5,color='C1',ls=':')
plt.plot(age,np.nanmedian(allmetcgma,axis=0),label='SFR-weighted XMD Analogs',lw=5,color='C2',ls='--')
plt.plot(age,np.nanmedian(allmetcgmx,axis=0),label='SFR-weighted XMDs',lw=5,color='C0')
plt.legend(loc='lower left',ncol=2,fontsize=30)
plt.yscale('log')
plt.ylim(1E-2,2)
plt.xlim(.5,13.9)
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.ylabel(r'$Z_{\rm CGM}/Z_\odot$')
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlim(.5,13.9)
ax.set_xlabel(r'$z$')
plt.savefig('allcgmmetvsz.png')


plt.clf()
plt.plot(age,np.nanmedian(allmetsla,axis=0),label='SFR-weighted LMD Analogs',lw=5,color='C3',ls='-.')
plt.plot(age,np.nanmedian(allmetsl,axis=0),label='SFR-weighted LMDs',lw=5,color='C1',ls='--')
plt.plot(age,np.nanmedian(allmetsa,axis=0),label='SFR-weighted XMD Analogs',lw=5,color='C2',ls=':')
plt.plot(age,np.nanmedian(allmetsx,axis=0),label='SFR-weighted XMDs',lw=5,color='C0')

plt.plot(age,np.nanmedian(allmassmetsla,axis=0),label='Mass-weighted LMD Analogs',lw=5,color='C3',ls=(0, (3, 1, 1, 1)))
plt.plot(age,np.nanmedian(allmassmetsl,axis=0),label='Mass-weighted LMDs',lw=5,color='C1',ls=(0, (5, 1)))
plt.plot(age,np.nanmedian(allmassmetsa,axis=0),label='Mass-weighted XMD Analogs',lw=5,color='C2',ls= (0, (1, 1)))
plt.plot(age,np.nanmedian(allmassmetsx,axis=0),label='Mass-weighted XMDs',lw=5,color='C0',ls=(5, (10, 3)))



#plt.fill_between(age,np.nanmedian(allmetsx,axis=0)-np.nanstd(allmetsx,axis=0),np.nanmedian(allmetsx,axis=0)+np.nanstd(allmetsx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsa,axis=0)-np.nanstd(allmetsa,axis=0),np.nanmedian(allmetsa,axis=0)+np.nanstd(allmetsa,axis=0),color='C1',alpha=.3)
plt.legend(loc='lower left',ncol=2,fontsize=25)
plt.yscale('log')
plt.ylim(1E-2,2)
plt.xlim(.5,13.9)
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.ylabel('Metal Fraction')
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlim(.5,13.9)
ax.set_xlabel(r'$z$')
plt.title('TNG100',y=.92)
plt.savefig('allsfrmetvsz_mvss.png')


plt.clf()
blnk,=plt.plot([],[],color='None')
sflaln,=plt.plot(age,np.nanmedian(allsfrla,axis=0),label='LMD Analogs',lw=5,color='C3',ls='-.')
sflln,=plt.plot(age,np.nanmedian(allsfrl,axis=0),label='LMDs',lw=5,color='C1',ls=':')
sfxaln,=plt.plot(age,np.nanmedian(allsfra,axis=0),label='XMD Analogs',lw=5,color='C2',ls='--')
sfxln,=plt.plot(age,np.nanmedian(allsfrx,axis=0),label='XMDs',lw=5,color='C0',ls='-')
plt.yscale('log')
plt.ylabel(r'SFR (M$_\odot$ yr$^{-1}$)')
plt.xlim(.5,13.9)
plt.xlabel(r'$t_{\rm age}$ (Gyr)')

ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlim(.5,13.9)
ax.set_xlabel(r'$z$')
plt.savefig('sfrhist.png')

allmfracyx=[]
allmfracyl=[]
allmfracyxa=[]
allmfracyla=[]

allmfracmx=[]
allmfracml=[]
allmfracmxa=[]
allmfracmla=[]

allmfrac3x=[]
allmfrac3l=[]
allmfrac3xa=[]
allmfrac3la=[]
for i in range(len(allmstx)):
    allmfracyx.append(allmstx[i][-3]/allmstx[i][-1])
    allmfracyl.append(allmstl[i][-3]/allmstl[i][-1])
    allmfracyxa.append(allmsta[i][-3]/allmsta[i][-1])
    allmfracyla.append(allmstla[i][-3]/allmstla[i][-1])

    allmfracmx.append(allmstx[i][-19]/allmstx[i][-1])
    allmfracml.append(allmstl[i][-19]/allmstl[i][-1])
    allmfracmxa.append(allmsta[i][-19]/allmsta[i][-1])
    allmfracmla.append(allmstla[i][-19]/allmstla[i][-1])

    allmfrac3x.append(allmstx[i][78]/allmstx[i][-1])
    allmfrac3l.append(allmstl[i][78]/allmstl[i][-1])
    allmfrac3xa.append(allmsta[i][78]/allmsta[i][-1])
    allmfrac3la.append(allmstla[i][78]/allmstla[i][-1])

print('young xmd',np.median(allmfracyx))
print('young lmd',np.median(allmfracyl))
print('young xmda',np.nanmedian(allmfracyxa))
print('young lmda',np.nanmedian(allmfracyla))

print('mid xmd',np.median(allmfracmx))
print('mid lmd',np.median(allmfracml))
print('mid xmda',np.nanmedian(allmfracmxa))
print('mid lmda',np.nanmedian(allmfracmla))

print('z0.3 xmd',np.median(allmfrac3x))
print('z0.3 lmd',np.median(allmfrac3l))
print('z0.3 xmda',np.nanmedian(allmfrac3xa))
print('z0.3 lmda',np.nanmedian(allmfrac3la))

plt.clf()
mstlaln,=plt.plot(age,np.nanmedian(allmstla,axis=0),label='LMD Analogs Stellar',lw=5,color='C3',ls='-.')
mstlln,=plt.plot(age,np.nanmedian(allmstl,axis=0),label='LMDs Stellar',lw=5,color='C1',ls='--')
mstxaln,=plt.plot(age,np.nanmedian(allmsta,axis=0),label='XMD Analogs Stellar',lw=5,color='C2',ls=':')
mstxln,=plt.plot(age,np.nanmedian(allmstx,axis=0),label='XMDs Stellar',lw=5,color='C0',ls='-')

plt.ylabel(r'$M_*$ (M$_\odot$)')
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.yscale('log')
plt.ylim(2E5,2E8)

axg=plt.twinx()
mglaln,=axg.plot(age,np.nanmedian(allmgla,axis=0),label='LMD Analogs Gas',lw=5,color='C3',ls=(0, (3, 10, 1, 10)))
mglln,=axg.plot(age,np.nanmedian(allmgl,axis=0),label='LMDs Gas',lw=5,color='C1',ls=(0, (5, 10)))
mgxaln,=axg.plot(age,np.nanmedian(allmga,axis=0),label='XMD Analogs Gas',lw=5,color='C2',ls=(0, (1, 3)))
mgxln,=axg.plot(age,np.nanmedian(allmgx,axis=0),label='XMDs Gas',lw=5,color='C0',ls=(5, (10, 5)))

mdlaln,=axg.plot(age,np.nanmedian(allmhla,axis=0),label='LMD Analogs DM',lw=5,color='C3',ls=(0, (3, 1, 1, 1)))
mdlln,=axg.plot(age,np.nanmedian(allmhl,axis=0),label='LMDs DM',lw=5,color='C1',ls=(0, (5, 1)))
mdxaln,=axg.plot(age,np.nanmedian(allmha,axis=0),label='XMD Analogs DM',lw=5,color='C2',ls=(0, (1, 1)))
mdxln,=axg.plot(age,np.nanmedian(allmhx,axis=0),label='XMDs DM',lw=5,color='C0',ls=(5, (10, 2)))
axg.set_yscale('log')
axg.set_ylim(5E7,5E10)
#plt.fill_between(age,np.nanmedian(allmetsx,axis=0)-np.nanstd(allmetsx,axis=0),np.nanmedian(allmetsx,axis=0)+np.nanstd(allmetsx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsa,axis=0)-np.nanstd(allmetsa,axis=0),np.nanmedian(allmetsa,axis=0)+np.nanstd(allmetsa,axis=0),color='C1',alpha=.3)
plt.xlim(.5,13.9)
axg.set_ylabel(r'$M_{g;DM}$ (M$_\odot$)')
plt.legend(handles=[mstxln,mstlln,mstxaln,mstlaln,mgxln,mglln,mgxaln,mglaln,mdxln,mdlln,mdxaln,mdlaln],labels=['XMDs Stellar','LMDs Stellar','XMD Analogs Stellar','LMD Analogs Stellar','XMDs Gas','LMDs Gas','XMD Analogs Gas','LMD Analogs Gas','XMDs DM','LMDs DM','XMD Analogs DM','LMD Analogs DM'],ncol=3,fontsize=20)
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlim(.5,13.9)
ax.set_xlabel(r'$z$')
plt.savefig('masshist.png')


plt.clf()

fig=plt.figure(figsize=(16,28))
mstlaln,=plt.plot(age,np.nanmedian(allmstla,axis=0),label='LMD Analogs Stellar',lw=5,color='C3',ls='-.')
mstlln,=plt.plot(age,np.nanmedian(allmstl,axis=0),label='LMDs Stellar',lw=5,color='C1',ls='--')
mstxaln,=plt.plot(age,np.nanmedian(allmsta,axis=0),label='XMD Analogs Stellar',lw=5,color='C2',ls=':')
mstxln,=plt.plot(age,np.nanmedian(allmstx,axis=0),label='XMDs Stellar',lw=5,color='C0',ls='-')

plt.ylabel(r'Mass (M$_\odot$)')
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.yscale('log')

mglaln,=plt.plot(age,np.nanmedian(allmgla,axis=0),label='LMD Analogs Gas',lw=5,color='C3',ls=(0, (3, 10, 1, 10)))
mglln,=plt.plot(age,np.nanmedian(allmgl,axis=0),label='LMDs Gas',lw=5,color='C1',ls=(0, (5, 10)))
mgxaln,=plt.plot(age,np.nanmedian(allmga,axis=0),label='XMD Analogs Gas',lw=5,color='C2',ls=(0, (1, 3)))
mgxln,=plt.plot(age,np.nanmedian(allmgx,axis=0),label='XMDs Gas',lw=5,color='C0',ls=(5, (10, 5)))

mdlaln,=plt.plot(age,np.nanmedian(allmhla,axis=0),label='LMD Analogs DM',lw=5,color='C3',ls=(0, (3, 1, 1, 1)))
mdlln,=plt.plot(age,np.nanmedian(allmhl,axis=0),label='LMDs DM',lw=5,color='C1',ls=(0, (5, 1)))
mdxaln,=plt.plot(age,np.nanmedian(allmha,axis=0),label='XMD Analogs DM',lw=5,color='C2',ls=(0, (1, 1)))
mdxln,=plt.plot(age,np.nanmedian(allmhx,axis=0),label='XMDs DM',lw=5,color='C0',ls=(5, (10, 2)))

#plt.fill_between(age,np.nanmedian(allmetsx,axis=0)-np.nanstd(allmetsx,axis=0),np.nanmedian(allmetsx,axis=0)+np.nanstd(allmetsx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsa,axis=0)-np.nanstd(allmetsa,axis=0),np.nanmedian(allmetsa,axis=0)+np.nanstd(allmetsa,axis=0),color='C1',alpha=.3)
plt.xlim(.5,13.9)
plt.ylim(1E6,5E10)
plt.yscale('log')
#plt.legend(handles=[mstxln,mstlln,mstxaln,mstlaln,mgxln,mglln,mgxaln,mglaln,mdxln,mdlln,mdxaln,mdlaln],labels=['XMDs Stellar','LMDs Stellar','XMD Analogs Stellar','LMD Analogs Stellar','XMDs Gas','LMDs Gas','XMD Analogs Gas','LMD Analogs Gas','XMDs DM','LMDs DM','XMD Analogs DM','LMD Analogs DM'],ncol=3,fontsize=20)
#plt.legend(handles=[mstxln,mstlln,mstxaln,mstlaln,mgxln,mglln,mgxaln,mglaln,mdxln,mdlln,mdxaln,mdlaln],labels=['XMDs','LMDs','XMD Analogs','LMD Analogs'],ncol=1,fontsize=20)
plt.legend(handles=[mstxln,mstlln,mstxaln,mstlaln,mgxln,mglln,mgxaln,mglaln,mdxln,mdlln,mdxaln,mdlaln],labels=['XMDs','LMDs','XMD Analogs','LMD Analogs'],ncol=1)

plt.text(9,2E7,'Stellar',fontsize=32)
plt.text(9,1.25E9,'Gas',fontsize=32)
plt.text(9,2E10,'Dark Matter',fontsize=32)
#plt.gca().xaxis.set_minor_locator(ticker.MultipleLocator(1))
plt.savefig('masshist2.png',dpi=300)

fig=plt.figure()
plt.clf()
for i in toplotx:
    plt.plot(age,allmetsx[i],color='C0',alpha=.3,ls='-')

for i in toplota:
    plt.plot(age,allmetsa[i],color='C2',alpha=.3,ls=':')

for i in toplotl:
    plt.plot(age,allmetsl[i],color='C1',alpha=.3,ls='--')

for i in toplotla:
    plt.plot(age,allmetsla[i],color='C3',alpha=.3,ls='-.')


plt.plot(age,np.nanmedian(allmetsla,axis=0),label='LMD Analogs',lw=5,color='C3',ls='-.')
plt.plot(age,np.nanmedian(allmetsl,axis=0),label='LMDs',lw=5,color='C1',ls='--')
plt.plot(age,np.nanmedian(allmetsa,axis=0),label='XMD Analogs',lw=5,color='C2',ls=':')
plt.plot(age,np.nanmedian(allmetsx,axis=0),label='XMDs',lw=5,color='C0',ls='-')
#plt.fill_between(age,np.nanmedian(allmetsx,axis=0)-np.nanstd(allmetsx,axis=0),np.nanmedian(allmetsx,axis=0)+np.nanstd(allmetsx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsa,axis=0)-np.nanstd(allmetsa,axis=0),np.nanmedian(allmetsa,axis=0)+np.nanstd(allmetsa,axis=0),color='C1',alpha=.3)
plt.legend()
plt.yscale('log')
plt.ylim(1E-2,2)
plt.xlim(.5,13.9)
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.ylabel(r'$Z_{\rm SFR}/Z_\odot$')
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlim(.5,13.9)
ax.set_xlabel(r'$z$')
plt.title('TNG100',y=0.92)
plt.savefig('allsfrmetvsz.png')

plt.clf()
for i in toplotx:
    plt.plot(age,allsizex[i],color='C0',alpha=.3)

for i in toplota:
    plt.plot(age,allsizea[i],color='C2',alpha=.3,ls=':')

for i in toplotl:
    plt.plot(age,allsizel[i],color='C1',alpha=.3,ls='--')

for i in toplotla:
    plt.plot(age,allsizela[i],color='C3',alpha=.3,ls='-.')

plt.plot(age,np.nanmedian(allsizela,axis=0),label='LMD Analogs',lw=5,color='C3',ls='-.')    
plt.plot(age,np.nanmedian(allsizel,axis=0),label='LMDs',lw=5,color='C1',ls='--')
plt.plot(age,np.nanmedian(allsizea,axis=0),label='XMD Analogs',lw=5,color='C2',ls=':')
plt.plot(age,np.nanmedian(allsizex,axis=0),label='XMDs',lw=5,color='C0')

plt.legend()
plt.yscale('log')
plt.xlim(.5,13.9)
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.ylabel('Gas Half mass radius (kpc)')
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlim(.5,13.9)
ax.set_xlabel(r'$z$')
plt.savefig('allsizeev.png')

plt.clf()
for i in toplotx:
    plt.plot(age,allmetsmx[i],color='C0',alpha=.3)

for i in toplota:
    plt.plot(age,allmetsma[i],color='C1',alpha=.3)

for i in toplotl:
    plt.plot(age,allmetsml[i],color='C2',alpha=.3)

for i in toplotla:
    plt.plot(age,allmetsmla[i],color='C3',alpha=.3)

    
plt.plot(age,np.nanmedian(allmetsmx,axis=0),label='XMDs',lw=5,color='C0')
plt.plot(age,np.nanmedian(allmetsml,axis=0),label='LMDs',lw=5,color='C1')
plt.plot(age,np.nanmedian(allmetsma,axis=0),label='XMD Analogs',lw=5,color='C2')
plt.plot(age,np.nanmedian(allmetsmla,axis=0),label='LMD Analogs',lw=5,color='C3')
#plt.fill_between(age,np.nanmedian(allmetsmx,axis=0)-np.nanstd(allmetsmx,axis=0),np.nanmedian(allmetsmx,axis=0)+np.nanstd(allmetsmx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsma,axis=0)-np.nanstd(allmetsma,axis=0),np.nanmedian(allmetsma,axis=0)+np.nanstd(allmetsma,axis=0),color='C1',alpha=.3)
plt.legend()
plt.yscale('log')
plt.ylim(5E-4,5E-3)
plt.xlim(.5,13.9)
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.ylabel('Mass-weighted Metal fraction')
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlim(.5,13.9)
ax.set_xlabel(r'$z$')
plt.savefig('allmassmetvsz.png')

plt.clf()
for i in toplotx:
    plt.plot(age,allmetsx[i]/allmetsx[i][-1],color='C0',alpha=.3,ls='-')

for i in toplota:
    plt.plot(age,allmetsa[i]/allmetsa[i][-1],color='C2',alpha=.3,ls=':')

for i in toplotl:
    plt.plot(age,allmetsl[i]/allmetsl[i][-1],color='C1',alpha=.3,ls='--')

for i in toplotla:
    plt.plot(age,allmetsla[i]/allmetsla[i][-1],color='C3',alpha=.3,ls='-.')

z0metsx=np.repeat(allmetsx[:,-1],100).reshape(len(allmetsx),100)
z0metsa=np.repeat(allmetsa[:,-1],100).reshape(len(allmetsa),100)
z0metsl=np.repeat(allmetsl[:,-1],100).reshape(len(allmetsl),100)
z0metsla=np.repeat(allmetsla[:,-1],100).reshape(len(allmetsla),100)

z0massmetsx=np.repeat(allmassmetsx[:,-1],100).reshape(len(allmetsx),100)
z0massmetsa=np.repeat(allmassmetsa[:,-1],100).reshape(len(allmetsa),100)
z0massmetsl=np.repeat(allmassmetsl[:,-1],100).reshape(len(allmetsl),100)
z0massmetsla=np.repeat(allmassmetsla[:,-1],100).reshape(len(allmetsla),100)
print(allmetsx.shape,z0metsx.shape)
plt.plot(age,np.nanmedian(allmetsla/z0metsla,axis=0),label='LMD Analaogs',lw=5,color='C3',ls='-.')
plt.plot(age,np.nanmedian(allmetsl/z0metsl,axis=0),label='LMDs',lw=5,color='C1',ls='--')
plt.plot(age,np.nanmedian(allmetsa/z0metsa,axis=0),label='XMD Analogs',lw=5,color='C2',ls=':')
plt.plot(age,np.nanmedian(allmetsx/z0metsx,axis=0),label='XMDs',lw=5,color='C0')
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
    
plt.savefig('allsfrmetfracgrowth.png')

plt.clf()
plt.plot(age,np.nanmedian(allmetsx/z0metsx,axis=0),label='XMDs',lw=5,color='C0')
plt.plot(age,np.nanmedian(allmetsa/z0metsa,axis=0),label='XMD Analogs',lw=5,color='C1')
plt.plot(age,np.nanmedian(allmetsl/z0metsl,axis=0),label='LMDs',lw=5,color='C2')
plt.plot(age,np.nanmedian(allmetsla/z0metsla,axis=0),label='LMD Analaogs',lw=5,color='C3')

plt.plot(age,np.nanmedian(allmassmetsx/z0massmetsx,axis=0),label='XMDs',lw=5,color='C0',ls='--')
plt.plot(age,np.nanmedian(allmassmetsa/z0massmetsa,axis=0),label='XMD Analogs',lw=5,color='C1',ls='--')
plt.plot(age,np.nanmedian(allmassmetsl/z0massmetsl,axis=0),label='LMDs',lw=5,color='C2',ls='--')
plt.plot(age,np.nanmedian(allmassmetsla/z0massmetsla,axis=0),label='LMD Analaogs',lw=5,color='C3',ls='--')
plt.legend()
plt.yscale('log')
plt.ylim(.1,10)
plt.xlim(.5,13.9)
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.ylabel(r'$Zs_{\rm SFR}/Z_{\rm SFR, z=0}$')
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlim(.5,13.9)
ax.set_xlabel(r'$z$')
    
plt.savefig('allsfrmetfracgrowth_vsmass.png')

plt.clf()
for i in toplotx:
    plt.plot(age,allsfrx[i],color='C0',alpha=.3)

for i in toplota:
    plt.plot(age,allsfra[i],color='C1',alpha=.3)

for i in toplotl:
    plt.plot(age,allsfrl[i],color='C2',alpha=.3)

for i in toplotla:
    plt.plot(age,allsfrla[i],color='C3',alpha=.3)

plt.plot(age,np.nanmedian(allsfrx,axis=0),label='XMDs',lw=5,color='C0')
plt.plot(age,np.nanmedian(allsfra,axis=0),label='XMD Analogs',lw=5,color='C1')
plt.plot(age,np.nanmedian(allsfrl,axis=0),label='LMDs',lw=5,color='C2')
plt.plot(age,np.nanmedian(allsfrla,axis=0),label='LMD Ananlogs',lw=5,color='C3')
#plt.fill_between(age,np.nanmedian(allmetsx,axis=0)-np.nanstd(allmetsx,axis=0),np.nanmedian(allmetsx,axis=0)+np.nanstd(allmetsx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsa,axis=0)-np.nanstd(allmetsa,axis=0),np.nanmedian(allmetsa,axis=0)+np.nanstd(allmetsa,axis=0),color='C1',alpha=.3)
plt.legend()
plt.yscale('log')
#plt.ylim(5E-4,5E-3)
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.xlim(.5,13.9)
plt.ylabel(r'SFR $({\rm M_\odot~yr^{-1}}$)')
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlabel(r'$z$')
ax.set_xlim(.5,13.9)
plt.savefig('allsfrvsz.png')

plt.clf()
for i in toplotx:
    plt.plot(age,allsfrx[i]/allsfrx[i][-1],color='C0',alpha=.3)

for i in toplota:
    plt.plot(age,allsfra[i]/allsfra[i][-1],color='C1',alpha=.3)

for i in toplotl:
    plt.plot(age,allsfrl[i]/allsfrl[i][-1],color='C2',alpha=.3)

for i in toplotla:
    plt.plot(age,allsfrla[i]/allsfrla[i][-1],color='C3',alpha=.3)

z0sfrx=np.repeat(allsfrx[:,-1],100).reshape(len(allsfrx),100)
z0sfra=np.repeat(allsfra[:,-1],100).reshape(len(allsfra),100)
z0sfrl=np.repeat(allsfrl[:,-1],100).reshape(len(allsfrl),100)
z0sfrla=np.repeat(allsfrla[:,-1],100).reshape(len(allsfrla),100)

plt.plot(age,np.nanmedian(allsfrx/z0sfrx,axis=0),label='XMDs',lw=5,color='C0')
plt.plot(age,np.nanmedian(allsfra/z0sfra,axis=0),label='XMD Analogs',lw=5,color='C1')
plt.plot(age,np.nanmedian(allsfrl/z0sfrl,axis=0),label='LMDs',lw=5,color='C2')
plt.plot(age,np.nanmedian(allsfrla/z0sfrla,axis=0),label='LMD Ananlogs',lw=5,color='C3')
#plt.fill_between(age,np.nanmedian(allmetsx,axis=0)-np.nanstd(allmetsx,axis=0),np.nanmedian(allmetsx,axis=0)+np.nanstd(allmetsx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsa,axis=0)-np.nanstd(allmetsa,axis=0),np.nanmedian(allmetsa,axis=0)+np.nanstd(allmetsa,axis=0),color='C1',alpha=.3)
plt.legend()
plt.yscale('log')
#plt.ylim(5E-4,5E-3)
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.xlim(.5,13.9)
plt.ylabel('SFR/SFR$(z=0)$')
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlabel(r'$z$')
ax.set_xlim(.5,13.9)
plt.savefig('allsfrnormvsz.png')

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
plt.xlim(.5,13.9)
#plt.ylim(5E-4,5E-3)
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlabel(r'$z$')
ax.set_xlim(.5,13.9)
plt.savefig('allmstz.png')

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
plt.xlim(.5,13.9)
#plt.ylim(5E-4,5E-3)
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlabel(r'$z$')
ax.set_xlim(.5,13.9)
plt.savefig('allmhz.png')
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
plt.savefig('allhmass.png')

plt.clf()
for i in toplotx:
    plt.plot(age,allmgx[i]/allsfrx[i],color='C0',alpha=.3)

for i in toplota:
    plt.plot(age,allmga[i]/allsfra[i],color='C1',alpha=.3)

for i in toplotl:
    plt.plot(age,allmgl[i]/allsfrl[i],color='C2',alpha=.3)

for i in toplotla:
    plt.plot(age,allmgla[i]/allsfrla[i],color='C3',alpha=.3)

plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.yscale('log')
plt.ylabel(r'$M_g/{\rm SFR} ({\rm Gyr})$')
plt.xlim(0,13.9)
#plt.ylim(5E-4,5E-3)
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlabel(r'$z$')
plt.savefig('ind_mgsfr.png')

plt.clf()

for i in toplotx[0:1]:
    plt.plot(age,allmgx[i],color='C0',alpha=.3)

for i in toplota[0:1]:
    plt.plot(age,allmga[i],color='C1',alpha=.3)

for i in toplotl[0:1]:
    plt.plot(age,allmgl[i],color='C2',alpha=.3)

for i in toplotla[0:1]:
    plt.plot(age,allmgla[i],color='C3',alpha=.3)


plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.yscale('log')
plt.ylabel(r'$M_g/{\rm SFR} ({\rm Gyr})$')
plt.xlim(0,13.9)
#plt.ylim(5E-4,5E-3)
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlabel(r'$z$')

ax2=plt.twinx()
ax2.set_yscale('log')
for i in toplotx[0:1]:
    ax2.plot(age,allsfrx[i],color='C0',ls='--',alpha=.3)

for i in toplota[0:1]:
    ax2.plot(age,allsfra[i],color='C1',alpha=.3,ls='--')

for i in toplotl[0:1]:
    ax2.plot(age,allsfrl[i],color='C2',alpha=.3,ls='--')

for i in toplotla[0:1]:
    ax2.plot(age,allsfrla[i],color='C3',alpha=.3,ls='--')


plt.savefig('ind_mgsfrb.png')


plt.clf()
plt.plot(age,np.nanmedian(allmgx/allsfrx,axis=0),label='XMDs',lw=5,color='C0')
plt.plot(age,np.nanmedian(allmga/allsfra,axis=0),label='XMD Analogs',lw=5,color='C1')
plt.plot(age,np.nanmedian(allmgl/allsfrl,axis=0),label='LMDs',lw=5,color='C2')
plt.plot(age,np.nanmedian(allmgla/allsfrla,axis=0),label='LMD Analogs',lw=5,color='C3')
#plt.fill_between(age,np.nanmedian(allmetsx,axis=0)-np.nanstd(allmetsx,axis=0),np.nanmedian(allmetsx,axis=0)+np.nanstd(allmetsx,axis=0),color='C0',alpha=.3)
#plt.fill_between(age,np.nanmedian(allmetsa,axis=0)-np.nanstd(allmetsa,axis=0),np.nanmedian(allmetsa,axis=0)+np.nanstd(allmetsa,axis=0),color='C1',alpha=.3)
plt.legend()
plt.xlabel(r'$t_{\rm age}$ (Gyr)')
plt.yscale('log')
plt.ylabel(r'$M_g/{\rm SFR} ({\rm Gyr})$')
plt.xlim(0,13.9)
#plt.ylim(5E-4,5E-3)
ax=plt.twiny()
ax.set_xticks(astrocosmo.age(zg[::-1]).value)
ax.set_xticklabels(zglabels[::-1])
ax.set_xlabel(r'$z$')
plt.savefig('allmgsfr.png')

wz1=np.where(age.value>6)[0]
wzp2=np.where(age.value>10)[0]
plt.clf()
print(np.array(allmgx[:,wz1]/allsfrx[:,wz1]).flatten())
plt.hist(np.log10(np.array(allmgx[:,wz1]/allsfrx[:,wz1]).flatten()),bins=30,density=1,alpha=.5,range=[9,14])
plt.hist(np.log10(np.array(allmga[:,wz1]/allsfra[:,wz1]).flatten()),bins=30,density=1,alpha=.5,range=[9,14])
plt.xlabel(r'SFE (Gyr)')
plt.ylabel(r'dN/dSFE')
plt.savefig('sfehistz1.png')

plt.clf()
plt.hist(np.log10(np.array(allmgx[:,wzp2]/allsfrx[:,wzp2]).flatten()),bins=30,density=1,alpha=.5,range=[9,14])
plt.hist(np.log10(np.array(allmga[:,wzp2]/allsfra[:,wzp2]).flatten()),bins=30,density=1,alpha=.5,range=[9,14])
plt.xlabel(r'SFE (Gyr)')
plt.ylabel(r'dN/dSFE')
plt.savefig('sfehistzp2.png')



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
plt.savefig('allhgrowth.png')

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
plt.savefig('allgasmass.png')

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
plt.savefig('allgasgrowth.png')


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
plt.savefig('allstarmass.png')

plt.clf()
z0mstx=np.repeat(allmstx[:,-1],100).reshape(len(allmstx),100)
z0msta=np.repeat(allmsta[:,-1],100).reshape(len(allmsta),100)
z0mstl=np.repeat(allmstl[:,-1],100).reshape(len(allmstl),100)
z0mstla=np.repeat(allmstla[:,-1],100).reshape(len(allmstla),100)

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
plt.savefig('allstargrowth.png')
