import numpy as np
import os
import h5py
import getbryan98
from astropy import cosmology,units,constants
from scipy.optimize import minimize_scalar
import glob
from astropy import table
from astropy.io import fits
import downloadcutout_snapid50
astrocosmo=cosmology.Planck15

hconst=0.6774

zs=np.array([20.05,14.99,11.98,10.976,10.0,9.389,9.002,8.449,8.012,7.595,7.236,7.005,6.491,6.011,5.847,5.530,5.228,4.996,4.665,4.428,4.177,4.008,3.709,3.491,3.283,3.008,2.896,2.733,2.578,2.444,2.316,2.208,2.103,2.002,1.904,1.823,1.743,1.667,1.604,1.531,1.496,1.414,1.358,1.302,1.248,1.206,1.155,1.114,1.074,1.036,0.997,0.951,0.923,0.887,0.851,0.817,0.791,0.757,0.733,0.700,0.676,0.644,0.621,0.598,0.576,0.546,0.525,0.503,0.482,0.461,0.440,0.420,0.40,0.380,0.361,0.348,0.329,0.310,0.298,0.273,0.261,0.244,0.226,0.214,0.1973,0.1804,0.1693,0.1527,0.1419,0.1258,0.1099,0.0994,0.0839,0.0737,0.0585,0.0485,0.0337,0.0240,0.0095,0])
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

obj=np.loadtxt('sampleselection_tng50_sfr5_notide_nomass.txt',delimiter=',',skiprows=1)

f=fits.open('tng50/alltng50gal_withr_withspin.fits')
def massr(rv,masses,dists):
    w=np.where(dists<rv)[0]
    return np.sum(masses[w])


def getoutflow(typei):
    wcat=np.where(obj[:,2]==typei)[0]
    #for i in range(10):
    #    for j in range(99):
    
    #for i in range(117,137):
    #for i in range(0,6):
    #for i in range(175,192):
    #for i in range(1722,1742):
    for i in wcat:
        if f[1].data.mstar[int(obj[i,0])]>3E8:
            continue
        print(i,obj[i,1])
        g=glob.glob('/Volumes/Elements/illustris/cutouts50/'+str(int(obj[i,1]))+'/*hdf5')
        g99=glob.glob('/Volumes/Elements/illustris/cutouts50/'+str(int(obj[i,1]))+'/*99.hdf5')
        g98=glob.glob('/Volumes/Elements/illustris/cutouts50/'+str(int(obj[i,1]))+'/*98.hdf5')
        if len(g)<2:
            os.system('mkdir /Volumes/Elements/illustris/cutouts50/'+str(int(obj[i,1]))+'/')
            if len(g99)==0:
                downloadcutout_snapid50.downloadcutout_snap(int(obj[i,1]),99)
                os.system('mv '+str(int(obj[i,1]))+'_99.hdf5 /Volumes/Elements/illustris/cutouts50/'+str(int(obj[i,1]))+'/')
            if len(g98)==0:
                downloadcutout_snapid50.downloadcutout_snap(int(obj[i,1]),98)
                os.system('mv '+str(int(obj[i,1]))+'_98.hdf5 /Volumes/Elements/illustris/cutouts50/'+str(int(obj[i,1]))+'/')
            g=glob.glob('/Volumes/Elements/illustris/cutouts50/'+str(int(obj[i,1]))+'/*hdf5')
            if len(g)==0:
                print('no files')
                continue

        snaps=[int(k.split('_')[1].split('.')[0]) for k in g]

        ssnaps=np.sort(snaps)

        prevpartidisf=[]
        prevpartidvir=[]
        prevpartid2rh=[]
        prevpartid5rh=[]
        
        moutflowvir=[]
        zoutflowvir=[]
        moutflow2rh=[]
        zoutflow2rh=[]
        moutflow5rh=[]
        zoutflow5rh=[]
    
        minflowvir=[]
        zinflowvir=[]
        minflow2rh=[]
        zinflow2rh=[]
        minflow5rh=[]
        zinflow5rh=[]

        mdiffusea=[]
        zdiffusea=[]
        mcondensea=[]
        zcondensea=[]

        usesnap=[]
        allmstar=[]
        allmvir=[]
        allrh=[]
        allsfr=[]
        allmetal=[]
        allmetalmass=[]
        allmsfgas=[]
        allmgasrh=[]
        allmgastot=[]
        alldensgas=[]
        print(ssnaps)
        for j in ssnaps:
            print(j)
            #ht=h5py.File('trees/'+str(obj[i,1])+'_tree.hdf5')


            print('/Volumes/Elements/illustris/cutouts50/'+str(int(obj[i,1]))+'/'+str(int(obj[i,1]))+'_'+str(j)+'.hdf5')
            h=h5py.File('/Volumes/Elements/illustris/cutouts50/'+str(int(obj[i,1]))+'/'+str(int(obj[i,1]))+'_'+str(j)+'.hdf5','r')

            if 'PartType1' not in h.keys() or 'PartType0' not in h.keys():
                continue
            usesnap.append(j)
            z=zs[j]

        
            centerx=np.mean(h['PartType1']['Coordinates'][:,0])
            centery=np.mean(h['PartType1']['Coordinates'][:,1])
            centerz=np.mean(h['PartType1']['Coordinates'][:,2])
            
            distgas=np.sqrt((h['PartType0']['Coordinates'][:,0]-centerx)**2+
                            (h['PartType0']['Coordinates'][:,1]-centery)**2+
                            (h['PartType0']['Coordinates'][:,2]-centerz)**2)

            distdm=np.sqrt((h['PartType1']['Coordinates'][:,0]-centerx)**2+
                           (h['PartType1']['Coordinates'][:,1]-centery)**2+
                           (h['PartType1']['Coordinates'][:,2]-centerz)**2)

            if 'PartType4' in h.keys():
                distst=np.sqrt((h['PartType4']['Coordinates'][:,0]-centerx)**2+
                               (h['PartType4']['Coordinates'][:,1]-centery)**2+
                               (h['PartType4']['Coordinates'][:,2]-centerz)**2)
            else:
                distst=[]
            
            alldist=np.hstack([distdm,distgas,distst])/hconst/(1+z)

            if 'PartType4' in h.keys():
                allmass=np.hstack([np.array([5.1E6]*len(h['PartType1']['ParticleIDs'])),h['PartType0']['Masses'][:]*1E10,h['PartType4']['Masses'][:]*1E10])
                allmstar.append(np.sum(h['PartType4']['Masses'][:]*1E10/hconst))
            else:
                allmass=np.hstack([np.array([5.1E6]*len(h['PartType1']['ParticleIDs'])),h['PartType0']['Masses'][:]*1E10])
                allmstar.append(0)
        
            dvir=getbryan98.getbryan98(astrocosmo,z)*astrocosmo.critical_density(z).to(units.M_sun/units.kpc**3).value

        
            #print(alldist,dvir)
            #print(massr(2,allmass,alldist),(4/3*np.pi*dvir*2**3),massr(50,allmass,alldist)/(4/3*np.pi*dvir*50**3),massr(100,allmass,alldist)/(4/3*np.pi*dvir*100**3))

            vir=minimize_scalar(lambda rv:(massr(rv,allmass,alldist)/hconst-(4/3*np.pi*dvir*rv**3))**2,bounds=[.5,200],method='bounded')

            rvir=vir['x']
            mvir=massr(vir['x'],allmass,alldist)
            #print(vir)
            rhalf=np.median(distst)
            allmvir.append(mvir)
            allrh.append(rhalf)
            wsfg=np.where(np.array(h['PartType0']['StarFormationRate'])>0)[0]
            wdensgas=np.where(np.array(h['PartType0']['Density'])>10**-3.25)[0]
            alldensgas.append(np.sum(np.array(h['PartType0']['Masses'])[wdensgas])/hconst*1E10)
            allmsfgas.append(np.sum(np.array(h['PartType0']['Masses'])[wsfg])/hconst*1E10)
            allmgastot.append(np.sum(np.array(h['PartType0']['Masses'][:])/hconst*1E10))
        
            allsfr.append(np.sum(h['PartType0']['StarFormationRate']))
            allmetal.append(np.mean(np.array(h['PartType0']['StarFormationRate'])*np.array(h['PartType0']['GFM_Metallicity'][:]))/np.mean(h['PartType0']['StarFormationRate']))
            allmetalmass.append(np.mean(np.array(h['PartType0']['GFM_Metallicity'])))
            wwithinrh=np.where(distgas<rhalf)[0]
            allmgasrh.append(np.sum(np.array(h['PartType0']['Masses'])[wwithinrh]))
        
            
            wwithinvir=np.where(distgas<rvir)[0]
            wwithin5rh=np.where(distgas<5*rhalf)[0]
            wwithin2rh=np.where(distgas<2*rhalf)[0]

            partidisf=np.array(h['PartType0']['ParticleIDs'][wsfg]).tolist()
            partidivir=np.array(h['PartType0']['ParticleIDs'][wwithinvir]).tolist()
            partidi5rh=np.array(h['PartType0']['ParticleIDs'][wwithin5rh]).tolist()
            partidi2rh=np.array(h['PartType0']['ParticleIDs'][wwithin2rh]).tolist()

            mdiffuse=0
            zdiffuse=0
            mcondense=0
            zcondense=0
            
            minflowi_vir=0
            zinflowi_vir=0
            moutflowi_vir=0
            zoutflowi_vir=0

            minflowi_5rh=0
            zinflowi_5rh=0
            moutflowi_5rh=0
            zoutflowi_5rh=0

            minflowi_2rh=0
            zinflowi_2rh=0
            moutflowi_2rh=0
            zoutflowi_2rh=0

            iddiffuse = set(prevpartidisf) - set(partidisf)
            idcondense = set(partidisf) - set(prevpartidisf)
            
            partidisf = list(partidisf)
            prevpartidisf = list(prevpartidisf)
            
            idoutflowvir=set(prevpartidvir) - set(partidivir)
            idinflowvir=set(partidivir) - set(prevpartidvir)
            
            partidivir=list(partidivir)
            prevpartidvir=list(prevpartidvir)
            
            idoutflow5rh=set(prevpartid5rh) - set(partidi5rh)
            idinflow5rh=set(partidi5rh) - set(prevpartid5rh)

            partidi5rh=list(partidi5rh)
            prevpartid5rh=list(prevpartid5rh)

            idoutflow2rh=set(prevpartid2rh) - set(partidi2rh)
            idinflow2rh=set(partidi2rh) - set(prevpartid2rh)

            partidi2rh=list(partidi2rh)
            prevpartid2rh=list(prevpartid2rh)

            for mp in iddiffuse:
                wdiffuse=prevpartidisf.index(mp)
                mdiffuse+=prevmasssf[wdiffuse]
                zdiffuse+=prevzsf[wdiffuse]

            for mp in idcondense:
                wcondense=partidisf.index(mp)
                mcondense+=h['PartType0']['Masses'][wsfg[wcondense]]/hconst*1E10
                zcondense+=h['PartType0']['GFM_Metallicity'][wsfg[wcondense]]*h['PartType0']['Masses'][wsfg[wcondense]]/hconst*1E10

            
            for mp in idoutflowvir:
                woutf=prevpartidvir.index(mp)
                moutflowi_vir+=prevmassvir[woutf]
                zoutflowi_vir+=prevzvir[woutf]
            
            for mp in idinflowvir:
                winf=partidivir.index(mp)
                minflowi_vir+=h['PartType0']['Masses'][wwithinvir[winf]]/hconst*1E10
                zinflowi_vir+=h['PartType0']['GFM_Metallicity'][wwithinvir[winf]]*h['PartType0']['Masses'][wwithinvir[winf]]/hconst*1E10

            for mp in idoutflow5rh:
                woutf=prevpartid5rh.index(mp)
                moutflowi_5rh+=prevmass5rh[woutf]
                
                zoutflowi_5rh+=prevz5rh[woutf]
            
            for mp in idinflow5rh:
                winf=partidi5rh.index(mp)
                minflowi_5rh+=h['PartType0']['Masses'][wwithin5rh[winf]]/hconst*1E10
                zinflowi_5rh+=h['PartType0']['GFM_Metallicity'][wwithin5rh[winf]]*h['PartType0']['Masses'][wwithin5rh[winf]]/hconst*1E10

            for mp in idoutflow2rh:
                woutf=prevpartid2rh.index(mp)
                moutflowi_2rh+=prevmass2rh[woutf]
                zoutflowi_2rh+=prevz2rh[woutf]
            
            for mp in idinflow2rh:
                winf=partidi2rh.index(mp)
                minflowi_2rh+=h['PartType0']['Masses'][wwithin2rh[winf]]/hconst*1E10
                zinflowi_2rh+=h['PartType0']['GFM_Metallicity'][wwithin2rh[winf]]*h['PartType0']['Masses'][wwithin2rh[winf]]/hconst*1E10

        
            mdiffusea.append(mdiffuse)
            zdiffusea.append(zdiffuse)
            mcondensea.append(mcondense)
            zcondensea.append(zcondense)

            moutflowvir.append(moutflowi_vir)
            zoutflowvir.append(zoutflowi_vir)
            minflowvir.append(minflowi_vir)
            zinflowvir.append(zinflowi_vir)

            moutflow5rh.append(moutflowi_5rh)
            zoutflow5rh.append(zoutflowi_5rh)
            minflow5rh.append(minflowi_5rh)
            zinflow5rh.append(zinflowi_5rh)

            moutflow2rh.append(moutflowi_2rh)
            zoutflow2rh.append(zoutflowi_2rh)
            minflow2rh.append(minflowi_2rh)
            zinflow2rh.append(zinflowi_2rh)

            prevpartidisf=np.array(h['PartType0']['ParticleIDs'][wsfg]).tolist()
            prevmasssf=np.array(h['PartType0']['Masses'][wsfg])/hconst*1E10
            prevzsf=np.array(h['PartType0']['GFM_Metallicity'][wsfg])*prevmasssf

            prevpartidvir=np.array(h['PartType0']['ParticleIDs'][wwithinvir]).tolist()
            prevmassvir=np.array(h['PartType0']['Masses'][wwithinvir])/hconst*1E10
            prevzvir=np.array(h['PartType0']['GFM_Metallicity'][wwithinvir])*prevmassvir

            prevpartid5rh=np.array(h['PartType0']['ParticleIDs'][wwithin5rh]).tolist()
            prevmass5rh=np.array(h['PartType0']['Masses'][wwithin5rh])/hconst*1E10
            prevz5rh=np.array(h['PartType0']['GFM_Metallicity'][wwithin5rh])*prevmass5rh

            prevpartid2rh=np.array(h['PartType0']['ParticleIDs'][wwithin2rh]).tolist()
            prevmass2rh=np.array(h['PartType0']['Masses'][wwithin2rh])/hconst*1E10
            prevz2rh=np.array(h['PartType0']['GFM_Metallicity'][wwithin2rh])*prevmass2rh

        for ar in [usesnap,zs[np.array(usesnap)],allmvir,allmstar,allmsfgas,allmgasrh,allrh,allsfr,allmetal,allmetalmass,minflowvir,zinflowvir,moutflowvir,zoutflowvir,minflow5rh,zinflow5rh,moutflow5rh,zoutflow5rh,minflow2rh,zinflow2rh,moutflow2rh,zoutflow2rh,allmgastot,alldensgas]:
            print('l',len(ar))

        tab=table.Table([usesnap,zs[np.array(usesnap)],allmvir,allmstar,allmsfgas,allmgastot,alldensgas,allmgasrh,allrh,allsfr,allmetal,allmetalmass,minflowvir,zinflowvir,moutflowvir,zoutflowvir,minflow5rh,zinflow5rh,moutflow5rh,zoutflow5rh,minflow2rh,zinflow2rh,moutflow2rh,zoutflow2rh,mdiffusea,zdiffusea,mcondensea,zcondensea],names=['snapshot','redshift','mvir','mstar','m_sf_gas','mgas_tot','m_dens_gas','m_gas_rh','rhstar','sfr','sfr_weight_metallicity','mass_weight_metallicity','minflow_vir','mzinflow_vir','moutflow_vir','mzoutflow_vir','minflow_5rh','mzinflow_5rh','moutflow_5rh','mzoutflow_5rh','minflow_2rh','mzinflow_2rh','moutflow_2rh','mzoutflow_2rh','mendsf','zendsf','mstartsf','zstartsf'])
        tab.write('flowinfo/flowinfo_50_'+str(typei)+'_'+str(int(obj[i,1]))+'.txt',format='csv',overwrite=True)

        
