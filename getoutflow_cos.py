import numpy as np
import h5py
import getbryan98
from astropy import cosmology,units,constants
from scipy.optimize import minimize_scalar
import glob
from astropy import table
from astropy.io import fits

astrocosmo=cosmology.Planck15

hconst=0.6774

zs=np.array([20.05,14.99,11.98,10.976,10.0,9.389,9.002,8.449,8.012,7.595,7.236,7.005,6.491,6.011,5.847,5.530,5.228,4.996,4.665,4.428,4.177,4.008,3.709,3.491,3.283,3.008,2.896,2.733,2.578,2.444,2.316,2.208,2.103,2.002,1.904,1.823,1.743,1.667,1.604,1.531,1.496,1.414,1.358,1.302,1.248,1.206,1.155,1.114,1.074,1.036,0.997,0.951,0.923,0.887,0.851,0.817,0.791,0.757,0.733,0.700,0.676,0.644,0.621,0.598,0.576,0.546,0.525,0.503,0.482,0.461,0.440,0.420,0.40,0.380,0.361,0.348,0.329,0.310,0.298,0.273,0.261,0.244,0.226,0.214,0.1973,0.1804,0.1693,0.1527,0.1419,0.1258,0.1099,0.0994,0.0839,0.0737,0.0585,0.0485,0.0337,0.0240,0.0095,0])


obj=np.loadtxt('sampleselection_tng50_sfr5_notide_nomass.txt',delimiter=',',skiprows=1)

f=fits.open('tng50/alltng50gal_withr_withspin.fits')

def massr(rv,masses,dists):
    w=np.where(dists<rv)[0]
    return np.sum(masses[w])

def getoutflowcos(typei):

    wcat=np.where(obj[:,2]==typei)[0]
    for i in wcat:

        print(i,obj[i,1])


        g=glob.glob('/Volumes/Elements/illustris/cutouts50/'+str(int(obj[i,1]))+'/*hdf5')
        snaps=[int(k.split('_')[1].split('.')[0]) for k in g]

        ssnaps=np.sort(snaps)
    
        prevpartidvir=[]
        prevpartid2rh=[]
        prevpartid5rh=[]
        
        prevmassvir=[]
        prevmass2rh=[]
        prevmass5rh=[]

        prevzvir=[]
        prevz2rh=[]
        prevz5rh=[]    

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

        usesnap=[]
        allmstar=[]
        allmvir=[]
        allrh=[]
        allsfr=[]
        allmetal=[]
        for j in ssnaps:
            #print(j)
            usesnap.append(j)
            z=zs[j]
            #ht=h5py.File('trees/'+str(obj[i,1])+'_tree.hdf5')
        

            h=h5py.File('/Volumes/Elements/illustris/cutouts50/'+str(int(obj[i,1]))+'/'+str(int(obj[i,1]))+'_'+str(j)+'.hdf5','r')

            if 'PartType1' not in h.keys() or 'PartType0' not in h.keys():
                allmvir.append(np.nan)
                allrh.append(np.nan)
                allsfr.append(np.nan)
                allmetal.append(np.nan)
                allmstar.append(np.nan)
                moutflowvir.append(np.nan)
                zoutflowvir.append(np.nan)
                minflowvir.append(np.nan)
                zinflowvir.append(np.nan)
                
                moutflow5rh.append(np.nan)
                zoutflow5rh.append(np.nan)
                minflow5rh.append(np.nan)
                zinflow5rh.append(np.nan)

                moutflow2rh.append(np.nan)
                zoutflow2rh.append(np.nan)
                minflow2rh.append(np.nan)
                zinflow2rh.append(np.nan)
                continue
        
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
            allsfr.append(np.sum(h['PartType0']['StarFormationRate']))
            allmetal.append(np.mean(np.array(h['PartType0']['StarFormationRate'])*np.array(h['PartType0']['GFM_Metallicity']))/np.mean(h['PartType0']['StarFormationRate']))
            
            wwithinvir=np.where(distgas<rvir)[0]
            wwithin5rh=np.where(distgas<5*rhalf)[0]
            wwithin2rh=np.where(distgas<2*rhalf)[0]

            partidivir=np.array(h['PartType0']['ParticleIDs'][wwithinvir]).tolist()
            partidi5rh=np.array(h['PartType0']['ParticleIDs'][wwithin5rh]).tolist()
            partidi2rh=np.array(h['PartType0']['ParticleIDs'][wwithin2rh]).tolist()

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
            
            for mp in idinflowvir:
                winf=partidivir.index(mp)
                minflowi_vir+=h['PartType0']['Masses'][wwithinvir[winf]]/hconst*1E10
                zinflowi_vir+=h['PartType0']['GFM_Metallicity'][wwithinvir[winf]]*h['PartType0']['Masses'][wwithinvir[winf]]/hconst*1E10

                
            for mp in idinflow5rh:
                winf=partidi5rh.index(mp)
                minflowi_5rh+=h['PartType0']['Masses'][wwithin5rh[winf]]/hconst*1E10
                zinflowi_5rh+=h['PartType0']['GFM_Metallicity'][wwithin5rh[winf]]*h['PartType0']['Masses'][wwithin5rh[winf]]/hconst*1E10

                
            for mp in idinflow2rh:
                winf=partidi2rh.index(mp)
                minflowi_2rh+=h['PartType0']['Masses'][wwithin2rh[winf]]/hconst*1E10
                zinflowi_2rh+=h['PartType0']['GFM_Metallicity'][wwithin2rh[winf]]*h['PartType0']['Masses'][wwithin2rh[winf]]/hconst*1E10

                
            #for mp in range(len(wwithinvir)):
            #    if h['PartType0']['ParticleIDs'][wwithinvir[mp]] not in prevpartid:
            #        minflowi+=h['PartType0']['Masses'][wwithinvir[mp]]/hconst*1E10
            #        zinflowi+=h['PartType0']['GFM_Metallicity'][wwithinvir[mp]]*h['PartType0']['Masses'][wwithinvir[mp]]/hconst*1E10

            #for mp in range(len(prevpartid)):
            #    if prevpartid[mp] not in h['PartType0']['ParticleIDs'][wwithinvir]:
            #        moutflowi+=prevmass[mp]
            #        zoutflowi+=prevz[mp]

        
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
        
            #vvir=np.sqrt(constants.G*mvir*units.M_sun/(rvir*units.kpc)).to(units.km/units.s)
            
            #comvelx=np.mean(h['PartType1']['Velocities'][:,0])
            #comvely=np.mean(h['PartType1']['Velocities'][:,1])
            #comvelz=np.mean(h['PartType1']['Velocities'][:,2])
            #gasvel=np.sqrt((h['PartType0']['Velocities'][:,0]-comvelx)**2+
            #               (h['PartType0']['Velocities'][:,1]-comvely)**2+
            #               (h['PartType0']['Velocities'][:,2]-comvelz)**2)
            
            #gasesc=[]
            #for gs in range(len(distgas)):
            
            #    gasescg=np.sqrt(constants.G*massr(distgas[gs],allmass,alldist)*units.M_sun/distgas[gs]/units.kpc).to(units.km/units.s)
            #    gasesc.append(gasvel[gs]/gasescg.value)

            prevpartidvir=prevpartidvir+np.array(h['PartType0']['ParticleIDs'][wwithinvir]).tolist()
            #prevmassvir=prevmassvir+list(np.array(h['PartType0']['Masses'][wwithinvir])/hconst*1E10)
            #prevzvir=prevzvir+list(np.array(h['PartType0']['GFM_Metallicity'][wwithinvir])*np.array(h['PartType0']['Masses'][wwithinvir])/hconst*1E10)

            prevpartid5rh=prevpartid5rh+np.array(h['PartType0']['ParticleIDs'][wwithin5rh]).tolist()
            #prevmass5rh=prevmass5rh+list(np.array(h['PartType0']['Masses'][wwithin5rh])/hconst*1E10)
            #prevz5rh=prevz5rh+list(np.array(h['PartType0']['GFM_Metallicity'][wwithin5rh])*np.array(h['PartType0']['Masses'][wwithin5rh])/hconst)
            
            prevpartid2rh=prevpartid2rh+np.array(h['PartType0']['ParticleIDs'][wwithin2rh]).tolist()
            #prevmass2rh=prevmass2rh+list(np.array(h['PartType0']['Masses'][wwithin2rh])/hconst*1E10)
            #prevz2rh=prevz2rh+list(np.array(h['PartType0']['GFM_Metallicity'][wwithin2rh])*np.array(h['PartType0']['Masses'][wwithin2rh])/hconst*1E10)

        print('l')
        for obb in [usesnap,zs[np.array(usesnap)],allmvir,allmstar,allrh,allsfr,allmetal,minflowvir,zinflowvir,minflow5rh,zinflow5rh,minflow2rh,zinflow2rh]:
            print(len(obb))
        tab=table.Table([usesnap,zs[np.array(usesnap)],allmvir,allmstar,allrh,allsfr,allmetal,minflowvir,zinflowvir,minflow5rh,zinflow5rh,minflow2rh,zinflow2rh],names=['snapshot','redshift','mvir','mstar','rhstar','sfr','sfr_weight_metallicity','mcosinflow_vir','mzcosinflow_vir','mcosinflow_5rh','mzcosinflow_5rh','mcosinflow_2rh','mzcosinflow_2rh'])
        tab.write('flowinfo/flowinfo_cos_'+str(int(obj[i,1]))+'.txt',format='csv',overwrite=True)
            #np.savetxt('flowinfo/flowinfo_'+str(int(obj[i,1]))+'.txt',np)
            #        mvir=np.sum(h['PartType0']['Masses'][:])*1E10+len(h['PartType1']['ParticleIDs'])*5.1E6+np.sum(h['PartType4']['Masses'][:])*1E10

        
        
            #       rvir=((mvir/(4/3*np.pi*dvir*astrocosmo.critical_density(z).to(units.M_sun/units.kpc**3)))**(1/3)).value
        
