from astropy.io import fits, ascii
import glob
import numpy as np
import gettemp
import h5py
from scipy.optimize import minimize_scalar

f=fits.open('alltnggal_withr_wmetals_wvmax.fits')
obj=np.loadtxt('sampleselection_sfr5_notide_nomass.txt',skiprows=1,delimiter=',')

wg=np.where((f[1].data.mstar<1E12) & (f[1].data.mstar>3E7) & (f[1].data.sfr>0) & (f[1].data.half_mass_rad_star>0))[0]

#scale factor for z=0.05, 0.1, z=0.5 -> ages of .7 Gyr, 1.3 Gyr, 5Gyr (ages not right cosmology, but should be fine)
agecut1=0.95
agecut2=0.91
agecut3=.67

allm1=[]
allm2=[]
allm3=[]
allm4=[]
allr1=[]
allr2=[]
allr3=[]
allr4=[]
alldist1=[]
alldist2=[]
alldist3=[]
alldist4=[]
alllum1=[]
alllum2=[]
alllum3=[]
alllum4=[]
allrall=[]
allmall=[]
alllumall=[]
allavgdist1=[]
allavgdist2=[]
allavgdist3=[]
allavgdist4=[]

h=0.6774
def massr(rv,masses,dists):
    w=np.where(dists<rv)[0]
    return np.sum(masses[w])

for i in range(len(obj)):
#for i in range(5):
    g=glob.glob('/Volumes/Elements/illustris/cutouts/'+str(int(obj[i,1]))+'/*hdf5')
    if len(g)==0:
        print('no files')
        allm1.append(0)
        allm2.append(0)
        allm3.append(0)
        allm4.append(0)
        allr1.append(0)
        allr2.append(0)
        allr3.append(0)
        allr4.append(0)
        alldist1.append(0)
        alldist2.append(0)
        alldist3.append(0)
        alldist4.append(0)
        alllum1.append(0)
        alllum2.append(0)
        alllum3.append(0)
        alllum4.append(0)
        allrall.append(0)
        allmall.append(0)
        allumall.append(0)
        allavgdist1.append(0)
        allavgdist2.append(0)
        allavgdist3.append(0)
        allavgdist4.append(0)
        continue

    hf=h5py.File('/Volumes/Elements/illustris/cutouts/'+str(int(obj[i,1]))+'/'+str(int(obj[i,1]))+'_99.hdf5')

    if 'PartType4' not in hf:
        allm1.append(0)
        allm2.append(0)
        allm3.append(0)
        allm4.append(0)
        allr1.append(0)
        allr2.append(0)
        allr3.append(0)
        allr4.append(0)
        alldist1.append(0)
        alldist2.append(0)
        alldist3.append(0)
        alldist4.append(0)
        alllum1.append(0)
        alllum2.append(0)
        alllum3.append(0)
        alllum4.append(0)
        allmall.append(0)
        allrall.append(0)
        alllumall.append(0)
        allavgdist1.append(0)
        allavgdist2.append(0)
        allavgdist3.append(0)
        allavgdist4.append(0)

        continue

    wst=np.where(hf['PartType4']['GFM_StellarFormationTime'][:]>0)[0]
    ages=hf['PartType4']['GFM_StellarFormationTime'][wst]

    wage1=wst[np.where(ages>agecut1)[0]]
    wage2=wst[np.where((ages<agecut1) & (ages>agecut2))[0]]
    wage3=wst[np.where((ages<agecut2) & (ages>agecut3))[0]]
    wage4=wst[np.where((ages<agecut3) & (ages>0))[0]]

    allmst=np.array(hf['PartType4']['Masses'][:])[wst]/.6774*1E10
    allmall.append(np.sum(allmst))

    allm1i=np.array(hf['PartType4']['Masses'][:])[wage1]/.6774*1E10
    allm2i=np.array(hf['PartType4']['Masses'][:])[wage2]/.6774*1E10
    allm3i=np.array(hf['PartType4']['Masses'][:])[wage3]/.6774*1E10
    allm4i=np.array(hf['PartType4']['Masses'][:])[wage4]/.6774*1E10

    glum=10**(-0.4*(hf['PartType4']['GFM_StellarPhotometrics'][:,4]-5.14))
    
    allm1.append(np.sum(allm1i))
    allm2.append(np.sum(allm2i))
    allm3.append(np.sum(allm3i))
    allm4.append(np.sum(allm4i))

    alllum1.append(np.sum(glum[wage1]))
    alllum2.append(np.sum(glum[wage2]))
    alllum3.append(np.sum(glum[wage3]))
    alllum4.append(np.sum(glum[wage4]))
    alllumall.append(np.sum(glum[wst]))

    allx=np.array(hf['PartType4']['Coordinates'][:,0])[wst]/.6774
    ally=np.array(hf['PartType4']['Coordinates'][:,1])[wst]/.6774
    allz=np.array(hf['PartType4']['Coordinates'][:,2])[wst]/.6774
    
    x1=np.array(hf['PartType4']['Coordinates'][:,0])[wage1]/.6774
    y1=np.array(hf['PartType4']['Coordinates'][:,1])[wage1]/.6774
    z1=np.array(hf['PartType4']['Coordinates'][:,2])[wage1]/.6774

    x2=np.array(hf['PartType4']['Coordinates'][:,0])[wage2]/.6774
    y2=np.array(hf['PartType4']['Coordinates'][:,1])[wage2]/.6774
    z2=np.array(hf['PartType4']['Coordinates'][:,2])[wage2]/.6774

    x3=np.array(hf['PartType4']['Coordinates'][:,0])[wage3]/.6774
    y3=np.array(hf['PartType4']['Coordinates'][:,1])[wage3]/.6774
    z3=np.array(hf['PartType4']['Coordinates'][:,2])[wage3]/.6774

    x4=np.array(hf['PartType4']['Coordinates'][:,0])[wage4]/.6774
    y4=np.array(hf['PartType4']['Coordinates'][:,1])[wage4]/.6774
    z4=np.array(hf['PartType4']['Coordinates'][:,2])[wage4]/.6774


    centerx1b=np.mean(x1*allm1i)/np.mean(allm1i)
    centery1b=np.mean(y1*allm1i)/np.mean(allm1i)
    centerz1b=np.mean(z1*allm1i)/np.mean(allm1i)
        
    centerx2b=np.mean(x2*allm2i)/np.mean(allm2i)
    centery2b=np.mean(y2*allm2i)/np.mean(allm2i)
    centerz2b=np.mean(z2*allm2i)/np.mean(allm2i)

    centerx3b=np.mean(x3*allm3i)/np.mean(allm3i)
    centery3b=np.mean(y3*allm3i)/np.mean(allm3i)
    centerz3b=np.mean(z3*allm3i)/np.mean(allm3i)
    
    centerx4b=np.mean(x4*allm4i)/np.mean(allm4i)
    centery4b=np.mean(y4*allm4i)/np.mean(allm4i)
    centerz4b=np.mean(z4*allm4i)/np.mean(allm4i)

    centerall=True
    if centerall:
        centerx1=np.mean(allx*allmst)/np.mean(allmst)
        centery1=np.mean(ally*allmst)/np.mean(allmst)
        centerz1=np.mean(allz*allmst)/np.mean(allmst)

        centerx2=np.mean(allx*allmst)/np.mean(allmst)
        centery2=np.mean(ally*allmst)/np.mean(allmst)
        centerz2=np.mean(allz*allmst)/np.mean(allmst)

        centerx3=np.mean(allx*allmst)/np.mean(allmst)
        centery3=np.mean(ally*allmst)/np.mean(allmst)
        centerz3=np.mean(allz*allmst)/np.mean(allmst)

        centerx4=np.mean(allx*allmst)/np.mean(allmst)
        centery4=np.mean(ally*allmst)/np.mean(allmst)
        centerz4=np.mean(allz*allmst)/np.mean(allmst)
        
    else:
        centerx1=centerx1b
        centery1=centery1b
        centerz1=centerz1b

        centerx2=centerx2b
        centery2=centery2b
        centerz2=centerz2b

        centerx3=centerx3b
        centery3=centery3b
        centerz3=centerz3b

        centerx4=centerx4b
        centery4=centery4b
        centerz4=centerz4b


    #print(centerx4,centerx4b)
    alldist1.append(np.sqrt((centerx1-centerx1b)**2+(centery1-centery1b)**2+(centerz1-centerz1b)**2))
    alldist2.append(np.sqrt((centerx2-centerx2b)**2+(centery2-centery2b)**2+(centerz2-centerz2b)**2))
    alldist3.append(np.sqrt((centerx3-centerx3b)**2+(centery3-centery3b)**2+(centerz3-centerz3b)**2))
    alldist4.append(np.sqrt((centerx4-centerx4b)**2+(centery4-centery4b)**2+(centerz4-centerz4b)**2))

    allavgdist1.append(np.mean(np.sqrt((centerx1-x1*allm1i)**2+(centery1-y1*allm1i)**2+(centerz1-z1*allm1i)**2))/np.mean(allm1i))
    allavgdist2.append(np.mean(np.sqrt((centerx2-x2*allm2i)**2+(centery2-y2*allm2i)**2+(centerz2-z2*allm2i)**2))/np.mean(allm2i))
    allavgdist3.append(np.mean(np.sqrt((centerx3-x3*allm3i)**2+(centery3-y3*allm3i)**2+(centerz3-z3*allm3i)**2))/np.mean(allm3i))
    allavgdist4.append(np.mean(np.sqrt((centerx4-x4*allm4i)**2+(centery4-y4*allm4i)**2+(centerz4-z4*allm4i)**2))/np.mean(allm4i))
    
    dista=np.sqrt((allx-centerx1)**2+(ally-centery1)**2+(allz-centerz1)**2)
    dist1=np.sqrt((x1-centerx1)**2+(y1-centery1)**2+(z1-centerz1)**2)
    dist2=np.sqrt((x2-centerx2)**2+(y2-centery2)**2+(z2-centerz2)**2)
    dist3=np.sqrt((x3-centerx3)**2+(y3-centery3)**2+(z3-centerz3)**2)
    dist4=np.sqrt((x4-centerx4)**2+(y4-centery4)**2+(z4-centerz4)**2)


    m12iall=minimize_scalar(lambda r12ia: abs(massr(r12ia,allmst,dista)-0.5*np.sum(allmst)),bounds=[0,2*np.median(dista)],method='bounded')
    m12i1=minimize_scalar(lambda r12i1: abs(massr(r12i1,allm1i,dist1)-0.5*np.sum(allm1i)),bounds=[0,2*np.median(dist1)],method='bounded')
    m12i2=minimize_scalar(lambda r12i2: abs(massr(r12i2,allm2i,dist2)-0.5*np.sum(allm2i)),bounds=[0,2*np.median(dist2)],method='bounded')
    m12i3=minimize_scalar(lambda r12i3: abs(massr(r12i3,allm3i,dist3)-0.5*np.sum(allm3i)),bounds=[0,2*np.median(dist3)],method='bounded')
    m12i4=minimize_scalar(lambda r12i4: abs(massr(r12i4,allm4i,dist4)-0.5*np.sum(allm4i)),bounds=[0,2*np.median(dist4)],method='bounded')
        
    
    allrall.append(m12iall['x'])
    allr1.append(m12i1['x'])
    allr2.append(m12i2['x'])
    allr3.append(m12i3['x'])
    allr4.append(m12i4['x'])

#print(allr1,allr2,allr3,allm1,allm2,allm3)
from astropy import table
#tab=table.Table([obj[:,0],obj[:,1],obj[:,2],allm1,allm2,allm3,allm4,allr1,allr2,allr3,allr4,alllum1,alllum2,alllum3,alllum4,alldist1,alldist2,alldist3,alldist4,allmall,allrall,alllumall],names=['allgalindex','haloid','type','myoung','mmid','mmid2','mold','ryoung','rmid','rmid2','rold','lumyoung','lummid','lummid2','lumold','dist1','dist2','dist3','dist4','mall','rall','lumall'])
tab=table.Table([obj[:,0],obj[:,1],obj[:,2],allm1,allm2,allm3,allm4,allr1,allr2,allr3,allr4,alllum1,alllum2,alllum3,alllum4,allavgdist1,allavgdist2,allavgdist3,allavgdist4,allmall,allrall,alllumall],names=['allgalindex','haloid','type','myoung','mmid','mmid2','mold','ryoung','rmid','rmid2','rold','lumyoung','lummid','lummid2','lumold','dist1','dist2','dist3','dist4','mall','rall','lumall'])
tab.write('mstagesizeinfo_ca.txt',format='csv',overwrite=True)
