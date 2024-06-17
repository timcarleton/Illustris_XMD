import numpy as np
import urllib.request
import shutil
import os
import glob
import h5py
from urllib.request import Request, urlopen

obj=np.loadtxt('sampleselection_tng50_sfr5_notide_nomass.txt',delimiter=',',skiprows=1)

wtry=np.where(obj[:,2]==4)[0]
#for i in range(len(obj)):
for i in wtry:
#for i in wtry[0:20]:
#for i in range(209,214,1):
    print(obj[i])
    g=glob.glob('/Volumes/Elements/illustris/cutouts50/'+str(int(obj[i,1]))+'/*')

    #print(g)
    #os.system('mkdir cutouts50/'+str(int(obj[i,1])))
    #g=glob.glob('cutouts50/'+str(int(obj[i,1]))+'/*')

    print(i)
    print('/Volumes/Elements/illustris/tng50trees/'+str(int(obj[i,1]))+'_tree.hdf5')
    h=h5py.File('/Volumes/Elements/illustris/tng50trees/'+str(int(obj[i,1]))+'_tree.hdf5','r')

    snap=100
    snapmin=np.max([60,np.min(h['SnapNum'][:])+1])
    j=0

    while h['SnapNum'][j]<snap and snap>snapmin:
        if '/Volumes/Elements/illustris/cutouts50/{:d}/{:d}_{:d}.hdf5'.format(int(obj[i,1]),int(obj[i,1]),h['SnapNum'][j]) in g:
            snap=h['SnapNum'][j]
            j=j+1
            continue
        url='http://www.tng-project.org/api/TNG50-1/snapshots/{:d}/subhalos/{:d}/cutout.hdf5'.format(int(h['SnapNum'][j]),int(h['SubfindID'][j]))

        print(url)
        savename='/Volumes/Elements/illustris/cutouts50/{:d}/{:d}_{:d}.hdf5'.format(int(obj[i,1]),int(obj[i,1]),h['SnapNum'][j])
        req=Request(url)
        print('a',i,j)
        req.add_header('API-key',os.environ['illustriskey'])
        try:
            with urlopen(req) as response:
                with open(savename, 'wb') as out_file:
                    shutil.copyfileobj(response, out_file)
        except:
            with urlopen(req) as response:
                with open(savename, 'wb') as out_file:
                    shutil.copyfileobj(response, out_file)

        snap=h['SnapNum'][j]
        j=j+1
        
