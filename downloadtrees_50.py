import numpy as np
import urllib.request
import shutil
import os
import glob
from urllib.request import Request, urlopen
from urllib.error import HTTPError

#obj=np.loadtxt('sampleselection_z55_sfr5.txt',delimiter=',',skiprows=1)
#obj=np.loadtxt('gcclust.txt')
obj=np.loadtxt('sampleselection_tng50_sfr5_notide_nomass.txt',delimiter=',',skiprows=1)
g=glob.glob('/Volumes/Elements/illustris/tng50trees/*')
for i in range(len(obj)):
#for i in [3089]:
#for i in [6]:
    
    if '/Volumes/Elements/illustris/tng50trees/{:d}_tree.hdf5'.format(int(obj[i,1])) not in g:

        url='http://www.tng-project.org/api/TNG50-1/snapshots/99/subhalos/{:d}/sublink/mpb.hdf5'.format(int(obj[i,1]))
        savename='{:d}_tree.hdf5'.format(int(obj[i,1]))
        req=Request(url)
        print(i,url)
        req.add_header('API-key',os.environ['illustriskey'])
        try:
            with urlopen(req) as response:
                with open(savename, 'wb') as out_file:
                    shutil.copyfileobj(response, out_file)
        except HTTPError as e:
            if e.code==404:
                print('skipping ' +str(obj[i,1]))
                continue
            else:
                print(e.reason, e.code)
                break

        os.system('sudo mv {:d}_tree.hdf5 /Volumes/Elements/illustris/tng50trees/'.format(int(obj[i,1])))

