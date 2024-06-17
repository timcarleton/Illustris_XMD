import numpy as np

def getperiodicdist(coord0,coords,boxsize=100):
    delta2=np.zeros(len(coords))
    cshape=coords.shape
    for i in range(cshape[1]):
        delta2i=np.min([(coords[:,i]-coord0[i])**2,(coord0[i]-coords[:,i]+boxsize)**2,(coord0[i]-coords[:,i]-boxsize)**2],axis=0)
        delta2+=delta2i
    return np.sqrt(delta2)
