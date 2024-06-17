import numpy as np

def getdeterminationcoef(x,y,**args):
    fit=np.polynomial.polynomial.Polynomial.fit(x,y,**args)
    fitval=fit(x)
    res2=np.sum((y-fitval)**2)
    ss2=np.sum((y-np.mean(y))**2)
    return 1-res2/ss2
