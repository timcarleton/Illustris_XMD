from astropy.io import ascii, fits
from astropy import table
import numpy as np


f=fits.open('alltnggal_withr_wmetals_wvmax.fits')
f50=fits.open('tng50/tng50_wvmax.fits')

x=ascii.read('sampleselection_sfr5_notide_nomass.txt')
x50=ascii.read('sampleselection_tng50_sfr5_notide_nomass.txt')

galtypes=np.zeros(len(x)).astype(np.str)

w1=np.where(x['type']==1)[0]
w2=np.where(x['type']==2)[0]
w3=np.where(x['type']==3)[0]
w4=np.where(x['type']==4)[0]

galtypes[w1]='XMD'
galtypes[w2]='LMD'
galtypes[w3]='XMD Analog'
galtypes[w4]='LMD Analog'

tab=table.Table([galtypes,x['haloid'],np.log10(f[1].data.gas_metal_sfr[x['allgalindex']])+8.69,np.log10(f[1].data.mstar[x['allgalindex']]),np.log10(f[1].data.sfr[x['allgalindex']]),np.log10(f[1].data.mgas[x['allgalindex']]),np.log10(f[1].data.mhalo[x['allgalindex']]),f[1].data.half_mass_rad_star[x['allgalindex']],np.log10(f[1].data.star_metal[x['allgalindex']])],names=['Type','IllustrisTNG100 Halo ID','12+log(O/H)','log M_*','log SFR','log M_gas','log M_halo','r12 star','Z_*'])
tab['12+log(O/H)'].format='.2f'
tab['log M_*'].format='.2f'
tab['log SFR'].format='.2f'
tab['log M_gas'].format='.2f'
tab['log M_halo'].format='.2f'
tab['r12 star'].format='.2f'
tab['Z_*'].format='.2f'
tab.write('papertable.txt',format='latex',overwrite=True)


galtypes50=np.zeros(len(x50)).astype(np.str)

w150=np.where(x50['type']==1)[0]
w250=np.where(x50['type']==2)[0]
w350=np.where(x50['type']==3)[0]
w450=np.where(x50['type']==4)[0]

galtypes50[w150]='XMD'
galtypes50[w250]='LMD'
galtypes50[w350]='XMD Analog'
galtypes50[w450]='LMD Analog'

tab50=table.Table([galtypes50,x50['haloid'],np.log10(f50[1].data.gas_metal_sfr[x50['allgalindex']])+8.69,np.log10(f50[1].data.mstar[x50['allgalindex']]),np.log10(f50[1].data.sfr[x50['allgalindex']]),np.log10(f50[1].data.mgas[x50['allgalindex']]),np.log10(f50[1].data.mhalo[x50['allgalindex']]),f50[1].data.half_mass_rad_star[x50['allgalindex']],np.log10(f50[1].data.star_metal[x50['allgalindex']])],names=['Type','IllustrisTNG50 Halo ID','12+log(O/H)','log M_*','log SFR','log M_gas','log M_halo','r12 star','Z_*'])
tab50['12+log(O/H)'].format='.2f'
tab50['log M_*'].format='.2f'
tab50['log SFR'].format='.2f'
tab50['log M_gas'].format='.2f'
tab50['log M_halo'].format='.2f'
tab50['r12 star'].format='.2f'
tab50['Z_*'].format='.2f'
tab50.write('papertable50.txt',format='latex',overwrite=True)
