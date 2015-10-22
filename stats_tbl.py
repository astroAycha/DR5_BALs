""" quick script to examine the properties of objects in each cluster: e.g., fraction of AlIII absorption in CIV clusters, the reddening properties for objects in each cluster, and whether objects are showing variablity or not. """

from astropy.table import Table, join


## table myBALs_red_var.fits contains results of cross-matching myBALCat.fits with krawczyk_reddening.fits (Krawczyk et al. 2015) and filiz2014.fits (Filiz Ak et al. 2014)
## table myBALs_red_var_xray.fits has the X-ray data from Sarah (and Robyn Smith).

## table myBALs_red_var_xray.fits includes masked arrays. Will not allow me to make histograms. Need to fix (fill masked cells).

t= Table.read('myBALs_var_red_xray_he2.csv')

t.keep_columns(['SDSSName', 'RA_1', 'DEC_1', 'M_i', 'MJD_spec', 'plate_1', 'fiberid', 'z_1', \
                'BI_SiIV', 'BIO_SiIV', 'EW_SiIV', 'Vmin_SiIV', 'Vmax_SiIV', \
                'BI_CIV', 'BIO_CIV', 'EW_CIV_1', 'Vmin_CIV', 'Vmax_CIV', \
                'BI_AlIII', 'BIO_AlIII', 'EW_AlIII', 'Vmin_AlIII', 'Vmax_AlIII', \
                'BI_MgII', 'BIO_MgII', 'EW_MgII', 'Vmin_MgII', 'Vmax_MgII', \
                'SN1700', 'logF1400', 'logF2500',\
                'E_B-V_1', 'E_B-V_2', 'R_6CM_2500A', 'DELTA_AOX', \
                'HeII_EW', 'alpha_UV', 'v_md', 'CF', 'FWHM', \
                'BI1', 'BI2', 'Delt', 'Separation'])

t.write('matched_tbl.csv')
t.write('matched_tbl.fits')


data= Table.read('myBALCat.fits')

c4= Table.read('./clusters/CIV3clstrs.fits')

si4= Table.read('./clusters/SiIV3clstrs.fits')

mg2= Table.read('./clusters/MgII3clstrs.fits')

c4data= join(c4, data, keys='SDSSName')
si4data= join(si4, data, keys= 'SDSSName')
mg2data= join(mg2, data, keys= 'SDSSName')

print len(c4data[(c4data['label'] ==2) & (c4data['BIO_AlIII'] >0)])*100./len(c4data[c4data['label'] ==2])

print len(si4data[(si4data['label'] ==0) & (si4data['BIO_AlIII'] >0)])*100./len(si4data[si4data['label'] ==0])

print len(mg2data[(mg2data['label'] ==2) & (mg2data['BIO_AlIII'] >0)])*100./len(mg2data[mg2data['label'] ==2])


hist(c4data['BIO_AlIII'][c4data['label'] ==1])


open file with cluster
open file with x-ray data
cross match with only cluster data
make histgrams
