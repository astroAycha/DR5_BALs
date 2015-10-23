""" quick script to examine the properties of objects in each cluster: e.g., fraction of AlIII absorption in CIV clusters, the reddening properties for objects in each cluster, and whether objects are showing variablity or not. """

from astropy.table import Table, join
from astropy import units as u




#cross match SDSS Quasar data with the DR5 BAL catalog. Then do some work to create a nice catalog file with the stuff i need.

#dr5qso= Table.read('dr5qso.fits') #full DR5 quasar catalog from SDSS
#dr5bal= Table.read('dr5BALCat.fits') #DR5 BAL quasars catalog, downloaded from VizieR

#dr5bal['SDSS'].name= "SDSSName" #change column name

#t= join(dr5qso, dr5bal, keys= 'SDSSName', join_type= 'right') #cross-match using SDSS name

#ended up doing the cross-matching with TOPCAT using the RA and DEC.

t= Table.read('SDSS_DR5BALs.fits')

t['z_2'].name= 'z'
t['plate_1'].name= 'plate'

t['BI_Si_'].name= 'BI_SiIV'
t['BIO_Si_'].name= 'BIO_SiIV'
t['EW_Si_'].name= 'EW_SiIV'
t['vmin_Si_'].name= 'Vmin_SiIV'
t['vmax_Si_'].name= 'Vmax_SiIV'
t['fd_Si_'].name= 'f_deep_SiIV'

t['BI_C_'].name= 'BI_CIV'
t['BIO_C_'].name= 'BIO_CIV'
t['EW_C_'].name= 'EW_CIV'
t['vmin_C_'].name= 'Vmin_CIV'
t['vmax_C_'].name= 'Vmax_CIV'
t['fd_C_'].name= 'f_deep_CIV'

t['BI_Al_'].name= 'BI_AlIII'
t['BIO_Al_'].name= 'BIO_AlIII'
t['EW_Al_'].name= 'EW_AlIII'
t['vmin_Al_'].name= 'Vmin_AlIII'
t['vmax_Al_'].name= 'Vmax_AlIII'
t['fd_Al_'].name= 'f_deep_AlIII'

t['BI_Mg_'].name= 'BI_MgII'
t['BIO_Mg_'].name= 'BIO_MgII'
t['EW_Mg_'].name= 'EW_MgII'
t['vmin_Mg_'].name= 'Vmin_MgII'
t['vmax_Mg_'].name= 'Vmax_MgII'
t['fd_Mg_'].name= 'f_deep_MgII'


#keep only columns I want. can check the current ones using t.colnames

t.keep_columns(['SDSSName', 'RA', 'DEC', 'z', 'M_i', \
                'plate', 'fiberid', 'MJD_spec', 'n_SDSS', \
                'BI_SiIV', 'BIO_SiIV', 'EW_SiIV', 'Vmin_SiIV', 'Vmax_SiIV', 'f_deep_SiIV', \
                'BI_CIV', 'BIO_CIV', 'EW_CIV', 'Vmin_CIV', 'Vmax_CIV', 'f_deep_CIV', \
                'BI_AlIII', 'BIO_AlIII', 'EW_AlIII', 'Vmin_AlIII', 'Vmax_AlIII', 'f_deep_AlIII', \
                'BI_MgII', 'BIO_MgII', 'EW_MgII', 'Vmin_MgII', 'Vmax_MgII', 'f_deep_MgII',\
                'flg', 'SN1700', 'logF1400', 'logF2500'])

t['EW_SiIV'].unit= t['EW_CIV'].unit = t['EW_AlIII'].unit= t['EW_MgII'].unit= 'Angstrom'

#remove row with problem in plate tarball. The download stopped in plate 0357. tarball file seems corrupt.
bad_rows=[]
f= open('bad_plates.txt', 'r') #this file contains a list of zombie plates. I found those plates using their sizes (314B) compared to (~78M) in the other "healty" plates.

badP= f.readlines()
f.close()

for i in range(0, len(t)):
    
    for j in badP:
        
        if t['plate'][i] == int(j.rstrip()):
            #print j.rstrip()
            bad_rows.append(i) #number of rows for the corrupt plate.

t.remove_rows(bad_rows)

print len(bad_rows), "rows were removed"
print "table now has ", len(t), "lines"

'''
## this is not what I really want. Will have to come back to it later.

#now separate flags: left: EmLost, middle: BALManyBadBins, right: BlueWingAbs, 1= MgII, 2= AlIII, 4= CIV, 8= SiIV

f1, f2, f3= [], [], []

for i in t['flg']:
    f1.append(int(i[0]))
    f2.append(int(i[1]))
    f3.append(int(i[2]))

c1= Column(f1, name= 'EmLostFlag')
c2= Column(f2, name= 'BMBBFlag')
c3= Column(f3, name= 'BlueWingFlag')

t.add_columns([c1, c2, c3], indexes=[30, 30, 30])
'''

t.write('myBALCat.fits')

## I used TOPCAT to cross-match this table with extinction_tab.fits and saved the new table as myBALCat.fits (same name)
##===========================================================================##

##Now cross-match this table with 5 other tables and do some cleaning:

#this table is a result of corrs-matching myBALCat.fits with: filiz2014.fits (Filiz Ak et al. 2015), krawczyk_reddning.fits (Krawczyk et al. 2015), tab1_Baskin15.txt (Baskin et al. 2015), and dr7_bh_May09_2011.fits (Shen et al. 2011 DR7 quasar catalog).
#xray_Robyn_final.fits (X-Ray data from Sarah) only has 120 of the objects so i did not include it. 'DELTA_AOX'

t= Table.read('myBALCat_var_red_he2_shen.fits')


#keep only columns I need:

t.keep_columns(['SDSSName', 'RA_1', 'DEC_1', 'M_i', 'MJD_spec', 'plate_1', 'fiberid', 'z_1', \
                'BI_SiIV', 'BIO_SiIV', 'EW_SiIV', 'Vmin_SiIV', 'Vmax_SiIV', 'f_deep_SiIV', \
                'BI_CIV', 'BIO_CIV', 'EW_CIV_1', 'Vmin_CIV', 'Vmax_CIV', 'f_deep_CIV', \
                'BI_AlIII', 'BIO_AlIII', 'EW_AlIII', 'Vmin_AlIII', 'Vmax_AlIII', 'f_deep_AlIII', \
                'BI_MgII', 'BIO_MgII', 'EW_MgII_1', 'Vmin_MgII', 'Vmax_MgII', 'f_deep_MgII', \
                'SN1700', 'logF1400', 'logF2500', 'flg', \
                'E_B-V_1', 'E_B-V_2', \
                'HeII_EW', 'alpha_UV', 'v_md', 'CF', 'FWHM', \
                'BI1', 'BI2', 'Delt', \
                'Z_HW', 'LOGLBOL', 'R_6CM_2500A',
                'LOGL_MGII', 'FWHM_MGII', 'EW_MGII_2', 'EW_FE_MGII', \
                'LOGL_CIV', 'FWHM_CIV', 'EW_CIV_2', 'VOFF_CIV_PEAK', 'LOGBH', 'LOGEDD_RATIO'])

#some renaming of columns to keep thigs tidy

t['RA_1'].name= 'RA'
t['DEC_1'].name= 'DEC'
t['plate_1'].name= 'plate'
t['z_1'].name= 'z'
t['EW_CIV_1'].name= 'EW_CIV'
t['EW_MgII_1'].name= 'EW_MgII'

t['HeII_EW'].name= 'HeII_EW_BLH'
t['alpha_UV'].name= 'alpha_UV_BLH'
t['v_md'].name= 'v_md_BLH'
t['CF'].name= 'CF_BLH'
t['FWHM'].name= 'FWHM_CIB_BAL_BLH'

t['LOGLBOL'].name= 'LOGLBOL_DR7'
t['R_6CM_2500A'].name= 'R_6CM_2500A_DR7'
t['LOGL_MGII'].name= 'LOGL_MGIIe_DR7'
t['FWHM_MGII'].name= 'FWHM_MGIIe_DR7'
t['EW_MGII_2'].name= 'EW_MGIIe_DR7'
t['EW_FE_MGII'].name= 'EW_FE_MGIIe_DR7'
t['LOGL_CIV'].name= 'LOGL_CIVe_DR7'
t['FWHM_CIV'].name= 'FWHM_CIVe_DR7'
t['EW_CIV_2'].name= 'EW_CIVe_DR7'
t['VOFF_CIV_PEAK'].name= 'VOFF_CIVe_PEAK_DR7'
t['LOGBH'].name= 'LOGBH_DR7'
t['LOGEDD_RATIO'].name= 'LOGEDD_RATIO_DR7'


#save as csv -some of the columns are masked (empty cells with nans that astropy table could not read for some reason). I hacked the file and replaced those nan cells with zeros.

t.write('myBALCat_xtra.csv')

## this table has many 'nan' values that astropy's table did not want to fill. I opend it as a spreadsheet and did a quick search and replace with 0 instead of nan.

##============================================================##


