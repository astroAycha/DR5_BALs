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


# flg column: left: EmLost, middle: BALManyBadBins, right: BlueWingAbs, 1= MgII, 2= AlIII, 4= CIV, 8= SiIV
# [0-9A-F] Flags EmLost (left), BALManyBadBins (middle) and BlueWingAbs (right) with
# 1=MgII, 2=AlIII, 4=CIV and 8=SiIV
# A= 10, B= 11, C= 12, D= 13, E= 14, F= 15
#example: C04 --> C=12= 4+8 = CIV+SiIV, 0=0, 4= CIV --> Emlost in CIV and SiIV, BALManyBadBins is 0, and BlueWingAbs in CIV

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
                'Dal1', 'E_B-V_1', 'Dal2', 'E_B-V_2', \
                'HeII_EW', 'alpha_UV', 'v_md', 'CF', 'FWHM', \
                'BI1', 'BI2', 'Delt', \
                'Z_HW', 'LOGBH', 'LOGLBOL', 'R_6CM_2500A',
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
t['FWHM'].name= 'FWHM_CIV_BAL_BLH'

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

#t.write('myBALCat_xtra.csv')

## filling empty or masked cells with -999

t['BIO_SiIV'] = t['BI_SiIV'].filled(-999)
t['BIO_CIV'] = t['BI_CIV'].filled(-999)
t['BIO_AlIII'] = t['BI_AlIII'].filled(-999)
t['BIO_MgII'] = t['BI_MgII'].filled(-999)
t['logF1400'] = t['logF1400'].filled(-999)
t['logF2500'] = t['logF2500'].filled(-999)


t['HeII_EW_BLH'] = t['HeII_EW_BLH'].filled(-999)
t['alpha_UV_BLH'] = t['alpha_UV_BLH'].filled(-999)
t['v_md_BLH'] = t['v_md_BLH'].filled(-999)
t['CF_BLH'] = t['CF_BLH'].filled(-999)
t['FWHM_CIV_BAL_BLH'] = t['FWHM_CIV_BAL_BLH'].filled(-999)

t['LOGLBOL_DR7'] = t['LOGLBOL_DR7'].filled(-999)
t['R_6CM_2500A_DR7'] = t['R_6CM_2500A_DR7'].filled(-999)
t['LOGL_MGIIe_DR7']= t['LOGL_MGIIe_DR7'].filled(-999)
t['FWHM_MGIIe_DR7'] = t['FWHM_MGIIe_DR7'].filled(-999)
t['EW_MGIIe_DR7'] = t['EW_MGIIe_DR7'].filled(-999)
t['EW_FE_MGIIe_DR7'] = t['EW_FE_MGIIe_DR7'].filled(-999)
t['LOGL_CIVe_DR7'] = t['LOGL_CIVe_DR7'].filled(-999)
t['FWHM_CIVe_DR7'] = t['FWHM_CIVe_DR7'].filled(-999)
t['EW_CIVe_DR7'] = t['EW_CIVe_DR7'].filled(-999)
t['VOFF_CIVe_PEAK_DR7'] = t['VOFF_CIVe_PEAK_DR7'].filled(-999)
t['LOGBH_DR7'] = t['LOGBH_DR7'].filled(-999)
t['LOGEDD_RATIO_DR7'] = t['LOGEDD_RATIO_DR7'].filled(-999)

t['Dal1'] = t['Dal1'].filled(-999)
t['E_B-V_1'] = t['E_B-V_1'].filled(-999)
t['Dal2'] = t['Dal2'].filled(-999)
t['E_B-V_2'] = t['E_B-V_2'].filled(-999)

t.write('myBALCat_xtra.csv')

## this table has many 'nan' values that astropy's table did not want to fill. I opend it as a spreadsheet and did a quick search and replace with -999 instead of nan.

##============================================================##


