
"""this is to record how i generated the data files i am using.
    """
from astropy.table import Table, Column, join
from astropy import units as u



## Feb 29 2016
## updating the tables (the old stuff are at the bottom of this file just for the record)

# starting with table 1 from Gibson et al. 2009 --> 5035 objects

# cross-match with the DR5 Schnider catalog to get plate-mjd-fiber used for download --> new_dr5_bals.fits

# used RA and DEC to get extinction data from http://irsa.ipac.caltech.edu/applications/DUST/ results are saved in extinction_tbl.txt

# cross match with the Shen et al. 2011 catalog to use the Hewett and Wild 2010 redhsifts --> found 5031 objects using RA and DEC with max error = 3 arcsec. file saved as tbl1.fits

## keep only relevant columns and rename some --> save as myBALs.fits

t= Table.read('tbl1.fits')

t.keep_columns(['Name', 'RA_1', 'DEC_1', 'z_1', 'SiIV-BI', 'SiIV-BIO', 'SiIV-EW', 'SiIV-vmin', 'SiIV-vmax', 'SiIV-fdeep', 'CIV-BI', 'CIV-BIO', 'CIV-EW', 'CIV-vmin', 'CIV-vmax', 'CIV-fdeep', 'AlIII-BI', 'AlIII-BIO', 'AlIII-EW', 'AlIII-vmin', 'AlIII-vmax', 'AlIII-fdeep', 'MgII-BI', 'MgII-BIO', 'MgII-EW', 'MgII-vmin', 'MgII-vmax', 'MgII-fdeep', 'SiIV-EmL', 'CIV-EmL', 'AlIII-EmL', 'MgII-EmL', 'SiIV-BMBB', 'CIV-BMBB', 'AlIII-BMBB', 'MgII-BMBB', 'CIV4-BWA', 'SN1700', 'logF1400', 'logF2500', 'M_i','MJD_spec', 'plate_1', 'fiberid', 'AV_SandF', 'Z_HW', 'LOGEDD_RATIO'])

t['Name'].name= 'SDSSName'
t['RA_1'].name= 'RA'
t['DEC_1'].name= 'DEC'
t['z_1'].name= 'z'
t['plate_1'].name= 'plate'

t.write("myBALs.fits")

# myBALs.fits is used to correct the spectra for Galactic extinction and redshift.



### big_tbl.fits is the result of corss-matching myBALs.fits with the Shen et al. 2010 data table, Filiz Ak et al. 2014 table 1, Baskin et al. 2015 table 1, and Krawczyk et al. 2015 table 1

tt= Table.read('big_tbl.fits')

##keep only the columns I need

tt.keep_columns(['SDSSName','RA_1','DEC_1','z_1','SiIV-BI','SiIV-BIO','SiIV-EW','SiIV-vmin','SiIV-vmax','SiIV-fdeep','CIV-BI','CIV-BIO','CIV-EW','CIV-vmin','CIV-vmax','CIV-fdeep','AlIII-BI','AlIII-BIO','AlIII-EW','AlIII-vmin','AlIII-vmax','AlIII-fdeep','MgII-BI','MgII-BIO','MgII-EW','MgII-vmin','MgII-vmax','MgII-fdeep','SiIV-EmL','CIV-EmL','AlIII-EmL','MgII-EmL','SiIV-BMBB','CIV-BMBB','AlIII-BMBB','MgII-BMBB','CIV4-BWA','SN1700','logF1400','logF2500','M_i','MJD_spec','plate_1','fiberid','AV_SandF','LOGEDD_RATIO_1','Z_HW_1','Dal1','E_Dal1','e_Dal2','E_B-V_1','E_E_B-V_1','e_E_B-V_2','Dal2','E_Dal3','e_Dal4','E_B-V_2','E_E_B-V_3','e_E_B-V_4','HeII_EW','alpha_UV','v_md','CF','FWHM','BI1','e_BI1','BI2','e_BI2','Delt'])

tt['RA_1'].name= 'RA'
tt['DEC_1'].name= 'DEC'
tt['z_1'].name= 'z'
tt['plate_1'].name= 'plate'
tt['LOGEDD_RATIO_1'].name= 'LOGEDD_RATIO'
tt['Z_HW_1'].name= 'Z_HW'

# calculate intrinsic extinction index and add column to table

cc= Column(name='int_alpha_nu', data= tt['Dal1']-0.28)
tt.add_column(cc)

tt.write('myBALsx.fits')

## replacing nan cells with -999 using astropy's fill_value did not work.
## saved as csv, opened as spreadsheet, replaced nan values with -999.
tt.write('myBALsx.csv')

############################

## the old stuff starts here:

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

## add a column

cc= Column(name='int_alpha_nu', data= t['Dal1']-0.28)
t.add_column(cc)

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
t['int_alpha_nu']= t['int_alpha_nu'].filled(-999)
t['E_B-V_1'] = t['E_B-V_1'].filled(-999)
t['Dal2'] = t['Dal2'].filled(-999)
t['E_B-V_2'] = t['E_B-V_2'].filled(-999)

t.write('myBALCat_xtra.csv')

## this table has many 'nan' values that astropy's table did not want to fill. I opend it as a spreadsheet and did a quick search and replace with -999 instead of nan.

##============================================================##


