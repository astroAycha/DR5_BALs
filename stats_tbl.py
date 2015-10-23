""" quick script to examine the properties of objects in each cluster: e.g., fraction of AlIII absorption in CIV clusters, the reddening properties for objects in each cluster, and whether objects are showing variablity or not. """

from astropy.table import Table, join
from astropy import units as u




#cross match SDSS Quasar data with the DR5 BAL catalog. Then do some work to create a nice catalog file with the stuff i need.

#dr5qso= Table.read('dr5qso.fits') #full DR5 quasar catalog from SDSS
#dr5bal= Table.read('dr5BALCat.fits') #DR5 BAL quasars catalog, downloaded from VizieR

#dr5bal['SDSS'].name= "SDSSName" #change column name

#t= join(dr5qso, dr5bal, keys= 'SDSSName', join_type= 'right') #cross-match using SDSS name

#ended up doing the cross-matching using TopCat and the RA and DEC.

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

t.write('myBALCat.fits')

## I used TopCat to cross-match this table with extinction_tab.fits and saved the new table as myBALCat.fits (same name)



## table myBALs_red_var.fits contains results of cross-matching myBALCat.fits with krawczyk_reddening.fits (Krawczyk et al. 2015) and filiz2014.fits (Filiz Ak et al. 2014)
## table myBALs_red_var_xray.fits has the X-ray data from Sarah (and Robyn Smith).

## table myBALs_red_var_xray.fits includes masked arrays. Will not allow me to make histograms. Need to fix (fill masked cells).

t= Table.read('myBALs_var_red_xray_he2_shen.csv')

#this table is a result of corrs-matching myBALCat.fits with: filiz2014.fits (Filiz Ak et al. 2015), krawczyk_reddning.fits (Krawczyk et al. 2015), xray_Robyn_final.fits (X-Ray data from Sarah), tab1_Baskin15.txt (Baskin et al. 2015), and dr7_bh_May09_2011.fits (Shen et al. 2011 DR7 quasar catalog).

#keep only interesting columns
t.keep_columns(['SDSSName', 'RA_1', 'DEC_1', 'M_i', 'MJD_spec', 'plate_1', 'fiberid', 'z_1', \
                'BI_SiIV', 'BIO_SiIV', 'EW_SiIV', 'Vmin_SiIV', 'Vmax_SiIV', \
                'BI_CIV', 'BIO_CIV', 'EW_CIV_1', 'Vmin_CIV', 'Vmax_CIV', \
                'BI_AlIII', 'BIO_AlIII', 'EW_AlIII', 'Vmin_AlIII', 'Vmax_AlIII', \
                'BI_MgII', 'BIO_MgII', 'EW_MgII', 'Vmin_MgII', 'Vmax_MgII', \
                'SN1700', 'logF1400', 'logF2500',\
                'E_B-V_1', 'E_B-V_2', 'R_6CM_2500A', 'DELTA_AOX', \
                'HeII_EW', 'alpha_UV', 'v_md', 'CF', 'FWHM', \
                'BI1', 'BI2', 'Delt', 'Separation'])

#some renaming of columns to keep thigs tidy

t['RA_1'].name= 'RA'
t['DEC_1'].name= 'DEC'
t['plate_1'].name= 'plate'
t['z_1'].name= 'z'
t['FWHM'].name= 'CIV_BAL_FWHM'

#save as csv -some of the columns are masked (empty cells with nans that astropy table could not read for some reason). I hacked the file and replaced those nan cells with zeros.

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
