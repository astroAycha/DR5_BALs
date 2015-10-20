""" Bulk download of spectra from the SDSS database
spectra are from the QSO database and have a redshift range of 1.6-2.1
To generate the links to the FITS files, I am getting the file name from the DR10 catalog of Paris et al. 2014
"""

import numpy as np
from astropy.table import Table, Column, join
from astropy import units as u
import subprocess
import wget
import tarfile
import shutil
import os


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
    
t['BI_C_'].name= 'BI_CIV'
t['BIO_C_'].name= 'BIO_CIV'
t['EW_C_'].name= 'EW_CIV'
t['vmin_C_'].name= 'Vmin_CIV'
t['vmax_C_'].name= 'Vmax_CIV'
    
t['BI_Al_'].name= 'BI_AlIII'
t['BIO_Al_'].name= 'BIO_AlIII'
t['EW_Al_'].name= 'EW_AlIII'
t['vmin_Al_'].name= 'Vmin_AlIII'
t['vmax_Al_'].name= 'Vmax_AlIII'
    
t['BI_Mg_'].name= 'BI_MgII'
t['BIO_Mg_'].name= 'BIO_MgII'
t['EW_Mg_'].name= 'EW_MgII'
t['vmin_Mg_'].name= 'Vmin_MgII'
t['vmax_Mg_'].name= 'Vmax_MgII'
    
    
#keep only columns I want. can check the current ones using t.colnames
    
t.keep_columns(['SDSSName', 'RA', 'DEC', 'z', 'M_i', \
                'plate', 'fiberid', 'MJD_spec', 'n_SDSS', \
                'BI_SiIV', 'BIO_SiIV', 'EW_SiIV', 'Vmin_SiIV', 'Vmax_SiIV', \
                'BI_CIV', 'BIO_CIV', 'EW_CIV', 'Vmin_CIV', 'Vmax_CIV', \
                'BI_AlIII', 'BIO_AlIII', 'EW_AlIII', 'Vmin_AlIII', 'Vmax_AlIII', \
                'BI_MgII', 'BIO_MgII', 'EW_MgII', 'Vmin_MgII', 'Vmax_MgII', \
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


def dr5_download(bals, plates_dir):

    """download a list of plates from SDSS DR5. Extract specific files with specific MJDs and FiberIDs 
        bals: file with BAL quasars from DR5 that also contains plate, mjd, and fiberid for each row.
        plates_dir: directory where you want the plates to be downloaded to
        balq_dir: directory where your data (FITS files) for the DR5 BAL catalog only will be saved
        
        these two directory will be created within the function
        
        """
    
    
    os.mkdir(plates_dir)

    catalog= Table.read(bals)
    #catalog['plate'].fill_value= 0000
    #all_plates= catalog['plate'].filled()
    all_plates= catalog['plate']
    plates= set(all_plates) #select only the unique values

    links_list= []
   
    for plate_num in plates:
        links_list.append('http://das.sdss.org/spectro/ss_tar_23/'+str(plate_num).zfill(4)+'.tar.gz')


    #now download. this will take a long time
    print "Go take a nap! This is gonna take a looong time"
    i= 0
    for line in links_list:
        
        wget.download(line, out=plates_dir)
        i+=1
        print "\n Downloaded", i, "of", len(links_list), "files"

    return
    
def spec_xtract(bals, plate_dir, balq_dir):

    """extract specific files from each plate"""

    os.mkdir(balq_dir)
    
    catalog= Table.read(bals)
    all_plates= catalog['plate']
    plates= set(all_plates) #select only the unique values
    
    ord_plates= sorted(plates)
    
    i= 0
    for p in ord_plates:
        i+=1
        print i
        mjd_ls= catalog['MJD_spec'][catalog['plate']== p]
        fib_ls= catalog['fiberid'][catalog['plate']== p]
        
        plate_num= str(p).zfill(4)
        print "plate", plate_num, "has", len(mjd_ls), "spectra in the DR5 BAL catalog"
        
        tarball= tarfile.open(plate_dir+"/"+plate_num+".tar.gz", "r")
        tarball.extractall()
        tarball.close()
        
        del tarball
        
        spec_file_names= []
        for f in range(len(mjd_ls)):

            spec_file_names.append(plate_num+"/spSpec/spSpec-"+str(mjd_ls[f])+"-"+plate_num+"-"+str(fib_ls[f]).zfill(3)+".fit")
            
        for spec_file in spec_file_names:
            shutil.copy2(spec_file, balq_dir)

        subprocess.call(["rm", "-R", plate_num])


    return

