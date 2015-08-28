""" 28 Aug 2015
read tables with clusters generated using the bal_cluster.py script.
cross match to find the spectra FITS files and create composite spectra
"""

import numpy as np

from astropy.table import Table, join
from astropy.io import fits

from scipy.stats import sigmaclip

def spec_compos(line, k):

    """create median composite spectra from a list of file names
    """

    clstr= Table.read(line+str(k)+"_tab.fits")
    colnames= clstr.colnames

    data= Table.read('myBALCat.fits')

    t= join(clstr, data, keys= 'SDSSName', join_type= 'left')


    compos_ls= [] #list of arrays, each array is a composite that represents one cluster
    spec_num= []  #number of objects in each composite (cluster) to be used in the plotting
    std_ls= [] # a list of stdev arrays for each composite

    for i in range(k):
        c= t[labels==i] #cluster with label= i
        clust_spec= np.arange(1100, 4000, 0.5) # wavelength -this is the base that spectra will be stacked to before they're medianed

        for q in range(len(cluster)):
            file_name= "./proc_data/spec-"+str(c['MJD'][q])+"-"+str(c['plate'][q])+"-"+str(c['fiber'][q])+"_proc.fits"
            spec=fits.open(file_name)
            flx= spec[0].data[1]
            wlen= spec[0].data[0]
            norm_flx= flx/np.median(flx[2360:2390]) # normalize spectra

            clust_spec= np.vstack((clust_spec, norm_flx)) # 2D array. 1st row: restframe wavelength, other rows have corrected fluxes of spectra from clusters (one for each row)
        del spec
    
    print "cluster", i+1, "has", len(clust_spec[1:]), "objects"
    
    spec_num.append(len(clust_spec[1:]))
    
    clipped_compo= [] # 3 sigma clipping before median combining the spectra
    std_array=[]
    
    for i in range(clust_spec.shape[1]):
        
        y= sigmaclip(clust_spec[1:,i], 3, 3)
        m=median(y[0])
        clipped_compo.append(m)
        s= stdev(y[0])
        std_array.append(s)
    
    
    compos_ls.append(clipped_compo) # list with the composites (compos[0] is composite from 1st cluster, compos[1] 2nd cluster,...)
    std_ls.append(std_array)

    #save each composite into a FITS file

    for m,n in zip(range(k), spec_num): #assumes there is a directory called composites in the working directory
    
        spec_name= "./composites/"+line+"_"+str(k)+"clstr"+str(m+1)+".fits"
        spec_file= np.vstack((wlen,compos_ls[m], std_ls[m]))
        
        hdu= fits.PrimaryHDU(spec_file)
        hdr= hdu.header
        hdr.set('SPEC_NUMBER', n)
        hdr.set('COMPOSITE', line+"-K"+str(k))
        hdr.set('PARAMETERS USED', 'EW, RHWHM, BHWHM')
        hdu.writeto(spec_name)


