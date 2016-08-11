""" 28 Aug 2015
read tables with clusters generated using the bal_cluster.py script.
cross match to find the spectra FITS files and create composite spectra
"""

import numpy as np

from astropy.table import Table, join
from astropy.io import fits

from scipy.stats import sigmaclip

def spec_compos(line, k, g):

    """create median composite spectra from a list of file names
    k: number of clusters
    g: directory where clusters are stored. e.g, g1/ , g2/ or g3/ ...
    """

    clstr= Table.read("./clusters/"+g+"/"+line+str(k)+"clstrs.fits")
    column_names= clstr.colnames
    
    data= Table.read('myBALs.fits')

    t= join(clstr, data, keys= 'SDSSName', join_type= 'left')
    print len(t)

    median_compos_ls= [] #list of arrays, each array is a medain composite that represents one cluster
    mean_compos_ls= [] #list of arrays, each array is a mean composite that represents one cluster
    spec_num= []  #number of objects in each composite (cluster) to be used in the plotting
    std_ls= [] # a list of stdev arrays for each composite

    for i in range(k):
        c= t[t['label']==i] #cluster with label= i
        clust_spec= np.arange(1100, 4000, 0.5) # wavelength -this is the base that spectra will be stacked to before they're medianed

        for q in range(len(c)):
            file_name= "./proc_data/spec-"+str(c['plate'][q])+"-"+str(c['MJD_spec'][q])+"-"+str(c['fiberid'][q]).zfill(4)+"_proc.fits"
            spec=fits.open(file_name)
            flx= spec[0].data[1]
            wlen= spec[0].data[0]
            # normalize spectra
            x1= (1975-1100)/2
            x2= (2000-1100)/2  # values used to normalize spectrum by. change 1975A and 2000A to the desired wavelength range
            norm_flx= flx/np.median(flx[x1:x2])

            clust_spec= np.vstack((clust_spec, norm_flx)) # 2D array. 1st row: restframe wavelength, other rows have corrected fluxes of spectra from clusters (one for each row)
            del spec
    
        print "cluster", i+1, "has", len(clust_spec[1:]), "objects"
    
        spec_num.append(len(clust_spec[1:]))
    
        clipped_median, clipped_mean= [], [] # 3 sigma clipping before median and mean combining the spectra
        std_array=[]
    
        for i in range(clust_spec.shape[1]):
        
            y= sigmaclip(clust_spec[1:,i], 3, 3)
            m= median(y[0])
            clipped_median.append(m)
            avg= mean(y[0])
            clipped_mean.append(avg)
            s= std(y[0])
            std_array.append(s)
    
    
        median_compos_ls.append(clipped_median) # list with the composites (compos[0] is composite from 1st cluster, compos[1] 2nd cluster,...)
        mean_compos_ls.append(clipped_mean)
        std_ls.append(std_array)

    #save each composite into a FITS file

    for m,n in zip(range(k), spec_num): #assumes there is a directory called composites in the working directory
    
        spec_name= "./composites/"+g+"/"+line+"_"+str(k)+"clstr"+str(m+1)+".fits"
        spec_file= np.vstack((wlen, median_compos_ls[m], mean_compos_ls[m], std_ls[m]))
        
        hdu= fits.PrimaryHDU(spec_file)
        hdr= hdu.header
        hdr.set('SPEC_NUMBER', n)
        hdr.set('COMPO', line+"-K"+str(k))
        param= '-'.join(column_names[:4])
        print param
        hdr.set('PARAM', param)
        hdu.writeto(spec_name)
    

    return

