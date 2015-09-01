""" unsupervised clustering on BALQ measurements form the Gibson et al. BALQ catalog from SDSS-DR5
"""

import numpy as np

from astropy.io import fits
from astropy.table import Table, join

from sklearn.cluster import KMeans
from sklearn import metrics

def bal_cluster(data_tab, line, k):

    """ unsupervised clustering using KMeans.
    param:
    data: data table
    line: which line to use for the clustering features
    k: number of clusters"""
    
    data= Table.read(data_tab)

    #selec the sample: line has an absorption trough: BI0 >0 , S/N >3, line-BMBB flag (BALManyBadBins) many bad bins in the spec as reported in the SDSS database

    s= data[(data['BIO_'+line] >0) & (data['SN1700'] >3)] # & (data['flg'] ==0)]
    
    print "sample has", len(s), "objects"

    #features
    s['z'].fill_value= 999
    redshift= s['z'].filled() #redshift
    s['BI_'+line].fill_value= 999
    bi= s['BI_'+line].filled() # balnicity: integration 3000-25000
    s['BIO_'+line].fill_value= 999
    bi0= s['BIO_'+line].filled() # modified balnicity: integration 0-25000
    s['EW_'+line].fill_value= 999
    ew= s['EW_'+line].filled() # restframe absorption EW
    s['Vmin_'+line].fill_value= 999
    vmin= s['Vmin_'+line].filled() # minimum velocity
    s['Vmax_'+line].fill_value= 999
    vmax= s['Vmax_'+line].filled() # maximum velocity
    
    if line== "MgII":
        lum= "logF2500"
    
    else:
        lum= "logF1400"
    
    s[lum].fill_value= 999
    cl= s[lum].filled() # Log of 1400 or 2500 monochromatic luminosity

    s['SDSSName'].fill_value= 999
    names= s['SDSSName'].filled() # SDSS name
    
    # list of features to be used in clustering
    f= [bi0, ew, vmin, vmax]# , cl]

    qs= np.column_stack(param for param in f) # 2D array to do clustering on

    #do the clustering

    kmeans= KMeans(init= 'k-means++', n_clusters= k, n_init= 10)
    kmeans.fit(qs)
    labels= kmeans.predict(qs)
    sc= metrics.silhouette_score(qs, labels)
    cntrs= kmeans.cluster_centers_
    print "Silhouette score= ", sc
    print "centroids: "+"\n", cntrs

    # save the results in a FITS table
    
    clstr_name= line+str(k)
    
    clstr_tab= Table([qs[:,0], qs[:,1], qs[:,2], qs[:,3], labels, names, redshift], \
                     names= ('BIO', 'EW', 'Vmin', 'Vmax', 'label', 'SDSSName', 'z'), \
                     dtype= ('float64', 'float64', 'float64', 'float64', 'int', 'S18', 'float64'))
    
    
    clstr_tab.write("./clusters/"+clstr_name+"clstrs.fits", format= 'fits')

    return
   

