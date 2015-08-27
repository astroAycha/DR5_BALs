""" unsupervised clustering on BALQ measurements form the Gibson et al. BALQ catalog from SDSS-DR5
"""

import numpy as np

from astropy.io import fits
from astropy.table import Table, join

from sklearn.cluster import KMeans
from sklearn import metrics

from scipy.stats import sigmaclip


def bal_cluster(data_tab, line, k):

    """ unsupervised clustering using KMeans.
    param:
    data: data table
    line: which line to use for the clustering features
    k: number of clusters"""
    
    data= Table.read(data_tab)

    #selec the sample: line has an absorption trough: BI0 >0 , S/N >3, line-BMBB flag (BALManyBadBins) many bad bins in the spec as reported in the SDSS database

    s= data[(data[line+'-BIO'] >0) & (data['SN1700'] >3) & (data[line+'-BMBB'] ==0)]
    
    print "sample has", len(s), "objects"

    #features
    s[line+'-BI'].fill_value= 999
    bi= s[line+'-BI'].filled() # balnicity: integration 3000-25000
    s[line+'-BIO'].fill_value= 999
    bi0= s[line+'-BIO'].filled() # modified balnicity: integration 0-25000
    s[line+'-EW'].fill_value= 999
    ew= s[line+'-EW'].filled() # restframe absorption EW
    s[line+'-vmin'].fill_value= 999
    vmin= s[line+'-vmin'].filled() # minimum velocity
    s[line+'-vmax'].fill_value= 999
    vmax= s[line+'-vmax'].filled() # maximum velocity
    s['logF1400'].fill_value= 999
    l1400= s['logF1400'].filled() # Log of 1400 monochromatic luminosity
    s['logF2500'].fill_value= 999
    l2500= s['logF2500'].filled() # Log of 2500 monochromatic luminosity
    s['Name'].fill_value= 999
    names= s['Name'].filled() # SDSS name

    # list of features to be used in clustering
    f= [bi0, ew, vmin, vmax, l1400, l2500]

    qs= np.column_stack(param for param in f) # 2D array to do clustering on

    #do the clustering

    kmeans= KMeans(init= 'k-means++', n_clusters= k, n_init= 10)
    kmeans.fit(qs)
    labels= kmeans.predict(qs)
    sc= metrics.silhouette_score(qs, labels)
    cntrs= kmeans.cluster_centers_
    print "Silhouette score= ", sc
    print "centroids: ", cntrs

    # save the results in a numpy array
    
    clstr_name= line+str(k)
    
    clstr_tab= Table([qs[:,0], qs[:,1], qs[:,2], qs[:,3], qs[:,4], qs[:,5], labels, names], \
                      names= ('BIO', 'EW', 'Vmin', 'Vmax', 'logF1400', 'logF2500', 'labels', 'SDSS_Names'), \
                      dtype= ('float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'int', 'S18'))
                      
    clstr_tab.write("./clusters/"+clstr_name+"tab.fits", format= 'fits')

    return
   

