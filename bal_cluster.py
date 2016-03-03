""" unsupervised clustering on BALQ measurements form the Gibson et al. BALQ catalog from SDSS-DR5
"""

import numpy as np

from astropy.io import fits
from astropy.table import Table, join

from sklearn.cluster import KMeans
from sklearn import metrics

def bal_cluster(line, k):

    """ unsupervised clustering using KMeans.
    param:
    line: which line to use for the clustering features
    k: number of clusters
    g1= normalized EW, Vmin, Vmax
    g2= normalized EW, Vmin, dV
    g3: normalized EW, Vmin, EW/dv
    g4: normalized EW, Vmax, EW/dv
    g5: normalized EW, Vmin, Vmax, dV
    g6: normalized EW, Vmin, Vmax, EW/dv
    
    """
    group = "g6/"
    colnames= ('EW', 'Vmin', 'Vmax', 'EW_dV', 'label', 'SDSSName', 'z')
    
    data= Table.read('myBALs.fits')
    #data= Table.read('myBALCat_xtra.csv', format= 'ascii.csv')

    #selec the sample: line has an absorption trough: BI0 >0 , S/N >3, and a redshift cutoff to restrict the bandwidth

    if line == "MgII":
        z1= 1.1
        z2= 2.2
        lum= "logF2500"
    
    else:
        z1= 1.79
        z2= 3.7
        lum= "logF1400"


    s= data[(data[line+'-BIO'] >0) & (data['SN1700'] >3) & (data['Z_HW'] >z1) & (data['Z_HW'] <z2)]# & (data[lum] !=-999)]
    
    print "sample has", len(s), "objects"


    #features
    """
    redshift= s['Z_HW']
    bi= s[line+'-BIO']
    ew= s[line+'-EW']
    vmin= s[line+'-vmin']
    vmax= s[line+'-vmax']
    fdeep= s[line+'-fdeep']

    """
    s['Z_HW'].fill_value= -999
    redshift= s['z'].filled() #redshift

    s[line+'-BI'].fill_value= -999
    bi= s[line+'-BI'].filled() # balnicity: integration 3000-25000

    s[line+'-BIO'].fill_value= -999
    bi0= s[line+'-BIO'].filled() # modified balnicity: integration 0-25000

    s[line+'-EW'].fill_value= -999
    ew= s[line+'-EW'].filled() # restframe absorption EW

    s[line+'-vmin'].fill_value= -999
    vmin= s[line+'-vmin'].filled() # minimum velocity

    s[line+'-vmax'].fill_value= -999
    vmax= s[line+'-vmax'].filled() # maximum velocity

    s[line+'-fdeep'].fill_value= -999
    fdeep= s[line+'-fdeep'].filled()
    
    s[lum].fill_value= -999
    cl= s[lum].filled() # Log of 1400 or 2500 monochromatic luminosity

    s['SDSSName'].fill_value= -999
    names= s['SDSSName'].filled() # SDSS name
    

    dv= vmax- vmin # delta v
    dd= ew/dv # some sort of an estimate of the trough depth
    
    #standardize (normalized) parameters before using them in clustering
    ew_n= (ew - mean(ew))/std(ew)
    vmin_n= (vmin - mean(vmin))/std(vmin)
    vmax_n= (vmax - mean(vmax))/std(vmax)
    dv_n = (dv - mean(dv))/std(dv)
    dd_n= (dd- mean(dd))/std(dd)
    
    #cl_n= (cl - mean(cl))/std(cl)
    
    # list of features to be used in clustering
    #f= [ew_n, vmin_n, dv_n]
    f= [ew_n, vmin_n, vmax_n, dd_n]

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
    
    
    # use with 3 features
    #clstr_tab= Table([qs[:,0], qs[:,1], qs[:,2], labels, names, redshift], \
                    # names= colnames, \
                    # dtype= ('float64', 'float64', 'float64', 'int', 'S18', 'float64'))
                     
    # use with 4 features
    clstr_tab= Table([qs[:,0], qs[:,1], qs[:,2], qs[:,3], labels, names, redshift], \
                                      names= colnames, \
                                      dtype= ('float64', 'float64', 'float64', 'float64', 'int', 'S18', 'float64'))


    clstr_tab.write("./clusters/"+group+clstr_name+"clstrs.fits", format= 'fits')
    
    
    #clstr_tab.write("./clusters/"+str(len(f))+"features/"+clstr_name+"clstrs.fits", format= 'fits')

    return
   

