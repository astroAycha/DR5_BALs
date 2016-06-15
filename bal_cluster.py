""" unsupervised clustering on BALQ measurements form the Gibson et al. BALQ catalog from SDSS-DR5
"""

import numpy as np

from astropy.io import fits
from astropy.table import Table, join, Column

from sklearn.cluster import KMeans
from sklearn import metrics

def bal_cluster(line, k, g):

    """ unsupervised clustering using KMeans.
    param:
    line: which line to use for the clustering features
    k: number of clusters
    g0: normalized EW, Vmin, Vmax
    g1: normalized EW, Vmin, dV
    g2: normalized EW, Vmax, dV
    g3: normalized EW, Vmin, EW/dv
    g4: normalized EW, Vmax, EW/dv
    g5: normalized EW, Vmin, Vmax, dV
    g6: normalized EW, Vmin, Vmax, EW/dv
    g7: normalized Vmin, Vmax, EW/dV
    g8: normalized EW, dV, EW/dV
    g9: normalized EW, dV
    
    g11: normalized EW, Vmax, BI0/dV
    g12: normalized EW, Vmin, Vmax, BI0/dV
    g13: normalized BI0, Vmax, BI0/dV
    g14: normalized BI0, Vmin, Vmax, BI0/dV
    g15: normalized BI0, dV, Vmax
    g16: normalized BI0, EW, Vmin, Vmax, dV, BI0/dV
    g17: normalized BI0, EW, Vmin, Vmax, BI0/dV
    
    """
    
    clstr_name= line+str(k)


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
    redshift= s['Z_HW'].filled() #redshift

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
    
    dp= bi0/dv # another dimensionless parameter to estimate trough depth
    
    #standardize (normalized) parameters before using them in clustering
    ew_n= (ew - mean(ew))/std(ew)
    vmin_n= (vmin - mean(vmin))/std(vmin)
    vmax_n= (vmax - mean(vmax))/std(vmax)
    dv_n = (dv - mean(dv))/std(dv)
    dd_n= (dd- mean(dd))/std(dd)
    bi_n= (bi0 - mean(bi0))/std(bi0)
    dp_n= (dp - mean(dp))/std(dp)
    
    
    #cl_n= (cl - mean(cl))/std(cl)


    if g==  'g0':
    
        f = [ew_n, vmin_n, vmax_n]
        colnames= ('EW', 'Vmin', 'Vmax')
        datatype= ('float64', 'float64', 'float64')

    elif g== 'g1':
        
        f = [ew_n, vmin_n, dv_n]
        colnames= ('EW', 'Vmin', 'dV')
        datatype= ('float64', 'float64', 'float64')

    elif g== 'g2':
    
        f = [ew_n, vmax_n, dv_n]
        colnames= ('EW', 'Vmax', 'dV')
        datatype= ('float64', 'float64', 'float64')

    elif g== 'g3':
        f= [ew_n, vmin_n, dd_n]
        colnames= ('EW', 'Vmin', 'EW_dV')
        datatype= ('float64', 'float64', 'float64')
    
    elif g== 'g4':
        f= [ew_n, vmax_n, dd_n]
        colnames= ('EW', 'Vmax', 'EW_dV')
        datatype= ('float64', 'float64', 'float64')

    elif g== 'g5':
        f = [ew_n, vmin_n, vmax_n, dv_n]
        colnames= ('EW', 'Vmin', 'Vmax', 'dV')
        datatype= ('float64', 'float64', 'float64', 'float64')
    
    elif g== 'g6':
        f= [ew_n, vmin_n, vmax_n, dd_n]
        colnames= ('EW', 'Vmin', 'Vmax', 'EW_dV')
        datatype= ('float64', 'float64', 'float64', 'float64')

    elif g== 'g7':
        f= [vmin_n, vmax_n, dd_n]
        colnames= ('Vmin', 'Vmax', 'EW_dV')
        datatype= ('float64', 'float64', 'float64')

    elif g== 'g8':
        f= [ew_n, dv_n, dd_n]
        colnames= ('EW', 'dV', 'EW_dV')
        datatype= ('float64', 'float64', 'float64')

    elif g== 'g9':
        f= [ew_n, dv_n]
        colnames= ('EW', 'dV')
        datatype= ('float64', 'float64')

    elif g== 'g11':
        f= [ew_n, vmax_n, dp_n]
        colnames= ('EW', 'Vmax', 'BI0_dV')
        datatype= ('float64', 'float64', 'float64')

    elif g== 'g12':
        f= [ew_n, vmin_n, vmax_n, dp_n]
        colnames= ('EW', 'Vmin', 'Vmax', 'BI0_dV')
        datatype= ('float64', 'float64', 'float64', 'float64')

    elif g== 'g13':
        f= [bi_n, vmax_n, dp_n]
        colnames= ('BI0', 'Vmax', 'BI0_dV')
        datatype= ('float64', 'float64', 'float64')

    elif g== 'g14':
        f= [bi_n, vmin_n, vmax_n, dp_n]
        colnames= ('BI0', 'Vmin', 'Vmax', 'BI0_dV')
        datatype= ('float64', 'float64', 'float64', 'float64')

    elif g== 'g15':
        f= [bi_n, vmax_n, dv_n]
        colnames= ('BI0', 'Vmax', 'dV')
        datatype= ('float64', 'float64', 'float64')

    elif g== 'g16':
        f= [bi_n, ew_n, vmin_n, vmax_n, dv_n, dp_n]
        colnames= ('BI0', 'EW', 'Vmin', 'Vmax', 'dV', 'BI0_dV')
        datatype= ('float64', 'float64', 'float64', 'float64', 'float64', 'float64')
    

    qs= np.column_stack(param for param in f) # 2D array to do clustering on

    #do the clustering

    kmeans= KMeans(init= 'k-means++', n_clusters= k, n_init= 10)
    kmeans.fit(qs)
    labels= kmeans.predict(qs)
    sc= metrics.silhouette_score(qs, labels)
    cntrs= kmeans.cluster_centers_

    ## file to save centroids and silhouette scores for every run
    param_f= open("./clusters/"+g+"/"+clstr_name+"param.txt", 'wr')

    param_f.write(g+", \t features: "+str(colnames)+ "\n")

    param_f.write("Silhouette score= "+ str(sc)+ "\n")

    param_f.write("centroids: "+"\n"+str(cntrs)+ "\n")

    param_f.close()

    print "Silhouette score= ", sc
    print "centroids: "+"\n", cntrs

    # save the results in a FITS table

    # table with clustering results
    clstr_tbl= Table([qs[:,b] for b in range(len(f))], names= colnames, dtype= datatype)

    print "table has", len(clstr_tbl), "objects"

    lab= Column(name= 'label', data= labels)

    clstr_tbl.add_columns([lab, names, redshift])

    clstr_tbl.write("./clusters/"+g+"/"+clstr_name+"clstrs.fits", format= 'fits')


    return
   

