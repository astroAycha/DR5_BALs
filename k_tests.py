""" calculate parameters to determine good values for K. One is sum of squared error and the other is the Sillhouette score. Do that for a range of K 3-7 for each line and plot the results
    """

import numpy as np

from astropy.io import fits
from astropy.table import Table, join

from sklearn.cluster import KMeans
from sklearn import metrics

def k_test(ff):
    
    '''
    param
    ff: number of features
    '''
    
    data= Table.read('myBALCat_xtra.csv')
    
    #selec the sample: line has an absorption trough: BI0 >0 , S/N >3, and a redshift cutoff to restrict the bandwidth
    line_z= [['CIV', 1.79, 3.7], ['SiIV', 1.79, 3.7], ['AlIII', 1.79, 3.7], ['MgII', 1.1, 2.2]]
    
    
    # lists to store silhouette scores
    sil=[]
    
    # lists of the sum of distances squared
    sos =[]


    for ss in line_z:
        line= ss[0]
        z1= ss[1]
        z2= ss[2]
        
        if line== "MgII":
            lum= "logF2500"
        
        else:
            lum= "logF1400"

        s= data[(data['BIO_'+line] >0) & (data['SN1700'] >3) & \
                (data['z'] >z1) & (data['z'] <z2) & (data[lum] !=-999)]


        #features
        s['z'].fill_value= -999
        redshift= s['z'].filled() #redshift
        s['BI_'+line].fill_value= -999
        bi0= s['BIO_'+line].filled() # modified balnicity: integration 0-25000
        s['EW_'+line].fill_value= -999
        ew= s['EW_'+line].filled() # restframe absorption EW
        s['Vmin_'+line].fill_value= -999
        vmin= s['Vmin_'+line].filled() # minimum velocity
        s['Vmax_'+line].fill_value= -999
        vmax= s['Vmax_'+line].filled() # maximum velocity
    
    
        s[lum].fill_value= -999
        cl= s[lum].filled() # Log of 1400 or 2500 monochromatic luminosity
    
        #standardize parameters before using them in clustering
        f1= (ew - mean(ew))/std(ew)
        f2= (vmin - mean(vmin))/std(vmin)
        f3= (vmax - mean(vmax))/std(vmax)
        f4= (cl - mean(cl))/std(cl)
    
        # list of features to be used in clustering
        
        if ff== 4:
            f= [f1, f2, f3, f4]
        
        if ff== 3:
            f= [f1, f2, f3] #use only vmin, vmax, ew for clustering
        
        qs= np.column_stack(param for param in f) # 2D array to do clustering on
    
    
        for k in range(3,8): #do the clustering on a range of K's
    
            kmeans= KMeans(init= 'k-means++', n_clusters= k, n_init= 10)
            kmeans.fit(qs)
            labels= kmeans.predict(qs)
            m1= metrics.silhouette_score(qs, labels)
            m2= kmeans.inertia_
            sil.append(m1)
            sos.append(m2)


    s1= array(sil)
    ss1= reshape(s1,(4,5))
    
    s2= array(sos)
    ss2= reshape(s2, (4,5))
    
    print ss1
    print ss2
    
    #plot results
    
    fig= figure(figsize=(8,12))
    subplots_adjust(hspace= 0.15)

    sns.set(font_scale= 1.4)
    sns.set_style("ticks", {'font.family': u'serif'})


    ax1= fig.add_subplot(211)

    #nullfmt   = NullFormatter()
    #ax1.xaxis.set_major_formatter(nullfmt)

    ax1.plot(range(3,8), ss2[0], marker= 'D', color= '0.1', ls='--', label= 'CIV, N= 2683')
    ax1.plot(range(3,8), ss2[1], marker= '*', color= '0.3', ls='-.', label= 'SiIV, N= 806')
    ax1.plot(range(3,8), ss2[2], marker= 'o', color= '0.5', ls=':', label= 'AlIII, N= 191')
    ax1.plot(range(3,8), ss2[3], marker= 'v', color= 'k', ls='-', lw=0.7, label= 'MgII, N= 108')

    ylabel(r'Sum of squares')
    xlim(2.9, 7.1)

    legend(numpoints=1)

    ax2= fig.add_subplot(212)

    ax2.plot(range(3,8), ss1[0], marker= 'D', color= '0.1', ls='--')
    ax2.plot(range(3,8), ss1[1], marker= '*', color= '0.3', ls='-.')
    ax2.plot(range(3,8), ss1[2], marker= 'o', color= '0.5', ls=':')
    ax2.plot(range(3,8), ss1[3], marker= 'v', color= 'k', ls='-', lw=0.7)

    xlabel(r'$K$')
    ylabel('Silhouette score')
    xlim(2.9, 7.1)
    
    return


