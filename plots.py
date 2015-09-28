""" plots similar to what is in lstr_plots.py but ordered according to the number of objects in each cluster 
May 8 2015"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from sklearn.cluster import KMeans
from mpl_toolkits.mplot3d import Axes3D
from operator import itemgetter
from glob import glob
import seaborn as sns


#sns.set_style('ticks')
#sns.set_context("paper", font_scale=2)




""" 2D scatter plots for clusters. read files from saved 2d numpy arrays
    the aray has the features in the "columns". the last column has the labels indicating which cluster each sample (row) belongs to.
    line: the emission line used in clustering: MgII, CIII], or CIV
    cluster: line_ew_hwhm -in this case. can refere to any other cluster name (eg, mg2_ew_hwhm)
    k: number of clusters used in creating the loaded numpy array. The labels column in the array goes from 0 to k-1
    feature 1: one of the features used in the clustering analysis. This will be plotted on the x-axis
    feature 2: to be plotted on the y-axis
    
    in the ew_hwhm clusters, the features are EW, BHWHM, and RHWHM given in the 0th, 1st and 2nd vectors on the array.
    
    """

def two_d_scatter(line, cluster, k, feature1, feature2, feature3):
    
    """ line: CIV, CIII, MGII
        cluster: e.g., c4_ew_hwhm or mg2_ew_hwhm
        k: number of clusters
        features: EW, BHWHM, RHWHM
        """


    #clstr_name= "./clusters/"+cluster+"_"+str(k)+"clstrs.npy"
    clstr= np.load(cluster)
    
    ew= clstr[:,0] #EW
    bhwhm= clstr[:,1] #BHWHM
    rhwhm= clstr[:,2] #RHWHM
    c_label= clstr[:,3] #label
    
    clstr_length=[]

    for c in range(0,k):
        clstr_length.append([c,len(clstr[c_label==c])])

    ordered_clstrs =sorted(clstr_length, key=itemgetter(1), reverse= True)
    print ordered_clstrs
    
    clr_ls= ['orange', 'navy', 'mediumvioletred','seagreen', 'khaki', 'cornflowerblue', 'brown' , 'olive', 'purple']

    fig= figure(figsize(8,8))
    ax= fig.add_subplot(111)

    text(0.05, 0.1, 'K= '+ str(k)+ ', Clustering features:'+' '+line+' ('+ feature1+','+feature2+','+ feature3+')',
         horizontalalignment='left', verticalalignment='center',
         transform=ax.transAxes, color= 'black')

    xlabel(line+" "+feature2+ " (km/s)")
    ylabel(line+" "+feature3+ " (km/s)")

    t=0
    for o in ordered_clstrs:
        print o[0], o[1]
        
        ax.scatter(bhwhm[c_label==o[0]], rhwhm[c_label==o[0]], c=clr_ls[t], label=str(o[1]), alpha= 0.7)

        t+=1

    legend(scatterpoints=1)

############

def twoD_cluster_kde(cluster_array, line):

    """ plot KDE for the clusters for the BHWHM, RHWHM
    still need to do a lot of work on it. but this should do for now
    cluster_array: a 2D numpy array with the clustering results
    line: a string with the line used for clustering. To be used as a title for the figure."""

    clstr= np.load(cluster_array)
    
  #  cmap_ls=['OrRd', 'PuBu', 'Purples', 'BuGn', 'RdPu', 'gray_r'] 'YlOrBr'
    cmap_ls= ['YlOrBr', 'Blues', 'RdPu', 'Greens', 'Greys', 'Reds']

    sns.set_style("ticks", {'font.family': u'sans-serif'})
   # sns.set(font_scale=1.5)
    
    fig= figure(figsize=(12, 12))
    ax= fig.add_subplot(111)
    
    xlabel('BHWHM (km/s)', fontsize=18)
    #xlabel(r'EW ($\AA$)', fontsize=18)
    ylabel('RHWHM (km/s)', fontsize=18)
    
    #xlim(0,50)
    xlim(0,8000)
    ylim(0,8000)
    
    x, y= [], []
    
    k_ls=[]
    
    clr_ls= ['orange', 'navy', 'mediumvioletred','seagreen', '0.5', 'red', 'cornflowerblue', 'brown' , 'olive', 'purple']
    
    for i in range(max(clstr[:,3].astype(int))+1):
    
        k_ls.append([i, (mean(clstr[:,2][clstr[:,3]==i]))])
        #k_ls.append([i, (mean(clstr[:,1][clstr[:,3]==i]) + mean(clstr[:,2][clstr[:,3]==i]))])
    
    #ord_k_ls= sorted(k_ls, key= itemgetter(1), reverse= True)
    
    ord_k_ls= sorted(k_ls, key= itemgetter(1))
    
    print ord_k_ls
    
    clstr_label= [ line+'-a', line+'-b', line+'-c', line+'-d', line+'-e', line+'-f']
    
    u=1
    for j in range(len(ord_k_ls)):
        k= ord_k_ls[j][0]
        print k
        
        u-=0.04
        x =mean(clstr[:,1][clstr[:,3]==k])
        y =mean(clstr[:,2][clstr[:,3]==k])
        n= len(clstr[:,2][clstr[:,3]==k])
        
        
        sns.kdeplot(clstr[:,1][clstr[:,3]==k], clstr[:,2][clstr[:,3]==k], shade=True, shade_lowest=False, alpha= 0.5, cmap= cmap_ls[j]) #n_levels= 5
        
        #scatter(x,y, marker= 'x', c='r', s=60)
    
        text(x, y, clstr_label[j], fontsize= 16, family= 'serif')
        text(0.05, u,  clstr_label[j]+", N="+str(n), transform=ax.transAxes, color= clr_ls[j], fontsize= 16, family= 'serif')

        ## now overplot BAL quasars

        b= Table.read('sample_myflags_BAL_only.fits')
        scatter(b['BHWHM_'+line], b['RHWHM_'+line], marker='o', s=3, color='0.7', alpha=0.5)
        #sns.kdeplot(b['BHWHM_'+line], b['RHWHM_'+line], shade=False, shade_lowest=False)

        #  plot 1:1 line
        z= arange(11000)
        plot(z,z,'k-')


#############################

def four_pan_cluster_kde(line, sample_name):
    
    """ plot KDE for the clusters for the BHWHM, RHWHM in 4 panels for k=3, 4, 5, and 6

        line: c3, c4, or mg2
        sample_name: main (no BALs), mixed (non BALs and BALs), bal (BALs only)
        """
    
    
    if line == "c3":
        line_name= "CIII"
        line_label= "CIII]"
        xlimits= (0, 7250)
        ylimits= (0, 10750)
    
    elif line== "c4":
        line_name= line_label= "CIV"
        xlimits= (0, 6250)
        ylimits= (0, 5750)

    elif line== "mg2":
        line_name= line_label = "MgII"
        xlimits= (0, 8750)
        ylimits= (0, 8250)

    if sample_name== "main":
        sample= "_ew_hwhm_"
        tab_name= "sample_myflags.fits"
    elif sample_name== "mixed":
        sample= "_ew_hwhm_mixed_"
        tab_name= "sample_mixed_myflags.fits"
    elif sample_name == "bal":
        sample="_ew_hwhm_bal_"
        tab_name= "sample_bal_myflags.fits"

    #t= Table.read(tab_name)

    
    sns.set_style("ticks", {'font.family': u'sans-serif'})
    
    fig= plt.figure(figsize=(12,10))
    fig.subplots_adjust(wspace=0.001, hspace=0.001)
    sns.set_style("ticks")
    fig1= fig.add_axes([0., 0., 1, 1])
    fig1.set_axis_off()
    fig1.set_xlim(0, 1)
    fig1.set_ylim(0, 1)
    fig1.text(.05, .5, line_label+ " RHWHM (km/s)", rotation='vertical' \
            , horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')
    fig1.text(.5, .03, line_label+ " BHWHM (km/s)", rotation='horizontal'\
              , horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')
    
    props= dict(boxstyle='round', alpha=0.5, color='w')
    
    cmap_ls= ['YlOrBr', 'Blues', 'RdPu', 'Greens', 'Greys', 'Reds']
    
    nullfmt   = NullFormatter()
    
    #ax1= fig.add_subplot(321)
    #sns.kdeplot(t['BHWHM_'+line_name], t['RHWHM_'+line_name] \
            ,cmap= 'summer', shade=True, shade_lowest=False, alpha=0.6)
    #xlim(xlimits)
    #ylim(ylimits)
    
    
    for (i,j) in zip(range(1,5), range(3,7)):
        
        ax= fig.add_subplot(3,2,i)
        
        if i ==1:
            ax.xaxis.set_major_formatter(nullfmt)
        
        elif i== 2:
            ax.xaxis.set_major_formatter(nullfmt)
            ax.yaxis.set_major_formatter(nullfmt)

        elif i ==4:
            ax.yaxis.set_major_formatter(nullfmt)
        
        text(.07, .9, "K="+str(j), transform=ax.transAxes \
            , horizontalalignment='center', verticalalignment='center', fontsize= 12, family= 'serif')
        
        z= arange(11000)
        plot(z,z,'k-', lw=.5) #plot 1:1 line
        
        clstr_name= "./clusters/"+line+sample+ str(j) +"clstrs.npy"
        
        print clstr_name
        
        clstr_array= np.load(clstr_name)
        
        xlim(xlimits)
        ylim(ylimits)
    
        clstr_num= []
    
        for k in range(max(clstr_array[:,3].astype(int))+1):
        
            #clstr_num.append([k, (len(clstr_array[clstr_array[:,3]== k]))])
            clstr_num.append([k, (mean(clstr_array[:,1][clstr_array[:,3]== k]) \
                               , (mean(clstr_array[:,2][clstr_array[:,3]== k])))])
    
        ordered_clstrs= sorted(clstr_num, key= itemgetter(1)) # reverse= True
        print ordered_clstrs
        
        ew=[]
        x, y= [], []
        n=[]

        cc= -1
        for c in ordered_clstrs:
        
            cc+=1
            
            if (min(clstr_array[:,1][clstr_array[:,3]==c[0]]) >0) & (min(clstr_array[:,2][clstr_array[:,3]==c[0]]) >0):
            
                sns.kdeplot(clstr_array[:,1][clstr_array[:,3]==c[0]], clstr_array[:,2][clstr_array[:,3]==c[0]] \
                ,cmap= cmap_ls[cc], shade=True, shade_lowest=False, alpha=0.6)
            
            ew.append(mean(clstr_array[:,0][clstr_array[:,3]==c[0]]))
            x.append(mean(clstr_array[:,1][clstr_array[:,3]==c[0]]))
            y.append(mean(clstr_array[:,2][clstr_array[:,3]==c[0]]))
            n.append(len(clstr_array[clstr_array[:,3]==c[0]]))
        
        clstr_label= ['a'+str(j),'b'+str(j),'c'+str(j),'d'+str(j),'e'+str(j),'f'+str(j)]
        clr_ls= ['orange', 'navy', 'mediumvioletred','seagreen', '0.5', 'red'] # 'cornflowerblue', 'brown' , 'olive', 'purple']
        
        #if j ==6:
            
        u=0.4
        for l in range(j):
            u-=0.06
                
            text(x[l]+150, y[l]-150, clstr_label[l], color= 'k', fontsize= 14, family= 'serif') #, bbox=props
            
            text(0.67, u,  line_label+"-"+clstr_label[l]+", N="+str(n[l]), transform=ax.transAxes \
            , color= clr_ls[l], fontsize= 12, family= 'serif')
            
        scatter(x,y, marker='D', color='w', s=[e for e in ew])

    return

        ## now overplot BAL quasars
        
        #b= Table.read('sample_myflags_BAL_only.fits')
        #scatter(b['BHWHM_'+line_name], b['RHWHM_'+line_name], marker='o', s=3, color='0.3', alpha=0.5)
        #sns.kdeplot(b['BHWHM_'+line_name], b['RHWHM_'+line_name], shade=False, shade_lowest=False)


##########################

## original script forked from Pauline
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
#from matplotlib.ticker import NullFormatter
import matplotlib.ticker as ticker
from seaborn import kdeplot

# based on http://www.astrobetter.com/blog/2014/02/10/visualization-fun-with-python-2d-histogram-with-1d-histograms-on-axes/

def kde_hist(line, sample_name,j):
    '''plot KDE for clusters for one value of K
    param:
    line: 'c3', 'c4', 'mg2'
    sample_name: 'main' for main sample with no BALs, 'mixed' for the sample with BALs and no BALs, 'bal' for BALs only
    j: number of clusters: 3, 4, 5, 6
       
        '''
    sns.set(font_scale= 1.5)
    sns.set_style("ticks", {'font.family': u'serif'})
    props= dict(boxstyle='round', alpha=0.9, color='w')
    
    clstr_name= "./clusters/"+line+sample+ str(j) +"clstrs.npy"
    
    # start with a rectangular Figure
    f = plt.figure(figsize=(12,12))
    
    # define where the axes go
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02
    
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.15, height]
    
    # add the axes to the figure
    ax2d = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    
    # no labels for the sidecar histograms, because the 2D plot has them
    nullfmt   = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    
    # the 2D plot:
    # note the all-important transpose!
    clstr_array= np.load(clstr_name)
    
    clstr_num= []
        
    for k in range(max(clstr_array[:,3].astype(int))+1):
        
        clstr_num.append([k, (mean(clstr_array[:,1][clstr_array[:,3]== k]), (mean(clstr_array[:,2][clstr_array[:,3]== k])))])
        
    ordered_clstrs= sorted(clstr_num, key= itemgetter(1)) # reverse= True
    print ordered_clstrs

    cmap_ls= ['YlOrBr', 'Blues', 'RdPu', 'Greens', 'Greys', 'Reds']

    ew, x, y, n, =[],[],[],[]
    cc= -1
    for c in ordered_clstrs:
        cc+=1
        
        if (min(clstr_array[:,1][clstr_array[:,3]==c[0]]) >0) & (min(clstr_array[:,2][clstr_array[:,3]==c[0]]) >0):
            sns.kdeplot(clstr_array[:,1][clstr_array[:,3]==c[0]], clstr_array[:,2][clstr_array[:,3]==c[0]] \
                        ,cmap= cmap_ls[cc], ax=ax2d, shade=True, shade_lowest=False, alpha=0.6)
        
        ew.append(mean(clstr_array[:,0][clstr_array[:,3]==c[0]]))
        x.append(mean(clstr_array[:,1][clstr_array[:,3]==c[0]]))
        y.append(mean(clstr_array[:,2][clstr_array[:,3]==c[0]]))
        n.append(len(clstr_array[clstr_array[:,3]==c[0]]))
    
    ax2d.scatter(x,y, marker='D', color='w', s=[e for e in ew])
        
    clstr_label= ['a'+str(j),'b'+str(j),'c'+str(j),'d'+str(j),'e'+str(j),'f'+str(j)]
    clr_ls= ['orange', 'navy', 'mediumvioletred','seagreen', '0.5', 'red'] # 'cornflowerblue', 'brown' , 'olive', 'purple']
    u=0.3
    for l in range(j):
        u-=0.04
            
        ax2d.text(x[l]+150, y[l]-150, clstr_label[l], color='r', fontsize= 14 ) #, bbox=props
        ax2d.text(0.67, u,  line_label+"-"+clstr_label[l]+", N="+str(n[l]) \
                  , transform=ax2d.transAxes, fontsize= 13, color= clr_ls[l])

    ax2d.set_xlabel(line_label+ " BHWHM (km/s)" )
    ax2d.set_ylabel(line_label+ " RHWHM (km/s)")
    ax2d.set_xlim(0,11000)
    ax2d.set_ylim(0,11000)

    b= Table.read("sample_bal_myflags.fits")
    ax2d.scatter(b['BHWHM_'+line_name], b['RHWHM_'+line_name], marker='o', s=1, color='0.5', label="BAL Quasars")
    z= arange(11000)
    ax2d.plot(z,z,'k-', lw=.5) #plot 1:1 line
    
    # the 1-D histograms: first the X-histogram
    sns.kdeplot(clstr_array[:,1], ax=axHistx, shade= False,  color='k', label= sample_label+" Sample"+"\n"+"N="+str(len(clstr_array)))
    sns.kdeplot(b['BHWHM_'+line_name], ax=axHistx, shade= False, lw= 2, color='c', label= "BALQ"+"\n"+"N="+str(len(b)))
    
    axHistx.set_xlim( ax2d.get_xlim()) # x-limits match the 2D plot
    axHistx.set_ylabel(line_label+' BHWHM')
    axHistx.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1.0e'))
    axHistx.yaxis.set_major_locator(ticker.MultipleLocator(base=0.0003)) # this locator puts ticks at regular intervals
    #gca().xaxis.set_major_locator(MaxNLocator(nbins=3, prune= 'both'))
    #axHistx.set_yticks([500, 1000, 1500])
    #axHistx.set_ylim(0,1000)
    axHistx.legend(prop={'size':12})
        
    # then the Y-histogram
    sns.kdeplot(clstr_array[:,2], ax=axHisty, vertical= True, shade= False, color='k', label= sample_label+" Sample"+"\n"+"N="+str(len(clstr_array)))
    sns.kdeplot(b['RHWHM_'+line_name], ax=axHisty, vertical= True, shade= False, lw= 2, color='c', label= "BALQ"+"\n"+"N="+str(len(b)))

    axHisty.set_ylim(ax2d.get_ylim()) # y-limits match the 2D plot
    axHisty.set_xlabel(line_label+' RHWHM')
    axHisty.xaxis.set_major_formatter(ticker.FormatStrFormatter('%1.0e'))
    axHisty.xaxis.set_major_locator(ticker.MultipleLocator(base=0.0003)) # this locator puts ticks at regular intervals
    #gca().yaxis.set_major_locator(MaxNLocator(nbins=3, prune= 'both'))
    #axHisty.set_xticks([500, 1000, 1500])
    #axHisty.set_xlim(0,1000)
    axHisty.legend(prop={'size':12})
    
    plt.show()
    return


########################

def plot_spec_three_pans(line, k):
    
    """ plot composite spectra in 4 panels:
        panel 1: Ly alpha, Si IV, C IV, (He II & OIII])
        panel 2: Al III, Si III], C III]
        panel 3: Mg II
        param:
        line: "SiIV", "CIV", "AlIII", or "MgII"
        k: number of clusters (3, 4, ...)
        """
    
    clstr_name= "./clusters/"+line+str(k)+"clstrs.fits"
    clstr= Table.read(clstr_name)
    
    clstr_num=[]
    for f in range(k):
        clstr_num.append([f, (mean(clstr['Vmin'][clstr['label']== f]), (mean(clstr['Vmax'][clstr['label']== f])))])
        
    ordered_clstrs= sorted(clstr_num, key= itemgetter(1)) #reverse= True
    print ordered_clstrs

    
    compo_list= []
    for r in ordered_clstrs:
        compo_name= "./composites/"+line+"_"+str(k)+"clstr"+str(r[0]+1)+".fits"
        spec= fits.open(compo_name)
        num_obj= spec[0].header['SPEC_NUMBER']
        compo_list.append([compo_name, num_obj])
    
    print compo_list


    fig= figure(figsize=(14,8))
    sns.set_style("ticks")
    fig1= fig.add_axes([0., 0., 1, 1])
    fig1.set_axis_off()
    fig1.set_xlim(0, 1)
    fig1.set_ylim(0, 1)
    fig1.text(.07, 0.5, r"Normalized Flux (erg s$^{-1}$ cm$^{-1}$ $\AA^{-1}$)", rotation='vertical', \
             horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')
    fig1.text(0.5, 0.01, r"Wavelength ($\AA$)", rotation='horizontal', \
              horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')

    
    line_mark= [[1215.7, 1240, 1305, 1335, 1396.8, 1549, 1640, 1663.5], \
                [1857, 1892, 1908], [2800]]
    line_labels= [[r'Ly$\alpha$', 'NV', 'OI + SiII', 'CII', 'SiIV', 'CIV', 'HeII', 'OIII]'], \
                  ['AlIII', 'SiIII]', 'CIII]'], ['MgII']]
    
    alphabet_list = ['a'+str(k), 'b'+str(k), 'c'+str(k), 'd'+str(k), 'e'+str(k), 'f'+str(k)]
    compo_labels= [line+"-"+ a for a in alphabet_list]
    
    clr_ls= ['orange', 'navy', 'mediumvioletred','seagreen', '0.5', 'red', 'cornflowerblue', 'brown' , 'olive', 'purple']
    
    splt_ls=[211, 223, 224]
    dx_ls= [(1150, 1700), (1800, 2000),  (2750, 2850)]
    dy_ls= [(0.3, 2), (0.75, 1.8),  (0.85, 1.8)]
    
    for s in range(3):
    
        ax= fig.add_subplot(splt_ls[s])
        xlim(dx_ls[s])
        ylim(dy_ls[s])
        
        for t in range(len(line_mark[s])):
            ax.axvline(line_mark[s][t], ls=':', c='k')
            ax.text(line_mark[s][t]-10, dy_ls[s][1]-(dy_ls[s][1]-dy_ls[s][0])/10, \
                    line_labels[s][t], rotation= 'vertical', fontsize= 14, family='serif')
        
        ii= dy_ls[s][1]
        for (sp, clr, clab) in zip(compo_list, clr_ls, compo_labels):
            n= sp[1]
            spec= fits.open(sp[0])
            wlen= spec[0].data[0]
            flx= spec[0].data[1]
            
            if n >25:
                plot(wlen, flx/flx[(dx_ls[s][0]-1100)*2], c= clr, lw= 2)
            
            ii-=0.1
            ax.text(1925, ii, clab+", N="+ str(n), color= clr, fontsize= 14, family= 'serif') #bbox=props

        #mean_compo= fits.open("./composites/mean_compo_"+sample_name+".fits")
        #mean_flx= mean_compo[0].data[1]
        
        #plot(wlen, mean_flx/mean_flx[(dx_ls[s][0]-1100)*2], c='k', lw=2, label= "Mean")
        
        #ax.text(2752, dy_ls[s][1]-0.15, sample_label+" Sample, "+"N="+ str(len(clstr_array))+"\n"+"Mean Composite ", color= 'k', fontsize= 14, family= 'serif')

    return


###########################
def plot_reprod(line, k):
    
    """ make plots with cluster centroids calculated 50 times to show reproducibility.
        read values from text files.
        """
    
    #lines = [('CIV', 3), ('CIV', 4), ('CIV', 5), ('CIII', 3), ('CIII', 4), ('CIII', 5), ('MGII', 3), ('MGII', 4), ('MGII', 5)]
    #for the BAL sample:lines = [('CIII', 3), ('CIII', 4), ('CIII', 5), ('CIII', 6), ('MGII', 3), ('MGII', 4), ('MGII', 5)]
    
    
    cntrs= loadtxt(line+str(k)+"_bal.txt") #sample with BALs
    #cntrs= loadtxt(line+str(k)+".txt") #use for the non-BAL sample
    
    fig= figure(figsize=(10,8))
    
    subplots_adjust(hspace = .05)
    ax1= fig.add_subplot(311)
    xlim(-4, 54)
    ylabel('EW '+line)
    ax1.set_xticks([10, 20, 30, 40, 50], [10, 20, 30, 40, 50])
    gca().yaxis.set_major_locator(MaxNLocator(nbins=7, prune= 'both'))
    
    for m in range(0, k*3, 3):
        ax1.scatter(range(50), cntrs[:, m], marker='s', edgecolor='k', facecolor='0.5')
    
    
    ax2= fig.add_subplot(312, sharex= ax1)
    xlim(-4, 54)
    ylabel("BHWHM "+line)
    gca().yaxis.set_major_locator(MaxNLocator(nbins=7, prune= 'both'))
    for m in range(1, k*3, 3):
        ax2.scatter(range(50), cntrs[:, m], marker='o', edgecolor='k', facecolor='0.5')
    
    ax3= fig.add_subplot(313, sharex= ax1)
    xlim(-4, 54)
    ylabel('RHWHM '+line)
    xlabel("Number of Repeats")
    gca().yaxis.set_major_locator(MaxNLocator(nbins=7, prune= 'both'))
    for m in range(2, k*3, 3):
        ax3.scatter(range(50), cntrs[:, m], marker='^', edgecolor='k', facecolor='0.5')


##########################

balt= Table.read('sample_bal_myflags.fits')



fig= figure()

ax1= fig.add_subplot(331)
ax1.hist(balt['BI_CIV'])
xlabel('BI_CIV')

ax2= fig.add_subplot(332)
ax2.hist(balt['VMIN_CIV_2000'])
xlabel('VMIN_CIV_2000')

ax3= fig.add_subplot(333)
ax3.hist(balt['VMAX_CIV_2000'])
xlabel('VMAX_CIV_2000')

ax4= fig.add_subplot(334)
ax4.hist(balt['AI_CIV'])
xlabel('AI_CIV')

ax5= fig.add_subplot(335)
ax5.hist(balt['VMIN_CIV_450'])
xlabel('VMIN_CIV_450')

ax6= fig.add_subplot(336)
ax6.hist(balt['VMAX_CIV_450'])
xlabel('VMAX_CIV_450')

ax7= fig.add_subplot(337)
ax7.hist(balt['REW_SIIV'])
xlabel('REW_SIIV')

ax8= fig.add_subplot(338)
ax8.hist(balt['REW_CIV'])
xlabel('REW_CIV')

ax9= fig.add_subplot(339)
ax9.hist(balt['REW_ALIII'])
xlabel('REW_ALIII')


###################
def bal_hist1(k):

    """make histograms for absorption properties in the BALQ sample
        k: number of clusters
        """
    sns.set_style('ticks', {'font.family': u'serif', 'xtick.direction': u'in'})

    t=Table.read('sample_bal_myflags.fits')

    cluster_array= np.load("./clusters/c3_ew_hwhm_bal_"+str(k)+"clstrs_name.npy")
    
    clstr_num= []
    for l in range(k):
        clstr_num.append([l, (mean(cluster_array[:,1].astype(float)[cluster_array[:,3].astype(int)== l]), \
                             (mean(cluster_array[:,2].astype(float)[cluster_array[:,3].astype(int)== l])))])

    ordered_clstrs= sorted(clstr_num, key= itemgetter(1))

    c_labels= []
    for ll in ordered_clstrs:
        c_labels.append(ll[0])
    print c_labels

    clr_ls= ['orange', 'navy', 'mediumvioletred','seagreen', '0.5', 'red', 'cornflowerblue', 'brown' , 'olive', 'purple']
    alphabet_list = ['a'+str(k), 'b'+str(k), 'c'+str(k), 'd'+str(k), 'e'+str(k), 'f'+str(k)]
    compo_labels= [line_label+"-"+ a for a in alphabet_list]


    bals= t[t['SDSS_NAME'] == cluster_array[:,4]]

    fig= figure(figsize=(14,10))

    ax1= fig.add_subplot(221)
    xlabel(r'EW SiIV abs trough ($\AA$)')
    sns.kdeplot(bals['REW_SIIV'], ax= ax1, color= 'k', lw=3, legend=False)

    for (i,cl, clr) in zip(c_labels, compo_labels, clr_ls):
        bals_c= t[(t['SDSS_NAME'] == cluster_array[:,4]) & (cluster_array[:,3].astype(int) == i)]
        sns.kdeplot(bals_c['REW_SIIV'], ax=ax1, color= clr, legend=False)


    ax2= fig.add_subplot(222)
    xlabel(r'EW CIV abs trough ($\AA$)')
    sns.kdeplot(bals['REW_CIV'], ax= ax2, label= "BALQ Sample, N="+str(len(t)), color= 'k', lw=3)
    
    for (i,cl, clr) in zip(c_labels, compo_labels, clr_ls):
        bals_c= t[(t['SDSS_NAME'] == cluster_array[:,4]) & (cluster_array[:,3].astype(int) == i)]
        sns.kdeplot(bals_c['REW_CIV'], ax= ax2, label= cl+", N= "+str(len(bals_c['REW_CIV'])), color= clr)


    ax3= fig.add_subplot(223)
    xlabel(r'EW AlIII abs trough ($\AA$)')
    sns.kdeplot(bals['REW_ALIII'], ax= ax3, color= 'k', lw=3, legend=False)
    
    for (i,cl, clr) in zip(c_labels, compo_labels, clr_ls):
        bals_c= t[(t['SDSS_NAME'] == cluster_array[:,4]) & (cluster_array[:,3].astype(int) == i)]
        sns.kdeplot(bals_c['REW_ALIII'], ax=ax3, color= clr, legend=False)


    return

##=====================

#plot redshift histograms to show dist for each line.

bals= Table.read('myBALCat.fits')

fig= figure(figsize=(10,8))

ax= fig.add_subplot(111)

hist(bals['z'][bals['BIO_CIV'] >0], histtype= 'step', normed= True, lw= 2, label= 'CIV')
hist(bals['z'][bals['BIO_SiIV'] >0], histtype= 'step', normed= True, lw= 2, label= 'SiIV')
hist(bals['z'][bals['BIO_AlIII'] >0], histtype= 'step', normed= True, lw= 2, label= 'AlIII')
hist(bals['z'][bals['BIO_MgII'] >0], histtype= 'step', normed= True, lw= 2, label= 'MgII')
ylim(0,1.3)

legend()







