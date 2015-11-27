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



## original script forked from Pauline
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
#from matplotlib.ticker import NullFormatter
import matplotlib.ticker as ticker
from seaborn import kdeplot

# based on http://www.astrobetter.com/blog/2014/02/10/visualization-fun-with-python-2d-histogram-with-1d-histograms-on-axes/

def kde_hist(line,k):
    '''plot KDE for clusters for one value of K
    param:
    line: 'SiIV', 'CIV', 'MgII'
    
    k: number of clusters: 3, 4, 5
       
        '''
    sns.set(font_scale= 1.5)
    sns.set_style("ticks", {'font.family': u'serif'})
    props= dict(boxstyle='round', alpha=0.9, color='w')
    
    clstr_tab= Table.read("./clusters/"+line+ str(k) +"clstrs.fits") #table with clustering parameters for this line
    
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
    clstr_num= []
        
    for j in range(k):
        
        clstr_num.append([j, (mean(clstr_tab['Vmin'][clstr_tab['label']== j]), (mean(clstr_tab['Vmax'][clstr_tab['label']== j])))])
        
    ordered_clstrs= sorted(clstr_num, key= itemgetter(1)) # reverse= True
    print ordered_clstrs

    cmap_ls= ['YlOrBr', 'Blues', 'RdPu', 'Greens', 'Greys', 'Reds']

    ew, x, y, n, =[],[],[],[]
    cc= -1
    for c in ordered_clstrs:
        cc+=1
        
        #if (min(clstr_tab['Vmin'][clstr_tab['label']==c[0]]) >0) & (min(clstr_tab['Vmax'][clstr_tab['label']==c[0]]) >0):
        sns.kdeplot(clstr_tab['Vmin'][clstr_tab['label']==c[0]], clstr_tab['Vmax'][clstr_tab['label']==c[0]] \
                        ,cmap= cmap_ls[cc], ax=ax2d, shade=True, shade_lowest=False, alpha=0.6)
        
        ew.append(mean(clstr_tab['EW'][clstr_tab['label']==c[0]]))
        x.append(mean(clstr_tab['Vmin'][clstr_tab['label']==c[0]]))
        y.append(mean(clstr_tab['Vmax'][clstr_tab['label']==c[0]]))
        n.append(len(clstr_tab[clstr_tab['label']==c[0]]))
    
    ax2d.scatter(x,y, marker='D', color='w', s=[e for e in ew])
        
    clstr_label= ['a'+str(k),'b'+str(k),'c'+str(k),'d'+str(k),'e'+str(k),'f'+str(k)]
    clr_ls= ['steelblue', 'olivedrab','orange', 'navy', 'mediumvioletred','seagreen', '0.5', 'red'] # 'cornflowerblue', 'brown' , 'olive', 'purple']
    u=0.3
    for l in range(k):
        u-=0.04
            
        ax2d.text(x[l]+150, y[l]-150, clstr_label[l], color='r', fontsize= 14 ) #, bbox=props
        ax2d.text(0.67, u,  line+"-"+clstr_label[l]+", N="+str(n[l]) \
                  , transform=ax2d.transAxes, fontsize= 13, color= clr_ls[l])

    ax2d.set_xlabel(line+ " Vmin (km/s)" )
    ax2d.set_ylabel(line+ " Vmax (km/s)")
    #ax2d.set_xlim(0,11000)
    #ax2d.set_ylim(0,11000)

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

###########

def clust_compos(line, k, f):

    '''plot composites in one panel
    k: number of clusters
    f: number of features
    '''
    
    #sns.palplot(sns.color_palette("Set2", 10)) draws the palette
    #sns.set_palette('Set2')
    
    if line== 'MgII':
        xlimit= (1700, 3100)
        r= arange(5,9)
    
    else:
        xlimit= (1300, 2165)
        r= arange(0,8)

    clstr_tbl= Table.read("./clusters/"+str(f)+"features/"+line+str(k)+"clstrs.fits")

    clstrs_ls=[]
    for o in range(k):
        clstrs_ls.append([o ,len(clstr_tbl[clstr_tbl['label'] ==o]),\
                          mean(clstr_tbl['Vmin'][clstr_tbl['label'] ==o]),\
                          mean(clstr_tbl['Vmax'][clstr_tbl['label'] ==o]), \
                          mean(clstr_tbl['EW'][clstr_tbl['label'] ==o]), \
                          mean(clstr_tbl['Lum'][clstr_tbl['label'] ==o])])

    ord_clstrs= sorted(clstrs_ls, key= itemgetter(2))

    print ord_clstrs


    #clr_ls= ['steelblue', 'olivedrab','orange', '0.4']
    #clr_ls= ["#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974", "#64B5CD"]
    #clr_ls= ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
    #nav= sns.light_palette("navy", as_cmap=True)
    #clrm_ls= ['Purples', 'Blues', 'Greys', 'Reds', nav, 'Greens']

    clr_ls = [sns.xkcd_rgb["windows blue"], sns.xkcd_rgb["amber"], sns.xkcd_rgb["greyish"], sns.xkcd_rgb["faded green"], sns.xkcd_rgb["pale red"], sns.xkcd_rgb["dusty purple"]]
    amb= sns.light_palette(sns.xkcd_rgb["amber"], as_cmap= True)
    
    clrm_ls= ['Blues', amb, 'Greys', 'Greens', 'Reds', 'Purples']

    
    clstr_name= ['a', 'b', 'c', 'd', 'e', 'f']
    

    fig= figure(figsize=(14,10))

    ax1= fig.add_subplot(221)
    xlabel(line + r'$V_{min}$ (km/s)')
    ylabel(line + r'$V_{max}$ (km/s)')
    
    i =1
    for c in ord_clstrs:
        l= c[0]
        print l
        if c[1] >30:
            sns.kdeplot(clstr_tbl['Vmin'][clstr_tbl['label'] == l], clstr_tbl['Vmax'][clstr_tbl['label'] == l], \
                    shade= True, shade_lowest= False, alpha= 0.5, cmap= clrm_ls[i-1], label= False)
        
        i+=1

    ax1.text(0.1, 0.8, line+" Sample"+"\n"+ "N= "+str(len(clstr_tbl)), color= 'k', fontsize= 18, transform=ax1.transAxes)

    vmin, vmax, ews, lum =[], [], [], []
    for c in ord_clstrs:
        if c[1] >30:
            vmin.append(c[2])
            vmax.append(c[3])
            ews.append(c[4])
            lum.append(c[5])

    ax1.scatter(vmin, vmax, s= [abs(e)*200 for e in ews], edgecolor= '#34495e', facecolor= 'w', marker= 'D')

    ax1.text(0.1, 0.8, line+" Sample"+"\n"+ "N= "+str(len(clstr_tbl)), color= 'k', fontsize= 18, transform=ax1.transAxes)


    for x in range(len(vmin)):
        text(vmin[x]-.1, vmax[x]-.1, clstr_name[x] , color= 'k', fontsize= 12)


    ax2= fig.add_subplot(222)
    xlabel(line + r'EW')
    ylabel(line + r'L$_{1400}$')
    
    i =1
    for c in ord_clstrs:
        l= c[0]
        print l
        if c[1] >30:
            sns.kdeplot(clstr_tbl['EW'][clstr_tbl['label'] == l], clstr_tbl['Vmax'][clstr_tbl['label'] == l], \
                        shade= True, shade_lowest= False, alpha= 0.5, cmap= clrm_ls[i-1], label= False)
        
        i+=1

    ax3= fig.add_subplot(212)
    xlim(xlimit)
    ylim(.1,3.2)
    xlabel(r'Restframe Wavelength ($\AA$)')
    ylabel(r'Normalized Flux (erg s$^{-1}$ cm$^{-1}$ $\AA^{-1}$)')

    line_mark= [1335, 1396.8, 1549, 1640, 1663.5, 1857, 1892, 1908, 2800]
    line_label= ['CII', 'SiIV', 'CIV', 'HeII', 'OIII]', 'AlIII', 'SiIII]', 'CIII]', 'MgII']

    #plot([1990,2065],[1,1], 'k-')
    ax3.arrow(2030, 1.3, -30, -.1, fc='k', ec='k')
    ax3.arrow(2030, 1.3, +30, -.1, fc='k', ec='k')
    text(2020, 1.35, r'FeIII', fontsize= 14, family='serif', color='k')
    
    for p in r:
        axvline(line_mark[p], ls= ':', color= '.5')
        text(line_mark[p], 3, line_label[p], rotation= 'vertical')
    

    i=1
    for c in ord_clstrs:
        l= c[0]
        compo_name= "./composites/"+str(f)+"features/"+line+"_"+str(k)+"clstr"+str(l+1)+".fits"
        spec= fits.open(compo_name)
        if c[1] > 20:
            plot(spec[0].data[0], spec[0].data[1]/spec[0].data[1][(2150-1100)*2], lw= 2, color= clr_ls[i-1])
        ax3.text(0.82, .9-i/15., line+"-"+clstr_name[i-1]+", N= "+str(len(clstr_tbl[clstr_tbl['label'] == l])), color= clr_ls[i-1], fontsize= 18, transform=ax3.transAxes)
        i+=1
    
    
    return

##considering adding pie chart here to show the fraction of objects with vatiablity
#pie([28,3,16,3,12], shadow= True, colors= clr_ls, startangle=30, radius= 0.6, labels=['a', 'b', 'c', 'd', 'e'])



    #prop_tbl= join(data, clstr_tbl, keys='SDSSName')
    #scatter(prop_tbl['Vmin_'+line][(prop_tbl['BI1']-prop_tbl['BI2']) !=-999], prop_tbl['Vmax_'+line][(prop_tbl['BI1']-prop_tbl['BI2']) !=-999], marker='o', s=5, color='k')



##############


def clstr_prop(line,k):

    """ 
    plot KDE of the clusters in the Vmin vs Vmax plane with the EW as a marker.
    lower panels: distributions of properties in each cluster: extinction, HeII, CF...
    
    Param: 
    line: 'CIV', 'SiIV', or 'MgII'
    k: number of clusters (3,4,5..)
    
    """

    data= Table.read('myBALCat_xtra.csv')

    clstr_tbl= Table.read("./clusters/4features/"+line+str(k)+"clstrs.fits")
    
    clstrs_ls=[]
    for o in range(k):
        clstrs_ls.append([o ,len(clstr_tbl[clstr_tbl['label'] ==o]),\
                           mean(clstr_tbl['Vmin'][clstr_tbl['label'] ==o]),\
                           mean(clstr_tbl['Vmax'][clstr_tbl['label'] ==o]), \
                           mean(clstr_tbl['EW'][clstr_tbl['label'] ==o]), \
                           mean(clstr_tbl['Lum'][clstr_tbl['label'] ==o])])
    
    
    ord_clstrs= sorted(clstrs_ls, key= itemgetter(2))
    
    print ord_clstrs

    clr_ls = [sns.xkcd_rgb["windows blue"], sns.xkcd_rgb["amber"], sns.xkcd_rgb["greyish"], sns.xkcd_rgb["faded green"], sns.xkcd_rgb["pale red"], sns.xkcd_rgb["dusty purple"]]
    amb= sns.light_palette(sns.xkcd_rgb["amber"], as_cmap= True)
    
    clrm_ls= ['Blues', amb, 'Greys', 'Greens', 'Reds', 'Purples']
    
    
    clstr_name= ['a', 'b', 'c', 'd', 'e', 'f']


    fig= figure(figsize=(16,10))

    fig1= fig.add_axes([0., 0., 1, 1])
    fig1.set_axis_off()
    fig1.set_xlim(0, 1)
    fig1.set_ylim(0, 1)
    fig1.text(.06, 0.54, r"Normalized Fraction", rotation='vertical', \
              horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')


    prop_tbl= join(data, clstr_tbl, keys='SDSSName')
    
    props= dict(boxstyle='round', facecolor='w', edgecolor='k')# , alpha=0.7)
    
    """
    
    ax1= fig.add_subplot(231)
    xlabel(line + r'$V_{min}$ (km/s)')
    xlabel(line + r'$V_{max}$ (km/s)')
    
    i =1
    for c in ord_clstrs:
        l= c[0]
        print l
        sns.kdeplot(clstr_tbl['Vmin'][clstr_tbl['label'] == l], clstr_tbl['Vmax'][clstr_tbl['label'] == l], \
                    shade= True, shade_lowest= False, alpha= 0.5, cmap= clrm_ls[i-1], label= False)
                    
        i+=1

    ax1.text(0.1, 0.8, line+" Sample"+"\n"+ "N= "+str(len(clstr_tbl)), color= 'k', fontsize= 16, transform=ax1.transAxes)

    vmin, vmax, ews=[], [], []
    for c in ord_clstrs:
        vmin.append(c[2])
        vmax.append(c[3])
        ews.append(c[4])

    scatter(vmin, vmax, s= [e*-5 for e in ews], edgecolor= '#34495e', facecolor= 'w', marker= 'D')
    
    for x in range(k):
        text(vmin[x]-1100, vmax[x]-1100, clstr_name[x] , color= 'k', fontsize= 12)

    """
    i= 0
    ax1= fig.add_subplot(321)
    for c in ord_clstrs:
        if c[1] >30:
            l= c[0]
            hist(prop_tbl['LOGEDD_RATIO_DR7'][(prop_tbl['label'] == l) & (prop_tbl['LOGEDD_RATIO_DR7'] !=-999) & (prop_tbl['LOGEDD_RATIO_DR7'] < 50000)], \
             bins=12, histtype= 'step', normed= True, color= clr_ls[i], lw= 2)
             
            i+=1
    
    ax1.text(0.65, 0.77,"log(E/Edd)", transform=ax1.transAxes, color= 'k', fontsize= 16)
    ax1.text(0.95, 0.85, "A", transform=ax1.transAxes, color= 'r', fontsize= 14, bbox= props)

    i =0
    ax2= fig.add_subplot(322)
    for c in ord_clstrs:
        if c[1] >30:
            l= c[0]
            hist(prop_tbl['E_B-V_1'][(prop_tbl['label'] == l) & (prop_tbl['E_B-V_1'] !=-999)], \
             bins=12, histtype= 'step', normed= True, color= clr_ls[i], lw= 2)
            ax2.text(0.62, .7-i/10., line+"-"+clstr_name[i]+", N= "+str(c[1]), color= clr_ls[i], fontsize= 16, transform=ax2.transAxes)
             
            i+=1
    
    ax2.text(0.7, 0.85,"E(B - V)", transform=ax2.transAxes, color= 'k', fontsize= 16)
    ax2.text(0.95, 0.85, "B", transform=ax2.transAxes, color= 'r', fontsize= 14, bbox= props)



    i =0
    ax3= fig.add_subplot(323)
    for c in ord_clstrs:
        if c[1] >30:
            l= c[0]
            hist(prop_tbl['alpha_UV_BLH'][(prop_tbl['label'] == l) & (prop_tbl['alpha_UV_BLH'] !=-999)], \
             bins=12, histtype= 'step', normed= True, color= clr_ls[i], lw= 2)
             
            i+=1

    ax3.text(0.8, 0.85,r"$\alpha_{UV}$", transform=ax3.transAxes, color= 'k', fontsize= 18)
    ax3.text(0.95, 0.85, "C", transform=ax3.transAxes, color= 'r', fontsize= 14, bbox= props)
    

    i =0
    ax4= fig.add_subplot(324)
    for c in ord_clstrs:
        if c[1] > 30:
            l= c[0]
            hist(prop_tbl['HeII_EW_BLH'][(prop_tbl['label'] == l) & (prop_tbl['HeII_EW_BLH'] !=-999) \
                                  & (prop_tbl['HeII_EW_BLH'] !=0)], bins=12, histtype= 'step', \
                                 normed= True, color= clr_ls[i], lw= 2)

            i+=1
    ax4.text(0.55, 0.82,r"EW(HeII) [$\AA$]", transform=ax4.transAxes, color= 'k', fontsize= 16)
    ax4.text(0.95, 0.85, "D", transform=ax4.transAxes, color= 'r', fontsize= 14, bbox= props)

    i =0
    ax5= fig.add_subplot(325)
    for c in ord_clstrs:
        if c[1] >30:
            l= c[0]
            hist(prop_tbl['v_md_BLH'][(prop_tbl['label'] == l) & (prop_tbl['v_md_BLH'] !=-999)], \
             bins=12, histtype= 'step', normed= True, color= clr_ls[i], lw= 2)
        
            i+=1

    ax5.text(0.65, 0.8,r"v$_{md}$ [km/s]", transform=ax5.transAxes, color= 'k', fontsize= 16)
    ax5.text(0.95, 0.85, "E", transform=ax5.transAxes, color= 'r', fontsize= 14, bbox= props)

    i =0
    ax6= fig.add_subplot(326)
    for c in ord_clstrs:
        if c[1] >30 :
            l= c[0]
            hist(prop_tbl['CF_BLH'][(prop_tbl['label'] == l) & (prop_tbl['CF_BLH'] !=-999)], \
             bins=12, histtype= 'step', normed= True, color= clr_ls[i], lw= 2)
        
            i+=1

    ax6.text(0.2, 0.85,"CF", transform=ax6.transAxes, color= 'k', fontsize= 16)
    ax6.text(0.05, 0.85, "F", transform=ax6.transAxes, color= 'r', fontsize= 14, bbox= props)


    return


#########

def clstr_2d(line,k):
    
    """
        plot KDE of the clusters in the Vmin vs Vmax plane with the EW as a marker.
        
        Param:
        line: 'CIV', 'SiIV', or 'MgII'
        k: number of clusters (3,4,5..)
        
        """
    
    data= Table.read('myBALCat_xtra.csv')
    
    clstr_tbl= Table.read("./clusters/3features/"+line+str(k)+"clstrs.fits")
    
    clstrs_ls=[]
    for o in range(k):
        clstrs_ls.append([o ,len(clstr_tbl[clstr_tbl['label'] ==o]),\
                          mean(clstr_tbl['Vmin'][clstr_tbl['label'] ==o]),\
                          mean(clstr_tbl['Vmax'][clstr_tbl['label'] ==o]), \
                          mean(clstr_tbl['EW'][clstr_tbl['label'] ==o])])
    
    ord_clstrs= sorted(clstrs_ls, key= itemgetter(2))
    
    print ord_clstrs

    clr_ls= ['steelblue', 'olivedrab','orange', '0.8']
    clrm_ls= ['Blues_r', 'Greens_r', 'copper_r', 'gray_r']
    
    
    fig= figure(figsize=(10,10))
    
    ax1= fig.add_subplot(111)
    
    i =0
    for c in ord_clstrs:
        l= c[0]
        print l
        sns.kdeplot(clstr_tbl['Vmin'][clstr_tbl['label'] == l], clstr_tbl['Vmax'][clstr_tbl['label'] == l], \
                    cmap= clrm_ls[i], n_level= 10)
        text(0.1, +0.1, str(len(clstr_tbl[clstr_tbl['label'] == l])), color= clr_ls[i], fontsize= 18, transform=ax1.transAxes)
                    
        i+=1
    
    prop_tbl= join(data, clstr_tbl, keys='SDSSName')
    
    #scatter(prop_tbl['Vmin_'+line][(prop_tbl['BI1']-prop_tbl['BI2']) !=-999], prop_tbl['Vmax_'+line][(prop_tbl['BI1']-prop_tbl['BI2']) !=-999], marker='o', s=5, color='k')

    return

