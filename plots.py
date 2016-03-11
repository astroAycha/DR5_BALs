
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table, join
from sklearn.cluster import KMeans
from mpl_toolkits.mplot3d import Axes3D
from operator import itemgetter
from glob import glob
import seaborn as sns
from scipy.stats import spearmanr


########################

def line_prof(line, k, g):
    
    """ plot composite spectra in 2 panels:
        panel 1: He II & OIII]
        panel 2: Al III, Si III], C III]
        
        param:
        line: "SiIV", "CIV", "AlIII", or "MgII"
        k: number of clusters (3, 4, 5 or 6)
        g: directory where composites (or clusters) are saved, e.g., g1/ or g2/
        """
    
    sns.set_style('ticks', {'font.family': u'serif', 'xtick.direction': u'in', 'ytick.direction': u'in'})
    
    tt= Table.read("./clusters/"+g+line+str(k)+"clstrs.fits")
    
    cutoff= 10
    
    props= dict(boxstyle='round', facecolor='w', edgecolor='k')# , alpha=0.7)
    
    clstrs_ls=[]
    for o in range(k):
        clstrs_ls.append([o ,len(tt[tt['label'] ==o]),\
                          mean(tt['Vmin'][tt['label'] ==o]),\
                          mean(tt['Vmax'][tt['label'] ==o]), \
                          mean(tt['EW'][tt['label'] ==o])])

    ord_clstrs= sorted(clstrs_ls, key= itemgetter(2))

    print ord_clstrs
    
    
    clr_ls = [sns.xkcd_rgb["windows blue"], sns.xkcd_rgb["dusty purple"], sns.xkcd_rgb["pale red"], \
              sns.xkcd_rgb["greyish"], sns.xkcd_rgb["faded green"], sns.xkcd_rgb["amber"], sns.xkcd_rgb["pale aqua"]]
        
    amb= sns.light_palette(sns.xkcd_rgb["amber"], as_cmap= True)
    aq= sns.light_palette(sns.xkcd_rgb["pale aqua"], as_cmap= True)
              
    clrm_ls= ['Blues', 'Purples' , 'Reds', 'Greys', 'Greens', amb, aq]
              
              
    clstr_name= ['a', 'b', 'c', 'd', 'e', 'f', 'g']
              
              
    fig= figure(figsize=(12,6))

    # define where the axes go
    leftp1, widthp1 = 0.1, 0.3
    leftp2, widthp2 = 0.4, 0.5
    bottom, height = 0.1, 0.8

    p1 = [leftp1, bottom, widthp1, height]
    p2 = [leftp2, bottom, widthp2, height]

    # add the axes to the figure
    ax1 = plt.axes(p1)
    ax2 = plt.axes(p2)


    # left panel for the HeII and OIII]
    #ax1= fig.add_subplot(121)
    ax1.set_xlim(1600,1700)
    ax1.set_ylim(0.8,1.35)
    xlabel(line + r'Rest Wavelength')
    ylabel(line + r' Arbitrary flux')
              
    i=1
    for c in ord_clstrs:
        l= c[0]
        compo_name= "./composites/"+g+line+"_"+str(k)+"clstr"+str(l+1)+".fits"
        spec= fits.open(compo_name)
        if c[1] > cutoff:
            ax1.plot(spec[0].data[0], spec[0].data[1]/spec[0].data[1][(1700-1100)*2], lw= 2, color= clr_ls[i-1])
            #ax1.text(0.82, .99-i/15., line+"-"+clstr_name[i-1]+", N= "+str(len(tt[tt['label'] == l])), color= clr_ls[i-1], fontsize= 18, transform=ax1.transAxes) # labels for each cluster with number of obejects. Text colors match the plot
        i+=1

    # right panel for the CIII] lend and FeIII
    #ax2= fig.add_subplot(122)
    ax2.set_xlim(1800,2150)
    ax2.set_ylim(1,2)
    xlabel(r'Restframe Wavelength ($\AA$)')
    ylabel(r'Normalized Flux (erg s$^{-1}$ cm$^{-1}$ $\AA^{-1}$)')
    
    line_mark= [1640, 1663.5, 1857, 1892, 1908]
    line_label= ['HeII', 'OIII]', 'AlIII', 'SiIII]', 'CIII]']
    
    # labels for FeIII
    ax2.arrow(2030, 1.4, -30, -.1, fc='k', ec='k')
    ax2.arrow(2030, 1.4, +30, -.1, fc='k', ec='k')
    text(2020, 1.45, r'FeIII', fontsize= 14, family='serif', color='k')
    
    #for p in r:
       # axvline(line_mark[p], ls= ':', color= '.5')
       # text(line_mark[p], 3, line_label[p], rotation= 'vertical')
    
    
    i=1
    for c in ord_clstrs:
        l= c[0]
        compo_name= "./composites/"+g+line+"_"+str(k)+"clstr"+str(l+1)+".fits"
        spec= fits.open(compo_name)
        if c[1] > cutoff:
            ax2.plot(spec[0].data[0], spec[0].data[1]/spec[0].data[1][(2150-1100)*2], lw= 2, color= clr_ls[i-1])
        #ax2.text(0.82, .99-i/15., line+"-"+clstr_name[i-1]+", N= "+str(len(tt[tt['label'] == l])), color= clr_ls[i-1], fontsize= 18, transform=ax2.transAxes) # labels for each cluster with number of obejects. Text colors match the plot
        i+=1

    
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

def clust_compos(line, k, g):

    '''plot composites in one panel, and other two panels to show clustering parameter space
    k: number of clusters
    f: number of features
    '''
    
    sns.set_style('ticks', {'font.family': u'serif', 'xtick.direction': u'in', 'ytick.direction': u'in'})
    
    if line== 'MgII':
        xlimit= (1700, 3100)
        r= arange(5,9)
        lum= "logF2500"
    
    else:
        xlimit= (1300, 2165)
        r= arange(0,8)
        lum = "logF1400"

    cutoff = 5 # change to 30 for CIV and SiIV, 20 for AlIII and 10 for MgII
    
    #clstr_tbl= Table.read("./clusters/4features/ew_vmin_vmax_deltv/"+line+str(k)+"clstrs.fits")

    clstr_tbl= Table.read("./clusters/"+g+"/"+line+str(k)+"clstrs.fits")

    #data= Table.read("myBALCat_xtra.csv", format= 'ascii.csv')
                          
    data = Table.read("myBALs.fits")
                          
    tt= join(clstr_tbl, data, keys= 'SDSSName')

    props= dict(boxstyle='round', facecolor='w', edgecolor='k')# , alpha=0.7)

    clstrs_ls=[]
    for o in range(k):
        clstrs_ls.append([o ,len(tt[tt['label'] ==o]),\
                          mean(tt[line+'-vmin'][tt['label'] ==o]),\
                          mean(tt[line+'-vmax'][tt['label'] ==o]), \
                          mean(tt[line+'-EW'][tt['label'] ==o]), \
                          mean(tt[lum][tt['label'] ==o])])

    ord_clstrs= sorted(clstrs_ls, key= itemgetter(2))

    print ord_clstrs


    clr_ls = [sns.xkcd_rgb["windows blue"], sns.xkcd_rgb["dusty purple"], sns.xkcd_rgb["pale red"], \
              sns.xkcd_rgb["greyish"], sns.xkcd_rgb["faded green"], sns.xkcd_rgb["amber"], sns.xkcd_rgb["pale aqua"]]
    
    amb= sns.light_palette(sns.xkcd_rgb["amber"], as_cmap= True)
    aq= sns.light_palette(sns.xkcd_rgb["pale aqua"], as_cmap= True)
    
    clrm_ls= ['Blues', 'Purples' , 'Reds', 'Greys', 'Greens', amb, aq]

    
    clstr_name= ['a', 'b', 'c', 'd', 'e', 'f', 'g']
    

    fig= figure(figsize=(14,10))

    ax1= fig.add_subplot(221)
    xlabel(line + r' V$_{min}$ (km/s)')
    ylabel(line + r' V$_{max}$ (km/s)')
    
    i =1
    for c in ord_clstrs:
        l= c[0]
        print l
        if c[1] > cutoff:
            sns.kdeplot(tt[line+'-vmin'][tt['label'] == l], tt[line+'-vmax'][tt['label'] == l], \
                    shade= True, shade_lowest= False, alpha= 0.5, cmap= clrm_ls[i-1], label= False, legend= False)

        
        i+=1

    #labels for full sample with number of clusters and number of objects
    ax1.text(0.1, 0.8, line+" Sample, K= "+str(k)+"\n"+ "N= "+str(len(clstr_tbl)), color= 'k', fontsize= 18, transform=ax1.transAxes)

    vmin, vmax, ews, lum =[], [], [], []
    for c in ord_clstrs:
        if c[1] > cutoff:
            vmin.append(c[2])
            vmax.append(c[3])
            ews.append(c[4])
            lum.append(c[5])

    #ax1.scatter(vmin, vmax, s= [abs(e)*10 for e in ews], edgecolor= '#34495e', facecolor= 'w', marker= 'D')

    for x in range(len(vmin)):
        ax1.text(vmin[x], vmax[x], clstr_name[x] , color= 'k', fontsize= 14 , multialignment= 'center', bbox= props)


    #panel for the Vmax vs EW space

    ax2= fig.add_subplot(222, sharey= ax1)
    xlabel(line + r' EW ($\AA$)')
    #ylabel(line + r' $V_{max}$ (km/s)')
    #ax2.axes.get_yaxis().set_visible(False)
    #ax2.set_xticklabels([])
    ax2.yaxis.tick_right()
    
    i =1
    for c in ord_clstrs:
        l= c[0]
        print l
        if c[1] > cutoff:
            sns.kdeplot(tt[line+'-EW'][tt['label'] == l], tt[line+'-vmax'][tt['label'] == l], \
                        shade= True, shade_lowest= False, alpha= 0.5, cmap= clrm_ls[i-1], label= False, legend= False)
        
        i+=1

    for x in range(len(vmin)):
        text(ews[x], vmax[x], clstr_name[x] , color= 'k', fontsize= 14, multialignment= 'center', bbox= props)

    subplots_adjust(wspace =0.01)

    #panel for the composites
    ax3= fig.add_subplot(212)
    xlim(xlimit)
    ylim(.1,3.2)
    xlabel(r'Restframe Wavelength ($\AA$)')
    ylabel(r'Normalized Flux (erg s$^{-1}$ cm$^{-1}$ $\AA^{-1}$)')

    line_mark= [1335, 1396.8, 1549, 1640, 1663.5, 1857, 1892, 1908, 2800]
    line_label= ['CII', 'SiIV', 'CIV', 'HeII', 'OIII]', 'AlIII', 'SiIII]', 'CIII]', 'MgII']

    # labels for FeIII
    ax3.arrow(2030, 1.4, -30, -.1, fc='k', ec='k')
    ax3.arrow(2030, 1.4, +30, -.1, fc='k', ec='k')
    text(2020, 1.45, r'FeIII', fontsize= 14, family='serif', color='k')
    
    for p in r:
        axvline(line_mark[p], ls= ':', color= '.5')
        text(line_mark[p], 3, line_label[p], rotation= 'vertical')
    

    i=1
    for c in ord_clstrs:
        l= c[0]
        compo_name= "./composites/"+g+"/"+line+"_"+str(k)+"clstr"+str(l+1)+".fits"
        spec= fits.open(compo_name)
        if c[1] > cutoff:
            plot(spec[0].data[0], spec[0].data[1]/spec[0].data[1][(2150-1100)*2], lw= 2, color= clr_ls[i-1])
        ax3.text(0.82, .99-i/15., line+"-"+clstr_name[i-1]+", N= "+str(len(clstr_tbl[clstr_tbl['label'] == l])), color= clr_ls[i-1], fontsize= 18, transform=ax3.transAxes) # labels for each cluster with number of obejects. Text colors match the plot
        i+=1
    
    
    return

##considering adding pie chart here to show the fraction of objects with vatiablity
#pie([28,3,16,3,12], shadow= True, colors= clr_ls, startangle=30, radius= 0.6, labels=['a', 'b', 'c', 'd', 'e'])



    #prop_tbl= join(data, clstr_tbl, keys='SDSSName')
    #scatter(prop_tbl['Vmin_'+line][(prop_tbl['BI1']-prop_tbl['BI2']) !=-999], prop_tbl['Vmax_'+line][(prop_tbl['BI1']-prop_tbl['BI2']) !=-999], marker='o', s=5, color='k')



##############


def clstr_prop(line,k):

    """ 
    histograms of the clusters in the Vmin vs Vmax plane with the EW as a marker.
    
    Param: 
    line: 'CIV', 'SiIV', 'AlIII', or 'MgII'
    k: number of clusters (3,4,5..)
    
    """

    data= Table.read('myBALCat_xtra.csv', format= 'ascii.csv')

    clstr_tbl= Table.read("./clusters/3features/"+line+str(k)+"clstrs.fits") # for 3 features
    
    clstrs_ls=[]
    for o in range(k):
        clstrs_ls.append([o ,len(clstr_tbl[clstr_tbl['label'] ==o]),\
                           mean(clstr_tbl['Vmin'][clstr_tbl['label'] ==o]),\
                           mean(clstr_tbl['Vmax'][clstr_tbl['label'] ==o]), \
                           mean(clstr_tbl['EW'][clstr_tbl['label'] ==o])])
    
    
    ord_clstrs= sorted(clstrs_ls, key= itemgetter(2)) #order the clusters according to their vmin (same as in the clusters/compo plots)
    
    print ord_clstrs

    clr_ls = [sns.xkcd_rgb["windows blue"], sns.xkcd_rgb["dusty purple"], sns.xkcd_rgb["pale red"], \
          sns.xkcd_rgb["greyish"], sns.xkcd_rgb["faded green"], sns.xkcd_rgb["amber"], sns.xkcd_rgb["pale aqua"]]

    clstr_name= ['a', 'b', 'c', 'd', 'e', 'f', 'g']

    cutoff = 5 # change to 30 for CIV and SiIV, 20 for AlIII and 10 for MgII


    fig= figure(figsize=(16,12))

    fig1= fig.add_axes([0., 0., 1, 1])
    fig1.set_axis_off()
    fig1.set_xlim(0, 1)
    fig1.set_ylim(0, 1)
    fig1.text(.06, 0.54, r"Normalized Dist", rotation='vertical', \
              horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')


    prop_tbl= join(data, clstr_tbl, keys='SDSSName')

    props= dict(boxstyle='round', facecolor='w', edgecolor='k')# , alpha=0.7)
    
    
    ax1= fig.add_subplot(421)
    i= 0
    j =0
    param= "LOGEDD_RATIO_DR7"
    hist_bins= arange(min(prop_tbl[param][prop_tbl[param] !=-999]), max(prop_tbl[param][prop_tbl[param] !=-999]), \
                      (max(prop_tbl[param][prop_tbl[param] !=-999])-min(prop_tbl[param][prop_tbl[param] !=-999]))/12)
    #print hist_bins
    
    print param, len(prop_tbl[prop_tbl[param] !=-999])
    
    for c in ord_clstrs:
        if c[1] > cutoff:
            l= c[0]
            hist(prop_tbl[param][(prop_tbl['label'] == l) & (prop_tbl[param] !=-999) & (prop_tbl[param] < 50000)], \
            bins= hist_bins, histtype= 'step', normed= True, color= clr_ls[i], lw= 2)
        
            ax1.text(0.05, .85-j/10., line+"-"+clstr_name[i], color= clr_ls[i], fontsize= 16, transform=ax1.transAxes)
            j+=1
        
        i+=1
    
    ax1.text(0.7, 0.77,r"log(L/L$_{\rm Edd}$)", transform=ax1.transAxes, color= 'k', fontsize= 16)
    ax1.text(0.95, 0.85, "A", transform=ax1.transAxes, color= 'r', fontsize= 14, bbox= props)


    ax2= fig.add_subplot(422)
    i =0
    param= "logF1400" # log monochromatic lum at 1400A from the Gibson catalog
    hist_bins= arange(min(prop_tbl[param][prop_tbl[param] !=-999]), max(prop_tbl[param][prop_tbl[param] !=-999]), \
                      (max(prop_tbl[param][prop_tbl[param] !=-999])-min(prop_tbl[param][prop_tbl[param] !=-999]))/12)

    print param, len(prop_tbl[prop_tbl[param] !=-999])
    
    for c in ord_clstrs:
        if c[1] >cutoff:
            l= c[0]
            hist(prop_tbl[param][(prop_tbl['label'] == l) & (prop_tbl[param] !=-999)], \
             bins= hist_bins, histtype= 'step', normed= True, color= clr_ls[i], lw= 2)
             
        i+=1
    
    ax2.text(0.05, 0.6, r"log L(1400$\AA$)"+"\n"+"[mW/m$^2$/Hz]", transform=ax2.transAxes, color= 'k', fontsize= 16)
    ax2.text(0.95, 0.85, "B", transform=ax2.transAxes, color= 'r', fontsize= 14, bbox= props)


    ax3= fig.add_subplot(423)
    i =0
    param= "E_B-V_1"
    hist_bins= arange(min(prop_tbl[param][prop_tbl[param] !=-999]), max(prop_tbl[param][prop_tbl[param] !=-999]), \
                      (max(prop_tbl[param][prop_tbl[param] !=-999])-min(prop_tbl[param][prop_tbl[param] !=-999]))/12)

    print param, len(prop_tbl[prop_tbl[param] !=-999])

    for c in ord_clstrs:
        if c[1] >cutoff:
            l= c[0]
            hist(prop_tbl[param][(prop_tbl['label'] == l) & (prop_tbl[param] !=-999)], \
             bins= hist_bins, histtype= 'step', normed= True, color= clr_ls[i], lw= 2)
             
        i+=1

    ax3.text(0.65, 0.8,r"E(B-V)", transform=ax3.transAxes, color= 'k', fontsize= 18)
    ax3.text(0.95, 0.85, "C", transform=ax3.transAxes, color= 'r', fontsize= 14, bbox= props)

    ax4= fig.add_subplot(424)
    i =0
    param= "int_alpha_nu"
    hist_bins= arange(min(prop_tbl[param][prop_tbl[param] !=-999]), max(prop_tbl[param][prop_tbl[param] !=-999]), \
                      (max(prop_tbl[param][prop_tbl[param] !=-999])-min(prop_tbl[param][prop_tbl[param] !=-999]))/12)
        
    print param, len(prop_tbl[prop_tbl[param] !=-999])
                      
    for c in ord_clstrs:
        if c[1] > cutoff:
            l= c[0]
            hist((prop_tbl[param][(prop_tbl['label'] == l) & (prop_tbl[param] !=-999) & (prop_tbl[param] !=0)]), \
                                       bins= hist_bins, histtype= 'step', normed= True, color= clr_ls[i], lw= 2)
                                      
        i+=1
    ax4.text(0.1, 0.8,r"Intrinsic $\alpha_\nu$", transform=ax4.transAxes, color= 'k', fontsize= 16)
    ax4.text(0.95, 0.85, "D", transform=ax4.transAxes, color= 'r', fontsize= 14, bbox= props)
    

    ax5= fig.add_subplot(425)
    i =0
    param= "HeII_EW_BLH"
    hist_bins= arange(min(prop_tbl[param][prop_tbl[param] !=-999]), max(prop_tbl[param][prop_tbl[param] !=-999]), \
                      (max(prop_tbl[param][prop_tbl[param] !=-999])-min(prop_tbl[param][prop_tbl[param] !=-999]))/12)

    print param, len(prop_tbl[prop_tbl[param] !=-999])

    for c in ord_clstrs:
        if c[1] > cutoff:
            l= c[0]
            hist(prop_tbl[param][(prop_tbl['label'] == l) & (prop_tbl[param] !=-999) & (prop_tbl[param] !=0)], \
                 bins= hist_bins, histtype= 'step', normed= True, color= clr_ls[i], lw= 2)

        i+=1
    ax5.text(0.55, 0.82,r"EW(HeII) [$\AA$]", transform=ax5.transAxes, color= 'k', fontsize= 16)
    ax5.text(0.95, 0.85, "E", transform=ax5.transAxes, color= 'r', fontsize= 14, bbox= props)


    ax6= fig.add_subplot(426)
    i =0
    param= "alpha_UV_BLH"
    hist_bins= arange(min(prop_tbl[param][prop_tbl[param] !=-999]), max(prop_tbl[param][prop_tbl[param] !=-999]), \
                      (max(prop_tbl[param][prop_tbl[param] !=-999])-min(prop_tbl[param][prop_tbl[param] !=-999]))/12)

    print param, len(prop_tbl[prop_tbl[param] !=-999])

    for c in ord_clstrs:
        if c[1] >cutoff:
            l= c[0]
            hist(prop_tbl[param][(prop_tbl['label'] == l) & (prop_tbl[param] !=-999)], \
             bins= hist_bins, histtype= 'step', normed= True, color= clr_ls[i], lw= 2)
        
        i+=1

    ax6.text(0.65, 0.8,r"$\alpha_{UV}$", transform=ax6.transAxes, color= 'k', fontsize= 16)
    ax6.text(0.95, 0.85, "F", transform=ax6.transAxes, color= 'r', fontsize= 14, bbox= props)


    ax7= fig.add_subplot(427)
    i =0
    param= "CF_BLH"
    hist_bins= arange(min(prop_tbl[param][prop_tbl[param] !=-999]), max(prop_tbl[param][prop_tbl[param] !=-999]), \
                      (max(prop_tbl[param][prop_tbl[param] !=-999])-min(prop_tbl[param][prop_tbl[param] !=-999]))/12)

    print param, len(prop_tbl[prop_tbl[param] !=-999])

    for c in ord_clstrs:
        if c[1] > cutoff:
            l= c[0]
            hist(prop_tbl[param][(prop_tbl['label'] == l) & (prop_tbl[param] !=-999)], \
             bins= hist_bins, histtype= 'step', normed= True, color= clr_ls[i], lw= 2)
        
        i+=1

    ax7.text(0.2, 0.85,"CF", transform=ax7.transAxes, color= 'k', fontsize= 16)
    ax7.text(0.05, 0.85, "G", transform=ax7.transAxes, color= 'r', fontsize= 14, bbox= props)


    ax8= fig.add_subplot(428)
    i =0
    param= "v_md_BLH"
    hist_bins= arange(min(prop_tbl[param][prop_tbl[param] !=-999]), max(prop_tbl[param][prop_tbl[param] !=-999]), \
                      (max(prop_tbl[param][prop_tbl[param] !=-999])-min(prop_tbl[param][prop_tbl[param] !=-999]))/12)
        
    print param, len(prop_tbl[prop_tbl[param] !=-999])
                      
    for c in ord_clstrs:
        if c[1] >cutoff:
            l= c[0]
            hist(prop_tbl[param][(prop_tbl['label'] == l) & (prop_tbl[param] !=-999)], \
                bins= hist_bins, histtype= 'step', normed= True, color= clr_ls[i], lw= 2)
                                      
        i+=1
            
    ax8.text(0.65, 0.8,r" v$_{md}$ [km/s]", transform=ax8.transAxes, color= 'k', fontsize= 16)
    ax8.text(0.95, 0.85, "H", transform=ax8.transAxes, color= 'r', fontsize= 14, bbox= props)

    return

#########

# compare alpha UV from Baskin et al 2014 and alpha lambda from Krawczyk et al 2014

#sns.set(font_scale= 1.5)
sns.set_style('ticks', {'font.family': u'serif', 'xtick.direction': u'in', 'ytick.direction': u'in'})


data= Table.read('myBALCat_xtra.csv', format= 'ascii.csv')

clstr= Table.read('./clusters/3features/CIV6clstrs.fits')

t= join(data, clstr, keys= 'SDSSName')

clstrs_ls=[]

k=6

for o in range(k):
    clstrs_ls.append([o ,len(clstr[clstr['label'] ==o]), \
                      mean(clstr['Vmin'][clstr['label'] ==o]),\
                      mean(clstr['Vmax'][clstr['label'] ==o]), \
                      mean(clstr['EW'][clstr['label'] ==o])])


ord_clstrs= sorted(clstrs_ls, key= itemgetter(2))

print ord_clstrs


clr_ls = [sns.xkcd_rgb["windows blue"], sns.xkcd_rgb["dusty purple"], sns.xkcd_rgb["pale red"], \
          sns.xkcd_rgb["greyish"], sns.xkcd_rgb["faded green"], sns.xkcd_rgb["amber"]]
    
clstr_name= ['a', 'b', 'c', 'd', 'e', 'f']

print len(t[(t['alpha_UV_BLH'] !=-999) & (t['Dal1']!=-999)])

fig= figure(figsize=(10,8))
xlabel(r"$\alpha_{UV}$", fontsize= 20)
ylabel(r"Intrinsic $\alpha_\nu$", fontsize= 20)
xlim(-3, 0.7)
ylim(-0.7,0.25)

text(-2.75, -0.6, "Red", color= 'r', fontsize= 20)
text(0.3, 0.15, "Blue", color= 'b', fontsize= 20)

# draw 1:1 line
x = y = arange(-3.,3.,.001)
plot(x,y, ls=':', c='k')

j= 0
for r in ord_clstrs:
    c= r[0]
    scatter(t['alpha_UV_BLH'][(t['alpha_UV_BLH'] !=-999) & (t['Dal1']!=-999) & (t['label']==c)],\
            -.28+t['Dal1'][(t['alpha_UV_BLH'] !=-999) & (t['Dal1']!=-999) & (t['label']== c)], c= clr_ls[j], s=100, alpha= 0.9)
    j+=1


## quick Mi vs z plot -for the CIV sample
#Tables are read as in the plot above

fig= figure(figsize=(10,9))
xlabel("Redshift", fontsize=20)
ylabel(r"$M_{\rm i}$", fontsize=20)
xlim(1.7,3.75)
ylim(-24.7, -30)

j= 0
for r in ord_clstrs:
    c= r[0]
    scatter(t['z_1'][t['label']==c], t['M_i'][t['label']== c], c= clr_ls[j], marker= 'o', s=10)
    j+=1

## same but for full sample and CIV sample overplotted

fig= figure(figsize=(10,9))
xlabel("Redshift", fontsize=20)
ylabel(r"$M_{\rm i}$", fontsize=20)

xlim(0,5.2)
ylim(-22.5, -29.8)
scatter(data['z'], data['M_i'], c= '0.8', marker= '.', label='Full BALQs sample')
sns.kdeplot(t['z_1'], t['M_i'], cmap='YlOrBr_r', shade_lowest= False, legend= False, n_levels= 10)
#scatter(t['z_1'], t['M_i'], c='r', marker= 'o', s=7, label='CIV sample')
#legend(loc= 4)


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

###############

## plot multiple panels with 2d scatter plots

sns.set_style('ticks', {'font.family': u'serif', 'xtick.direction': u'in', 'ytick.direction': u'in'})

data= Table.read('myBALCat_xtra.csv', format= 'ascii.csv')

clstr= Table.read('./clusters/3features/CIV6clstrs.fits')

t= join(data, clstr, keys= 'SDSSName')

fig= figure(figsize=(16,12))

##subplots
nl= 15 #number of contours
cm= "OrRd" #color map
#cm= sns.cubehelix_palette(light=1, as_cmap=True)

ax1= fig.add_subplot(4,3,1)
px = "Vmin_CIV"
py = "LOGEDD_RATIO_DR7"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
#text(0.1, 0.85, str(len(t[t[py] != -999]))+" objects", color='k', fontsize=18, transform=ax1.transAxes)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.1, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
     color='k', fontsize=14, transform=ax1.transAxes)

ax2= fig.add_subplot(4,3,2)
px = "Vmax_CIV"
py = "LOGEDD_RATIO_DR7"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.1, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
     color='k', fontsize=14, transform=ax2.transAxes)

ax3= fig.add_subplot(4,3,3)
px = "EW_CIV"
py = "LOGEDD_RATIO_DR7"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
xlim(-60,10)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.1, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
     color='k', fontsize=14, transform=ax3.transAxes)

ax4= fig.add_subplot(4,3,4)
px = "Vmin_CIV"
py = "logF1400"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
#text(0.1, 0.85, str(len(t[t[py] != -999]))+" objects", color='k', fontsize=18, transform=ax4.transAxes)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.1, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
    color='k', fontsize=14, transform=ax4.transAxes)
     
ax5= fig.add_subplot(4,3,5)
px = "Vmax_CIV"
py = "logF1400"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.1, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
    color='k', fontsize=14, transform=ax5.transAxes)
     
ax6= fig.add_subplot(4,3,6)
px = "EW_CIV"
py = "logF1400"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
xlim(-60,10)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.1, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
    color='k', fontsize=14, transform=ax6.transAxes)
     
ax7= fig.add_subplot(4,3,7)
px = "Vmin_CIV"
py = "int_alpha_nu"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
#text(0.1, 0.85, str(len(t[t[py] != -999]))+" objects", color='k', fontsize=18, transform=ax7.transAxes)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.1, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
    color='k', fontsize=14, transform=ax7.transAxes)
     
ax8= fig.add_subplot(4,3,8)
px = "Vmax_CIV"
py = "int_alpha_nu"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.1, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
    color='k', fontsize=14, transform=ax8.transAxes)
     
ax9= fig.add_subplot(4,3,9)
px = "EW_CIV"
py = "int_alpha_nu"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
xlim(-60,10)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.1, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
    color='k', fontsize=14, transform=ax9.transAxes)
     
ax10= fig.add_subplot(4,3,10)
px = "Vmin_CIV"
py = "HeII_EW_BLH"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
#text(0.1, 0.85, str(len(t[t[py] != -999]))+" objects", color='k', fontsize=18, transform=ax10.transAxes)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.1, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
    color='k', fontsize=14, transform=ax10.transAxes)
     
ax11= fig.add_subplot(4,3,11)
px = "Vmax_CIV"
py = "HeII_EW_BLH"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.1, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
    color='k', fontsize=14, transform=ax11.transAxes)
     
ax12= fig.add_subplot(4,3,12)
px = "EW_CIV"
py = "HeII_EW_BLH"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
xlim(-60,10)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.1, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
    color='k', fontsize=14, transform=ax12.transAxes)
     
#axes labels
fig1= fig.add_axes([0., 0., 1, 1])
fig1.set_axis_off()
fig1.set_xlim(0, 1)
fig1.set_ylim(0, 1)

fig1.text(.23, 0.93, r" CIV Vmin (km/s)", rotation='horizontal', horizontalalignment='center',verticalalignment='center', fontsize= 18, family= 'serif')

fig1.text(0.5, 0.93, r"CIV Vmax (km/s)", rotation='horizontal', horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')

fig1.text(0.78, 0.925, r"CIV EW ($\AA$)", rotation='horizontal', horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')

fig1.text(0.082, 0.81, r"log L/L$_{\rm Edd}$", rotation='vertical', horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')

fig1.text(0.082, 0.6, r"log L(1400)", rotation='vertical', horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')

fig1.text(0.082, 0.4, r"Intrin $\alpha_\nu$", rotation='vertical', horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')

fig1.text(0.082, 0.19, r"EW(HeII) ($\AA$)", rotation='vertical', horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')
###############################################

## keep only a few of the previous subplots

## plot multiple panels with 2d scatter plots

sns.set_style('ticks', {'font.family': u'serif', 'xtick.direction': u'in', 'ytick.direction': u'in'})

data= Table.read('myBALCat_xtra.csv', format= 'ascii.csv')

clstr= Table.read('./clusters/3features/CIV6clstrs.fits')

t= join(data, clstr, keys= 'SDSSName')

fig= figure(figsize=(15,9))

#subplots
nl= 15 #number of contours
cm= "PuBuGn_r" #color map

ax1= fig.add_subplot(231)
xlim(5500, -24100)
px = "Vmin_CIV"
py = "int_alpha_nu"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.1, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
     color='k', fontsize=14, transform=ax1.transAxes)

ax2= fig.add_subplot(232)
xlim(2500, -31200)
px = "Vmax_CIV"
py = "int_alpha_nu"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.1, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
     color='k', fontsize=14, transform=ax2.transAxes)

ax3= fig.add_subplot(233)
xlim(5,-52)
px = "EW_CIV"
py = "int_alpha_nu"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.6, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
     color='k', fontsize=14, transform=ax3.transAxes)

ax4= fig.add_subplot(234)
xlim(5500, -24100)
ylim(-5,15)
px = "Vmin_CIV"
py = "HeII_EW_BLH"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.6, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
     color='k', fontsize=14, transform=ax4.transAxes)

ax5= fig.add_subplot(235)
xlim(2500, -31200)
ylim(-5,15)
px = "Vmax_CIV"
py = "HeII_EW_BLH"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.6, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
     color='k', fontsize=14, transform=ax5.transAxes)

ax6= fig.add_subplot(236)
xlim(5,-52)
ylim(-5,15)
px = "EW_CIV"
py = "HeII_EW_BLH"
sns.kdeplot(t[px][t[py] != -999], t[py][t[py] != -999], cmap=cm, n_levels= nl, shade_lowest= False, legend= False)
s= spearmanr(t[px][t[py] != -999], t[py][t[py] != -999])
text(0.6, 0.75 ,"r= "+"{:03.2f}".format(s[0])+"\n"+"p= "+"{:3.2f}".format(s[1]), \
     color='k', fontsize=14, transform=ax6.transAxes)

#axes labels
fig1= fig.add_axes([0., 0., 1, 1])
fig1.set_axis_off()
fig1.set_xlim(0, 1)
fig1.set_ylim(0, 1)

fig1.text(.23, 0.05, r" CIV Vmin (km/s)", rotation='horizontal', horizontalalignment='center',verticalalignment='center', fontsize= 18, family= 'serif')

fig1.text(0.5, 0.05, r"CIV Vmax (km/s)", rotation='horizontal', horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')

fig1.text(0.78, 0.04, r"CIV EW ($\AA$)", rotation='horizontal', horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')

fig1.text(0.09, 0.71, r"Intrin $\alpha_\nu$", rotation='vertical', horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')

fig1.text(0.09, 0.3, r"EW(HeII) ($\AA$)", rotation='vertical', horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')



#####
### this is similar to the previous plot (kde of params vs vmin, vmax and EW) but using mean of param in each cluster

sns.set_style('ticks', {'font.family': u'serif', 'xtick.direction': u'in', 'ytick.direction': u'in'})

data= Table.read('myBALCat_xtra.csv', format= 'ascii.csv')

clstr= Table.read('./clusters/3features/CIV6clstrs.fits')

t= join(data, clstr, keys= 'SDSSName')

clstrs_ls=[]

k=6

for o in range(k):
    clstrs_ls.append([o ,len(clstr[clstr['label'] ==o]), \
                      mean(clstr['Vmin'][clstr['label'] ==o]),\
                      mean(clstr['Vmax'][clstr['label'] ==o]), \
                      mean(clstr['EW'][clstr['label'] ==o])])


oc= sorted(clstrs_ls, key= itemgetter(2)) #ordered clusters

clr_ls = [sns.xkcd_rgb["windows blue"], sns.xkcd_rgb["dusty purple"], sns.xkcd_rgb["pale red"], \
          sns.xkcd_rgb["greyish"], sns.xkcd_rgb["faded green"], sns.xkcd_rgb["amber"], sns.xkcd_rgb["pale aqua"]]
    
clstr_name= ['a', 'b', 'c', 'd', 'e', 'f']
var_ls= [21*20, 21*20, 42*20, 40*20, 30*20, 10*20]

fig= figure(figsize=(16,12))


#subplots

ax1= fig.add_subplot(4,3,1, sharex= ax10)
ax1.set_xticklabels('',visible=False)
xlim(100,-18500)
px = "Vmin_CIV"
py = "LOGEDD_RATIO_DR7"
m1, m2=[], []
for i in oc:
    c= i[0]
    m1.append(mean(t[px][t['label'] == c]))
    m2.append(mean(t[py][(t[py]!=-999) & (t['label']== c)]))
scatter(m1,m2, color= clr_ls, marker= 'o', s= var_ls, alpha=0.8)
for j in range(k):
    text(m1[j], m2[j], clstr_name[j], color='k')

ax2= fig.add_subplot(4,3,2, sharex= ax11)
ax2.set_xticklabels('',visible=False)
ax2.set_yticklabels('',visible=False)
xlim(-3500, -23500)
px = "Vmax_CIV"
py = "LOGEDD_RATIO_DR7"
m1, m2=[], []
for i in oc:
    c= i[0]
    m1.append(mean(t[px][t['label'] == c]))
    m2.append(mean(t[py][(t[py]!=-999) & (t['label']== c)]))
scatter(m1,m2, color= clr_ls, marker= 'o', s= var_ls, alpha=0.8)
for j in range(k):
    text(m1[j], m2[j], clstr_name[j], color='k')

ax3= fig.add_subplot(4,3,3, sharex= ax12)
ax3.set_xticklabels('',visible=False)
ax3.set_yticklabels('',visible=False)
xlim(-1,-52)
px = "EW_CIV"
py = "LOGEDD_RATIO_DR7"
m1, m2=[], []
for i in oc:
    c= i[0]
    m1.append(mean(t[px][t['label'] == c]))
    m2.append(mean(t[py][(t[py]!=-999) & (t['label']== c)]))
scatter(m1,m2, color= clr_ls, marker= 'o', s= var_ls, alpha=0.8)
for j in range(k):
    text(m1[j], m2[j], clstr_name[j], color='k')

ax4= fig.add_subplot(4,3,4, sharex= ax10)
ax4.set_xticklabels('',visible=False)
xlim(100,-18500)
px = "Vmin_CIV"
py = "logF1400"
m1, m2=[], []
for i in oc:
    c= i[0]
    m1.append(mean(t[px][t['label'] == c]))
    m2.append(mean(t[py][(t[py]!=-999) & (t['label']== c)]))
scatter(m1,m2, color= clr_ls, marker= 'o', s= var_ls, alpha=0.8)
for j in range(k):
    text(m1[j], m2[j], clstr_name[j], color='k')


ax5= fig.add_subplot(4,3,5, sharex= ax11)
ax5.set_xticklabels('',visible=False)
ax5.set_yticklabels('',visible=False)
xlim(-3500, -23500)
px = "Vmax_CIV"
py = "logF1400"
m1, m2=[], []
for i in oc:
    c= i[0]
    m1.append(mean(t[px][t['label'] == c]))
    m2.append(mean(t[py][(t[py]!=-999) & (t['label']== c)]))
scatter(m1,m2, color= clr_ls, marker= 'o', s= var_ls, alpha=0.8)
for j in range(k):
    text(m1[j], m2[j], clstr_name[j], color='k')

ax6= fig.add_subplot(4,3,6, sharex= ax12)
ax6.set_xticklabels('',visible=False)
ax6.set_yticklabels('',visible=False)
xlim(-1,-52)
px = "EW_CIV"
py = "logF1400"
m1, m2=[], []
for i in oc:
    c= i[0]
    m1.append(mean(t[px][t['label'] == c]))
    m2.append(mean(t[py][(t[py]!=-999) & (t['label']== c)]))
scatter(m1,m2, color= clr_ls, marker= 'o', s= var_ls, alpha=0.8)
for j in range(k):
    text(m1[j], m2[j], clstr_name[j], color='k')

ax7= fig.add_subplot(4,3,7, sharex= ax10)
ax7.set_xticklabels('',visible=False)
xlim(100,-18500)
px = "Vmin_CIV"
py = "int_alpha_nu"
m1, m2=[], []
for i in oc:
    c= i[0]
    m1.append(mean(t[px][t['label'] == c]))
    m2.append(mean(t[py][(t[py]!=-999) & (t['label']== c)]))
scatter(m1,m2, color= clr_ls, marker= 'o', s= var_ls, alpha=0.8)
for j in range(k):
    text(m1[j], m2[j], clstr_name[j], color='k')

ax8= fig.add_subplot(4,3,8, sharex= ax11)
ax8.set_xticklabels('',visible=False)
ax8.set_yticklabels('',visible=False)
xlim(-3500, -23500)
px = "Vmax_CIV"
py = "int_alpha_nu"
m1, m2=[], []
for i in oc:
    c= i[0]
    m1.append(mean(t[px][t['label'] == c]))
    m2.append(mean(t[py][(t[py]!=-999) & (t['label']== c)]))
scatter(m1,m2, color= clr_ls, marker= 'o', s= var_ls, alpha=0.8)
for j in range(k):
    text(m1[j], m2[j], clstr_name[j], color='k')

ax9= fig.add_subplot(4,3,9, sharex= ax12)
ax9.set_xticklabels('',visible=False)
ax9.set_yticklabels('',visible=False)
xlim(-1,-52)
px = "EW_CIV"
py = "int_alpha_nu"
m1, m2=[], []
for i in oc:
    c= i[0]
    m1.append(mean(t[px][t['label'] == c]))
    m2.append(mean(t[py][(t[py]!=-999) & (t['label']== c)]))
scatter(m1,m2, color= clr_ls, marker= 'o', s= var_ls, alpha=0.8)
for j in range(k):
    text(m1[j], m2[j], clstr_name[j], color='k')

ax10= fig.add_subplot(4,3,10)
xlim(100,-18500)
px = "Vmin_CIV"
py = "HeII_EW_BLH"
m1, m2=[], []
for i in oc:
    c= i[0]
    m1.append(mean(t[px][t['label'] == c]))
    m2.append(mean(t[py][(t[py]!=-999) & (t['label']== c)]))
scatter(m1,m2, color= clr_ls, marker= 'o', s= var_ls, alpha=0.8)
for j in range(k):
    text(m1[j], m2[j], clstr_name[j], color='k')

ax11= fig.add_subplot(4,3,11)
ax11.set_yticklabels('',visible=False)
xlim(-3500, -23500)
px = "Vmax_CIV"
py = "HeII_EW_BLH"
m1, m2=[], []
for i in oc:
    c= i[0]
    m1.append(mean(t[px][t['label'] == c]))
    m2.append(mean(t[py][(t[py]!=-999) & (t['label']== c)]))
scatter(m1,m2, color= clr_ls, marker= 'o', s= var_ls, alpha=0.8)
for j in range(k):
    text(m1[j], m2[j], clstr_name[j], color='k')

ax12= fig.add_subplot(4,3,12)
ax12.set_yticklabels('',visible=False)
xlim(-1,-52)
px = "EW_CIV"
py = "HeII_EW_BLH"
m1, m2=[], []
for i in oc:
    c= i[0]
    m1.append(mean(t[px][t['label'] == c]))
    m2.append(mean(t[py][(t[py]!=-999) & (t['label']== c)]))
scatter(m1,m2, color= clr_ls, marker= 'o', s= var_ls, alpha=0.8)
for j in range(k):
    text(m1[j], m2[j], clstr_name[j], color='k')

subplots_adjust(wspace =0.04, hspace= 0.04)
#ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')

#axes labels
fig1= fig.add_axes([0., 0., 1, 1])
fig1.set_axis_off()
fig1.set_xlim(0, 1)
fig1.set_ylim(0, 1)

#x-axis bottom
fig1.text(.23, 0.05, r" CIV Vmin (km/s)", rotation='horizontal', horizontalalignment='center',verticalalignment='center', fontsize= 18, family= 'serif')

fig1.text(0.5, 0.05, r"CIV Vmax (km/s)", rotation='horizontal', horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')

fig1.text(0.78, 0.05, r"CIV EW ($\AA$)", rotation='horizontal', horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')

#y-axis
fig1.text(0.07, 0.81, r"log L/L$_{\rm Edd}$", rotation='vertical', horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')

fig1.text(0.07, 0.6, r"log L(1400)"+"\n"+" [mW/m$^2$/Hz]", rotation='vertical', horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')

fig1.text(0.07, 0.4, r"Intrin $\alpha_\nu$", rotation='vertical', horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')

fig1.text(0.07, 0.19, r"EW(HeII) ($\AA$)", rotation='vertical', horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')

###############

## quick 3D plots to show clusters

#from astropy.table import Table
#from mpl_toolkits.mplot3d import Axes3D
#import seaborn as sns

c4= Table.read('./clusters/3features/ew_vmin_vmax/CIV6clstrs.fits')

clstrs_ls=[]

k=6

for o in range(k):
    clstrs_ls.append([o ,len(c4[c4['label'] ==o]), \
                      mean(c4['Vmin'][c4['label'] ==o]),\
                      mean(c4['Vmax'][c4['label'] ==o]), \
                      mean(c4['EW'][c4['label'] ==o])])


oc= sorted(clstrs_ls, key= itemgetter(2)) #ordered clusters (with Vmin)

clr_ls = [sns.xkcd_rgb["windows blue"], sns.xkcd_rgb["dusty purple"], sns.xkcd_rgb["pale red"], \
          sns.xkcd_rgb["greyish"], sns.xkcd_rgb["faded green"], sns.xkcd_rgb["amber"], sns.xkcd_rgb["pale aqua"]]

fig= figure(figsize=(10,8))

ax= plt.axes(projection='3d')

#plot projections
#ax.scatter(c4['EW'], c4['Vmin'], marker= '+', color='k', zdir='z', zs=-2)
#ax.scatter(c4['Vmin'], c4['Vmax'], marker= '+', color='k', zdir='x', zs=0)
#ax.scatter(c4['EW'], c4['Vmax'], marker= '+', color='k', zdir='y', zs=-2.5)

j= 0
for c in oc:
    i= c[0]

    ax.scatter(xs= c4['EW'][c4['label']== i], ys= c4['Vmin'][c4['label']== i], zs= c4['Vmax'][c4['label']== i], c= clr_ls[j])

    j+= 1


ax.set_xlabel('CIV EW')
ax.set_ylabel('CIV Vmin')
ax.set_zlabel('CIV Vmax')


