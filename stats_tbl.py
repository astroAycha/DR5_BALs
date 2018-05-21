""" this will soon (!) include code to figure out the properties of each cluster with regards to the extra data added to the Gibson catalog (mainly reddening and variablity. and maybe BH masses...)
"""

import numpy as np
from operator import itemgetter
from astropy.table import Table, join


## used in paper
def overlap(g):
    
    """ table to show overlap between our sample from Gibson et al. 2015 and other datasets: Filiz Ak 2014, Baskin et al 2015, Krawczyk et al. 2015, and Shen et al. 2011
        param:
        g: group number: refers to the clustering runs in bal_cluster.py each with a different combination of features
        """
    
    #line_ls= ["CIV", "SiIV", "AlIII", "MgII"]
    line_ls= ["CIV"]
    
    tbl= open("cntrs_num_tbl"+str(g)+"x.txt", 'wrb') # new file to recored numbers in a latex table format
    
    data= Table.read('myBALsx.fits')
    
    alph= ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    
    for l in line_ls:
        
        if l == 'MgII':
            lum = 'logF2500'
        else:
            lum = 'logF1400'
    
        for k in range(3,9):
            
            tt= Table.read("./clusters/"+g+"/"+l+str(k)+"clstrs.fits")
            
            t= join(data, tt, keys= 'SDSSName')
            
            tbl.write(l+ "&" + str(len(t)) + "& "+ "K="+ str(k)+ "\n")
            
            clstrs_ls=[]
      
            for o in range(k):
                clstrs_ls.append([l, k, o ,len(t[t['label'] ==o]),\
                                  mean(t[l+'-EW'][t['label'] ==o]),\
                                  mean(t[l+'-vmax'][t['label'] ==o]), \
                                  mean(t[l+'-vmin'][t['label'] ==o]),\
                                  mean(t['EW_dV'][t['label']== o])])
            
            oc= sorted(clstrs_ls, key= itemgetter(4)) #ordered clusters
            
            nur=[] # nubmber of objects from Filiz Ak
            baskin= [] # number of objects from Baskin
            krawczyk= [] #number of objects from krawczyk
            shen= [] # number of objects from Shen
            var= [] # variability
            si4= [] # SiIV only. no AlIII
            sial=[] # SiIV with or without AlIII
            
            for x in oc:
                q= x[2] #cluster label
                
                nur.append(len(t[(t['label'] ==q) & (t['Delt'] != -999)]))
                
                baskin.append(len(t[(t['label'] ==q) & (t['HeII_EW'] != -999)]))
                
                krawczyk.append(len(t[(t['label'] ==q) & (t['E_B-V_1'] != -999)]))
                
                shen.append(len(t[(t['label'] ==q) & (t['LOGEDD_RATIO'] != -999)]))
                
                var.append(len(t[(t['label'] ==q) & (abs(t['BI1']-t['BI2']) >0)])*100./len(t[t['label']==q]))
                
                si4.append(len(t[(t['label'] ==q) & (t['SiIV-BIO'] >0) & (t['AlIII-BIO'] ==0)])*100./len(t[t['label']==q]))
                
                sial.append(len(t[(t['label'] ==q) & (t['SiIV-BIO'] >0) & (t['AlIII-BIO'] >0)])*100./len(t[t['label']==q]))
            
            
            for (c,j) in zip(oc, range(k)):
                a= alph[j]
                tbl.write("& {} & ${:d}$ & ${:06.2f}$ & ${:06.2f}$ & ${:06.2f}$ & ${:02.1f}$ & ${:02.1f}$ & ${:02.1f}$ & ${:02.1f}$ & ${}$ & ${}$ & ${}$ & ${}$ \n".format(l+"-"+a, c[3], c[4], c[5], c[6], c[7], si4[j], sial[j], var[j], nur[j], baskin[j], krawczyk[j], shen[j]))
                
                #tbl.write("{},{:d},{:06.2f},{:06.2f},{:06.2f},{:02.1f},{:02.1f},{:02.1f},{},{},{},{} \n".format(l+"-"+a, c[3], c[4], c[5], c[6], si4[j], sial[j], var[j], nur[j], baskin[j], krawczyk[j], shen[j]))

    tbl.close()
    
    return


###################

def clstr_cntrs(g):
    
    """ create a table with the centroids for the clusters for k= 3,4 and 5. formated for the manuscript
        param:
        g: group number: refers to the clustering runs in bal_cluster.py each with a different combination of features
        """
    
    #line_ls= ["CIV", "SiIV", "AlIII", "MgII"]
    line_ls= ["CIV"]
    
    tbl= open("cntrs_tbl"+str(g)+".txt", 'wrb')
    
    data= Table.read('myBALsx.fits')
    
    alph= ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    
    for l in line_ls:
        
        if l == 'MgII':
            lum = 'logF2500'
        else:
            lum = 'logF1400'
    
        for k in range(3,9):
            
            tt= Table.read("./clusters/"+g+"/"+l+str(k)+"clstrs.fits")
            
            t= join(data, tt, keys= 'SDSSName')
            
            tbl.write(l+ "&" + str(len(t)) + "& "+ "K="+ str(k)+ "\n")
            
            clstrs_ls=[]
            
            for o in range(k):
                clstrs_ls.append([l, k, o ,len(t[t['label'] ==o]),\
                                  mean(t[l+'-vmin'][t['label'] ==o]),\
                                  mean(t[l+'-vmax'][t['label'] ==o]), \
                                  mean(t[l+'-EW'][t['label'] ==o])])
        
            oc= sorted(clstrs_ls, key= itemgetter(4)) #ordered clusters
            
            
            var= [] #variability
            si4= [] # SiIV only. no AlIII
            sial=[] # SiIV with or without AlIII
            
            for x in oc:
                q= x[2] #cluster label
                
                var.append(len(t[(t['label'] ==q) & (abs(t['BI1']-t['BI2']) >0)])*100./len(t[t['label']==q]))
                
                si4.append(len(t[(t['label'] ==q) & (t['SiIV-BIO'] >0) & (t['AlIII-BIO'] ==0)])*100./len(t[t['label']==q]))
                
                sial.append(len(t[(t['label'] ==q) & (t['SiIV-BIO'] >0) & (t['AlIII-BIO'] >0)])*100./len(t[t['label']==q]))
            
            
            for (c,j) in zip(oc, range(k)):
                a= alph[j]
                tbl.write("{} & {:d} & {:06.2f} & {:06.2f} & {:06.2f} & {:02.1f} & {:02.1f} & {:02.1f} \n".format(l+"-"+a, c[3], c[4], c[5], c[6], var[j], si4[j], sial[j]))
                
    tbl.close()
    
    return

###################


#check fraction of objects with BI0(AlIII) >0 in each of the CIV and SiIV clusters:

def frac(line, k):

    """ print %'s of objects from each cluster which has BI0_Alll >0
    
    """

    data= Table.read('myBALCat_xtra.csv')

    c= Table.read("./clusters/3features/"+line+str(k)+"clstrs.fits")

    t= join(data, c, keys='SDSSName')
    
    q= 'BIO_AlIII' # change to BIO_SiIV
    #q= 'BIO_SiIV'
    
    clstrs_ls=[]
    for o in range(k):
        clstrs_ls.append([o ,len(t[t['label'] ==o]),\
                          mean(t['Vmin'][t['label'] ==o]),\
                          mean(t['Vmax'][t['label'] ==o]), \
                          mean(t['EW'][t['label'] ==o])])
    
    oc= sorted(clstrs_ls, key= itemgetter(2))

    #print oc

    for x in oc:
        l= x[0]
        print "N= "+str(len(t[t['label'] ==l]))+", N("+q+")= "+str(len(t[(t['label'] ==l) & (t[q] >0)]))+ \
            " = "+str(len(t[(t['label'] ==l) & (t[q] >0)])*100./len(t[t['label']==l]))+"%"

    return

#################

def var_frac(line, k):

    """fraction of objects with variablity measurements (i.e., have different BIs as measured in Filiz Ak et al. 2014)
    """

    data= Table.read('myBALCat_xtra.csv')
    
    c= Table.read("./clusters/3features/"+line+str(k)+"clstrs.fits")
    
    t= join(data, c, keys='SDSSName')
    
    clstrs_ls=[]
    for o in range(k):
        clstrs_ls.append([o ,len(t[t['label'] ==o]),\
                          mean(t['Vmin'][t['label'] ==o]),\
                          mean(t['Vmax'][t['label'] ==o]), \
                          mean(t['EW'][t['label'] ==o])])
    
    oc= sorted(clstrs_ls, key= itemgetter(2))

    #print oc
    
    for x in oc:
        l= x[0]
        print "N= "+str(len(t[t['label'] ==l]))+", "+ \
            "N(var)= "+str(len(t[(t['label'] ==l) & (abs(t['BI1']-t['BI2']) >0)]))+"= " \
            +str(len(t[(t['label'] ==l) & (abs(t['BI1']-t['BI2']) >0)])*100./len(t[t['label']==l]))+"%"
    
    return


#################

## pie charts to show the fraction of objects with vatiablity
#

def pies(line, k):
    
    """pie charts to show the fraction of objects with other absorption troughs and variablity
        """
    
    data= Table.read('myBALCat_xtra.csv')
    
    c= Table.read("./clusters/3features/"+line+str(k)+"clstrs.fits")
    
    t= join(data, c, keys='SDSSName')
    
    clstrs_ls=[]
    for o in range(k):
        clstrs_ls.append([o ,len(t[t['label'] ==o]),\
                          mean(t['Vmin'][t['label'] ==o]),\
                          mean(t['Vmax'][t['label'] ==o]), \
                          mean(t['EW'][t['label'] ==o])])

    oc= sorted(clstrs_ls, key= itemgetter(2))

    var= [] #variability
    si4= [] # SiIV only. no AlIII
    sial=[] # SiIV with or without AlIII
    al3= [] # AlIII abs present

    for x in oc:
        l= x[0]
        '''print "N= "+str(len(t[t['label'] ==l]))+", "+ \
            "N(var)= "+str(len(t[(t['label'] ==l) & (abs(t['BI1']-t['BI2']) >0)]))+"= " \
            +str(len(t[(t['label'] ==l) & (abs(t['BI1']-t['BI2']) >0)])*100./len(t[t['label']==l]))+"%"
            '''
        
        var.append(len(t[(t['label'] ==l) & (abs(t['BI1']-t['BI2']) >0)])*100./len(t[t['label']==l]))
        
        si4.append(len(t[(t['label'] ==l) & (t['BIO_SiIV'] >0) & (t['BIO_AlIII'] ==0)])*100./len(t[t['label']==l]))
        
        sial.append(len(t[(t['label'] ==l) & (t['BIO_SiIV'] >0)])*100./len(t[t['label']==l]))
        
        al3.append(len(t[(t['label'] ==l) & (t['BIO_AlIII'] >0)])*100./len(t[t['label']==l]))

    print var

    ex1, ex2, ex3, ex4= [], [], [], []
    
    for i in var:
        if i == max(var):
            ex1.append(0.1)
        else:
            ex1.append(0)

    print var
    print ex1

    for i in si4:
        if i == max(si4):
            ex2.append(0.1)
        else:
            ex2.append(0)

    print si4
    print ex2

    for i in sial:
        if i == max(sial):
            ex3.append(0.1)
        else:
            ex3.append(0)

    for i in al3:
        if i == max(al3):
            ex4.append(0.1)
        else:
            ex4.append(0)


    clr_ls = [sns.xkcd_rgb["windows blue"], sns.xkcd_rgb["dusty purple"], sns.xkcd_rgb["pale red"], \
          sns.xkcd_rgb["greyish"], sns.xkcd_rgb["faded green"], sns.xkcd_rgb["amber"]]

    labels = ['a', 'b', 'c', 'd', 'e', 'f']


    fig= figure(figsize=(12,12))

    ax1= fig.add_subplot(221)

    pie(var, shadow= True, colors= clr_ls, startangle=30, radius= 1, autopct='%1.1f%%', \
        labels= labels[:k], explode= ex2)

    text(0.1, 0.001, "Var", transform= ax1.transAxes, size= 14, color='k', family= 'serif')


    ax2= fig.add_subplot(222)

    pie(si4, shadow= True, colors= clr_ls, startangle=30, radius= 1, autopct='%1.1f%%', \
        labels= labels[:k], explode= ex2)

    text(0.1, 0.001, "SiIV only abs", transform= ax2.transAxes, size= 14, color='k', family= 'serif')

    text(0.3, 0.99, line+", K= "+str(k), transform= ax2.transAxes, size= 16, color='red', family= 'serif')

    ax3= fig.add_subplot(223)
    
    pie(sial, shadow= True, colors= clr_ls, startangle=30, radius= 1, autopct='%1.1f%%', \
        labels= labels[:k], explode= ex3)
        
    text(0.1, 0.001, "SiIV abs", transform= ax3.transAxes, size= 14, color='k', family= 'serif')
        
    text(0.3, 0.99, line+", K= "+str(k), transform= ax3.transAxes, size= 16, color='red', family= 'serif')


    ax4= fig.add_subplot(224)

    pie(al3, shadow= True, colors= clr_ls, startangle=30, radius= 1, autopct='%1.1f%%',\
        labels= labels[:k], explode= ex4)

    text(0.1, 0.001, "AlIII abs", transform= ax4.transAxes, size= 14, color='k', family= 'serif')


    return



