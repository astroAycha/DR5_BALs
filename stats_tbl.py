""" this will soon (!) include code to figure out the properties of each cluster with regards to the extra data added to the Gibson catalog (mainly reddening and variablity. and maybe BH masses...)
"""

import numpy as np
from operator import itemgetter
from astropy.table import Table

def clstr_cntrs():

    """ create a table with the centroids for the clusters for k= 3,4 and 5. formated for the manuscript
    """

    line_ls= ["CIV", "SiIV", "AlIII", "MgII"]

    tbl= open('cntrs_tbl.txt', 'wrb')
    
    #vmin, vmax, ew, num=[], [], [], []

    for l in line_ls:
    
        for k in range(3,6):
        
            t= Table.read("./clusters/3features/"+l+str(k)+"clstrs.fits")
            
            tbl.write(l+ "&" + str(len(t)) + "& "+ "K="+ str(k)+ "\n")
            
            clstrs_ls=[]
            
            for o in range(k):
                clstrs_ls.append([l, k, o ,len(t[t['label'] ==o]),mean(t['Vmin'][t['label'] ==o]),\
                                  mean(t['Vmax'][t['label'] ==o]), mean(t['EW'][t['label'] ==o])])
    
            oc= sorted(clstrs_ls, key= itemgetter(4)) #ordered clusters
            
            for c in oc:
                tbl.write("{:d}, ({:06.2f},{:06.2f},{:06.2f}) \n".format(c[3], c[4], c[5], c[6]))
            
           
    tbl.close()
    
    return

###################

#check fraction of objects with BI0(AlIII) >0 in each of the CIV and SiIV clusters:

def al_frac(line, k, q):

    """ print %'s of objects from each cluster which has a quantity >0
    
    """

    data= Table.read('myBALCat_xtra.csv')

    c= Table.read("./clusters/3features/"+line+str(k)+"clstrs.fits")

    t= join(data, c, keys='SDSSName')

    for l in range(k):
        print len(t[t['label'] ==l]), len(t[(t['label'] ==l) & (t[q] >0)])*100./len(t['label']==l)

    return

#################
def var_frac(line, k):

    """fraction of objects with variablity measurements (i.e., have different BIs as measured in Filiz Ak et al. 2014)
    """

    data= Table.read('myBALCat_xtra.csv')
    
    c= Table.read("./clusters/3features/"+line+str(k)+"clstrs.fits")
    
    t= join(data, c, keys='SDSSName')
    print len(t)

    for l in range(k):
        print len(t[t['label'] ==l]), len(t[(t['label'] ==l) & (abs(t['BI1']-t['BI2']) >0)]), len(t[(t['label'] ==l) & (abs(t['BI1']-t['BI2']) >0)])*100./len(t['label']==l)
    
    return



