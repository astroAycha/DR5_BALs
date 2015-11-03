""" this will soon (!) include code to figure out the properties of each cluster with regards to the extra data added to the Gibson catalog (mainly reddening and variablity. and maybe BH masses...)
"""

import numpy as np
from astropy.table import Table

def clstr_cntrs():

    """ create a table with the centroids for the clusters for k= 3,4 and 5. formated for the manuscript
    """

    line_ls= ["CIV", "SiIV", "AlIII", "MgII"]

    tbl= open('cntrs_tbl.txt', 'wrb')
    
    #vmin, vmax, ew, num=[], [], [], []

    for k in range(3,6):
    
        for l in line_ls:
        
            t= Table.read("./clusters/3features/"+l+str(k)+"clstrs.fits")
            tbl.write(l+"-"+str(len(t)) + "\t"+ "K="+ str(k)+ "\n")
            
            for c in range(k):
                tbl.write(str(mean(t['Vmin'][t['label'] ==c]))+"\t"+ str(mean(t['Vmax'][t['label'] ==c]))+"\t"+ str(mean(t['EW'][t['label'] ==c]))+"\t"+ str(len(t[t['label'] ==c]))+"\n")
            
    tbl.close()
            
    return