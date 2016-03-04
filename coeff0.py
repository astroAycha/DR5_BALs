
""" quick script to read fits files and compare the COEFF0 keyword in their headers.
    Compare SDSS-DR5 spectra from an old repo to the ones in a new repo
"""

from astropy.io import fits
import random
from glob import glob

def coeff():
    """ compare COEFF0 for a 100 objects selected randomly
    """
    
    all_old_spec= glob('./old_data/*fit')
    old_spec_ls= random.sample(all_old_spec, 100)
    
    d=[]
    
    for s in old_spec_ls:

        mjd= s[18:23]
        p= s[24:28]
        f= s[29:32]
        
        #print mjd, p, f
    
        new_name= "./data/spec-"+p+"-"+mjd+"-"+f.zfill(4)+".fits"
        
        new_spec= fits.open(new_name)
        old_spec= fits.open(s)
        
        #n_coeff= 10**new_spec[1].data.field(1)[0]
        n_coeff= 10**new_spec[0].header['COEFF0']
        o_coeff= 10**old_spec[0].header['COEFF0']
    
        d.append(n_coeff-o_coeff)

    hist(d, bins=20, histtype='step', lw= 2)

    xlabel('new COEFF0 - old COEFF0')
    ylabel('Num')