
import numpy as np
from astropy.io import fits

from astropy.table import Table, join


def spec_compo_compare(line, k, g, c):

    """ plot comosite spectrum and overplot individual spectra for comparison.
        param: line: used in the clustering CIV, SiIV, AlIII, or MgII
                k: number of clusters: 3-8
                g: group name (each group contains clusters generated using different set of clustering features (see bal_clusters.py)
                 g0,... g9
                c: cluster number: 0 up to k-1
    """

    compo_name= "./composites/"+g+"/"+line+"_"+str(k)+"clstr"+str(c+1)+".fits"
    
    compo= fits.open(compo_name)

    data = Table.read("myBALs.fits")
    
    clstr_data= Table.read("./clusters/"+g+"/"+line+str(k)+"clstrs.fits")
    
    clstr= clstr_data[clstr_data['label'] == c]
    
    t= join(clstr, data, keys= "SDSSName")
    
    print len(t)
    
    fig= figure(figsize=(16,10))
    
    for i in range(len(t)):
        try:
            plot(compo[0].data[0], compo[0].data[1], color= 'r', lw= 2)
            spec_name= "./proc_data/spec-"+str(t['plate'][i])+"-"+str(t['MJD_spec'][i])+"-"+str(t['fiberid'][i]).zfill(4)+"_proc.fits"
        
            spec= fits.open(spec_name)
            
            plot(spec[0].data[0], spec[0].data[1]/spec[0].data[1][(2150-1100)*2], color= '0.4')
            xlim(1300,2200)
            ylim(-1,5)

            resume = input("Press Enter to plot next spectrum on list.")
        
        except SyntaxError:
            pass
            clf()



##############
def spec_display(spec_ls, n1, n2):
    
    """ read a list of spectra and display them. Read input and use as flag (for either low SNR or BAL quasar).
        
        spec_ls: numpy array with the quasar sample as selected in quasar_cluster.py
        flag= 0 keep
        flag= 1 reject
        
        n1: start at line number n1
        n2: stop at line number n2
        
        """
    
    data= np.load(spec_ls)
    
    sample= data[n1:n2+1]
    print "Looking at lines", n1, "to", n2
    
    flag_ls=[]
    names=[]
    
    wavelen= np.arange(1100, 4000, 0.1)  #wavelength array
    fig= figure(figsize(20,8))
    
    for i in range(len(sample)):
        print "Looking at spectrum number", i+1
        try:
            spectrum_name= "./proc_data/spec-"+str(sample['PLATE'][i])+"-"+str(sample['MJD'][i])+"-"+str(sample['FIBERID'][i]).zfill(4)+"_proc.fits"
            spec= fits.open(spectrum_name)
            flx= spec[0].data[1]
            
            plot(wavelen, flx)
            xlim(1200, 3100)
            ylim(-1,4)
            axvline(1397, ls=':')
            axvline(1549, ls=':')
            axvline(1908, ls=':')
            axvline(2800, ls=':')
            text(2500, 3, sample['SDSS_NAME'][i])
            print "Flags: 0= keep, 1= reject"
            flag= input()
            flag_ls.append(flag)
            names.append(sample['SDSS_NAME'][i])
            resume = input("Press Enter to plot next spectrum on list.")
        
        except SyntaxError:
            pass
            clf()

new_array= np.column_stack((names, flag_ls))
    save("myflags_"+str(n1)+"_to_"+str(n2)+".npy", new_array)
