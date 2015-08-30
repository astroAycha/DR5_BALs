"""This script was written for the clustering project on the SDSS-DR10.
Modified on Aug 27th 2015 to use with the BAL quasars in the BALQ catalog on DR5 (Gibson et al. 2009).

Read a raw spectrum from the SDSS QSO DR5 catalog.
    correct for Galactic extinction and redshift and save as a new fits file
    """

import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy import units as u
from specutils import extinction


def spec_proc(data_tab):

    s= Table.read(data_tab)


    for i in range(len(s)):
        #print i
        z= s['z'][i]
        Av= s['AV_SandF'][i] #Av (Schlafly and Finkbeiner 2011) see extinction_tab.txt
        spec_name= './data/spSpec-'+str(s['MJD_spec'][i])+'-'+str(s['plate'][i]).zfill(4)+'-'+str(s['fiberid'][i]).zfill(3)+'.fit' # fits file name
        spec= fits.open(spec_name) # read file
        flx= spec[0].data[0] # flux
        wlen= 10.**(spec[0].header['coeff0']+ spec[0].header['coeff1'] * np.arange(len(spec[0].data[0]))) # wavelength, coeff0: starting wavelength, coeff1: dispersion
        ext_lambda= extinction.extinction_ccm89(wlen * u.angstrom, a_v= Av, r_v= 3.1)
        tau_lambda= ext_lambda/((1 * u.angstrom) * 1.086) # for some reason the extinction has a unit of angstrom (should be magnitude i.e. unitless)
        dered_flx= flx * np.exp(tau_lambda) # extinction corrected flux
        rest_wlen= wlen /(1+z) # restframe wavelength
        x= np.arange(1100, 4000, 0.5) # the wavelength values for the rebinned spectrum with bin width = 0.5 A
        rebin_dered_flx= np.interp(x, rest_wlen, dered_flx) # resampled flux with delta lambda= 0.5 A
        proc_spec= np.vstack((x, rebin_dered_flx)) # 2D array. 1st row: restframe wavelength, 2nd row: corrected flux
        hdu= fits.PrimaryHDU(proc_spec) # HDU with flux and wavelngth to be saved as a fits file
        new_spec_name= './proc_data/spec-'+str(s['plate'][i])+'-'+str(s['MJD_spec'][i])+'-'+str(s['fiberid'][i]).zfill(4)+'_proc.fits' # name of new processed spectrum
        hdu.writeto(new_spec_name)

    return