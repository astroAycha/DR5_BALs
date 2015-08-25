""" unsupervised clustering on BALQ measurements form the Gibson et al. BALQ catalog from SDSS-DR5
"""

import numpy as np

from astropy.io import fits
from astropy.table import Table

from sklean.cluster import KMeans

from scipy.stats import sigmaclip


data= Table.read()