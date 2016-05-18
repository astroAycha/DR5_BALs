## 18 May 2016
## testing using manifold learning on spectra to reduce dimentionality and apply clustering

import numpy as np

import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter

from sklearn import manifold

from astropy.io import fits
from astropy.table import Table

Axes3D

## read the spectra and create a 2D
data= Table.read('myBALs.fits')

Z= np.arange(1100, 4000, 0.5)

for s in range(1000):
    spec_name= "./proc_data/spec-"+str(data['plate'][s])+"-"+str(data['MJD_spec'][s])+"-"+str(data['fiberid'][s]).zfill(4)+"_proc.fits"
    spec= fits.open(spec_name)
    wlen= spec[0].data[0]
    flx= spec[0].data[1]
    norm_flx= flx/np.median(flx[2360:2390]) # normalize spectra
    Z= np.vstack((Z, norm_flx))



X= Z[:,1820:1960]
n_neighbors= 10
n_components= 3

fig = plt.figure(figsize=(15, 8))
plt.suptitle("Manifold Learning with %i points, %i neighbors"
             % (1000, n_neighbors), fontsize=14)

try:
    # compatibility matplotlib < 1.0
    ax = fig.add_subplot(251, projection='3d')
    ax.scatter(X[:, 0], X[:, 1], X[:, 2])
    ax.view_init(4, -72)
except:
    ax = fig.add_subplot(251, projection='3d')
    plt.scatter(X[:, 0], X[:, 2])

methods= ['standard', 'ltsa', 'hessian', 'modified']
labels= ['LLE', 'LTSA', 'Hessian LLE', 'Modified LLE']


for i, method in enumerate(methods):

    Y= manifold.LocallyLinearEmbedding(n_neighbors, n_components,
                                       eigen_solver= 'auto',
                                       method= method).fit_transform(X)
                                       
    print Y.shape

    ax = fig.add_subplot(252 + i)
    plt.scatter(Y[:, 0], Y[:, 1])
    plt.title("%s" % (labels[i]))
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())
    plt.axis('tight')
