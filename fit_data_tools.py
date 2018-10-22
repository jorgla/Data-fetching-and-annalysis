import numpy as np
import random
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import h5py
from os import walk
from mpl_toolkits.axes_grid1 import make_axes_locatable

def ploy(x, a, b, c):
    return a*x**b+c

def poly2(x, a, b, c):
    return a*x**2 + b*x + c

def gaussian(x,a,b,c,d):
    return np.abs(a)*np.exp(-4*np.log(2)*(x-b)**2./(c**2))+d
# where c is FWHM

def gaussian0(x,a,b,c):
    return np.abs(a)*np.exp(-(x-b)**2/(2*c**2))+d
# where c is standerd devation 

def fit_poly(function,x,y,p0=None,sigma=None,bounds=None):
    popt,pcov = curve_fit(function,x,y,p0=p0,sigma=sigma, bounds=bounds)
    x = np.arange(0,1e4)
    curve = function(x,*popt)
    perr = np.sqrt(np.diag(pcov))
    return popt,x,curve,perr

def fit(function,x,y,p0=None,sigma=None,bounds=(-np.inf, np.inf)):
    popt,pcov = curve_fit(function,x,y,p0=p0,sigma=sigma, bounds=bounds)
    x = np.arange(0,3e3,1e-2)
    curve = function(x,*popt)
    perr = np.sqrt(np.diag(pcov))
    return popt,x,curve,perr

def discover_datas(path):
    for (dirpath, dirnames, datanames) in walk(path):
        break
        
    if len(datanames) != 0 :
        mask = np.ones(len(datanames), dtype=bool) 
        for i in np.arange(len(datanames)) :
            if datanames[i][0] == '.' : 
                mask[i] = 0
        datanames = np.array(datanames)
        datanames = datanames[mask]
        
    return datanames