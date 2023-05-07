import warnings
warnings.filterwarnings('ignore')
# import sys
# import os

import numpy as np
import ccdprocx
from astropy.time import Time
from astropy.table import QTable, Table
import astroalign as aa
from astropy.wcs import WCS
# from astropy.utils.console import ProgressBar, ProgressBarOrSpinner
from tqdm import tqdm
# from astropy.coordinates import SkyCoord as acoord
# import astropy.units as un
from astropy.io import fits

from astropy.nddata.ccddata import CCDData
CCDData._config_ccd_requires_unit = False
# photometry imports
# from photutils.psf import IntegratedGaussianPRF, DAOGroup
# from photutils.background import MMMBackground, MADStdBackgroundRMS
# from astropy.modeling.fitting import LevMarLSQFitter
# from astropy.stats import gaussian_sigma_to_fwhm
# from photutils.psf import IterativelySubtractedPSFPhotometry
from photutils.aperture import CircularAperture, aperture_photometry, CircularAnnulus
# from photutils import DAOStarFinder
# from astropy.stats import mad_std
# from astroquery.simbad import Simbad

from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground



'''
Ceres is the handler of observational image series in Dorado.
'''

__all__ = ['Ceres']

class Ceres:
    '''
        The Ceres class encapsulates a set of astronomical data from a single night of observation.
        Ceres can handle multiple stacks of data in different filters and perform a variety of
        actions on them in an orgainized fashion.

        Attributes
        ----------

        filters: dictionary

        data: Stack array

        bias: CCDdata

        time: 'astropy.Time' 

        datestr: str

        
    '''
    
    def __init__(self, filters = {}, data = [], bias = None, time = None, datestr = None):
        # metadata
        self.filters = filters
        self.data = data
        self.time = time
        self.datestr = datestr

        self.bias = bias
        
        # location, weather, timezone, camera, observer array

        # TODO figure out how to clean this up
        try:
            self.date = Time(int(self.time.mjd), format = 'mjd')
        except:
            self.date = None
        
        if datestr == None:
            try:
                epoch = self.date.ymdhms
                if (epoch['hour'] + Dorado.UTCoffset) < 0:
                    day = str(epoch['day'] - 1)
                    day2 = str(epoch['day'])
                    month = str(epoch['month'])

                    if (epoch['day'] - 1) < 10:
                        day = '0' + str(epoch['day'] - 1)
                    if epoch['day'] < 10:
                        day2 = '0' + str(epoch['day'])

                    if epoch['month'] < 10:
                        month = '0' + str(epoch['month'])

                else: 
                    day = str(epoch['day'])
                    day2 = str(epoch['day'] + 1)
                    month = str(epoch['month'])

                    if epoch['day'] < 10:
                        day = '0' + str(epoch['day'])
                        if epoch['day'] < 9:
                            day2 = '0' + str(epoch['day'] + 1)

                    if epoch['month'] < 10:
                        month = '0' + str(epoch['month'])

                self.datestr = str(epoch['year']) + '-' + month + '-' + day + '+' + day2
            except:
                self.datestr = datestr
        
    def add_stack(self, stack):
        '''
        Add a Dorado stack to the current Ceres instance (self).
        
        Parameters
        ----------
        stack: dorado.stack
            Instance of dorado.stack class to add to self.
        
        '''
        
        # eventually stacks themelves should have some metadata 
        # to denote stuff like calibration status
        # TODO do I need to add times?
        self.filters[stack.filter] = len(self.data)
        self.data.append(stack)
        
    def rem_stack(self, filter):
        '''
        Remove a Dorado stack to the current Ceres instance (self).
        
        Parameters
        ----------
        filter: str
            String representation of the relevent filter to remove from self.
        
        '''
        # TODO delete time strings
        del self.data[self.filters[filter]]
        

class Bias_Ceres:
    def __init__(self, file_list = [], data = [], datestr = ''):
        self.time = time
        try:
            self.date = Time(int(self.time.mjd), format = 'mjd')
        except:
            self.date = None
        if datestr == '':
            try:
                epoch = self.date.ymdhms
                if (epoch['hour'] + Dorado.UTCoffset) < 0:
                    day = str(epoch['day'] - 1)
                    day2 = str(epoch['day'])
                    month = str(epoch['month'])

                    if (epoch['day'] - 1) < 10:
                        day = '0' + str(epoch['day'] - 1)
                    if epoch['day'] < 10:
                        day2 = '0' + str(epoch['day'])

                    if epoch['month'] < 10:
                        month = '0' + str(epoch['month'])

                else: 
                    day = str(epoch['day'])
                    day2 = str(epoch['day'] + 1)
                    month = str(epoch['month'])

                    if epoch['day'] < 10:
                        day = '0' + str(epoch['day'])
                        if epoch['day'] < 9:
                            day2 = '0' + str(epoch['day'] + 1)

                    if epoch['month'] < 10:
                        month = '0' + str(epoch['month'])

                self.datestr = str(epoch['year']) + '-' + month + '-' + day + '+' + day2
            except:
                self.datestr = datestr

        # if np.any(file_list)