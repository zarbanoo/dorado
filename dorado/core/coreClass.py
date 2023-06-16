import warnings
warnings.filterwarnings('ignore')

from ..ceres import Ceres
from ..stack import Stack
from ..core.readerClass import *
from ..dorphot.dorphotClass import *
from ..target.targetClass import *
from ..core.cacheClass import *

import ccdprocx

from astropy.nddata.ccddata import CCDData
CCDData._config_ccd_requires_unit = False
import numpy as np
from astropy import config as _config
from astropy.utils.misc import isiterable
import astropy.units as un
from astropy.time import Time


from astroquery.astrometry_net import AstrometryNet
from astroquery.exceptions  import TimeoutError
from astropy.wcs import WCS

import os
import datetime
from pathlib import Path


ast = AstrometryNet()

'''
Filer is the handler of the Dorado system,
'''

__all__ = ['Dorado', 'Dorado_core']


def get_builtins():
    """Due to the way Python works, ``__builtins__`` can strangely be either a module or a dictionary,
    depending on whether the file is executed directly or as an import. I couldn’t care less about this
    detail, so here is a method that simply returns the namespace as a dictionary."""
    return getattr( __builtins__, '__dict__', __builtins__ )


class Dorado_core:
    '''
    
    '''
    def __init__(self):
        # open and use logger
        # make function to create data class from hardware, processed, or raw data folder
        # needs function to import data into raw
        # fix get_night()
        # add new calibration frames to folders
        # find most recent calibration frame if none
        # -> figure out if none
        # clear cache function
        self.stardir = os.getcwd()
        self.rootname = 'dorado'
        self.config_dir = _config.get_config_dir(self.rootname)
        # read in config from config_dir
        
        self.dordir = Path(self.config_dir).parent
        self.init_dir()
        self.unit = un.adu

        ## class storage zone -----------------------
        # redo this so theres an easy storage and readout
        # reader
        self.reader = aico_reader() # ussually we'd want to grab this from a configuration file, should name be aicoReader?
        self.cache = cache_core()   # add cache handler to Dorado core instance
        # dorphot
        self.dorphot = aicoPhot()
        # ceres & stacks
        self.ceres_keys = {}
        self.ceres = [] # maybe call this series?
        # self.filters[stack.filter] = len(self.data)
        # self.data.append(stack)

        # targets
        self.target_keys = {}
        self.targets = [] # TODO choose target class
        # all early ts instances will belong to a target;
        # timeseries
        self.timeseries = [] # exclude this for the time being as 

        ### read from config
        # site information (hardcoded atm)
        self.UTCoffset = -5

    def target(self, key):
        return self.targets[self.target_keys[key]]
    
    def cere(self, key):
        return self.ceres[self.ceres_keys[key]]
    
    def init_dir(self, tess = False):
        """
        init_dir initializes the dorado home directory '$user/.dorado/'.
        
        Parameters
        ----------
        tess: boolean
            Whether to initialize a TESS data directory in './data/'
        Returns
        -------
        
        """
        self.enter_dordir()
        os.makedirs('./data/wrk',            exist_ok = True)
        os.makedirs('./data/flats',          exist_ok = True)
        os.makedirs('./data/bias',           exist_ok = True)
        os.makedirs('./data/darks',          exist_ok = True)
        os.makedirs('./data/raw',            exist_ok = True)
        os.makedirs('./data/graphical',      exist_ok = True)
        os.makedirs('./data/projects',       exist_ok = True)
        os.makedirs('./data/targets',        exist_ok = True)
        os.makedirs('./logs',                exist_ok = True)
        os.makedirs('./cache',               exist_ok = True)
        os.makedirs('./cache/astrometryNet', exist_ok = True)
        if tess:
            os.makedirs('./data/tess',       exist_ok = True)
        self.exit_dordir()
        
    def enter_dordir(self):
        """
        enter_dordir changes the workin directory to self.dordir which by default is
        '$user/.dorado/'.
        
        Parameters
        ----------

        Returns
        -------
        
        """
        os.chdir(self.dordir)

    def exit_dordir(self):
        """
        exit_dordir changes the workin directory to the starting directory at runtime.
        
        Parameters
        ----------

        Returns
        -------
        
        """
        os.chdir(self.stardir)
        
    def mkceres(self, date, name = None, sub = 'raw', target = None, calibrated = False, aligned = False):
        """
        mkceres is a wrapper function for self.reader.mkceres with an additional arguement for setting the keyname 
        for the ceres instance in self.ceres_keys. The default keyname is the datestring for the ceres instance.
        
        Parameters
        ----------
        name: string
            keyname for the ceres instance. Default is the datestring for the ceres instance.
        
        self.reader.mkceres args** 
        
        Returns
        -------
        
        """
        if target == None:
            ctemp = self.reader.mkceres(date, sub = sub, target = target, calibrated = calibrated, aligned = aligned)
        else:
            ctemp = self.reader.mkceres(date, sub = sub, target = self.targets[self.target_keys[target]], calibrated = calibrated, aligned = aligned)
        # add logic to print out a summary

        if name == None:
            name = ctemp.datestr
            # TODO handle if no date found either
            print('No series nickname given, defaulting to date string: ', name)
        else:
            print('Call series as: ', name, '\n')
        self.ceres_keys[name] = len(self.ceres)
        self.ceres.append(ctemp)
    
    def mktrgt(self, name, constructor = Target):
        """
        mktrgt takes a target name and constructs a dorado.target object that is returned to self.targets.
        
        Parameters
        ----------
        name: string
            target name as found in SIMBAD
        Returns
        -------
        
        """
        # TODO allow for nickname in keys
        # add dictionary
        # , coordinates = None
        self.target_keys[name] = len(self.targets)
        self.targets.append(constructor(name))
        
    def savewrk(self, cr, filters = None, save_series = True, save_base = True):
        """
        savewrk is a wrapper function for self.reader.savewrk.
        Parameters
        ----------

        self.reader.savewrk args**
        Returns
        -------
        
        """
        self.reader.savewrk(cr, filters, save_series, save_base)

    # logic to set reader and dorphot
        
        
        
    # move to ceres
    def getDateString(self, cr, utc = None):
        """
        getDateString computes the date string for the given ceres object. ceres time is represented
        in UTC while datestring is represented in local time. A UTC offset must be provided either by self.UTCoffest
        or by the utc arguement.
        
        Parameters
        ----------
        cr : str
            the ceres object key to get the date string for.
        utc: float
            utc offset to apply to the date string. Default is None.

        Returns
        -------
        
        """

        ## TODO :: look into UTC wrecking stuff
        # if the hour is less than the utc offset of the site then the utc date is one ahead of local time
        epoch = self.ceres[self.ceres_keys[cr]].date.ymdhms

        # epoch = date.ymdhms
        if (utc):
            
            if (epoch['hour'] + utc) < 0:
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
        else:
            if (epoch['hour'] + self.UTCoffset) < 0:
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


        datestr = str(epoch['year']) + '-' + month + '-' + day + '+' + day2
        self.ceres[self.ceres_keys[cr]].datestr = datestr
        
    # move to reader
    def mkwrk(self, cr):
        """
        mkwrk produces a working folder for saving a ceres object to the dorado data working directory
        for later use.
        
        Parameters
        ----------
        cr : str
            ceres key for ceres instance to create working directory for.
        Returns
        -------
        
        """
        # TODO this should be a wrapper for reader
        if self.ceres[self.ceres_keys[cr]].datestr == None:
            self.getDateString(cr) 
        datestr = self.ceres[self.ceres_keys[cr]].datestr
        wrkdir = self.dordir / 'data' / 'wrk'
        os.makedirs(wrkdir / datestr,                  exist_ok = True)
        os.makedirs(wrkdir / datestr / 'aligned',      exist_ok = True)
        os.makedirs(wrkdir / datestr / 'calibrated',   exist_ok = True)
        os.makedirs(wrkdir / datestr / 'uncalibrated', exist_ok = True)
        os.makedirs(wrkdir / datestr / 'WCS',          exist_ok = True)
        # figures, targets, log, omitted images, observation metadata

    # move to utils or reader
    def saveWCS(self, cr, filters = None):
        """
        saveWCS takes a ceres object and an optional filters argument and saves the contained WCS
        data to the ceres objects working directory.
        
        Parameters
        ----------
        cr : str
            ceres object key string for the relevent ceres object.
        filters : str or list of str
            filters of WCS data to be saved.

        Returns
        -------
        
        """
        wrkdir = self.dordir / 'data' / 'wrk'
        if self.ceres[self.ceres_keys[cr]].datestr == None:
            self.getDateString(cr)
        datestr = self.ceres[self.ceres_keys[cr]].datestr
        self.mkwrk(cr)
        if filters == None:
            filters = self.ceres[self.ceres_keys[cr]].filters.keys()
        
        for filter in filters:
            stack = self.ceres[self.ceres_keys[cr]].data[self.ceres[self.ceres_keys[cr]].filters[filter]]
            if stack.wcs == None:
                self.dorphot.getWCS(filter, self) # TODO this is now in dorphot, should it be moved?
            fname = str(filter) + '-solved.fits'
            solved = self.ceres[self.ceres_keys[cr]].data[self.ceres[self.ceres_keys[cr]].filters[filter]].solved
            solved.data = solved.data.astype('uint16')
            solved.mask = None
            solved.uncertainty = None
            solved.write(wrkdir / datestr / 'WCS' / fname, overwrite = True)
        # not done might be redundant atm also shouldnt this belong to reader?
    
    # move to ceres or reader
    def saveBase(self, cr, filters = None):
        """
        saveBase save the base frame of a ceres object if any.
        Parameters
        ----------
        cr: string
            ceres key string for relevent ceres object.
        filters: str or list of str
            filters to save base frame of a ceres object.
        Returns
        -------
        
        """
        wrkdir = self.dordir / 'data' / 'wrk'
        if self.ceres[self.ceres_keys[cr]].datestr == None:
            self.getDateString(cr)
        datestr = self.ceres[self.ceres_keys[cr]].datestr

        print('Saved to data/wrk/', datestr)
        self.mkwrk(cr)

        if filters == None:
            filters = self.ceres[self.ceres_keys[cr]].filters.keys()
        
        for filter in filters:
            imname = filter + '_base.fits'
            fname = wrkdir / datestr / imname
            base = self.ceres[self.ceres_keys[cr]].data[self.ceres[self.ceres_keys[cr]].filters[filter]].base
            base.data        = base.data.astype('uint16') 
            # base.mask        = base.mask.astype('uint16')
            # base.uncertainty = base.uncertainty.astype('uint16')
            base.mask = None
            base.uncertainty = None
            base.write(fname, overwrite = True)
        # not done might be redundant atm also shouldnt this belong to reader?
        
        # merge header

Dorado = Dorado_core()

# this gone pollute the namespace but ay lmao
G = get_builtins()
G[ 'G' ] = G
G['Dorado'] = Dorado



# we need to pass the calibration files(if none grab the most recent out of 
# dorado by closest mjd)
# if there is no calibration files in dorado send a warning and set 
# calibration files to none
# handle multi filter

