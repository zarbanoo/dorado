import warnings
warnings.filterwarnings('ignore')

from ..ceres import Ceres
from ..stack import Stack
from ..core.coreClass import *
from ..core.utils import *

# import utils
import ccdprocx

# from astropy.nddata.ccddata import CCDData
# CCDData._config_ccd_requires_unit = False
import astropy.units as un
from astropy.time import Time
from astropy.io import fits
from astropy.utils.misc import isiterable

import lightkurve as lk
from lightkurve.lightcurve import TessLightCurve as tlc

import os
from tqdm import tqdm
import numpy as np

__all__ = ['reader', 'aico_reader'] #, 'tess_reader'

class reader:
    '''
        readers need to have a dirscan method and a savewrok method. Readers should also probably have a mkceres method.
    '''
    def __init__(self):
        # should report to logger and read from config here
        # will this carry over to inherited classes?
        self.desired_datatype = 'uint16' # datatype to force onto images
        
    def diread(self, dirarray):
        """
        diread intakes a filepath array and catalogues the contents by either 'is file' or 'is directory'
        and returns the resulting lists.
        
        Parameters
        ----------
        dirarray: str array
            array containing the path to the desired directory.
        Returns
        -------
        
        """
        if isiterable(dirarray):
            path = Dorado.dordir
            for dir in dirarray:
                path = path / dir
        else:
            path = dirarray
        # path = self.dordir / 'data' / 'raw' / date
        contents = os.scandir(path = path)
        files = []
        directories = []
        for entry in contents:
            if not entry.name.startswith('.'):
                if entry.is_file():
                    files.append(entry)
                if entry.is_dir():
                    directories.append(entry)
        return files, directories
    
    def newdat(self):
        """
        newdat is an nfinished function that will scan the '$user/.dorado/data/raw' directory
        for new data to be processed.
        
        Parameters
        ----------

        Returns
        -------
        
        """
        # find data that hasn't been processed yet
        print('searching for unprocessed data...')
    
    def force16(self, hdu):
        """
        An unfinished convinience function for forcing the datatype of CCDData to be 16-bit instead of 32-bit
        for consistency and filesize optimization.
        
        Parameters
        ----------
        hdu: CCDData
            hdu to force 16-bit datatype.
        Returns
        -------
        
        NOTE:: in the future this will handle arrays of data instead of single HDU's and will allow for specifying the bit datatype.
        """
        print('This function is not implemented yet. See mkBias() for example functionality.')
    



# TODO Make it so that each class isnt loaded unless chosen (core class takes string and then imports correct class)
class aico_reader(reader):
    '''
    
    '''
    def __init__(self):
        self.biasstr = ['Bias', 'Bias', 'bias', 'BIAS']
        self.flatsstr = ['FLAT', 'FlatField', 'flat', 'Flat', 'Flats', 'flats', 'FLATS', 'FlatFields']
        self.lightsstr = ['lights', 'Lights', 'LIGHTS']
        self.stack_types = {'bias':self.biasstr, 'flats':self.flatsstr, 'lights':self.lightsstr}
        self.unique_lights = True # whether to find lights by non calibration files (True) or via self.lightstr (False)
        self.require_cal = False  # whether calibration frames are essential, this should be dynamically set by the read calibration level
        self.calibration_window = 7
        
    def _single_level(self, files):
        bias    = self._read_stack(files, 'bias')
        flats  = self._read_stack(files, 'flats')
        lights = self._read_stack(files, 'lights')
        return bias, flats, lights
        
    def _multi_level_top(self, files, directories):
        biasdir = [s for s in directories if s.name in self.biasstr]
        # There should never be multiple bias directories; unless bias were taken with multiple binning settings or taken multiple times.
        # Consult chief data officer about this. TODO
        if len(biasdir) > 0:
            bias_files, _ = self.diread(biasdir[0])
            bias    = self._read_stack(bias_files, 'bias')
        else:
            bias = []
        # It is possible to have multiple flats directories for different filters, they should all be contained in a single base flats directory tho
        flatsdir = [s for s in directories if s.name in self.flatsstr]
        if len(flatsdir) > 0:
            flats = self._multi_level_sub(flatsdir[0], 'flats')
        else:
            flats = []
        if self.unique_lights:
            lightsdir = [s for s in directories if (s.name not in self.flatsstr) and (s.name not in self.biasstr)] # may hand back more than one directory
        else:
            lightsdir = [s for s in directories if s.name in self.lightsstr]  # may hand back more than one directory
        if len(lightsdir) > 0: # should also handle more than one rogue dir
            lights = self._multi_level_sub(lightsdir[0], 'lights')
        else:
            lights = []
        return bias, flats, lights

    def _multi_level_sub(self, subdir, stack_type):
        files, directories = self.diread(subdir)
        sub_files = []
        if len(directories) == 0 & len(files) > 0:
            print('Single directory ', stack_type, ' organization format detected.')
            sub_files  = self._read_stack(files, stack_type)
        elif len(files) == 0 & len(directories) > 0:
            print('Multi directory ', stack_type, ' organization format detected.')
            for d in directories:
                # should not be any further sub folders
                f, _ = self.diread(d)
                sub_files += self._read_stack(f, stack_type) # should we continue to hand stack_type if only the desired frame type should be present this deep in the mix?
        elif len(files) == 0 & len(directories) == 0:
            print('No', stack_type, '  detected.')
        else:
            print('Mix directory ', stack_type, ' organization format detected.')
            sub_files  = self._read_stack(files, stack_type)
            for d in directories:
                # should not be any further sub folders
                f, _ = self.diread(d)
                sub_files += self._read_stack(f, stack_type)
        return sub_files
            
    def _read_stack(self, files, stack_type):
        if (stack_type == 'lights') & (self.unique_lights):
            stackstr = self.stack_types['bias'] + self.stack_types['flats'] # NOTE:: in the future this should be modded to include all non light stack types.
        else:
            try:
                stackstr = self.stack_types[stack_type]
            except:
                stackstr = stack_type
            
        stack_list = []
        for s in files:
            good = True # janky way of making sure we dont double up when a few str dont match
            for st in stackstr:
                if (stack_type == 'lights') & (self.unique_lights):
                    if (str(st)) in str(s.name):
                        good = False
                        
                else:
                    if (str(st)) in str(s.name):
                        stack_list.append(s)
            if (stack_type == 'lights') & (self.unique_lights):
                if good:
                    stack_list.append(s)
        stack_files = []
        if len(stack_list) > 50:
            n = 0 # this is a very janky way of writing this, enjoy
            for f in (pbar := tqdm(stack_list, colour = 'green')):
                n += 1
                pbar.set_description('Reading ' + stack_type + ' : ' + str(n))
                pbar.refresh()
                hdu = CCDDatax.read(f.path)
                stack_files.append(hdu)
        elif len(stack_list) > 0:
            for f in stack_list:
                hdu = CCDDatax.read(f.path)
                stack_files.append(hdu)
        else:
            raise Exception('Error locating files: ', stack_type, ' not found.')
        return stack_files

    def dirscan(self, dirarray):
        '''
        Dirscan takes a filepath array and scans the resulting directory for Bias, Flats, and Lights
        based on a simple filename pattern. Data can be either all in this single directory or in appropriately 
        named sub directories(i.e. '/bias', '/flats', and '/lights'). 
        
        **This currently only supports single filter data** 
        
        Parameters
        ----------
        dirarray: str array
            array containing the path to the desired directory.
        '''
        # TODO handle multi filter searching 
        path = Dorado.dordir
        for dir in dirarray:
            path = path / dir
        files, directories = self.diread(path)

        # TODO check if files are actually fits files before we break the fits reader
        if len(directories) == 0:
            if len(files) == 0:
                raise Exception('No viable data found')

            else:
                print('Single directory level organization format detected.')
                return self._single_level(files)


        elif len(files) == 0:
            print('Multi directory level organization format detected.')
            # People need to adopt a standard way of saving stuff or I will go nuts trying to come up with all the different ways they can organize(or disorganize) their data.
            return self._multi_level_top(files, directories) # what if there is usable data in single level and we only search multi-level
        
    def _chkPhi(self, stack_files):
        # sequence through headers saving every filter keyword
        phi_list = []
        phi_stacks = {}
        for im in stack_files:
            phi = im.header['FILTER'] # hardcoded filter keyword for now
            if phi in phi_list: # check if filter is already accounted for 
                phi_stacks[phi].append(im)
            else:
                phi_list.append(phi)
                phi_stacks[phi] = [im]
        for p in phi_list:
            print('Found ', len(phi_stacks[p]), ' ', p, ' frames.')
        # Should we go ahead and construct stacks first before returning?
        return phi_list, phi_stacks

    def _chkCal(self):
        pass
    
    def mkceres(self,  date, sub = 'raw', target = None, calibrated = False, aligned = False):
            '''
            mkceres creates a ceres object from a given datestring for an observation and an optional 
            dorado.target instance. mkceres can be pointed to calibrated or aligned data via the
            calibrated and aligned boolean flags (Note, in this case 'sub' should be set to 'wrk').
            
            Parameters
            ----------
            date: date string
                Instance of dorado.stack class to add to self.
            sub: string
                pointer string denoting subfolder to '$user/.dorado/data' (dorado directory/data/).
                Default is 'raw'
            target: dorado.target
                dorado.target instance to use as the target for the created dorado.ceres
                object. Default is None

            calibrated: boolean
                sets whether the imported images are calibrated and within an 'calibrated'
                folder. Default is False.

            aligned: boolean
                sets whether the imported images are aligned and within an 'aligned'
                folder. Default is False.
            '''
            
            # TODO needs ability to handle multifilter directories
            # TODO auto handle setting 'wrk' as subfolder
            # TODO ignore looking for calibration files if calibrated or aligned 
            if aligned:
                dirarray = ['data', sub, date, 'aligned']
            elif calibrated:
                dirarray = ['data', sub, date, 'calibrated']
            else:
                dirarray = ['data', sub, date]

            biases, flats, lights = self.dirscan(dirarray)

            print(len(flats),  ' flats found.')
            print(len(biases), ' bias frames found.')
            print(len(lights), ' lights found')
            
            # lets split the flats into filters
            flat = {}
            if len(flats) > 0:
                phi_flat_list, phi_flat_stacks = self._chkPhi(flats)
                print('Flat filters found : ', phi_flat_list)
                for f in phi_flat_list:
                    print('Making ', f, ' flat.')
                    flat[f] = self.mkFlat(phi_flat_stacks[f])
            else:
                phi_flat_list = [] # set empty array to check if lights filters are found in flat filters.
            # Make the  Base Bias if possile
            if len(biases) > 0:
                bias = self.mkBias(biases)
                
            if len(lights) > 0:
                epoch = lights[0].header['DATE-OBS']
                datestr = getDateString(epoch) # TODO 
                mjdstr  = str(int(np.round(Time(epoch, format='fits').mjd, 0)))
            else:
                raise Exception("No usable lights found") # NOTE user may want to run calibration files without lights, or have lights without calibration files

            # Get Base Bias if none in raw
            if len(biases) == 0:
                bias = self.getBias(mjdstr) # find bias that are near in time
                
            # lets split the lights into filters
            phi_list, phi_stacks = self._chkPhi(lights)
            print('Light filters found : ', phi_list)
        
            for phi in phi_list: # This is not a very stable section, needs some love later to handle niche cases
                # find missing flats near in time to data
                if phi not in phi_flat_list:
                    # print('Retrieving ', phi, ' flat.')
                    phi_flat = self.getFlat(mjdstr, phi)
                    flat[phi] = phi_flat
            
            # Construct cere class
            cere = Ceres(time = Time(epoch, format='fits'), datestr=datestr)
            cere.bias = bias
                    
            # Construct stacks and hand them to Ceres
            for phi in phi_list:
                cere.add_stack(Stack(phi_stacks[phi], flat = flat[phi],  calibrated = calibrated, aligned = aligned, target = target))
                
            return cere
        
    def mkFlat(self, flats):
            """
            mkFlat takes  a list of flats to construct a calibrated flatfield image.
            
            Parameters
            ----------
            flats: array[CCDdata]
                    array of raw flatfields. ------------> this needs to be corrected for the new image storage format

            Returns
            -------
            flat: CCDdata
                    The combined calibrated flatfield image.
            """
            date = Time(flats[0].header['DATE-OBS'], format='fits').mjd
            filt = flats[0].header['filter']
            ## TODO :: standardize filter names and include telescope profiles
            fname = str(int(date)) + '_' + str(filt) + '_flat.fits'
            flatdir = Dorado.dordir / 'data' / 'flats' 
            contents = os.scandir(path = flatdir)
            save = True
            for entry in contents:
                if fname in entry.name:
                    save = False
                    flat = CCDDatax.read(flatdir / fname) #, unit = Dorado.unit) ## NOTE edited
            # TODO Allow RGB data
            
            
            if save:
                import timeit
                start = timeit.default_timer()
                # from datetime import datetime
                # start=datetime.now()
                flat = ccdprocx.combine(flats, method='median', mem_limit=2e9, sigma_clip=True, dtype=np.int16)
                stop = timeit.default_timer()
                print('Time to stack ', filt, ' flats: ', stop - start)
                # print datetime.now()-start
                # c = ccdprocx.Combiner(flats)
                # c.sigma_clipping()
                # flat = c.median_combine()
                # , method = 'average',
                #                     sigma_clip = True, sigma_clip_low_thresh = 5, sigma_clip_high_thresh = 5,
                #                     sigma_clip_func = np.ma.median, sigma_clip_dev_func = mad_std, unit = self.unit)
                flat.header['stacked'] = True
                flat.header['numsubs'] = len(flats)
                flat.header['DATE-OBS'] = flats[0].header['DATE-OBS']
                flat.header['filter'] = flats[0].header['filter']
                ## TODO :: There is probably more missing keywords in the combined header, where's Waldo...

                
                flat.data = flat.data.astype('uint16') 

                print('Saving Flat for later use')
                flat.write(flatdir / fname)
            else:
                print('Flat for filter and date already saved.')

            return flat
        
    def mkBias(self, biasIFC):
            """
            mkBias takes a list of bias images to construct 
            a combined bias image. ------------> this needs to be corrected for the new image storage format
            
            Parameters
            ----------
            biasIFC: array[CCDdata]
                    array of raw bias images.

            Returns
            -------
            bias: CCDdata
                    The combined bias image.
            """
            # Allow specification of median or mean
            # Allow RGB data

            date = Time(biasIFC[0].header['DATE-OBS'], format='fits').mjd
            fname = str(int(date)) + '_Bias.fits'
            biasdir = Dorado.dordir / 'data' / 'bias' 
            contents = os.scandir(path = biasdir)
            save = True
            for entry in contents:
                if fname in entry.name:
                    save = False
                    bias = CCDDatax.read(biasdir / fname) #, unit = Dorado.unit) ## NOTE edited
            if save:
                bias = ccdprocx.combine(biasIFC, method = 'average', mem_limit=2e9, dtype=np.int16) #, unit = Dorado.unit) ## NOTE edited
                bias.meta['stacked'] = True
                bias.header['numsubs'] = len(biasIFC)
                # date = Time(bias.header['DATE-OBS'], format='fits').mjd
                bias.data = bias.data.astype('uint16')
                print('Saving Bias for later use')
                bias.write(biasdir / fname)
            else:
                print('Bias for date already saved.')
            return bias
    
    def savewrk(self, cr, filters = None):
        """
        savewrk takes a ceres key and a list of filters to save and saves the desired data
        to the dorado working directory '$user/.dorado/data/wrk/'. 
        
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to save.
        filter: str
            String representation of the relevent filter. Default is all filters in ceres instance.

        Returns
        -------
        
        """
        # TODO mod fplate to accept cr name
        if Dorado.ceres[Dorado.ceres_keys[cr]].datestr == None:
            Dorado.getDateString(cr) 
        # wrkdir = Dorado.dordir / 'data' / 'wrk'

        datestr = Dorado.ceres[Dorado.ceres_keys[cr]].datestr
        print('Save to data/wrk/', datestr)
        Dorado.mkwrk(cr)

        if filters == None:
            filters = Dorado.ceres[Dorado.ceres_keys[cr]].filters.keys()
        
        if type(filters) == type([]):
            for phi in filters:
                self._save_phi(cr, phi)
        elif type(filters) == str:
            self._save_phi(cr, filters)
        else:
            raise Exception('Cannot parse filter arguement. Did you enter the string of a single filter or a list of filter strings?')
        
    def _save_phi(self, cr, phi):
        wrkdir  = Dorado.dordir / 'data' / 'wrk'
        datestr = Dorado.ceres[Dorado.ceres_keys[cr]].datestr
        fildat  = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[phi]]
        
        if (fildat.target == None):
            fplate = str(int(Dorado.ceres[Dorado.ceres_keys[cr]].date.mjd)) + '-' + phi + '_'
        else:
            fplate = str(fildat.target.name) + '-' + str(int(Dorado.ceres[Dorado.ceres_keys[cr]].date.mjd)) + '-' + phi + '_' 
            
        if (fildat.calibrated == True) and (fildat.aligned == True): 
            wrdir = wrkdir / datestr / 'aligned'
            os.makedirs(wrkdir / datestr / 'aligned' / phi, exist_ok = True)
            wrdir = wrdir / phi
            fsub = '_ca'
        elif (fildat.aligned == True):
            wrdir = wrkdir / datestr / 'aligned'
            os.makedirs(wrkdir / datestr / 'aligned' / phi, exist_ok = True)
            wrdir = wrdir / phi
            fsub = '_a'
        elif (fildat.calibrated == True):
            wrdir = wrkdir / datestr / 'calibrated'
            os.makedirs(wrkdir / datestr / 'calibrated' / phi, exist_ok = True)
            wrdir = wrdir / phi
            fsub = '_c'
        else: 
            wrdir = wrkdir / datestr / 'uncalibrated'
            os.makedirs(wrkdir / datestr / 'uncalibrated' / phi, exist_ok = True)
            wrdir = wrdir / phi
            fsub = ''
        
        if fildat.base != None:
            fname = fplate + '_base.fits'
            fildat.base.write(wrkdir / datestr / fname, overwrite = True)
            
        if fildat.solved != None: 
            fname = fplate + '_solved.fits'
            fildat.solved.write(wrkdir / datestr / fname, overwrite = True)
            
        for p in range(len(fildat.data)):
            image = fildat.data[p]
            tstr = str(Time(image.header['DATE-OBS'], format='fits', out_subfmt='date_hms').value)
            image.data = image.data.astype('uint16') 
            # image.mask = image.mask.astype('uint16') 
            image.mask = None
            image.uncertainty = None
            # image.uncertainty = image.uncertainty.astype('uint16') 
            fname = fplate + str(p) + '_' + tstr + fsub + '.fits'
            image.write(wrdir / fname, overwrite = True)

    def getBias(self, mjdstr):
        # NOTE This function does not account for binning/image size
        target_mjd = int(mjdstr)
        biasdir = Dorado.dordir / 'data' / 'bias' 
        contents = os.scandir(path = biasdir)
        candidates = []
        # print(target_mjd)
        for entry in contents:
            entry_mjd = int(entry.name[0:5]) # grab first five characters of filename. Should be the mjd number
            # print(entry_mjd)
            if np.abs(target_mjd - entry_mjd) <= self.calibration_window:
                candidates.append(entry_mjd)
        if len(candidates) == 0:
            print('No viable bias frames found. Try increasing your calibration window.')
            bias = None
        else:
            mjd_to_use = candidates[np.argmin(np.array(candidates) - entry_mjd)] #locate the mjd of the closest file in time
            print('Retrieved bias from mjd: ', mjd_to_use)
            fname = str(mjd_to_use) + '_Bias.fits'
            bias = CCDDatax.read(biasdir / fname) #, unit = Dorado.unit) ## NOTE edited
        return bias
    
    def getFlat(self, mjdstr, phi):
        target_mjd = int(mjdstr)
        flatdir = Dorado.dordir / 'data' / 'flats' 
        contents = os.scandir(path = flatdir)
        candidates = []
        # print(target_mjd)
        for entry in contents:
            if phi in entry.name:
                entry_mjd = int(entry.name[0:5]) # grab first five characters of filename. Should be the mjd number
                # print(entry_mjd)
                if np.abs(target_mjd - entry_mjd) <= self.calibration_window:
                    candidates.append(entry_mjd)
        if len(candidates) == 0:
            print('No viable flat frames found for ', phi, '. Try increasing your calibration window or using a different filter.')
            flat = None
        else:
            mjd_to_use = candidates[np.argmin(np.array(candidates) - entry_mjd)]
            print('Retrieving ', phi, ' flat from mjd: ', mjd_to_use)
            fname = str(mjd_to_use) + '_' + str(phi) + '_flat.fits'
            flat = CCDDatax.read(flatdir / fname) #, unit = Dorado.unit) ## NOTE edited
        return flat

# class tess_reader:
#     def __init__(self):
#         Dorado.init_dir(tess = True)
#         # check if a name query is needed to retrieve data
#         # initialize 

#     def mkceres(self, tid, sub = 'raw', target = None):
#         fits_file = self.dirscan()
#         tpf = lk.open(fits_file)

#         with fits.open(fits_file, mode="readonly") as hdulist:
#             tess_bjds = hdulist[1].data['TIME'] # Should I use 'TIMECORR'?
#             tess_bjds_corr = hdulist[1].data['TIMECORR'] 
#             raw_counts = hdulist[1].data['RAW_CNTS']
#             calibrated_fluxes = hdulist[1].data['FLUX'] # Should I use raw or calibrated?
#             flux_err = hdulist[1].data['FLUX_ERR']
#             tid = hdulist[0].header['TICID']
#             TSTART = hdulist[0].header['TSTART']
#             TSTOP = hdulist[0].header['TSTOP']
#             tidstr = 'TIC ' + str(tid)
        
#         cere = Ceres()
#         # need to make a stack, and reconfigure the ceres class to be more generic in data storage. maybe make the dorphot function
#         # a wrapper, or make the ability to produce a target pixel file of aico data.
#         return cere


#     def dirscan(self, Dorado, dirarray):
#         path = Dorado.dordir
#         for dir in dirarray:
#             path = path / dir

#         files, directories = Dorado.diread(path)

#         if len(files) == 0:
#             print('No files found...')
#         elif len(files) > 1:
#             for file in files:
#                 if '_tp.fits' in file:
#                     tp_file = file    
#         elif len(files) == 1:
#             tp_file = files[0]

#         return tp_file

