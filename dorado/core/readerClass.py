import warnings
warnings.filterwarnings('ignore')

from ..ceres import Ceres
from ..stack import Stack
from ..core.coreClass import *

import ccdprocx

from astropy.nddata.ccddata import CCDData
CCDData._config_ccd_requires_unit = False
from astropy.time import Time
from astropy.io import fits

import lightkurve as lk
from lightkurve.lightcurve import TessLightCurve as tlc

import os
from tqdm import tqdm

__all__ = ['aico_reader'] #, 'tess_reader'

# TODO Make it so that each class isnt loaded unless chosen (core class takes string and then imports correct class)
class aico_reader:
    '''
    
    '''
    def __init__(self):
        self.temp = None

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
        files, directories = Dorado.diread(path)

        biasstr = ['Bias', 'Bias', 'bias', 'BIAS']

        flatsstr = ['FLAT', 'FlatField', 'flat', 'Flat', 'Flats', 'flats', 'FLATS', 'FlatFields']
        
        lightsstr = ['lights', 'Lights', 'LIGHTS']

        if len(directories) == 0:
            if len(files) == 0:
                raise Exception('No viable data found')

            else:
                print('Single directory level organization format detected.')
                # print('Reading files.')
                # compile these into master frame and pass them to Filer
                # path = self.dordir / 'data' / 'raw' / date
                biasl = []
                for strbias in biasstr:
                    for s in files:
                        if (str(strbias)) in str(s.name):
                            biasl.append(s)
                bias = []
                for i in range(len(biasl)):
                    hdu = CCDData.read(biasl[i].path) #, unit = Dorado.unit) ## NOTE edited
                    bias.append(hdu)
                # print('Bias searched.')
                flatsl = []
                for strflat in flatsstr:
                    for s in files:
                        if strflat in s.name:
                            flatsl.append(s)
                flats = []
                for i in flatsl:
                    hdu = CCDData.read(i.path) #, unit = Dorado.unit) ## NOTE edited
                    flats.append(hdu)
                # print('flats searched.')
                # strip into ceres (check if multi-filter)
                lightsl = []
                nonlight = []
                for b in biasl:
                    nonlight.append(b.name)
                for f in flatsl:
                    nonlight.append(f.name)

                for s in files:
                    if (str(s.name) not in nonlight):
                        lightsl.append(s)

                # lightsl = [s for s in files if (s.name not in biasl) and (s.name not in flatsl)]
                lights = []
                print('Reading lights...')
                for i in tqdm(lightsl, colour = 'green'):
                    hdu = CCDData.read(i.path) #, unit = Dorado.unit) ## NOTE edited
                    ## Trying to enforce 16bit data instead of 32 or 64
                    hdu.data = hdu.data.astype('uint16')
                    lights.append(hdu)
                # print('lights searched.')
                return bias, flats, lights


        elif len(files) == 0:
            print('Multi directory level organization format detected.')
            # print('Reading directories.')
            # read into this and compile into master bias, pass to Filer
            biasdir = [s for s in directories if s.name in biasstr]
            # check if multifilter, compile, pass to Filer
            flatsdir = [s for s in directories if s.name in flatsstr]
            # check if multifilter (or subdirectories) and pass to ceres
            lightsdir = [s for s in directories if (s.name not in flatsstr) and (s.name not in biasstr)]
            biasl, _ = Dorado.diread(biasdir[0])
            bias = []
            for i in biasl:
                hdu = CCDData.read(i.path) #, unit = Dorado.unit) ## NOTE edited
                bias.append(hdu)


            for ldir in lightsdir:
                files, directories = Dorado.diread(ldir)
                if len(directories) == 0:
                    if len(files) == 0:
                        raise Exception('No viable light data found')
                    else:
                        print('Single directory lights organization format detected.')
                        # print('Reading files.')
                        lights = []
                        print('Reading lights...')
                        for i in tqdm(files, colour = 'green'):
                            hdu = CCDData.read(i.path) #, unit = Dorado.unit) ## NOTE edited
                            ## Trying to enforce 16bit data instead of 32 or 64
                            hdu.data = hdu.data.astype('uint16')
                            lights.append(hdu)
                        # filter = ImageFileCollection(ldir).values('filter', unique = True)
                        # lights = [filter, lightsarr]

                        
                elif len(files) == 0:
                    print('Multi directory lights organization format detected.')
                    # print('Reading directories.')
                    lights = []



            for fdir in flatsdir:
                files, directories = Dorado.diread(fdir)
                if len(directories) == 0:
                    if len(files) == 0:
                        raise Exception('No viable light data found')
                    else:
                        print('Single directory flats organization format detected.')
                        # print('Reading files.')
                        flats = []
                        for i in files:
                            hdu = CCDData.read(i.path) #, unit = Dorado.unit) ## NOTE edited
                            flats.append(hdu)

                elif len(files) == 0:
                    print('Multi directory flats organization format detected.')
                    # print('Reading directories.')
                    flats = []

            return bias, flats, lights
        
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
            if aligned:
                dirarray = ['data', sub, date, 'aligned']
            elif calibrated:
                dirarray = ['data', sub, date, 'calibrated']
            else:
                dirarray = ['data', sub, date]

            biasIFC, flats, lights = self.dirscan(dirarray)

            print(len(flats), ' flats found.')
            print(len(biasIFC), ' bias frames found.')
            print(len(lights), ' lights found')
            
            

            # save these frames
            ## TODO :: get missing calibraation files
            ## if data aligned or calibrated, cal frames arent needed

            ## TODO :: look into UTC wrecking stuff
            if len(biasIFC) == 0:
                cere = Ceres(time = Time(lights[0].header['DATE-OBS'], format='fits'))
                ## TODO get date string is in filer
                # self.getDateString(cere)
            else:
                bias = self.mkBias(biasIFC)
                cere = Ceres(bias = bias, time = Time(lights[0].header['DATE-OBS'], format='fits'))
                # self.getDateString(cere)

            # for f in filter_data:
            # cere.flats[flat.header['filter']] = flat

            ## TODO :: multifilter fun
            if len(flats) == 0:
                cere.add_stack(Stack(lights, calibrated = calibrated, aligned = aligned, target = target))
            else:
                flat = self.mkFlat(flats)
                cere.add_stack(Stack(lights, flat = flat, calibrated = calibrated, aligned = aligned, target = target))
            

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
                    flat = CCDData.read(flatdir / fname) #, unit = Dorado.unit) ## NOTE edited
            # TODO Allow RGB data
            
            
            if save:
                c = ccdprocx.Combiner(flats)
                c.sigma_clipping()
                flat = c.median_combine()
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
                    bias = CCDData.read(biasdir / fname) #, unit = Dorado.unit) ## NOTE edited
            if save:
                bias = ccdprocx.combine(biasIFC, method = 'average') #, unit = Dorado.unit) ## NOTE edited
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
        wrkdir = Dorado.dordir / 'data' / 'wrk'

        datestr = Dorado.ceres[Dorado.ceres_keys[cr]].datestr
        print('Saved to data/wrk/', datestr)
        Dorado.mkwrk(cr)

        if filters == None:
            filters = Dorado.ceres[Dorado.ceres_keys[cr]].filters.keys()
        
        for filter in filters:
            fildat = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
            
            if (fildat.target == None):
                fplate = str(int(Dorado.ceres[Dorado.ceres_keys[cr]].date.mjd)) + '-' + filter + '_'
            else:
                fplate = str(fildat.target.name) + '-' + filter + '_' 


            if (fildat.calibrated == True) and (fildat.aligned == True): 
                wrdir = wrkdir / datestr / 'aligned'
                if len(filters) != 1:
                    os.makedirs(wrkdir / datestr / 'aligned' / filter, exist_ok = True)
                    wrdir = wrdir / filter
                fsub = '_ca'
            elif (fildat.calibrated == True):
                wrdir = wrkdir / datestr / 'calibrated'
                if len(filters) != 1:
                    os.makedirs(wrkdir / datestr / 'calibrated' / filter, exist_ok = True)
                    wrdir = wrdir / filter
                fsub = '_c'
            else: 
                wrdir = wrkdir / datestr / 'uncalibrated'
                if len(filters) != 1:
                    os.makedirs(wrkdir / datestr / 'uncalibrated' / filter, exist_ok = True)
                    wrdir = wrdir / filter
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
                fname = fplate + str(p) + '_' + tstr + fsub + '.fits'
                image.write(wrdir / fname, overwrite = True)



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

