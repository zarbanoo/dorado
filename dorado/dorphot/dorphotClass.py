from tqdm import tqdm
import numpy as np
import os

from ..core.coreClass import *
from ..timeseries.timeseriesClass import timeSeries
from ..graf.grafClass import star_chart
import ccdprocx

from astropy.time import Time
from astropy.nddata.ccddata import CCDData
CCDData._config_ccd_requires_unit = False
import astroalign as aa
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales as ppps


from astropy.visualization import MinMaxInterval, ZScaleInterval, SquaredStretch, SqrtStretch, AsinhStretch, LinearStretch, LogStretch, HistEqStretch, PowerStretch, SinhStretch

from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import SigmaClip, gaussian_sigma_to_fwhm
from photutils import Background2D, MedianBackground

from photutils.aperture import CircularAperture, aperture_photometry, CircularAnnulus

from astroquery.gaia import Gaia
import astropy.units as u
from astropy.table import Table
import time
from astropy.coordinates import SkyCoord
from scipy.ndimage import median_filter

import skimage.exposure as skie
from skimage.feature import blob_dog, blob_log, blob_doh
# https://scikit-image.org/docs/stable/api/skimage.feature.html#skimage.feature.blob_log
# https://scikit-image.org/docs/stable/auto_examples/features_detection/plot_blob.html
# StarSeeker is a method inspired by the above links and developed
# For the precursor to DORADO, DRACO circa 2019.

__all__ = ['aicoPhot', 'dracoPhot']

class aicoPhot:
    def __init__(self):
        # TODO make observatory class
        self.temp = None
        # list of calibration frames from disk
    def calibrate(self, cr,  filter, use_med_cr_removal = False, rb = 0, use_lac_cr_removal = False, scln = False, scln_xy = [1,-1, 1,-1]):
        '''
        calibrate performs 'pre-processing' CCDData calibration to a data stack within a ceres object. This 
        involves Bias and Flatfield correction and optional median value based cosmic ray removal based on the
        ccdprocx.cosmicray_median() method.
        
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to calibrate.
        filter: str
            String representation of the relevent filter.
        use_med_cr_removal: boolean
            controls whether to use cosmic ray removal in calibration, may affect runtime. Default is False.
        rb: int
            the rbox value for ccdprocx.cosmicray_median(). Default is 0.
        use_lac_cr_removal : boolean
            controls whether to use lacosmic ray removal in calibration, may affect runtime. Default is False
            ::NOTE:: the package ccdproc is required.
        scln: boolean
            controls whether to crop out scanlines. Default is False
        scln_xy: array
            array of position values for line cropping. [xmin, xmax, ymin, ymax]. Default is [1, -1, 1, -1]
        '''
        # for bla in series: add bias corrected = True to header
        stack = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        flat = stack.flat
        bias = Dorado.ceres[Dorado.ceres_keys[cr]].bias
        c_series = []

        print('Calibrating')
        if use_lac_cr_removal:
                try:
                    import ccdproc
                except Exception:
                    print('Failed to open ccdproc. Is it installed?')
        for im in tqdm(stack.data, colour = 'green'):
            # bar.update()
            im.data = im.data.astype('uint16') 
            flat.data = flat.data.astype('uint16') 
            bias.data = bias.data.astype('uint16') 
            im = ccdprocx.ccd_process(im, master_bias = bias, master_flat = flat)
            if scln:
                im.data = im.data[scln_xy[0]:scln_xy[1], scln_xy[2]:scln_xy[3]]
                im.mask = im.mask[scln_xy[0]:scln_xy[1], scln_xy[2]:scln_xy[3]]
                im.uncertainty = im.uncertainty[scln_xy[0]:scln_xy[1], scln_xy[2]:scln_xy[3]]
            if use_med_cr_removal:
                im = ccdprocx.cosmicray_median(im, rbox = rb)
            if use_lac_cr_removal:
                im == ccdproc.cosmicray_lacosmic(im)
            im.data = im.data.astype('uint16') 
            c_series.append(im)
        # add more flags for header
        Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].data = c_series
        Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].calibrated = True
        
    def imarith(self, cr, filter, operator, operand):
        '''
        imarith is a basic replication of the IRAF image arithmatic tool imarith.
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to perform image arithmatic on.
        filter: str
            String representation of the relevent filter.
        operator: string
            String of operator to be used for arithmatic. Supported operations are '+', '-', '/', and '*'
        '''
        # TODO should this be in stack? like have a wrapper here?
        # mod to check datatype using type()
        # mod to remove im_count and make possible to use single image NOTE might already be done
        # mod to accomodate CCDdata object NOTE might already be done
        series = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        for i in range(len(series)):
            if (operator == '+'):
                series[i].data = series[i].data  + operand
            elif (operator == '-'):
                series[i].data = series[i].data - operand
            elif (operator == '/'):
                series[i].data = series[i].data  / operand
            elif (operator == '*'):
                series[i].data = series[i].data  * operand
        
        Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]] = series
        

    def getWCS(self, cr, filter, alignto = None, cache = True):
        '''
        getWCS obtains WCS information for an image either via previously solved data in the cache or
        by passing the image to astrometryNet via dorado.plate_solve().
        
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to get WCS information for.
        filter: str
            String representation of the relevent filter.
        alignto: int
            Index of image to use as reference. Default is stack.alignto
        cache: boolean
            Controls whether to call astrometryNet for solve or use solved result stored in cache 
            from previous run. Default is True.
        '''
        # TODO mod so cache results are target specific
        series = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        if alignto == None:
            alignto = series.alignTo
        if cache:
            hdulist = fits.open(Dorado.dordir / 'cache' / 'astrometryNet' / 'solved.fits') 
            Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].wcs = WCS(hdulist[0].header, hdulist)
            Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].solved = CCDData.read(Dorado.dordir / 'cache' / 'astrometryNet' / 'solved.fits')
            hdulist.close()
        else:
            toalign = series.data[alignto]
            fname, cachedir = Dorado.mkcacheObj(toalign, 'astrometryNet')
            path = [cachedir, fname]
            writearray = [cachedir, 'solved.fits']
            solved, wcs_header = Dorado.plate_solve(path, writearray = writearray)
            Dorado.delcacheObj( fname, 'astrometryNet')
            Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].wcs = WCS(wcs_header)
            Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].solved = solved
    
    def align(self, cr, filter, alignto = None, getWCS = True, cache = False, ds = 2, ma = 5):
        '''
        align aligns a filter stack within a specified ceres object to a specified frame (default is first
        frame). align can optionally retrieve the corresponding WCS data for the aligned stack.
        
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to align.
        filter: str
            String representation of the relevent filter.
        alignto: int
            Index of image to use as reference. Default is stack.alignto
        getWCS: boolean
            Controls whether to obtain WCS information for stack.
        cache: boolean
            Controls whether to call astrometryNet for solve or use solved result stored in cache 
            from previous run. Default is True.
        ds: int or float
            Sets the detection sigma value for astroalign.register. Default is 2.
        ma: int or float
            Sets the minimum area value for astroalign.register. Default is 5.
        '''
        series = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        
        if alignto == None:
            alignto = series.alignTo
        else: 
            series.alignTo = alignto
        toalign = series.data[alignto]
        ## TODO :: make this use ceres.getWCS()
        if getWCS:
            if cache:
                toalign =  CCDData.read(Dorado.dordir / 'cache' / 'astrometryNet' / 'solved.fits') #, unit = Dorado.unit)
                hdulist = fits.open(Dorado.dordir / 'cache' / 'astrometryNet' / 'solved.fits') 
                Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].wcs = WCS(hdulist[0].header, hdulist)
                hdulist.close()
                Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].solved = toalign
            else:
                fname, cachedir = Dorado.mkcacheObj(toalign, 'astrometryNet')
                path = [cachedir, fname]
                writearray = [cachedir, 'solved.fits']
                solved, wcs_header = Dorado.plate_solve(path, writearray = writearray)
                toalign = solved
                Dorado.delcacheObj( fname, 'astrometryNet')
                Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].wcs = WCS(wcs_header)
                Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].solved = solved
                # delete cache object
                # save solved to target

        aa_series = []
        skipped = []

        print('Aligning')
        for image in tqdm(series.data, colour = 'green'):
            # bar.update()
            try:
                img, _ = aa.register(image.data, toalign.data, detection_sigma = ds, min_area = ma)
                aaim = image
                aaim.data = img
                aa_series.append(aaim)
            except:
                skipped.append(image)
                # print('Image skipped')
        if len(skipped) != 0:
            print(len(skipped), ' images skipped.')
            ## TODO :: need to redo times and such for less ims

        Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].data = aa_series
        Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].aligned = True
    
    def differential_magnitude(self, flux1, flux2, mag2):
        '''
        Probably could use target class inputs and uncertainties
        '''
        mag1 = -2.5 * np.log10(flux1/flux2) + mag2
        return mag1
        
    def apPhot(self, cr, filter, toid, control_toid = None, shape = 21, unc = 0.1):
        '''
        apPhot performs basic aperture photometry based on photutils.aperture_photometry on a target
        within stack. Target photometry can optionally be compared to a control target within the stack 
        via the differential photometry method.
        
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to perform aperture photometry on.
        filter: str
            String representation of the relevent filter.
        toid: string

        control_toid: string

        shape: int

        unc: float

        '''
        # TODO this needs a better name
        # TODO get seeing from PSF
        stack = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        # TODO if no wcs, complain alot
        w = stack.wcs
        toi = Dorado.targets[Dorado.target_keys[toid]]
        xy = w.wcs_world2pix(toi.coords.ra.deg, toi.coords.dec.deg, 1)
        ra = toi.coords.ra.deg
        dec = toi.coords.dec.deg
        pos = [(float(xy[0]), float(xy[1]))]
        aperture = CircularAperture(pos, r = shape)
        annulus_aperture = CircularAnnulus(pos, r_in = shape + 2, r_out = shape + 5)
        apers = [aperture, annulus_aperture]

        # TODO can a control be chosen? can more than one be chosen?
        if control_toid != None:
            control_toi = Dorado.targets[Dorado.target_keys[control_toid]]
            xyc = w.wcs_world2pix(control_toi.coords.ra.deg, control_toi.coords.dec.deg, 1)
            posc = [(float(xyc[0]), float(xyc[1]))]
            aperturec = CircularAperture(posc, r = shape)
            annulus_aperturec = CircularAnnulus(posc, r_in = shape + 2, r_out = shape + 5)
            apersc = [aperturec, annulus_aperturec]


        timestr = []
        exptimes = []
        ray = []
        decx = []
        x = []
        y = []
        flux = [] # NOTE this is exptime corrected
        fluxunc = []
        apsum = [] # NOTE this is raw aperture sum
        apsum_unc = []


        print('Performing photometry')
        for image in tqdm(stack.data, colour = 'green'):
            error = unc * image.data
            results = aperture_photometry(image, apers) #, error = error)
            bkg_mean = results['aperture_sum_1'] / annulus_aperture.area
            bkg_sum = bkg_mean * aperture.area
            results['flux_fit'] = results['aperture_sum_0'] - bkg_sum
            
            timestr.append(image.header['DATE-OBS']) # TODO This is most likely needed, but verify 
            exptimes.append(image.header['EXPTIME']) # TODO is this also needed
            ray.append(ra) # TODO is this really needed?
            decx.append(dec)
            x.append(results['xcenter'][0]) # TODO is this needed either?
            y.append(results['ycenter'][0])

            if control_toid != None:
                resultsc = aperture_photometry(image, apersc) # , error = error)
                bkg_meanc = resultsc['aperture_sum_1'] / annulus_aperturec.area
                bkg_sumc = bkg_meanc * aperturec.area
                resultsc['flux_fit'] = resultsc['aperture_sum_0'] - bkg_sumc

                apsum.append(results['flux_fit'][0] - resultsc['flux_fit'][0])
                flux.append((results['flux_fit'][0] - resultsc['flux_fit'][0])/image.header['EXPTIME'])
            else:
                apsum.append(results['flux_fit'][0])
                flux.append(results['flux_fit'][0]/image.header['EXPTIME'])

            fluxunc.append(1) ## TODO:: modify this to account for exposure time and control
            apsum_unc.append(1)
        times = Time(timestr)
        ts = timeSeries(times = times, flux = flux, exptimes = exptimes, x = x, y = y, 
        ra = ray, dec = decx, flux_unc = fluxunc, apsum = apsum, apsum_unc = apsum_unc)

        Dorado.targets[Dorado.target_keys[toid]].filters[filter] = len(toi.ts)
        Dorado.targets[Dorado.target_keys[toid]].ts.append(ts) 
        # TODO accomodate targets embedded in core (list of targets)
        # TODO the name for this function needs updating
        
    def mkBase(self, cr, filter, sigClip = False, minmax = False):
        '''
        mkBase stacks filter data within a ceres object to produce a base or 'average' stacked image.
        
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to create a base stacked image for.
        filter: str
            String representation of the relevent filter.
        sigClip: boolean

        minmax: boolean

        '''
        ## TODO :: add the option to change the combination method. Right now default is 
        # sigma clipped median combination.
        series = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        # toalign = series.data[series.alignTo]

        c = ccdprocx.Combiner(series.data)
        if minmax:
            c.minmax_clipping(min_clip = 0.1)
        if sigClip:
            c.sigma_clipping()
        base = c.median_combine()
        base.header['stacked'] = True
        base.header['numsubs'] = len(series.data)
        base.header['DATE-OBS'] = series.data[0].header['DATE-OBS']
        base.header['filter'] = series.filter

        Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].base = base
        ## TODO :: sort out what is in the header of this base file.

        
    def calBase(self, cr, filter):
        '''
        calBase calibrates a base image for a stack by computing a 2D background for the image
        and removing it from the base image.
        
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to perform base image calibration on.
        filter: str
            String representation of the relevent filter.
        '''
        # TODO this needs hella optimization and direction
        # I am 99% sure this is to remove the background gradient from the stacked base image
        img = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].base
        norm = ImageNormalize(stretch=SqrtStretch())
        sigma_clip = SigmaClip(sigma=3.)
        bkg_estimator = MedianBackground()
        bkg = Background2D(img, (50, 50), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        base = img
        base.data = img.data / bkg.background.value
        base.data[np.isnan(base.data)] = 0
        Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].base = base



class dracoPhot:
    def __init__(self,  limit_Mag = 16, search_bounds = [30, 30]):
        #TODO make observatory class
        self.limit_Mag = limit_Mag
        self.search_bounds = search_bounds
        self.scale = 'mm'
        self.transform = 'sqrt'
        
    def getWCS(self, cr, filter, alignto = None, cache = True):
        '''
        getWCS obtains WCS information for an image either via previously solved data in the cache or
        by passing the image to astrometryNet via dorado.plate_solve().
        
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to get WCS information for.
        filter: str
            String representation of the relevent filter.
        alignto: int
            Index of image to use as reference. Default is stack.alignto
        cache: boolean
            Controls whether to call astrometryNet for solve or use solved result stored in cache 
            from previous run. Default is True.
        '''
        # TODO mod so cache results are target specific
        series = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        if alignto == None:
            alignto = series.alignTo
        
        if cache:
            hdulist = fits.open(Dorado.dordir / 'cache' / 'astrometryNet' / 'solved.fits') 
            Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].wcs = WCS(hdulist[0].header, hdulist)
            Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].solved = CCDData.read(Dorado.dordir / 'cache' / 'astrometryNet' / 'solved.fits')
            hdulist.close()
        else:
            if alignto == 'base':
                toalign = series.base
            else:
                toalign = series.data[alignto]
            fname, cachedir = Dorado.mkcacheObj(toalign, 'astrometryNet')
            path = [cachedir, fname]
            writearray = [cachedir, 'solved.fits']
            solved, wcs_header = Dorado.plate_solve(path, writearray = writearray)
            Dorado.delcacheObj( fname, 'astrometryNet')
            Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].wcs = WCS(wcs_header)
            Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].solved = solved
    
    def imarith(self, cr, filter, operator, operand):
        '''
        imarith is a basic replication of the IRAF image arithmatic tool imarith.
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to perform image arithmatic on.
        filter: str
            String representation of the relevent filter.
        operator: string
            String of operator to be used for arithmatic. Supported operations are '+', '-', '/', and '*'
        '''
        # TODO should this be in stack? like have a wrapper here?
        # mod to check datatype using type()
        # mod to remove im_count and make possible to use single image NOTE might already be done
        # mod to accomodate CCDdata object NOTE might already be done
        series = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        for i in range(len(series)):
            if (operator == '+'):
                series[i].data = series[i].data  + operand
            elif (operator == '-'):
                series[i].data = series[i].data - operand
            elif (operator == '/'):
                series[i].data = series[i].data  / operand
            elif (operator == '*'):
                series[i].data = series[i].data  * operand
        
        Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]] = series
    
    def differential_magnitude(self, flux1, flux2, mag2):
        '''
        Probably could use target class inputs and uncertainties
        '''
        mag1 = -2.5 * np.log10(flux1/flux2) + mag2
        return mag1
        
    def apPhot(self, cr, filter, toid, search = True, read_date = None) :
        '''
        apPhot performs basic aperture photometry based on photutils.aperture_photometry on a target
        within stack. Target photometry can optionally be compared to a control target within the stack 
        via the differential photometry method.
        
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to perform aperture photometry on.
        filter: str
            String representation of the relevent filter.
        toid: string


        '''
        # TODO this needs a better name
        # TODO get seeing from PSF
        stack = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        # TODO if no wcs, complain alot
        w = stack.wcs
        self.projectdir = Dorado.dordir / 'data' / 'projects' / toid 
        self.projectdatedir = Dorado.dordir / 'data' / 'projects' / toid / Dorado.ceres[Dorado.ceres_keys[cr]].datestr
        os.makedirs(self.projectdatedir, exist_ok = True)
        out_filename_prefix  = toid + '_'
        self.get_stars(cr, filter, toid, search, read_date)
        print('Initial FWHM: ', np.mean(self.stars['FWHM']), '+/-', np.std(self.stars['FWHM']), 'px')
        sxy =  w.proj_plane_pixel_scales() #ppps(w)
        print('Pixel scales: ', sxy)
        print('Initial FWHM: ', (np.mean(self.stars['FWHM']) * sxy[0]).to(u.arcsec), '+/-', (np.std(self.stars['FWHM']) * sxy[0]).to(u.arcsec), 'px')
        print('Performing Photometry...')
        # The run table is most likely superseeded by the log table 
        # side note, the log table is less of a log and more of results
        # since it doesnt log the procedure, and instead contains results
        run = Table(names = ('time', '', 'sky',))
        # TODO:: Maybe add instrument temp, stuff like airmass, alt/az, airtemp, other params
        self.log = Table(names = ('time', 'image', 'exptime', 'zp_m', 'zp_b', 'sky', 'FWHM', 'FWHM_std', 'seeing'))
        
        for i in tqdm(range(len(stack.data)), colour = 'green'):
            im = stack.data[i]
            tstr = str(Time(im.header['DATE-OBS'], format='fits').mjd)
            imname = tstr + '_' + str(i)
            imPhot = photo(im, self.stars, w, im_index = i)
            imPhot.apPhot_step()
            # Should the PSF photometry step go here?
            zpv = imPhot.get_zero_point()
            imPhot.mag_calibrate()
            imPhot.set_time()
            outstr = out_filename_prefix + imname + '.fits'
            self.log.add_row(imPhot.write(self.projectdatedir / outstr))
        self.log.write(self.projectdatedir / (out_filename_prefix + 'log.fits'), overwrite = True)
        print('Photometry completed.')
        # TODO:: write out a summary table to terminal??

    def inIm(self, tab, cr, filter, border = 0):
        stack = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        height, width = stack.data[stack.alignTo].data.shape
        num = len(tab)
        tab  = tab[tab['y']     >=      0 + border]
        tab  = tab[tab['y']     <= height - border]
        tab  = tab[tab['x']     >=      0 + border]
        tab  = tab[tab['x']     <=  width - border]

        print('Went from ', num , ' stars in input, to ', len(tab), ' stars in field area.')
        return tab

    def get_field(self, cr, filter, toid):
        startTime = time.time()
        stack = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        w = stack.wcs
        im = stack.data[stack.alignTo]
        width = u.Quantity(self.search_bounds[0], u.arcmin)
        height = u.Quantity(self.search_bounds[1], u.arcmin)
        coords = SkyCoord.from_name(toid) 
        # Make the Query
        # Columns we want to keep
        columnse = [ 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'teff_val', 'dist']
        # and what we want to name those columns
        columns = ['g_mag', 'bp_mag', 'rp_mag', 'teff', 'dist']
        # Result return limit, default was 50
        Gaia.ROW_LIMIT = 10000
        re = Gaia.query_object(coordinate=coords, width=width, height=height) 
        # Parse the results
        de = [str(d) for d in re['DESIGNATION']] # This fixes the bad return from Gaia
        # Construct the star object table
        r = Table((de, re['ra'], re['dec']), names = ('DESIGNATION', 'ra', 'dec'))
        # Fill in rest of columns
        for (s, se )in zip(columns, columnse):
            r[s] = re[se]
        # Limit the results to a set magnitude
        r = r[r['rp_mag'] <= self.limit_Mag]
        # Get stellar object pixel position
        r['x'], r['y'] = w.world_to_pixel(SkyCoord(r['ra'], r['dec']))
        # This is a hardcoded scale size for plotting
        r['r'] = 28 - r['rp_mag'].value # Do we keep it?
        # Figure out which results are in the image
        r = self.inIm(r, cr, filter, 100)
        # Round of the sig figs 
        r.round(4)
        # output runtime
        executionTime = np.round((time.time() - startTime), 3)
        print('Execution time in seconds for Gaia lookup: ' + str(executionTime))
        self.star_chart(r, im, 'Gaia', w)
        return r

    def starSeeker(self, cr, filter):
        if self.scale == 'mm':
            interval = MinMaxInterval()
        elif self.scale == 'z':
            interval = ZScaleInterval()
        if self.transform == 'sqrt':
            stretch = SqrtStretch()
        elif self.transform == 'squared':
            stretch = SquaredStretch()
        elif self.transform == 'asinh':
            stretch = AsinhStretch()
        elif self.transform == 'sinh':
            stretch = SinhStretch()
        elif self.transform == 'linear':
            stretch = LinearStretch()
        elif self.transform == 'log':
            stretch = LogStretch()
        elif self.transform == 'power':
            stretch = PowerStretch()
        elif self.transform == 'hist':
            stretch = HistEqStretch()
        startTime = time.time()
        stack = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        w = stack.wcs
        im = stack.data[stack.alignTo]
        data = im.data
        # mf = median_filter(data, size= 5)
        # datamf = data - mf
        # limg = np.arcsinh(datamf)
        norm = ImageNormalize(data, interval=interval, stretch=stretch)
        limg = norm(data)
        limg = limg / limg.max()
        low = np.percentile(limg, 1)
        high = np.percentile(limg, 99.5)
        opt_img  = skie.exposure.rescale_intensity(limg, in_range=(low,high))
        # Laplacian blob detection
        stars =  blob_log(opt_img, max_sigma=25, min_sigma = 5, num_sigma=10, threshold=.2)
        # full width at half max should be calculated from original unscaled image
        # how can we do this per image?
        fwhm = gaussian_sigma_to_fwhm * stars[:,2]
        # Convert from sigma to radii in the 3rd column.
        stars[:, 2] = stars[:, 2] * np.sqrt(2)
        y, x, r= stars[:, 0], stars[:, 1], stars[:, 2]
        results = Table((x, y, r, fwhm), names = ('x', 'y', 'r', 'FWHM'))
        co = w.pixel_to_world(stars[:,1], stars[:,0])
        results['ra']  = co.ra
        results['dec'] = co.dec
        results = self.inIm(results, cr, filter, 100)
        print('Found ', len(results), ' stars via starseeker.')
        executionTime = np.round((time.time() - startTime), 3)
        print('Execution time in seconds for starSeeker: ' + str(executionTime))
        img = im.copy()
        img.data = opt_img
        self.star_chart(results, img, 'Starseeker', w)
        return results, img
    
    def get_stars(self, cr, filter, toid, search, read_date):
        # add these back if needed :: , limit_Mag = 16, search_bounds = [30, 20]
        # ask if stars should be saved with  flag
        if search:
            gaia_stars = self.get_field(cr, filter, toid)
            s3, opt_img = self.starSeeker(cr, filter)
            gaia_stars['detection_separation'] = np.zeros(len(gaia_stars))
            gaia_stars['detection_x'] = np.zeros(len(gaia_stars))
            gaia_stars['detection_y'] = np.zeros(len(gaia_stars))
            gaia_stars['detection_r'] = np.zeros(len(gaia_stars))
            gaia_stars['FWHM'] = np.zeros(len(gaia_stars))
            matched = Table(names = gaia_stars.colnames, dtype = gaia_stars.dtype)
            print('Matching stars..')
            self.unmatched = 0
            for s in (pbar := tqdm(gaia_stars, colour = 'green')):
                pbar.set_description('Matching star : ' + str(s['DESIGNATION']))
                pbar.refresh()
                sm = self.match_star(s, s3)
                if sm != False:
                    matched.add_row(sm)
            self.stars = matched # should this really be an internal list
            # Or should it belong to anoter class
            # projectdir = Dorado.dordir / 'data' / 'projects' / toid 
            os.makedirs(self.projectdir / 'stars', exist_ok = True)
            self.stars.write(self.projectdir / 'stars' / 'stars.fits', overwrite = True)
            print(self.unmatched, ' stars were not matched.')
            # crop out unmatched stars
        else:
            if read_date != None:
                fname = 'stars_' + str(read_date) + '.fits'
                self.stars = Table.read(self.projectdir / 'stars' / fname)
            else:
                print('No read date given for saved star catalogue.')

    def match_star(self, star, s3):
        sx, sy = star['x'], star['y']
        sr = star['r']
        separation = np.sqrt((sx - s3['x'])**2 + (sy - s3['y'])**2)
        candidate = s3[separation <= sr]
        sep = separation[separation <= sr]
        if len(sep) > 1:
            print(len(sep))
            print(candidate)
            candidate = candidate[sep == np.min(sep)][0]
            sep = sep[sep == np.min(sep)][0]
            star['detection_separation'] = sep
            star['detection_x'] = candidate['x']
            star['detection_y'] = candidate['y']
            star['detection_r'] = candidate['r']
            star['FWHM']        = candidate['FWHM']
            return star
        elif len(sep) == 0:
            self.unmatched += 1
            return False
        else:
            star['detection_separation'] = sep
            star['detection_x'] = candidate['x']
            star['detection_y'] = candidate['y']
            star['detection_r'] = candidate['r']
            star['FWHM']        = candidate['FWHM']
            return star

    def star_chart(self, stars, im, toid, w):
        sc = star_chart(im, title = toid, wcs = w)
        sc.plt_stars(stars, label = 'Stars')
        sc.add_compass()
        sc.add_scale()
        sc.plot()




## TODO :: figure out how to handle aligning and processing planetary images and image with low star counts



def mag(flux, exp = 1, unc = None):
    mag_inst = -2.5 * np.log10(flux / exp)
    if unc:
        mag_inst_unc = unc / (flux * 2.30258509)
        return mag_inst, mag_inst_unc
    else:
        return mag_inst

class photo:
    '''
    Modify this to add PSF stuff
    '''
    def __init__(self, image, stars, wcs, time = None, im_index = 0):
        self.image = image
        self.im_index = im_index
        self.stars = stars
        self.wcs = wcs
        # currently not in use, grabbing from im header becaaus Im a degenerate
        self.time = time
        self.set_time()
        try:
            self.exp = image.header['EXPTIME']
        except:
            try:
                self.exp = image.header['EXPOSURE']
            except:
                print('ERROR: No image exposure length info found in fits header using default keywords.')
        # this needs a more rigorous treatment, but for now its a good aproximation.
        self.sky = np.mean((np.mean(self.image.data), np.median(self.image.data)))
        self.fwhm = 0
        self.seeing = 0
    
    def apPhot_step(self):
        self.stars['aperture_sum'] = np.zeros(len(self.stars))
        self.stars['inst_mag'] = np.zeros(len(self.stars))
        self.stars['aperture_sum_raw'] = np.zeros(len(self.stars))
        self.stars['inst_mag_raw'] = np.zeros(len(self.stars))
        for i in range(len(self.stars)):
            pos = (self.stars[i]['x'], self.stars[i]['y'])
            # TODO :: why is there no annulus? seriously, this is basic photometry
            # and I haven't even made an annulus aperture. pathetic
            shape = float(1.2 * self.stars[i]['detection_r']) # its possible this was handing back a row instead of a float
            if shape <= 0:
                print('negative or zero shape encountered', shape)
            aperture = CircularAperture(pos, r=shape)
            annulus_aperture = CircularAnnulus(positions = pos, r_in = (shape + 2), r_out = (shape + 5))
            apers = [aperture, annulus_aperture]
            phot_table = aperture_photometry(self.image.data, apers)
            bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area
            bkg_sum = bkg_mean * aperture.area
            self.stars[i]['aperture_sum_raw'] = phot_table['aperture_sum_0']
            self.stars[i]['aperture_sum'] = phot_table['aperture_sum_0'] - bkg_sum
            self.stars[i]['inst_mag_raw'] = mag(phot_table['aperture_sum_0'], self.exp)
            self.stars[i]['inst_mag'] = mag(phot_table['aperture_sum_0'] - bkg_sum, self.exp)
    
    def get_zero_point(self):
        stars = self.stars[self.stars['inst_mag'] <= -2.3] # Heres a place to mod, this cutoff value
        self.zero_point_vals = np.polyfit(stars['inst_mag'], stars['rp_mag'], 1) 
        self.zero_point = np.poly1d(self.zero_point_vals)
        return self.zero_point_vals
    
    def mag_calibrate(self):
        # Well this is elegent, its a single line of code, well, it was, then I ruined it with this comment
        self.stars['fit_mag'] = self.zero_point(self.stars['inst_mag'])
    
    def set_time(self):
        self.stars['time'] = [float(Time(self.image.header['DATE-OBS'], format='fits').mjd) for i in range(len(self.stars))]
    
    def write(self, filename):
        self.stars.write(filename, overwrite = True)
        # why am I doing this in the same action as writing the table? who knows really
        # I could pretend that it was to add the filename to the log table
        # TODO:: add filename to log table :)
        mean_fwhm, std_fwhm = np.mean(self.stars['FWHM']), np.std(self.stars['FWHM'])
        self.fwhm = mean_fwhm
        return [float(self.stars['time'][0]), self.im_index, self.exp,  self.zero_point_vals[0], self.zero_point_vals[1], self.sky, self.fwhm, std_fwhm, self.seeing]





# class tessPhot:
    # def __init__(self):
    #     self.temp = None

    # def dorphot(self, cere, filter, aperture_shape, annulus_shape):
    #     stack = cere.data[cere.filters[filter]]
        
    #     if stack.ts == []:
    #         #tess_bjds, make sure to read them in as time objects for cere
    #         stack.ts = timeSeries()
    #     ts = stack.ts

    #     perture = RectangularAperture((5,4), 3,3) 
    #     annulus_aperture = RectangularAnnulus((5,4), 4, 5, 5) 
    #     # Combine the aperture
    #     aps = [aperture , annulus_aperture]

    #     for s in tqdm(range(stack.length), colour = 'green'):
    #         # Place the aperature on the current stamp and perform aperture photometry
    #         phot_table = aperture_photometry(calibrated_fluxes[s,:,:], aps)
    #         # Compute the background values from the annulus photometry
    #         bkg_sum = (phot_table['aperture_sum_1'] / annulus_aperture.area) * aperture.area
    #         # Subtract the background from main aperture
    #         phot_table['residual_aperture_sum'] = phot_table['aperture_sum_0'] - bkg_sum
    #         # Store the aperture results 
    #         aperture_sums.append(float(phot_table['aperture_sum_0']) - float(bkg_sum))
    #         # Store the flux error
    #         err.append(np.mean(flux_err[s,:,:]))
