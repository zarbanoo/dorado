import warnings
warnings.filterwarnings('ignore')
import os
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord as acoord
import numpy as np
import astropy.units as un
from scipy.interpolate import InterpolatedUnivariateSpline, LSQUnivariateSpline

__all__ = ['Target', 'Fournax']

class Target:
    '''
        The TOI class represents an astronomical target of interest (TOI) and handles the targets relevent attributes.
        Unpassed target parameters will be gathered via astroquery (SIMBAD).


        Attributes
        ----------

        name: str
            name of target in string format.
    '''
    def __init__(self, name):
        self.name = name
        # Call Simbad for relevent data
        try:
            s = Simbad()
            r = s.query_object(self.name)
            self.coords = acoord(ra = r['RA'], dec = r['DEC'], unit = (un.hourangle, un.deg), frame = 'icrs')
        except:
            raise Exception('Error initializing Target.')
        
        self.filters = {}
        self.ts = []
        # TODO also include size angular size and other info like magnitude and type
        self.x = None
        self.y = None
        
    def calcmag(self, filter):
        '''
        calcmag converts the targets flux and associated uncertainty into an 
        instrumental magnitude and uncertainty.
        WARNING:: This function is currently not complete and no garentee is given
        on its compatability or reliability.

        Parameters
        ----------
        filter: str
            String representation of the relevent filter.

        Returns
        -------
        None

        Sets
        ----
        self.ts[filter]['mag'], self.ts[filter]['mag_unc']

        '''
        flux  = self.ts[self.filters[filter]]['flux']
        flux_unc = self.ts[self.filters[filter]]['flux_unc']
        magnitudes = -2.5 * np.log10(flux / 15)
        mag_unc = flux_unc / ((flux / 15) * 2.30258509)
        self.ts[self.filters[filter]]['mag'] = magnitudes
        self.ts[self.filters[filter]]['mag_unc'] = mag_unc

    def record(self, cr, saveType = 'fits'):
        '''
        record writes each filters timeseries to the dorado working data directory
        for the relevent date and target.

        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance for which the timeseries was derived. The save location 
            will be the working directory for this instance.
            
        saveType: str
            String representation of the file extension to save to. Default is 'fits'. See 
            'astropy.table - write()' for acceptable values.

        Returns
        -------
        None
        '''
        wrkdir = Dorado.dordir / 'data' / 'wrk'
        if Dorado.ceres[Dorado.ceres_keys[cr]].datestr == None:
            Dorado.ceres[Dorado.ceres_keys[cr]].datestr = Dorado.getDateString(cr)
        datestr = Dorado.ceres[Dorado.ceres_keys[cr]].datestr
        wrdir = wrkdir / datestr
        Dorado.mkwrk(cr)
        for fi in self.filters.keys():
            wrts = self.ts[self.filters[fi]]
            fname = str(self.name) + '_' + str(fi) + '-' + str(int(Dorado.ceres[Dorado.ceres_keys[cr]].date.mjd)) + '.' + saveType
            wrts.toTable(self.name)
            wrts.table.write(wrdir / fname, overwrite = True)
        
    def export(self, objectClass = None):
        '''
        export will record the TOI object into the dorado targets directory for future use 
        and reference. This function is not implemented yet.

        Parameters
        ----------
        objectClass: str
            Class of object to save the target under. 

        Notes
        -----
        Examples of object classes may be: 'star', 'galaxy', 'exoplanet', 'minor planet', 'satellite', 
        'white dwarf', 'nebula', 'messier_object', 'O_star', 'binary', 'globular_cluster', 'open_cluster', 
        'galaxy_cluster', 'quasar', 'AGN'. 
        
        Users can craft their own object naming schemes.
        '''
        tardir = Dorado.dordir / 'data/targets'
        if objectClass:
            os.makedirs(tardir / objectClass, exist_ok = True)
            tardir = tardir / objectClass
        



'''
Fournax is an abbreviation of Fourier numerical astronomy extension, its name is a backronym styled to match the constellation 'fornax'. 

'''

from scipy.signal import find_peaks

class Fournax(Target):

    '''
        The Fournax class extends the TOI target class to provide a consistent simple, yet robust interface to targets with regular or semi-regular photometric variability for the purposes of lightcurve/timeseries fourier analysis.  

        Fournax is an abbreviation of Fourier numerical astronomy extension, its name is a backronym styled to match the constellation 'fornax'. 

        Attributes
        ----------

        name: str
            name of target in string format.

        epoch: float or astropy.Time
            Epoch for the targets ephemeris of which the theoretical times of extrema will be extrapolated from.

        period: float
            The period between extrema of interest corresponding to the ephemeris epoch given. 

    '''
    def __init__(self, name, epoch = None, period = None): 
        ## Inherit from Ceres object (date, ts, etc.) 
        Target.__init__(self, name)

        self.freq = [] # Array of observed frequencies (Raw) --> convert to table with amplitudes
        ## NOTE:: Should frequencies be cleaned for aliasing, can the spread of aliasing peak
        #  values provide a measure of uncertainty on the fundamental frequency?

        self.Operiod = None # Observed period
        self.Operiod_unc = None # Uncertainty on period
        self.Ofreq = None # observed frequencies


        ## verify ephemeris TODO:: make into a dummy function
        if (epoch == None) and (period == None):
            print('No ephemeris data given for target')

        elif (epoch != None) and (period != None):
            self.epoch = epoch
            self.period = period # Must have units (TODO:: consider falling back to default unit if none)
            ## should these values be combined into an ephemeris dictionary to accomodate for an observational period to be determined later?

        elif (epoch != None):
            self.epoch = epoch 
            print('Period data not given.')
        elif (period != None):
            self.period = period
            print('Epoch data not given.')

        else:
            print('Unknown ephemeris error encountered. Please report via Dorado Github issues page.')
            print('Period given: ', period)
            print('Epoch given: ', epoch)
    
    def OMinusC(self, fi):
        '''
            OMinusC takes observed time(s) of max light for a repeating variable star and ephemeris data and returns O-C values as well as the corresponding cycle.



            Returns
            -------

            cycle: float, list, or array
                    The cycle corresponding to the time(s) of max light

            OmC: float, list, or array
                    The O-C value for the time(s) of max light
        '''

        ## TODO:: accomodate an array of toml values. This will shift this function from returning values to a wrapper function to setting values.
        
        cycle_ref = (self.ts[self.filters[fi]].toml - self.epoch) / self.period
        cycle = np.round(cycle_ref)
        OmC = self.ts[self.filters[fi]].toml - (self.epoch + cycle * self.period)

        self.ts[self.filters[fi]].OmC = OmC
        self.ts[self.filters[fi]].cycle = cycle

    def superfit(self, fi, terms, s):
        '''
            superfit takes a raw timeseries and performs a multistep smoothed fit of the data
            which includes a spline fit tailored with the curvature of knots from a spline
            fit of a Fourier fit. The result is then expanded and then convoluted with a
            'blackman' window function.
            Parameters
            ----------
            x, y: array
                    timeseries arrays to be fit
            terms: int
                    terms to retain in the fourier fit.
            s: int
                    new data array size 
            Returns
            -------
            X, Y: array
                    The superfit timeseries arrays.
        '''

        ## NOTE:: The name is tacky.
        ## TODO:: accomodate flux uncertainties
        ## TODO:: Zero mean of signal.

        y = self.ts[self.filters[fi]].flux
        x = self.ts[self.filters[fi]].times.value ## TODO:: convert to float friendly format like mjd

        # Fourier fit the data to model curvature
        f = np.fft.rfft(y)
        # Null or zero coefficients above ammount of series "terms"
        # This corresponds to undesired high-frequency terms
        f[terms+1:] = 0
        # Collapse back into function space, result is smoothed Fourier curve
        F = np.fft.irfft(f)
        # Create a spline fit of the fourier fit to extract knots
        tispl = InterpolatedUnivariateSpline(x[:len(F)], F, k = 5)
        # feed knots into spline of raw data
        LSQspl = LSQUnivariateSpline(x, y, tispl.get_knots()[1:-1]) 

        X = np.linspace(np.min(x), np.max(x), s)
        Y = self.smooth(LSQspl(X), window = 'blackman')[5:-5]

        self.ts[self.filters[fi]].fit_flux = Y
        self.ts[self.filters[fi]].fit_times = X

    def smooth(self, x, window_len=11, window='hanning'):
        """smooth the data using a window with requested size.
        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal 
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.
        input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
                flat window will produce a moving average smoothing.
        output:
        the smoothed signal
        example:
        t=linspace(-2,2,0.1)
        x=sin(t)+randn(len(t))*0.1
        y=smooth(x)
        see also: 
        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
        scipy.signal.lfilter
        TODO: the window parameter could be the window itself if an array instead of a string
        NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
        https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
        2017-07-13 (last modified), 2006-10-31 (created)
        Section author: Unknown[1], GaelVaroquaux, Unknown[142], Unknown[143], Unknown[144], Unknown[145], Unknown[146], Unknown[147], WesTurner, Christian Gagnon, clecocel
        """
        ## Should this be removed from the class and instead be relegated to being a utility? Can this utiity be located in Dorado dependencies?
        if x.ndim != 1:
                raise ValueError("smooth only accepts 1 dimension arrays.")

        if x.size < window_len:
                raise ValueError("Input vector needs to be bigger than window size.")


        if window_len<3:
                return x


        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError("Window is none of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


        s = np.r_[x[window_len-1:0:-1], x, x[-2:-window_len-1:-1]]
        # print(len(s))
        if window == 'flat': # moving average
                w = np.ones(window_len,'d')
        else:
                w = eval('np.' + window +'(window_len)')

        y = np.convolve(w/w.sum(), s, mode='valid')

        return y

    def analyze(self, fi, graphical = True, terms = None, s = None):
        '''
        analyze is a wrapper function that runs a pipeline of analysis functions in the 
        Fournax target class, self.superfit, self.tomlFind, self.OMinusC, and self.fourFind
        Parameters
        ----------
        filter: string
                filter string of filter to perform analysis on.
        
        graphical: boolean
                default is True. Controls whether or not graphical results are output 
                (NOTE:: parameter currently not in use)
        '''
        if terms == None:
            terms = int(np.floor(len(self.ts[self.filters[fi]].flux)/3))
        if s == None:
            s = len(self.ts[self.filters[fi]].flux)
            
        self.superfit(filter = fi, terms = terms, s = s)
        # Find times of max light
        self.tomlFind(fi = fi)
        # calculate O-C
        self.OMinusC(fi = fi)
        # frequency analysis
        self.fourFind(fi = fi)
    
    def fourFind(self, fi, fitted = True):
        '''
        fourFind locates distinct peak structures in a Fourier power spectra 
        generated via a photometric timeseries stored within self.ts.
        Parameters
        ----------
        filter: string
            filter string of filter to perform analysis on.
        fitted: boolean
                default is True. Controls whether fitted timeseries or raw
                timeseries is used for analysis.
        '''
        if fitted and (self.ts[self.filters[fi]].fit_flux != []):
            Y = self.ts[self.filters[fi]].fit_flux
            X = self.ts[self.filters[fi]].fit_times 
        else:
            Y = self.ts[self.filters[fi]].flux
            X = self.ts[self.filters[fi]].times 
        # Set up a Fourier power spectra from photometric amplitude values
        # Compute real valued Fourier transform
        f = np.fft.fft(Y)
        p = np.square(np.abs(f))
        # timestep currently defaults to units of days, whereas exposure time is in seconds
        timestep = (np.mean(self.ts[self.filters[fi]].exptimes) * un.s).to(un.day).value

        # Build an array of frequencies to plot against
        freq_vec = np.fft.fftfreq(len(Y), d = timestep)

        # find frequency peaks
        peaks, _ = find_peaks(p, height = np.mean(p))

        # Return frequency vector for plotting?
        self.freq = peaks
        self.freq_vec = freq_vec
        self.power_vec = p

    def tomlFind(self, fi, fitted = True):
        '''
        tomlFind locates the distinct times of max light in a photometric timeseries
        stored in self.ts
        Parameters
        ----------
        filter: string
            filter string of filter to perform analysis on.
        fitted: boolean
                default is True. Controls whether fitted timeseries or raw
                timeseries is used for analysis.
        '''
        if fitted and (self.ts[self.filters[fi]].fit_flux != []):
            Y = self.ts[self.filters[fi]].fit_flux
            X = self.ts[self.filters[fi]].fit_times 
        else:
            Y = self.ts[self.filters[fi]].flux
            X = self.ts[self.filters[fi]].times 
        
        peaks, _ = find_peaks(Y, height = 1.1 * np.mean(Y))

        toml = X[peaks]

        self.ts[self.filters[fi]].toml = toml




class TESSeract(Fournax):
    '''
    This class is currently a work in progress. In the future it will handle targets
    obtained from TESS(the Transiting Exoplanet Survey Satellite). Stay tuned!
    '''
    def __init__(self, tid = None, sector = None, type = None):
        self.tid = tid # find this if possible
        self.sector = sector # set this later if not given, with all possibilities
        self.type = type

        
