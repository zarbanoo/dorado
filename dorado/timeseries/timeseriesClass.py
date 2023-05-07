#ts = QTable([times, exptimes, x, y, ray, decx, flux, fluxunc], names=('time', 'exptime', 'x', 'y', 'ra', 'dec', 'flux', 'flux_unc'), meta={'name': filter})

import numpy as np
from astropy.time import Time
from astropy.table import QTable, Table
from astropy.io import fits



__all__ = ['timeSeries']



class timeSeries:
    '''
        The timeSeries

        Attributes
        ----------

        name: str
            name of target in string format.

    '''
    def __init__(self, times, flux, name = 'timeseries', exptimes = [], x = [], y = [], ra = [], dec = [], flux_unc = [], 
        apsum = [], apsum_unc = [], fit_times = [], fit_flux = [], toml = [], OmC = [], cycle = [], table = None): 

        # TODO append info on time system
        self.times = times
        self.flux = flux
        self.name = name
        self.exptimes = exptimes

        self.x = x
        self.y = y
        self.ra = ra
        self.dec = dec

        self.flux_unc = flux_unc
        self.apsum = apsum
        self.apsum_unc = apsum_unc

        self.fit_times = fit_times
        self.fit_flux = fit_flux

        self.toml = toml
        self.OmC = OmC
        self.cycle = cycle
        self.table = table

        # self.symbo = None # make a symbolic expression to represent the curve analytically.

    def toTable(self, name = None):
        '''
        toTable is a convinience function that produces an atropy.Table object from
        the data stored within a dorado.timeseries instance.
        Parameters
        ----------
        name: string
            The name to be used for the table. If none is given, self.name (default is 'timeseries') is used.
        '''
        if (name == None):
            name = self.name
        # name is table name
        # TODO verify the mjd part is exporting properly
        self.table = QTable([self.times.mjd, self.flux], names=('time','flux'), meta={'name': name}) 

        colnom = ['flux_unc', 'exptime', 'apsum', 'apsum_unc', 'x', 'y', 'ra', 'dec']
        cols = [self.flux_unc, self.exptimes, self.apsum, self.apsum_unc, self.x, self.y, self.ra, self.dec]

        for col in range(len(colnom)):
            if cols[col] != []:
                try:
                    self.table[colnom[col]] = cols[col]
                except:
                    print('Error merging', colnom[col], ' with table.')
    
    def graph(self, c_scale = False, err = True):
        '''
        The graph function produces a graph of the photometric timeseries data as
        a lightcurve.
        '''
        # TODO will also need graph of frequencies
        # TODO check if there is a fit involved
        import matplotlib.pyplot as plt # NOTE where is the best place to put this
        if self.table == None:
            self.toTable()
        
        # set up plot
        # TODO set up plotting style
        fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
        ax.set_xlabel('Time')
        ax.set_ylabel('Flux')
        try:  
            ax.set_title(str(self.name) + ' Lightcurve') # TODO decide on title
        except:
            ax.set_title('Lightcurve') # TODO decide on title
        if c_scale:
            print('colour scaling not implemented yet, sorry bud...')
            z = 'yellow' # self.flux - np.min(self.flux)
        else:
            z = 'gold'
        if err:
            ax.errorbar(self.times.mjd, self.flux, xerr = self.exptimes / 2, yerr = self.flux_unc, c = z)
        else:
            ax.plot(self.times.mjd, self.flux, c = z)
        ax.grid()
        plt.tight_layout()
        plt.show()
        
