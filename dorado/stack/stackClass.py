import warnings
warnings.filterwarnings('ignore')
from astropy.time import Time

__all__ = ['Stack']

''' 
'dorado.stack' holds the base Stack class for handling correlated stacks of image objects. 
The Stack class can be used alone, but is designed to work best from within other Dorado 
components such as Ceres.
'''


class Stack:
    '''
    The Stack class encapsulates correlated stacks of image objects and their associated data.
    Stack exceeds simply storing stacks of images in an array or list by expanding upon the array
    object to hold additional data about the image stack such as a WCS object, the enclosed target,
    and the correlated calibration frames.

    Attributes
    ----------

    data: CCDdata array
        Array of CCDdata images.

    flat: CCDdata
        Flatfield calibration frame for data stack. Optional

    filter: str
        Name of filter data was collected in. Optional.

    times: array or list-like
        Time for each image as an 'astropy.Time' object. Optional.

    calibrated: Boolean
        Whether the data is calibrated or not. Default is 'None'. Optional.

    aligned: Boolean
        Whether the data is aligned or not. Default is 'None'. Optional.

    target: TOI object
        Instance of the TOI astronomical target class containing the target of interest in the stack. Optional.

    alignTo: int
        Index of the image which all other stack images should be aligned to. Default is  0. Optional.

    '''

    ## TODO :: auto identify targets in stack

    def __init__(self, data, flat = None, filter = '', times = [], calibrated = None,
        aligned = None, target = None, alignTo = 0):
        self.data = data
        self.flat = flat
        self.filter = filter
        self.length = len(data)
        self.calibrated = calibrated
        self.aligned = aligned

        self.target = target
        self.target_info = {} # dictionary of values, call simbad?

        self.times = times 
        self.wcs = None
        self.alignTo = alignTo
        self.solved = None
        self.base = None
        # include things like flux uncertainty etc.
        # save data
        # if flat is none, call clippy for help

        if self.filter == '':
            try:
                self.filter = data[0].header['filter']
            except:
                self.filter = ''
        
        if self.times == []:
            try: 
                self.get_times()
                # print('Times set')
            except:
                self.times = []
        
    def get_times(self):
        '''
        get_times sequences through the stack data and passes the FITS header timestamp
        for 'DATE-OBS' into 'astropy.time' and sets the resulting array of times as
        'self.times'.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Sets
        ----
        self.times: array or list-like
            array of 'astropy.time' objects for each image timestamp.

        '''
        times = []
        for im in self.data:
            times.append(Time(im.header['DATE-OBS'], format='fits'))
        self.times = times
        

    def get_target_info(self, target = None):
        '''
        get_target_info is a convinience function for setting an instance of TOI
        as the Stack target object.
        
        Parameters
        ----------
        target: TOI object
            Instance of the TOI astronomical target class containing the target of interest in the stack.

        Returns
        -------
        None

        Sets
        ----
        self.target: TOI object
            Target of interest in the stack.
        '''
        # TODO does this need to be a wrapper or need more here
        if target != None:
            self.target = target
        
   