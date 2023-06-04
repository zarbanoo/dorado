from astroquery.astrometry_net import AstrometryNet
from astroquery.exceptions  import TimeoutError
from astropy.wcs import WCS

ast = AstrometryNet()

__all__ =  ['plate_solve', 'get_night', 'getDateString', 'CCDDatax']

from astropy.nddata.ccddata import CCDData
CCDData._config_ccd_requires_unit = False
import astropy.units as un

CCDDatax = CCDData
CCDDatax.unit = un.adu # this is what the kids call hardcoding and I'm doing it because astropy refuses to fix their godforsaken unit issue with CCDData
# why did CCDData decide to randomly break? who knows
# why did CCDData start refusing to achknowledge me telling it explicitely multiple times that it doesn't need a unit
# because when I hand it a unit it ignores it and instead just reads it from the header? I don't know, does astropy team know?
# probably not. But its a thing and I really wish it wasn't at this hour of the night.


def plate_solve(dirarray, data = None, writearray = False):
    """
    plate_solve takes a dirarray pointing to a fits file to plate solve using astrometryNet. It returns
    the solved image with the WCS header and the WCS header itself if the solve succeeded. If 'writearray'
    is set with a filepath array, a copy of the solved image will be saved.
    
    Parameters
    ----------
    dirarray: str array
        filepath array to image to be solved.
    data: CCDData
        optional image data to combine with the WCS header on return.save. Default is data to be solved.
    Returns
    -------
    
    """
    # TODO better handle cache from here and check if API key
    # can we print a summary of the solve?
    # path = Dorado.dordir
    from pathlib import Path
    path = Path()
    for dir in dirarray:
        path = path / dir

    if data == None:
        data = CCDDatax.read(path) #, unit = Dorado.unit) ## NOTE edited

    trying = True
    submission_id = None
    num = 0

    while trying:
            try:
                if not submission_id:
                    wcs_header = ast.solve_from_image(path, force_image_upload=True, submission_id=submission_id, solve_timeout=300)
                else:
                    print('Monitoring: try #', num)
                    wcs_header = ast.monitor_submission(submission_id, solve_timeout=300)
            except TimeoutError as e:
                print(TimeoutError)
                num = num + 1
                print('Timed out: try #', num)
                submission_id = e.args[1]

            if wcs_header != None:
                # got a result, so terminate while loop
                trying = False
    if wcs_header:
        # Code to execute when solve succeeds
        print('Solve succeeded! :)')
        wcs_hdu = data
        wcs_hdu.header = wcs_header
        if writearray:
            path = Path()
            for dir in writearray:
                path = path / dir
            wcs_hdu.write(path, overwrite = True)

        return wcs_hdu, wcs_header
    else:
        # Code to execute when solve fails
        print('Solve failed! :(')
        return 


import datetime

def get_night():
    """
    get_night obtains a timestring for the most recent(previous) night based on local/hardware
    time provided by datetime. The format follows yyyy-mm-(dd-1)+dd where (dd-1) is last nights 
    day of the month. This currently does not support dates at which the start of the observing
    night was the last day of the month.

    Parameters
    ----------
    None

    Returns
    -------
    night: str
            Timestring for the most recent night.
    """
    # currently does not support first/last of the month
    # option for last night or tonight
    # UTC or local time?
    date = datetime.date.today()
    year = date.year
    month = date.month
    day = date.day
    # date1 = date2 - 1
    # night = str(year) + '-' + str(month) + '-' + str(date1) + '+' + str(date2)


    if day < 10:
        daystr = '0' + str(day)
    else:
        daystr = str(day)
    if day < 9:
        day2str = '0' + str(day + 1)
    else:
        day2str = str(day + 1)
    if month < 10:
        monthstr = '0' + str(month)
    else:
        monthstr = str(month)


    night = str(year) + '-' + monthstr + '-' + daystr + '+' + day2str




    return night
        
from astropy.time import Time

def getDateString(epoch, utc = None):
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
    if type(epoch) == str:
        epoch = Time(epoch, format='fits')
    elif type(epoch) != type(Time):
        raise Exception("Could not parse input")
    # epoch = date.ymdhms <- desired time instance format
    # Fall back to Dorado UTC offset if UTC offset is not provided
    if not (utc):
        utc = Dorado.UTCoffset
    epdt = epoch.datetime
    # if the hour is less than the utc offset of the site then the utc date is one ahead of local time
    if (epdt.hour + utc) < 0:
        day = str(epdt.day - 1)
        day2 = str(epdt.day)
        month = str(epdt.month)
        if (epdt.day - 1) < 10:
            day = '0' + str(epdt.day - 1)
        if epdt.day < 10:
            day2 = '0' + str(epdt.day)
        if epdt.month < 10:
            month = '0' + str(epdt.month)

    else: 
        day = str(epdt.day)
        day2 = str(epdt.day + 1)
        month = str(epdt.month)
        if epdt.day < 10:
            day = '0' + str(epdt.day)
        if epdt.day < 9:
            day2 = '0' + str(epdt.day + 1)
        if epdt.month < 10:
            month = '0' + str(epdt.month)

    return str(epdt.year) + '-' + month + '-' + day + '+' + day2
    