import os

__all__ = ['cache_core']

class cache_core:
    '''
    '''
    def __init__(self):
        pass

    def mkcacheObj(self, object, subcache = False):
        """
        mkcacheObj is a convenience function that creates a cache object for a .fits compatable
        object in the self.dordir/cache directory.
        
        Parameters
        ----------
        object: CCDData
            object to be saved in the cache directory for access.
        subcache: str
            optional subcache file for cache organization. Optional
        Returns
        -------
        
        """
        if subcache:
            cachedir = Dorado.dordir / 'cache' / subcache
            dirarray = ['cache', subcache]
        else:
            cachedir = Dorado.dordir / 'cache' 
            dirarray = ['cache']
        # if isiterable(object):
        #     # print('This function does not currently support iterable objects.')
        #     return warnings.WarningMessage('This function does not currently support iterable objects.')
        # else:
        #     files, _ = self.diread(dirarray)
        #     fname = 'cache_object_' + str(len(files) + 1) + '.fits'
        #     object.write(fname)
        #     return(fname, cachedir)

        # check if iterable
        # NOTE moved to reader
        files, _ = Dorado.reader.diread(dirarray)
        # TODO will the cache object always be a fits file? what about JSON?
        fname = 'cache_object_' + str(len(files) + 1) + '.fits'
        # using object name seems fishy, is it a protected word?
        object.write(cachedir / fname, overwrite = True)
        return(fname, cachedir)
        
    def delcacheObj(self, fname, subcache = False):
        """
        delcacheObj is a convenience method that deletes a selected cache object given is file name.
        
        Parameters
        ----------
        fname: string
            filestring of cache object to delete.
        subcache: string
            filestring of subcache containing the cache object to delete.
        Returns
        -------
        
        """
        # TODO extend to clearing cache by accepting None type as fname and reading the dir for all files
        if subcache:
            cachedir = Dorado.dordir / 'cache' / subcache
        else:
            cachedir = Dorado.dordir / 'cache'

        os.remove(cachedir / fname)
        