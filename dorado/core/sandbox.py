def _single_level(self, files):
    bias    = self._read_stack(files, 'bias')
    flats  = self._read_stack(files, 'flats')
    lights = self._read_stack(files, 'lights')
    return bias, flats, lights
    
def _multi_level_top(self, files, directories):
    biasdir = [s for s in directories if s.name in self.biasstr]
    # There should never be multiple bias directories; unless bias were taken with multiple binning settings or taken multiple times.
    # Consult chief data officer about this. TODO
    bias_files, _ = self.diread(biasdir[0])
    bias    = self._read_stack(bias_files, 'bias')
    # It is possible to have multiple flats directories for different filters, they should all be contained in a single base flats directory tho
    flatsdir = [s for s in directories if s.name in self.flatsstr]
    flats = self._multi_level_sub(flatsdir[0], 'flats')
    if self.unique_lights:
        lightsdir = [s for s in directories if (s.name not in self.flatsstr) and (s.name not in self.biasstr)] # may hand back more than one directory
    else:
        lightsdir = [s for s in directories if s.name in self.lightsstr]  # may hand back more than one directory
    lights = self.multi_level_sub(lightsdir[0], 'lights')
    return bias, flats, lights

def _multi_level_sub(self, subdir, stack_type):
    files, directories = self.diread(subdir)
    sub_files = []
    if len(directories) == 0:
        print('Single directory ', stack_type, ' organization format detected.')
        sub_files  = self._read_stack(files, stack_type)
    elif len(files) == 0:
        print('Multi directory ', stack_type, ' organization format detected.')
        for d in directories:
            # should not be any further sub folders
            f, _ = self.diread(d)
            sub_files += _read_stack(f, stack_type) # should we continue to hand stack_type if only the desired frame type should be present this deep in the mix?
    return sub_files
        
    

def _chkPhi(self, stack_files):
    # sequence through headers saving every filter keyword
    phi_list = []
    phi_stacks = {}
    for im in stack_files:
        phi = im.header['FILTER'] # hardcoded filter keyword for now
        if phi in phi_list: # check if filter is already accounted for 
            phi_stacks[phi].append(im)
        else:
            phi_stacks[phi] = [im]
    for p in phi_list:
        print('Found ', len(phi_stacks[p]), ' ', p, ' frames.')
    
    return phi_list, phi_stacks

def chkCal(self, frames):
    
    pass

def _read_stack(self, files, stack_type):
    if (stack_type == 'lights') & (self.unique_lights):
        stackstr = self.stack_types['bias'] + self.stack_types['flats'] # NOTE:: in the future this should be modded to include all non light stack types.
    else:
        try:
            stackstr = self.stack_types[stack_type]
        except:
            stackstr = stack_type
        
    stack_list = []
    for st in stackstr:
        for s in files:
            if (stack_type == 'lights') & (self.unique_lights):
                if (str(st)) not in str(s.name):
                    stack_list.append(s)
            else:
                if (str(st)) in str(s.name):
                    stack_list.append(s)
    stack_files = []
    if len(stack_list) > 50:
        n = 0 # this is a very janky way of writing this, enjoy
        for f in (pbar := tqdm(stack_list, colour = 'green')):
            n += 1
            pbar.set_description('Reading ' + stack_type + ' : ' + str(n))
            pbar.refresh()
            hdu = CCDData.read(f.path)
            stack_files.append(hdu)
    else:
        for f in stack_list:
            hdu = CCDData.read(f.path)
            stack_files.append(hdu)
    return stack_files

def getBias(mjdstr):
    target_mjd = int(mjdstr)
    biasdir = Dorado.dordir / 'data' / 'bias' 
    contents = os.scandir(path = biasdir)
    candidates = []
    for entry in contents:
        entry_mjd = int(entry.name[0:4]) # grab first five characters of filename. Should be the mjd number
        if np.abs(target_mjd - entry_mjd) <= self.calibration_window:
            candidates.append(entry_mjd)
    if len(candidates) == 0:
        print('No viable bias frames found. Try increasing your calibration window.')
    else:
        mjd_to_use = candidates[np.argmin(candidates- entry_mjd)]
    
    fname = str(mjd_to_use) + '_Bias.fits'
    bias = CCDData.read(biasdir / fname) #, unit = Dorado.unit) ## NOTE edited
    return bias

def getFlat(self, mjdstr, phi):
    target_mjd = int(mjdstr)
    flatdir = Dorado.dordir / 'data' / 'flats' 
    contents = os.scandir(path = flatdir)
    candidates = []
    for entry in contents:
        if phi in entry.name:
            entry_mjd = int(entry.name[0:4]) # grab first five characters of filename. Should be the mjd number
            if np.abs(target_mjd - entry_mjd) <= self.calibration_window:
                candidates.append(entry_mjd)
    if len(candidates) == 0:
        print('No viable flat frames found. Try increasing your calibration window or using a different filter.')
    else:
        mjd_to_use = candidates[np.argmin(candidates- entry_mjd)]
    
    fname = str(mjd_to_use) + '_' + str(phi) + '_flat.fits'
    flat = CCDData.read(flatdir / fname) #, unit = Dorado.unit) ## NOTE edited
    return flat