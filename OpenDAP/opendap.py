import netCDF4
import urllib2
import logging
import numpy as np
import pandas as pd
from collections import OrderedDict
from datetime import datetime, timedelta
import xml.etree.ElementTree as ElementTree


class Dataset:
    '''OpenDAP dataset object
    
    Similar to netCDF4.Dataset but approaches OpenDAP directory with
    successing netCDF files as if it were a single netCDF file. It
    parses the catalog.xml file to objtain an index of all available
    netCDF files. If the first dimension of a variable is a time
    dimension, the entire dataset can be indexed using datatime
    objects rather than indices.

    Attributes
    ----------
    variables : OrderedDict
        Dictionary with variable objects
    dimensions : OrderedDict
        Dictionary with dimension objects

    Examples
    --------
    >> dt1 = datetime(2013,6,1)
    >> dt1 = datetime(2013,6,31)
    >> with Dataset('http://opendap.deltares.nl/thredds/catalog/opendap/tudelft/thermal_ir/kijkduin/') as nc:
    >>     df = nc.variables['temperature'][dt1:dt2:6,:,:]
    >> plt.imshow(df.ix[0])

    '''
    
    variables = OrderedDict()
    dimensions = OrderedDict()
    
    time_var = None
    time_units = 'seconds since 1970-01-01 00:00:00'
    time = []
    
    ncfiles = []
    file2idx = []
    idx2file = []
    
    def __init__(self, url, time_var='time', max_files=np.inf, *args, **kwargs):
        '''Initialization method

        Parameters
        ----------
        url : str
            URL to OpenDAP directory
        time_var : str, optional
            Name of time variable
        max_files : int, optional
            Maximum number of files to include in dataset
        '''
        
        self.url = url
        self.time_var = time_var
        self.max_files = max_files
        
        self.read_catalog()
        self.read_files()
        
        
    def __enter__(self, *args, **kwargs):
        return self
    
    
    def __exit__(self, *args, **kwargs):
        pass
    
    
    def __repr__(self):
        s = ''
        s += 'OpenDAP Dataset:\n'
        s += '  url:        %s\n' % self.url
        s += '  files:      %d\n' % len(self.ncfiles)
        s += '  length:     %d\n' % len(self.time)
        s += '  dimensions: %d\n' % len(self.dimensions)
        s += '  variables:  %d\n\n' % len(self.variables)
        s += repr(self.variables)
        return s
        
        
    def read_catalog(self):
        '''Read hosted ncfiles from catalog XML'''
        
        fp = urllib2.urlopen('%scatalog.xml' % self.url)
        root = ElementTree.parse(fp)
        fp.close()

        self.ncfiles = []
        for item in root.iter():
            if item.tag.endswith('dataset'):
                if item.attrib['name'].endswith('.nc'):
                    self.ncfiles.append(item.attrib['name'])
            if len(self.ncfiles) >= self.max_files:
                break
        self.ncfiles.sort()


    def read_files(self):
        '''Read number of time steps, variables and dimensions from individual netCDF files'''
        
        for i, nc in self.iterate_files():
            if nc:
                nt = len(nc.variables[self.time_var])
                self.time.extend([None] * nt)
                self.idx2file.extend([i] * nt)
                self.file2idx.extend(range(nt))
                
                self.variables.update(nc.variables)
                self.dimensions.update(nc.dimensions)

        for k in self.variables.iterkeys():
            self.variables[k] = Variable(self, k)
            
        for k in self.dimensions.iterkeys():
            self.dimensions[k] = Dimension(self, k)
            
        
    def iterate_files(self, fids=None):
        '''Iterator over all files returning file id's and netCDF4 objects

        Parameters
        ----------
        fids : list
            List of file ids that should be included in iteration

        Example
        -------
        >> # iterate over first 5 files
        >> for i, nc in iterate_files(fid=range(5)):
        >>     print nc.variables
        '''
        
        for fid, ncfile in enumerate(self.ncfiles):
            if fids is not None and fid not in fids: continue
                
            url = '%s/%s' % (self.url, ncfile)
            
            try:
                with netCDF4.Dataset(url) as nc:
                    yield (fid, nc)
            except:
                logging.warn('Cannot read %s' % url)
                yield (fid, None)
    
    
    def ncattrs(self):
        '''Retrieve netCDF attributes

        Returns
        -------
        list
            List with attributes for each netCDF file
        '''
        
        attrs = []
        for fid, nc in self.iterate_files():
            if nc:
                attrs.append(nc.ncattrs())
        return attrs
    
    
    def getncattr(self, name):
        '''Retrieve netCDF attribute values

        Parameters
        ----------
        name : str
            Name of attribute

        Returns
        -------
        list
            List with attribute values for each netCDF file
        '''
        
        attrs = []
        for fid, nc in self.iterate_files():
            if nc:
                attrs.append(nc.getncattr(name))
        return attrs
    
    
class Variable:
    '''OpenDAP variable object

    Implements a getter function to retrieve data over multiple file
    in an OpenDAP.Dataset dataset object.

    Attributes
    ----------
    var : str
        Name of variable to be referenced
    dimensions : OrderedDict
        Dictionary with dimension objects
    '''
    
    dimensions = OrderedDict()
    
    var = None
    
    slice_time = slice(None, None, None)
    slice_u = slice(None, None, None)
    slice_v = slice(None, None, None)
    
    def __init__(self, dataset, var):
        '''Initialization method

        Parameters
        ----------
        dataset : OpenDAP.Dataset
            Parent dataset object where variable belongs to
        var : str
            Name of variable to be referenced
        '''
        
        self.dataset = dataset
        self.var = var
        
        self.dimensions = Dimension(self, var)
    
    
    def __getitem__(self, s):
        '''Getter function with datetime slicing

        Parameters
        ----------
        s : tuple
            3-tuple with slice for each dimension

        Returns
        -------
        pandas.DataFrame or pandas.Panel
            Pandas object with concatenated data from netCDF files
        '''
        
        if type(s) is not tuple:
            raise IndexError('Invalid number of dimensions [1]')
        elif len(s) != 3:
            raise IndexError('Invalid number of dimensions [%d]' % len(s))
        else:
            self.slice_time = self.__convert_datetime_slice(s[0])
            self.slice_u = s[1]
            self.slice_v = s[2]

            return self.__get_netcdf_data()
            

    def ncattrs(self):
        '''Retrieve netCDF attributes

        Returns
        -------
        list
            List with attributes for each netCDF file
        '''
        
        attrs = []
        for fid, nc in self.dataset.iterate_files():
            if nc:
                attrs.append(nc.variables[self.var].ncattrs())
        return attrs
    
    
    def getncattr(self, name):
        '''Retrieve netCDF attribute values

        Parameters
        ----------
        name : str
            Name of attribute

        Returns
        -------
        list
            List with attribute values for each netCDF file
        '''
        
        attrs = []
        for fid, nc in self.dataset.iterate_files():
            if nc:
                attrs.append(nc.variables[self.var].getncattr(name))
        return attrs

            
    def __get_netcdf_data(self):
        '''Retrieve and concatenate data from netCDF files

        Uses set slices and bisection search algorithm to retrieve
        data from multiple netCDF files in dataset.

        Returns
        -------
        pandas.DataFrame or pandas.Panel
            Pandas object with concatenated data from netCDF files
        '''
        
        fids = np.asarray(self.dataset.idx2file[self.slice_time])
        idxs = np.asarray(self.dataset.file2idx[self.slice_time])

        time = []
        data = []
        for i, nc in self.dataset.iterate_files(fids):
            if nc:
                idx = idxs[fids == i]
                
                units = self.__getunits(nc.variables[self.dataset.time_var])
                t = netCDF4.num2date(nc.variables[self.dataset.time_var][idx], units)
                
                d = nc.variables[self.var][idx,
                                           self.slice_u,
                                           self.slice_v]
                
                time.extend(t)
                data.extend(d)
                
        if len(time) <= 1:
            return pd.DataFrame(np.squeeze(data))
        else:
            return pd.Panel(data, items=time)
    
    
    def __convert_datetime_slice(self, s):
        '''Convert datetimes in slice to indices

        Parameters
        ----------
        s : slice or index or datetime
            Index to be converted
        
        Returns
        -------
        slice or index or datetime
            Converted index
        '''
        
        if type(s) is slice:
            start = s.start
            stop = s.stop
            step = s.step

            if type(start) is datetime:
                start = self.__datetime2index(start, method='higher')
            if type(stop) is datetime:
                stop = self.__datetime2index(stop, method='lower')
            
            return slice(start, stop, step)
        else:
            if type(s) is datetime:
                s = self.__datetime2index(s, method='nearest')
            return s

    
    def __datetime2index(self, dt, method='nearest'):
        '''Search routine time axis

        Time axes are not read initially. If data is requested this
        routine searches for the index that corresponds to a given
        datetime. The routine uses a bisection search. If the exact
        datetime is found, the corresponding index is
        returned. Otherwise three matching methods are distinguished:
        nearest, lower or higher, that return the nearest index, the
        last index before the given datetime or the first index after
        the given datetime.

        Parameters
        ----------
        dt : datetime
            Datetime object searched
        method : str, optional
            Matching method (nearest, lower or higher)

        Returns
        -------
        int
            Time index corresponding to datetime
        '''
        
        ix_low, ix_high = self.__getbounds(dt)
        
        n = 0
        while ix_high - ix_low > 1:
            
            # run bisection search
            idx = int(round(ix_low + (ix_high - ix_low) / 2.))
            fid = self.dataset.idx2file[idx]
            
            for i, nc in self.dataset.iterate_files(fids=[fid]):
                n = n + 1
                
                if n > len(self.dataset.ncfiles):
                    raise RuntimeError('Cannot read netCDF files')
                
                if nc:
                    ii = np.where(np.asarray(self.dataset.idx2file) == i)[0]
                    i1, i2 = min(ii), max(ii)+1
                    if None in self.dataset.time[i1:i2]:
                        var = nc.variables[self.dataset.time_var]
                        units = self.__getunits(var)
                        self.dataset.time[i1:i2] = netCDF4.num2date(var[:],
                                                                    units)

                        ix_low, ix_high = self.__getbounds(dt)
                    elif fid > 0:
                        fid = fid - 1
                    else:
                        fid = len(self.dataset.ncfiles) - 1
                elif fid > 0:
                    fid = fid - 1
                else:
                    fid = len(self.dataset.ncfiles) - 1
            
        if method == 'lower':
            return ix_low
        elif method == 'higher':
            return ix_high
        elif dt - self.dataset.time[ix_low] <= self.dataset.time[ix_high] - dt:
            return ix_low
        else:
            return ix_high
                
                
    def __getbounds(self, dt):
        '''Find time bounds along time axes

        Parameters
        ----------
        dt : datetime
            Datetime object searched

        Returns
        -------
        2-tuple
            Last known index before given datatime
            and first knwon index after given datetime
        '''
        
        ix_low = -1
        ix_high = len(self.dataset.time)
        for i, t in enumerate(self.dataset.time):
            if t is None:
                continue
            elif t == dt: # exact match
                ix_low, ix_high = i
                break
            elif t < dt:
                ix_low = i
            elif t > dt:
                ix_high = i
                break
               
        return ix_low, ix_high
    
    
    def __getunits(self, var):
        '''Get unit description for netCDF variable'''
        
#        try:
#            units = var.units
#        except:
        units = self.dataset.time_units
        return units
        
        
class Dimension:
    '''OpenDAP dimension object

    Not implemented.
    '''
    
    def __init__(self, dataset, var):
        '''Initialization method

        Parameters
        ----------
        dataset : OpenDAP.Dataset
            Parent dataset object where dimension belongs to
        var : str
            Name of variable to be referenced
        '''
        
        self.dataset = dataset
        self.var = var
        
        
    def __getitem__(self, s):
        raise NotImplementedError
