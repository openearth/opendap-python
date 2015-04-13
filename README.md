opendap-python
==============

Similar to netCDF4.Dataset but approaches OpenDAP directory with
successing netCDF files as if it were a single netCDF file. It
parses the catalog.xml file to obtain an index of all available
netCDF files. If the first dimension of a variable is a time
dimension, the entire dataset can be indexed using datatime
objects rather than indices.

Example
-------

```python
url = 'http://opendap.deltares.nl/thredds/catalog/opendap/tudelft/thermal_ir/kijkduin/'
with OpenDAP.Dataset(url) as nc:
  dt1 = datetime.datetime(2013,6,1)
  dt1 = datetime.datetime(2013,6,5)
  df = nc.variables['temperature'][dt1:dt2:10,:,192]
  
  # plot timestack
  fig, axs = subplots()
  axs.matshow(df, cmap='coolwarm')
```

Known issues
------------

* Dimension object is not yet implemented
* Datetime slicing with timedelta steps is not yet implemented
* NetCDF4.Dataset functions need to be piped
