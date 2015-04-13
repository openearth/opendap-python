# opendap-python

Similar to netCDF4.Dataset but approaches OpenDAP directory with
successing netCDF files as if it were a single netCDF file. It
parses the catalog.xml file to obtain an index of all available
netCDF files. If the first dimension of a variable is a time
dimension, the entire dataset can be indexed using datatime
objects rather than indices.
