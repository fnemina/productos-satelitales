import numpy as np
from netCDF4 import Dataset
from osgeo import gdal
from osgeo import osr

# Define KM_PER_DEGREE
KM_PER_DEGREE = 111.32

# GOES-16 Extent (satellite projection) [llx, lly, urx, ury]

# GOES16_EXTENT = [-5434894.885056, -5434894.885056, 5434894.885056, 5434894.885056]

# GOES-16 Spatial Reference System
sourcePrj = osr.SpatialReference()
sourcePrj.ImportFromProj4('+proj=geos +h=35786023.0 +a=6378137.0 +b=6356752.31414 +f=0.00335281068119356027 +lat_0=0.0'
                          ' +lon_0=-75.0 +sweep=x +no_defs')

# Lat/lon WSG84 Spatial Reference System
targetPrj = osr.SpatialReference()
targetPrj.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')


def export_image(image, path):
    driver = gdal.GetDriverByName('netCDF')
    return driver.CreateCopy(path, image, 0)


def get_geo_t(extent, nlines, ncols):
    # Compute resolution based on data dimension
    resx = (extent[2] - extent[0]) / ncols
    resy = (extent[3] - extent[1]) / nlines
    return [extent[0], resx, 0, extent[3], 0, -resy]


def get_scale_offset(path, var):
    nc = Dataset(path)
    scale = nc.variables[var].scale_factor
    offset = nc.variables[var].add_offset
    nc.close()
    return scale, offset


def get_extent(path, _var=None):
    nc = Dataset(path)
    h = nc.variables['goes_imager_projection'].perspective_point_height
    x1 = nc.variables['x_image_bounds'][0] * h
    x2 = nc.variables['x_image_bounds'][1] * h
    y1 = nc.variables['y_image_bounds'][1] * h
    y2 = nc.variables['y_image_bounds'][0] * h
    goes_extent = [x1, y1, x2, y2]
    nc.close()
    return goes_extent


# noinspection DuplicatedCode
def remap(path, extent, resolution, driver, var):
    """Read scale/offset from file"""
    goes_extent = get_extent(path)
    # Build connection info based on given driver name
    if driver == 'NETCDF':
        connection_info = 'NETCDF:\"' + path + '\":' + var
    else:  # HDF5
        connection_info = 'HDF5:\"' + path + '\"://' + var
    # Open NetCDF file (GOES-16 data)  
    raw = gdal.Open(connection_info, gdal.GA_ReadOnly)
    # Setup projection and geo-transformation
    raw.SetProjection(sourcePrj.ExportToWkt())
    raw.SetGeoTransform(get_geo_t(goes_extent, raw.RasterYSize, raw.RasterXSize))
    # Compute grid dimension
    if resolution is not None:
        sizex = int(((extent[2] - extent[0]) * KM_PER_DEGREE) / resolution) - 1
        sizey = int(((extent[3] - extent[1]) * KM_PER_DEGREE) / resolution) - 1
    else:
        sizex = raw.RasterXSize
        sizey = raw.RasterYSize
    # Get memory driver
    mem_driver = gdal.GetDriverByName('MEM')
    # Create grid
    grid = mem_driver.Create('grid', sizex, sizey, 1, gdal.GDT_Float32)
    # Setup projection and geo-transformation
    grid.SetProjection(targetPrj.ExportToWkt())
    grid.SetGeoTransform(get_geo_t(extent, grid.RasterYSize, grid.RasterXSize))
    # Perform the projection/resampling
    gdal.ReprojectImage(raw, grid, sourcePrj.ExportToWkt(), targetPrj.ExportToWkt(), gdal.GRA_NearestNeighbour,
                        options=['NUM_THREADS=ALL_CPUS'])
    # Close file
    raw = None
    # Read grid data
    array = grid.ReadAsArray()
    # Mask fill values (i.e. invalid values)
    np.ma.masked_where(array, array == -1, False)
    return grid, array
