#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""
import logging
import time as t
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
from netCDF4 import Dataset
from osgeo import osr, gdal
from plotters_lib.cpt_convert import load_cpt
from PIL import Image

from config.constants import EXTENT_WEBMET
from config.logging_conf import GOES_LOGGER_NAME
from wrf_api.goes_api_ingest import post_img_to_api

from plotters_lib.plot import get_geo_t

logger = logging.getLogger(GOES_LOGGER_NAME)

matplotlib.use("agg")

# Define KM_PER_DEGREE
KM_PER_DEGREE = 111.32

# GOES-16 Extent (satellite projection) [llx, lly, urx, ury]
GOES16_EXTENT = [-5434894.885056,
                 -5434894.885056,
                 5434894.885056,
                 5434894.885056]

RESOLUTION = 1.

def convertir_negro_transparente(ruta):
    """
    Convierte el color negro de una imagen a transparente.
    """
    img = Image.open(ruta)
    img = img.convert("RGBA")
    pixdata = img.load()
    width, height = img.size
    for y in range(height):
        for x in range(width):
            if pixdata[x, y] == (0, 0, 0, 255):
                pixdata[x, y] = (0, 0, 0, 0)
    img.save(ruta, "PNG")

def generar_imagen_webmet(nombre, _metadatos, path_imagenes, canal, extent=EXTENT_WEBMET):
    start = t.time()

    icanal = _metadatos['icanal']

    # Parametross de calibracion
    offset = _metadatos['offset']
    scale = _metadatos['scale']

    # Parametros de proyeccion
    lat_0 = str(_metadatos['lat_0'])
    lon_0 = str(_metadatos['lon_0'])
    h = str(_metadatos['h'])
    a = str(_metadatos['a'])
    b = str(_metadatos['b'])
    f = str(_metadatos['f'])

    # %% lectura y extraccion de informacion de la pasada
    connection_info = 'HDF5:\"' + nombre + '\"://' + canal.variable

    raw = gdal.Open(connection_info, gdal.GA_ReadOnly)

    # driver = raw.GetDriver().LongName

    band = raw.GetRasterBand(1)
    bandtype = gdal.GetDataTypeName(band.DataType)
    print(bandtype)

    # %% Proyecciones

    # GOES-16 Spatial Reference System
    source_prj = osr.SpatialReference()
    proj_str = f"+proj=geos +h={h} +a={a} +b={b} +f={f} lat_0={lat_0} +lon_0={lon_0} +sweep=x +no_defs"
    source_prj.ImportFromProj4(proj_str)
    # Lat/lon WSG84 Spatial Reference System https://epsg.io/3857
    target_prj = osr.SpatialReference()
    target_prj.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

    # Setup projection and geo-transformation
    raw.SetProjection(source_prj.ExportToWkt())
    raw.SetGeoTransform(get_geo_t(GOES16_EXTENT,
                                  raw.RasterYSize,
                                  raw.RasterXSize))

    # Compute grid dimension
    sizex = int(((extent[2] - extent[0]) * KM_PER_DEGREE) / RESOLUTION)
    sizey = int(((extent[3] - extent[1]) * KM_PER_DEGREE) / RESOLUTION)

    # Get memory driver
    mem_driver = gdal.GetDriverByName('MEM')

    # Create grid
    grid = mem_driver.Create('grid', sizex, sizey, 1, gdal.GDT_Float32)

    # Setup projection and geo-transformation
    grid.SetProjection(target_prj.ExportToWkt())
    grid.SetGeoTransform(get_geo_t(extent, grid.RasterYSize, grid.RasterXSize))

    # Perform the projection/resampling

    gdal.ReprojectImage(
        raw,
        grid,
        source_prj.ExportToWkt(),
        target_prj.ExportToWkt(),
        gdal.GRA_NearestNeighbour,
        options=['NUM_THREADS=ALL_CPUS']
    )

    # Read grid data
    array1 = grid.ReadAsArray()

    # Mask fill values (i.e. invalid values)
    np.ma.masked_where(array1, array1 == -1, False)

    # %% Calibracion
    array = array1 * scale + offset

    grid.GetRasterBand(1).SetNoDataValue(-1)
    grid.GetRasterBand(1).WriteArray(array)

    # %% Plot the Data ========================================
    # Create the basemap reference for the Rectangular Projection
    plt.clf()
    fig = plt.figure(frameon=False)
    fig.set_size_inches(25. * array.shape[1] / array.shape[0], 25, forward=False)
    proj3857 = ccrs.epsg(3857)
    ax = plt.axes([0., 0., 1., 1.], projection=proj3857)
    ax.set_axis_off()
    fig.add_axes(ax)
    # Get map extent
    # Reconfigura el extent del mapa para cartopy
    extent_map = [extent[0], extent[2], 
                    extent[1], extent[3]] # [west, east, south, north]
    # 3857 es WGS84 (LatLOn)
    ax.set_extent(extent_map, crs=proj3857)

    # Converts a CPT file to be used in Python
    cpt = load_cpt(canal.cptfile2)

    # Makes a linear interpolation
    cpt_convert = LinearSegmentedColormap('cpt', cpt)
    # Plot the GOES-16 channel with the converted CPT colors
    # (you may alter the min and max to match your preference)
    if canal.visible:
        ax.imshow(
            array,
            origin='upper',
            cmap='gray',
            vmin=0.,
            vmax=1.,
            extent=extent_map, 
            transform=proj3857
        )
    else:
        temp = array - 273.15
        temp[temp > -5] = 50
        ax.imshow(
            temp,
            origin='upper',
            cmap=cpt_convert,
            vmin=-90,
            vmax=50,
            extent=extent_map, 
            transform=proj3857
        )

    # grabar a png
    path_imagen = path_imagenes + '/GOES16' + '_\
' + _metadatos['fecha_img'].strftime('%Y%m%dT%H%M%S') + 'Z_C' + str(_metadatos['icanal']) + '.png'
    plt.savefig(path_imagen, transparent=True)
    plt.clf()
    plt.close()
    # Close file
    raw = None
    convertir_negro_transparente(path_imagen)
    print('- finished! Time:', t.time() - start, 'seconds')