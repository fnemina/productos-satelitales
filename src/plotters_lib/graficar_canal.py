#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""
import os
import logging
from pathlib import Path
import datetime
import time as t
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as image
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
# from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator)
import fiona
import shapely.geometry as sgeom
from netCDF4 import Dataset
from osgeo import osr, gdal
from PIL import Image
from plotters_lib.cpt_convert import load_cpt
from plotters_lib.plot import get_geo_t

from config.constants import LOGO, SHAPEFILES, EXTENT
from config.logging_conf import GOES_LOGGER_NAME
from wrf_api.goes_api_ingest import post_img_to_api

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


def generar_imagen(nombre: str, _metadatos: dict, path_imagenes: str, canal, extent=EXTENT):
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
    # print(bandtype)

    # %% Proyecciones

    # GOES-16 Spatial Reference System
    source_prj = osr.SpatialReference()
    proj_str = f"+proj=geos +h={h} +a={a} +b={b} +f={f} lat_0={lat_0} +lon_0={lon_0} +sweep=x +no_defs"
    source_prj.ImportFromProj4(proj_str)
    # Lat/lon WSG84 Spatial Reference System
    target_prj = osr.SpatialReference()
    target_prj.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

    # Setup projection and geo-transformation
    raw.SetProjection(source_prj.ExportToWkt())
    raw.SetGeoTransform(get_geo_t(GOES16_EXTENT, raw.RasterYSize, raw.RasterXSize))

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
    plt.figure(figsize=(10, 9.5))

    # 4326 es WGS84 (LatLOn)
    # Reconfigura el extent del mapa para cartopy
    extent_map = [extent[0], extent[2], 
                    extent[1], extent[3]] # [west, east, south, north]
    # 4326 es WGS84 (LatLOn)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent_map, crs=ccrs.PlateCarree())

    # Draw the shapefiles
    # Agrego líneas de costa
    ax.coastlines(resolution='10m', color='black', linewidth=0.8)
    ax.stock_img()

    # Agrego grilla en la figura
    gl = ax.gridlines(draw_labels=False, color='white', alpha=0.5, linestyle=':', linewidth=0.6)
    
    # Abro archivo shapefile con departamentos
    shpfile = f"{SHAPEFILES}/departamentos.shp"
    
    with fiona.open(shpfile) as records:
        geometries = [sgeom.shape(shp['geometry']) for shp in records]
        # Agrego a la figura cada uno de los departamentos
    ax.add_geometries(geometries, ccrs.PlateCarree(), edgecolor='#808080', facecolor='none', linewidth=0.2,
                        alpha=1)
    
    shpfile = f"{SHAPEFILES}/008_límites_provinciales.shp"
    with fiona.open(shpfile) as records:
        geometries = [sgeom.shape(shp['geometry']) for shp in records]
        # Agrego a la figura cada uno de los departamentos
    
    ax.add_geometries(geometries, ccrs.PlateCarree(), edgecolor='#FFFFFF', facecolor='none', linewidth=0.6,
                        alpha=1)
    crs = ccrs.PlateCarree()

    plt.subplots_adjust(left=0.02, right=0.98, top=1, bottom=0.02)

    # Converts a CPT file to be used in Python
    cpt = load_cpt(canal.cptfile)

    # Makes a linear interpolation
    cpt_convert = LinearSegmentedColormap('cpt', cpt)

    # Plot the GOES-16 channel with the converted CPT colors
    # (you may alter the min and max to match your preference)
    if canal.visible:
        map = ax.imshow(
            array,
            origin='upper',
            cmap='gray',
            vmin=0.,
            vmax=1.,
            extent=extent_map, 
            transform=crs
        )
    else:
        temp = array - 273.15
        map = ax.imshow(
            temp,
            origin='upper',
            cmap=cpt_convert,
            vmin=-90,
            vmax=50,
            extent=extent_map, 
            transform=crs
        )
    # Redefino la grilla
    parallels=list(np.arange(extent_map[2]//5*5, extent_map[3], 5))
    meridians=list(np.arange(extent_map[0]//5*5, extent_map[1], 5))
    gl.xlocator = mticker.FixedLocator(meridians)
    gl.ylocator = mticker.FixedLocator(parallels)
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter() 

    # Add a black rectangle in the bottom to insert the image description
    # Max Lon - Min Lon
    lon_difference = (extent[2] - extent[0])
    current_axis = plt.gca()
    current_axis.add_patch(Rectangle(
        (extent[0], extent[1]),
        lon_difference,
        lon_difference * 0.020,
        alpha=1,
        zorder=3,
        facecolor='black'
    ))

    titulo_negro = " GOES-16 ABI Canal %02d %s UTC" % (icanal, _metadatos['fecha_img'].strftime('%Y-%m-%d %H:%M'))
    institucion = "CONAE-Argentina"
    # Add the image description inside the black rectangle
    # Max lat - Min lat
    lat_difference = (extent[3] - extent[1])
    plt.text(extent[0], extent[1] + lat_difference * 0.005, titulo_negro,
             horizontalalignment='left', color='white', size=7)
    plt.text(extent[2], extent[1] + lat_difference * 0.005, institucion,
             horizontalalignment='right', color='white', size=7)

    
    # Insert the colorbar at the right
    im_ratio = array.shape[0]/array.shape[1]
    print(im_ratio)
    cb = plt.colorbar(map, ax=ax, location='bottom', fraction=0.06*im_ratio, pad=0.01*im_ratio, aspect=60)
    # # Remove the colorbar outline
    cb.outline.set_visible(True)
    # # Remove the colorbar ticks
    cb.ax.tick_params(width=0)
    # # Put the colobar labels inside the colorbar
    cb.ax.xaxis.set_tick_params(pad=0)
    # # Change the color and size of the colorbar labels
    cb.ax.tick_params(axis='x', colors='black', labelsize=8)
    cb.set_label(canal.unidad)

    img = image.imread(LOGO)
    plt.figimage(
        img,
        25,
        100,
        # ax.figure.bbox.xmax - 160,
        # ax.figure.bbox.ymax - 70,
        zorder=1
    )
    # ax.text(0,
    #         1.10,
    #         canal.nombre,
    #         verticalalignment='top',
    #         transform=ax.transAxes,
    #         fontsize=20
    #         )
    # ax.text(0,
    #         1.03,
    #         metadatos['fecha_img'].strftime('%Y-%m-%d %H:%M') + ' UTC',
    #         verticalalignment='top',
    #         transform=ax.transAxes,
    #         fontsize=12
    #         )
    # ax.text(1,
    #         1.03,
    #         'GOES-16 ABI Canal %02d' % icanal,
    #         horizontalalignment='right',
    #         verticalalignment='top',
    #         transform=ax.transAxes,
    #         fontsize=10
    #         )

    # grabar a png
    fecha = _metadatos['fecha_img']
    fecha_str = fecha.strftime('%Y-%m-%d_%H_%M')
    path_fecha = fecha.strftime('%Y_%m/%d')
    path_imagen = f"{path_imagenes}/C{_metadatos['icanal']}_ARG{fecha_str}_WGS84.png"
    plt.savefig(path_imagen)
    plt.clf()
    plt.close()
    print('- finished! Time:', t.time() - start, 'seconds')
    img_api_path = f"GOES/{path_fecha}/{canal.codigo}/C{_metadatos['icanal']}_ARG{fecha_str}_WGS84.png"
    post_img_to_api(img_api_path, fecha, producto=canal.codigo, campo_prod='short_name')