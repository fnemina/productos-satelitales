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
from plotters_lib.plot import get_geo_t, get_metadatos, convertir_negro_transparente

from config.constants import LOGO, SHAPEFILES, EXTENT, EXTENT_WEBMET, NOAA_GOES16_PATH
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


def graficar_geocolor(fecha: datetime.datetime, _canales: dict, datos: str, mapas_out: str, extent=EXTENT):
    
    c1 = _canales['C01']
    c2 = _canales['C02']
    c3 = _canales['C03']
    path_fecha = fecha.strftime('%Y_%m/%d/')
    files_1 = sorted(os.listdir(datos + '/' + path_fecha + c1.codigo))
    files_1 = [x for x in files_1 if x[-3:] == '.nc']
    files_2 = sorted(os.listdir(datos + '/' + path_fecha + c2.codigo))
    files_2 = [x for x in files_2 if x[-3:] == '.nc']
    files_3 = sorted(os.listdir(datos + '/' + path_fecha + c3.codigo))
    files_3 = [x for x in files_3 if x[-3:] == '.nc']
    files = zip(files_1, files_2, files_3)
    path_imagenes = f"{mapas_out}/{path_fecha}GeoColor"
    os.system('mkdir -p ' + path_imagenes)
    for file in files:
        
        metadatos_1 = get_metadatos(f"{NOAA_GOES16_PATH}/{path_fecha}{c1.codigo}/{file[0]}")
        metadatos_2 = get_metadatos(f"{NOAA_GOES16_PATH}/{path_fecha}{c2.codigo}/{file[1]}")
        metadatos_3 = get_metadatos(f"{NOAA_GOES16_PATH}/{path_fecha}{c3.codigo}/{file[2]}")
        if metadatos_1['fecha_img'] != metadatos_2['fecha_img'] or metadatos_1['fecha_img'] != metadatos_3['fecha_img']:
            logger.warning('ERROR EN GEOCOLOR, ARCHIVOS NO SINCRONIZADOS')
            continue
        img_base_path = f"GeoColor_ARG{metadatos_1['fecha_img'].strftime('%Y-%m-%d_%H_%M')}_WGS84.png"
        path_imagen = Path(f"{path_imagenes}/{img_base_path}")
        if path_imagen.is_file():
            continue
        print('Generando imagen para %s' % file[0])

        nombre_1 = f"{NOAA_GOES16_PATH}/{path_fecha}{c1.codigo}/{file[0]}"
        nombre_2 = f"{NOAA_GOES16_PATH}/{path_fecha}{c2.codigo}/{file[1]}"
        nombre_3 = f"{NOAA_GOES16_PATH}/{path_fecha}{c3.codigo}/{file[2]}"

        start = t.time()

        def gdal_array(nombre, canal, _metadatos):
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
            _proj_str = f"+proj=geos +h={h} +a={a} +b={b} +f={f} lat_0={lat_0} +lon_0={lon_0} +sweep=x +no_defs"
            source_prj.ImportFromProj4(_proj_str)
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
            array = array1 * _metadatos['scale'] + _metadatos['offset']

            grid.GetRasterBand(1).SetNoDataValue(-1)
            grid.GetRasterBand(1).WriteArray(array)

            # print(array.min(),array.max())
            # minimo = array.min()
            # maximo = array.max()
            # array = (array - minimo) / (maximo - minimo)
            # print(array.min(),array.max())

            array = np.maximum(array, 0.0)
            array = np.minimum(array, 1.0)
            array = np.sqrt(array)
            array = np.maximum(array, 0.07)
            # Close file
            raw = None
            return array

        array_1 = gdal_array(nombre_1, c1, metadatos_1)
        array_2 = gdal_array(nombre_2, c2, metadatos_2)
        array_3 = gdal_array(nombre_3, c3, metadatos_3)

        print(array_1.min(), array_1.max())
        print(array_2.min(), array_2.max())
        print(array_3.min(), array_3.max())

        def rebin(a, shape):
            sh = shape[0], a.shape[0] // shape[0], shape[1], a.shape[1] // shape[1]
            return a.reshape(sh).mean(-1).mean(1)

        array_2 = rebin(array_2, array_1.shape)

        minimo = array_1.min()
        maximo = array_1.max()
        array_1 = (array_1 - minimo) / np.maximum((maximo - minimo) - 0.05, 0.5)
        array_1 = np.maximum(array_1, 0.0)
        array_1 = np.minimum(array_1, 1.0)
        minimo = array_2.min()
        maximo = array_2.max()
        array_2 = (array_2 - minimo) / np.maximum((maximo - minimo) - 0.05, 0.6)
        array_2 = np.maximum(array_2, 0.0)
        array_2 = np.minimum(array_2, 1.0)
        minimo = array_3.min()
        maximo = array_3.max()
        array_3 = (array_3 - minimo) / np.maximum((maximo - minimo) - 0.05, 0.7)
        array_3 = np.maximum(array_3, 0.0)
        array_3 = np.minimum(array_3, 1.0)
        true_green = 0.48358168 * array_2 + 0.45706946 * array_1 + 0.06038137 * array_3
        _geocolor = np.stack([array_2, true_green, array_1], axis=2)
        
        # %% Plot the Data ========================================
        # Create the basemap reference for the Rectangular Projection
        plt.clf()
        plt.figure(figsize=(10, 9.5))
        # Reconfigura el extent del mapa para cartopy
        extent_map = [extent[0], extent[2], 
                      extent[1], extent[3]] # [west, east, south, north]
        # 4326 es WGS84 (LatLOn)
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(extent_map, crs=ccrs.PlateCarree())
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
        # Ploteo la imagen
        ax.imshow(_geocolor, origin='upper', extent=extent_map, transform=crs)
        
        # Redefino la grilla
        parallels=list(np.arange(extent_map[2]//5*5, extent_map[3], 5))
        meridians=list(np.arange(extent_map[0]//5*5, extent_map[1], 5))
        gl.xlocator = mticker.FixedLocator(meridians)
        gl.ylocator = mticker.FixedLocator(parallels)
        gl.xformatter = LongitudeFormatter()
        gl.yformatter = LatitudeFormatter() 

        plt.subplots_adjust(left=0.02, right=0.98, top=1, bottom=0.02)

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

        titulo_negro = " GOES-16 ABI GeoColor %s UTC" % metadatos_1['fecha_img'].strftime('%Y-%m-%d %H:%M')
        institucion = "CONAE-Argentina"
        # Add the image description inside the black rectangle
        # Max lat - Min lat
        lat_difference = (extent[3] - extent[1])
        plt.text(extent[0], extent[1] + lat_difference * 0.005, titulo_negro,
                 horizontalalignment='left', color='white', size=7)
        plt.text(extent[2], extent[1] + lat_difference * 0.005, institucion,
                 horizontalalignment='right', color='white', size=7)

        ax = plt.gca()
        img = image.imread(LOGO)
        plt.figimage(img, 25, 100, zorder=1)
        # grabar a png
        plt.savefig(path_imagen)
        plt.clf()
        plt.close()
        logger.info(f"- finished! Time: {t.time() - start} seconds")
        api_base_path = f"GOES/{path_fecha}GeoColor/{img_base_path}"
        post_img_to_api(api_base_path, metadatos_1['fecha_img'], producto="GeoColor", campo_prod='short_name')
        #exit()
