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
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
# from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from osgeo import osr, gdal
from PIL import Image
from plotters_lib.cpt_convert import load_cpt

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


def get_geo_t(extent, nlines, ncols):
    # Compute resolution based on data dimension
    resx = (extent[2] - extent[0]) / ncols
    resy = (extent[3] - extent[1]) / nlines
    return [extent[0], resx, 0, extent[3], 0, -resy]


def get_metadatos(file):
    nc = Dataset(file)
    _metadatos = {'icanal': nc.variables['band_id'][:][0], 'cols': nc.variables['x'].shape[0],
                  'rows': nc.variables['y'].shape[0], 't_0': float(nc.variables['t'][0]),
                  't_start': nc.variables['t'].units[14:33], 'offset': nc.variables['CMI'].add_offset,
                  'scale': nc.variables['CMI'].scale_factor,
                  'proj': nc.variables['goes_imager_projection'].grid_mapping_name,
                  'lat_0': nc.variables['geospatial_lat_lon_extent'].getncattr('geospatial_lat_center'),
                  'lon_0': nc.variables['geospatial_lat_lon_extent'].getncattr('geospatial_lon_center'),
                  'h': nc.variables['goes_imager_projection'].getncattr('perspective_point_height'),
                  'a': nc.variables['goes_imager_projection'].getncattr('semi_major_axis'),
                  'b': nc.variables['goes_imager_projection'].getncattr('semi_minor_axis'),
                  'f': 1 / float(nc.variables['goes_imager_projection'].getncattr('inverse_flattening')),
                  'fecha_img': datetime.datetime.strptime(nc.time_coverage_start[:16], '%Y-%m-%dT%H:%M')}
    nc.close()
    return _metadatos


def print_shapes(
        m, #: Basemap, #TODO: remove basemap
        parallels=np.arange(-90, 90, 5),
        meridians=np.arange(0, 360, 5),
        prov=True,
        dep=True
):
    """
    Dibuja las linas de los mapas.

    Parametros
    ----------
    m : basemap
        basemap al cual se le van a dibujar los shapes
    parallels : np.arrange, opcional
        Arreglo con las coordenadas donde van paralelos
        Por defecto: np.arange(-90,90,5).
    meridians : np.arrange, opcional
        Arreglo con las coordenadas donde van meridianos
        Por defecto: np.arange(0,360,5).
    prov : Bool, opcional
        Si prov=True, dibuja los contornos de las provincias
        Por defecto: True.
    dep : Bool, opcional
        Si dep=True, dibuja los contornos de los departamentos
        Por defecto: True.
    """
    m.drawparallels(parallels, labels=[1, 0, 0, 0], color='#FFFFFF', fontsize=0)
    m.drawmeridians(meridians, labels=[0, 0, 0, 1], color='#FFFFFF', fontsize=0)
    if prov:
        m.readshapefile(
            shapefile=SHAPEFILES + '/provincias',
            name='prov',
            drawbounds=True,
            zorder=None,
            linewidth=0.60,
            color='#FFFFFF',
            antialiased=1
        )
    if dep:
        m.readshapefile(
            shapefile=SHAPEFILES + '/departamentos',
            name='dep',
            drawbounds=True,
            zorder=None,
            linewidth=0.20,
            color='#808080',
            antialiased=1
        )


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


# def generar_imagen(nombre: str, _metadatos: dict, path_imagenes: str, canal, extent=EXTENT):
#     start = t.time()

#     icanal = _metadatos['icanal']

#     # Parametross de calibracion
#     offset = _metadatos['offset']
#     scale = _metadatos['scale']

#     # Parametros de proyeccion
#     lat_0 = str(_metadatos['lat_0'])
#     lon_0 = str(_metadatos['lon_0'])
#     h = str(_metadatos['h'])
#     a = str(_metadatos['a'])
#     b = str(_metadatos['b'])
#     f = str(_metadatos['f'])

#     # %% lectura y extraccion de informacion de la pasada
#     connection_info = 'HDF5:\"' + nombre + '\"://' + canal.variable

#     raw = gdal.Open(connection_info, gdal.GA_ReadOnly)

#     # driver = raw.GetDriver().LongName

#     band = raw.GetRasterBand(1)
#     bandtype = gdal.GetDataTypeName(band.DataType)
#     # print(bandtype)

#     # %% Proyecciones

#     # GOES-16 Spatial Reference System
#     source_prj = osr.SpatialReference()
#     proj_str = f"+proj=geos +h={h} +a={a} +b={b} +f={f} lat_0={lat_0} +lon_0={lon_0} +sweep=x +no_defs"
#     source_prj.ImportFromProj4(proj_str)
#     # Lat/lon WSG84 Spatial Reference System
#     target_prj = osr.SpatialReference()
#     target_prj.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

#     # Setup projection and geo-transformation
#     raw.SetProjection(source_prj.ExportToWkt())
#     raw.SetGeoTransform(get_geo_t(GOES16_EXTENT, raw.RasterYSize, raw.RasterXSize))

#     # Compute grid dimension
#     sizex = int(((extent[2] - extent[0]) * KM_PER_DEGREE) / RESOLUTION)
#     sizey = int(((extent[3] - extent[1]) * KM_PER_DEGREE) / RESOLUTION)

#     # Get memory driver
#     mem_driver = gdal.GetDriverByName('MEM')

#     # Create grid
#     grid = mem_driver.Create('grid', sizex, sizey, 1, gdal.GDT_Float32)

#     # Setup projection and geo-transformation
#     grid.SetProjection(target_prj.ExportToWkt())
#     grid.SetGeoTransform(get_geo_t(extent, grid.RasterYSize, grid.RasterXSize))

#     # Perform the projection/resampling

#     gdal.ReprojectImage(
#         raw,
#         grid,
#         source_prj.ExportToWkt(),
#         target_prj.ExportToWkt(),
#         gdal.GRA_NearestNeighbour,
#         options=['NUM_THREADS=ALL_CPUS']
#     )

#     # Read grid data
#     array1 = grid.ReadAsArray()

#     # Mask fill values (i.e. invalid values)
#     np.ma.masked_where(array1, array1 == -1, False)

#     # %% Calibracion
#     array = array1 * scale + offset

#     grid.GetRasterBand(1).SetNoDataValue(-1)
#     grid.GetRasterBand(1).WriteArray(array)

#     # %% Plot the Data ========================================
#     # Create the basemap reference for the Rectangular Projection
#     plt.clf()
#     plt.figure(figsize=(10, 9.5))

#     # 4326 es WGS84 (LatLOn)
#     bmap = Basemap(resolution='h', llcrnrlon=extent[0], llcrnrlat=extent[1],
#                    urcrnrlon=extent[2], urcrnrlat=extent[3], epsg=4326)

#     # Draw the shapefiles
#     print_shapes(bmap)

#     plt.subplots_adjust(left=0.02, right=0.98, top=1, bottom=0.02)

#     # Converts a CPT file to be used in Python
#     cpt = load_cpt(canal.cptfile)

#     # Makes a linear interpolation
#     cpt_convert = LinearSegmentedColormap('cpt', cpt)

#     # Plot the GOES-16 channel with the converted CPT colors
#     # (you may alter the min and max to match your preference)
#     if canal.visible:
#         bmap.imshow(
#             array,
#             origin='upper',
#             cmap='gray',
#             vmin=0.,
#             vmax=1.
#         )
#     else:
#         temp = array - 273.15
#         bmap.imshow(
#             temp,
#             origin='upper',
#             cmap=cpt_convert,
#             vmin=-90,
#             vmax=50
#         )

#     # Add a black rectangle in the bottom to insert the image description
#     # Max Lon - Min Lon
#     lon_difference = (extent[2] - extent[0])
#     current_axis = plt.gca()
#     current_axis.add_patch(Rectangle(
#         (extent[0], extent[1]),
#         lon_difference,
#         lon_difference * 0.020,
#         alpha=1,
#         zorder=3,
#         facecolor='black'
#     ))

#     titulo_negro = " GOES-16 ABI Canal %02d %s UTC" % (icanal, _metadatos['fecha_img'].strftime('%Y-%m-%d %H:%M'))
#     institucion = "CONAE-Argentina"
#     # Add the image description inside the black rectangle
#     # Max lat - Min lat
#     lat_difference = (extent[3] - extent[1])
#     plt.text(extent[0], extent[1] + lat_difference * 0.005, titulo_negro,
#              horizontalalignment='left', color='white', size=7)
#     plt.text(extent[2], extent[1] + lat_difference * 0.005, institucion,
#              horizontalalignment='right', color='white', size=7)

#     # Insert the colorbar at the right
#     cb = bmap.colorbar(location='bottom', size='2%', pad='1%')
#     # Remove the colorbar outline
#     cb.outline.set_visible(True)
#     # Remove the colorbar ticks
#     cb.ax.tick_params(width=0)
#     # Put the colobar labels inside the colorbar
#     cb.ax.xaxis.set_tick_params(pad=0)
#     # Change the color and size of the colorbar labels
#     cb.ax.tick_params(axis='x', colors='black', labelsize=8)
#     cb.set_label(canal.unidad)

#     ax = plt.gca()
#     img = image.imread(LOGO)
#     plt.figimage(
#         img,
#         25,
#         100,
#         # ax.figure.bbox.xmax - 160,
#         # ax.figure.bbox.ymax - 70,
#         zorder=1
#     )
#     # ax.text(0,
#     #         1.10,
#     #         canal.nombre,
#     #         verticalalignment='top',
#     #         transform=ax.transAxes,
#     #         fontsize=20
#     #         )
#     # ax.text(0,
#     #         1.03,
#     #         metadatos['fecha_img'].strftime('%Y-%m-%d %H:%M') + ' UTC',
#     #         verticalalignment='top',
#     #         transform=ax.transAxes,
#     #         fontsize=12
#     #         )
#     # ax.text(1,
#     #         1.03,
#     #         'GOES-16 ABI Canal %02d' % icanal,
#     #         horizontalalignment='right',
#     #         verticalalignment='top',
#     #         transform=ax.transAxes,
#     #         fontsize=10
#     #         )

#     # grabar a png
#     fecha = _metadatos['fecha_img']
#     fecha_str = fecha.strftime('%Y-%m-%d_%H_%M')
#     path_fecha = fecha.strftime('%Y_%m/%d')
#     path_imagen = f"{path_imagenes}/C{_metadatos['icanal']}_ARG{fecha_str}_WGS84.png"
#     plt.savefig(path_imagen)
#     plt.clf()
#     plt.close()
#     print('- finished! Time:', t.time() - start, 'seconds')
#     img_api_path = f"GOES/{path_fecha}/{canal.codigo}/C{_metadatos['icanal']}_ARG{fecha_str}_WGS84.png"
#     post_img_to_api(img_api_path, fecha, producto=canal.codigo, campo_prod='short_name')


# def generar_imagen_webmet(nombre, _metadatos, path_imagenes, canal, extent=EXTENT_WEBMET):
#     start = t.time()

#     icanal = _metadatos['icanal']

#     # Parametross de calibracion
#     offset = _metadatos['offset']
#     scale = _metadatos['scale']

#     # Parametros de proyeccion
#     lat_0 = str(_metadatos['lat_0'])
#     lon_0 = str(_metadatos['lon_0'])
#     h = str(_metadatos['h'])
#     a = str(_metadatos['a'])
#     b = str(_metadatos['b'])
#     f = str(_metadatos['f'])

#     # %% lectura y extraccion de informacion de la pasada
#     connection_info = 'HDF5:\"' + nombre + '\"://' + canal.variable

#     raw = gdal.Open(connection_info, gdal.GA_ReadOnly)

#     # driver = raw.GetDriver().LongName

#     band = raw.GetRasterBand(1)
#     bandtype = gdal.GetDataTypeName(band.DataType)
#     print(bandtype)

#     # %% Proyecciones

#     # GOES-16 Spatial Reference System
#     source_prj = osr.SpatialReference()
#     proj_str = f"+proj=geos +h={h} +a={a} +b={b} +f={f} lat_0={lat_0} +lon_0={lon_0} +sweep=x +no_defs"
#     source_prj.ImportFromProj4(proj_str)
#     # Lat/lon WSG84 Spatial Reference System https://epsg.io/3857
#     target_prj = osr.SpatialReference()
#     target_prj.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

#     # Setup projection and geo-transformation
#     raw.SetProjection(source_prj.ExportToWkt())
#     raw.SetGeoTransform(get_geo_t(GOES16_EXTENT,
#                                   raw.RasterYSize,
#                                   raw.RasterXSize))

#     # Compute grid dimension
#     sizex = int(((extent[2] - extent[0]) * KM_PER_DEGREE) / RESOLUTION)
#     sizey = int(((extent[3] - extent[1]) * KM_PER_DEGREE) / RESOLUTION)

#     # Get memory driver
#     mem_driver = gdal.GetDriverByName('MEM')

#     # Create grid
#     grid = mem_driver.Create('grid', sizex, sizey, 1, gdal.GDT_Float32)

#     # Setup projection and geo-transformation
#     grid.SetProjection(target_prj.ExportToWkt())
#     grid.SetGeoTransform(get_geo_t(extent, grid.RasterYSize, grid.RasterXSize))

#     # Perform the projection/resampling

#     gdal.ReprojectImage(
#         raw,
#         grid,
#         source_prj.ExportToWkt(),
#         target_prj.ExportToWkt(),
#         gdal.GRA_NearestNeighbour,
#         options=['NUM_THREADS=ALL_CPUS']
#     )

#     # Read grid data
#     array1 = grid.ReadAsArray()

#     # Mask fill values (i.e. invalid values)
#     np.ma.masked_where(array1, array1 == -1, False)

#     # %% Calibracion
#     array = array1 * scale + offset

#     grid.GetRasterBand(1).SetNoDataValue(-1)
#     grid.GetRasterBand(1).WriteArray(array)

#     # %% Plot the Data ========================================
#     # Create the basemap reference for the Rectangular Projection
#     plt.clf()
#     fig = plt.figure(frameon=False)
#     fig.set_size_inches(25. * array.shape[1] / array.shape[0], 25, forward=False)
#     ax = plt.Axes(fig, [0., 0., 1., 1.])
#     ax.set_axis_off()
#     fig.add_axes(ax)

#     # 3857 es WGS84 (LatLOn)
#     bmap = Basemap(resolution='h', llcrnrlon=extent[0], llcrnrlat=extent[1],
#                    urcrnrlon=extent[2], urcrnrlat=extent[3], epsg=3857)

#     # Converts a CPT file to be used in Python
#     cpt = load_cpt(canal.cptfile2)

#     # Makes a linear interpolation
#     cpt_convert = LinearSegmentedColormap('cpt', cpt)

#     # Plot the GOES-16 channel with the converted CPT colors
#     # (you may alter the min and max to match your preference)
#     if canal.visible:
#         bmap.imshow(
#             array,
#             origin='upper',
#             cmap='gray',
#             vmin=0.,
#             vmax=1.
#         )
#     else:
#         temp = array - 273.15
#         temp[temp > -5] = 50
#         bmap.imshow(
#             temp,
#             origin='upper',
#             cmap=cpt_convert,
#             vmin=-90,
#             vmax=50
#         )

#     # grabar a png
#     path_imagen = path_imagenes + '/GOES16' + '_\
# ' + _metadatos['fecha_img'].strftime('%Y%m%dT%H%M%S') + 'Z_C' + str(_metadatos['icanal']) + '.png'
#     plt.savefig(path_imagen, transparent=True)
#     plt.clf()
#     plt.close()
#     # Close file
#     raw = None
#     convertir_negro_transparente(path_imagen)
#     print('- finished! Time:', t.time() - start, 'seconds')


# def geocolor(fecha: datetime.datetime, _canales: dict, datos: str, mapas_out: str, extent=EXTENT):
    
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

        # 4326 es WGS84 (LatLOn)
        bmap = Basemap(resolution='h', llcrnrlon=extent[0], llcrnrlat=extent[1],
                       urcrnrlon=extent[2], urcrnrlat=extent[3], epsg=4326)

        # Draw the shapefiles
        print_shapes(bmap)

        plt.subplots_adjust(left=0.02, right=0.98, top=1, bottom=0.02)
        bmap.imshow(_geocolor, origin='upper')

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


def fire_temperature(fecha, _canales, datos, mapas_out, extent=EXTENT):
    c5 = _canales['C05']
    c6 = _canales['C06']
    c7 = _canales['C07']
    path_fecha = fecha.strftime('%Y_%m/%d/')
    files_1 = sorted(os.listdir(datos + '/GOES/' + path_fecha + c5.codigo))
    files_1 = [x for x in files_1 if x[-3:] == '.nc']
    files_2 = sorted(os.listdir(datos + '/GOES/' + path_fecha + c6.codigo))
    files_2 = [x for x in files_2 if x[-3:] == '.nc']
    files_3 = sorted(os.listdir(datos + '/GOES/' + path_fecha + c7.codigo))
    files_3 = [x for x in files_3 if x[-3:] == '.nc']
    files = zip(files_1, files_2, files_3)
    path_imagenes = mapas_out + '/' + path_fecha + 'FireTemperature'
    os.system('mkdir -p ' + path_imagenes)
    for file in files:
        metadatos_1 = get_metadatos(datos + '/GOES/' + path_fecha + c5.codigo + '/' + file[0])
        metadatos_2 = get_metadatos(datos + '/GOES/' + path_fecha + c6.codigo + '/' + file[1])
        metadatos_3 = get_metadatos(datos + '/GOES/' + path_fecha + c7.codigo + '/' + file[2])
        if metadatos_1['fecha_img'] != metadatos_2['fecha_img'] or metadatos_1['fecha_img'] != metadatos_3['fecha_img']:
            print('ERROR EN FIRE TEMPERATURE, ARCHIVOS NO SINCRONIZADOS')
            continue
        if os.path.isfile(path_imagenes + '/FireTemperature_ARG' + metadatos_1['fecha_img'].strftime('%Y-%m-%d_%H_%M')
                          + '_WGS84.png'):
            continue
        print('Generando imagen para %s' % file[0])

        nombre_1 = datos + '/GOES/' + path_fecha + c5.codigo + '/' + file[0]
        nombre_2 = datos + '/GOES/' + path_fecha + c6.codigo + '/' + file[1]
        nombre_3 = datos + '/GOES/' + path_fecha + c7.codigo + '/' + file[2]

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

            # BORRAR
            print("SizeX = %d, SizeY %d [km]" % (sizex, sizey))
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

            array = np.maximum(array, 0.0)
            array = np.minimum(array, 1.0)
            array = np.sqrt(array)
            array = np.maximum(array, 0.07)
            # Close file
            raw = None
            return array

        array_1 = gdal_array(nombre_1, c5, metadatos_1)
        array_2 = gdal_array(nombre_2, c6, metadatos_2)
        array_3 = gdal_array(nombre_3, c7, metadatos_3)

        print(array_1.min(), array_1.max())
        print(array_2.min(), array_2.max())
        print(array_3.min(), array_3.max())

        def rebin(a, shape):
            sh = shape[0], a.shape[0] // shape[0], shape[1], a.shape[1] // shape[1]
            return a.reshape(sh).mean(-1).mean(1)

        array_2 = rebin(array_2, array_1.shape)
        # Normalize each channel by the appropriate range of values  e.g. R = (R-minimum)/(maximum-minimum)
        R = (array_3 - 273) / (333 - 273)
        G = (array_2 - 0) / (1 - 0)
        B = (array_1 - 0) / (0.75 - 0)

        # Apply range limits for each channel. RGB values must be between 0 and 1
        R = np.clip(R, 0, 1)
        G = np.clip(G, 0, 1)
        B = np.clip(B, 0, 1)

        # Apply the gamma correction to Red channel.
        #   corrected_value = value^(1/gamma)
        gamma = 0.4
        R = np.power(R, 1 / gamma)

        # The final RGB array :)
        RGB = np.dstack([R, G, B])

        plt.clf()
        plt.figure(figsize=(20, 20))
        # plt.figure(figsize=(10, 9.5))

        # 4326 es WGS84 (LatLOn)
        bmap = Basemap(resolution='h', llcrnrlon=extent[0], llcrnrlat=extent[1],
                       urcrnrlon=extent[2], urcrnrlat=extent[3], epsg=4326)

        # Draw the shapefiles
        print_shapes(bmap)

        plt.subplots_adjust(left=0.02, right=0.98, top=1, bottom=0.02)
        bmap.imshow(RGB, origin='upper')

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

        titulo_negro = " GOES-16 ABI Fire %s UTC" % metadatos_1['fecha_img'].strftime('%Y-%m-%d %H:%M')
        institucion = "OHMC - CONAE-Argentina"
        # Add the image description inside the black rectangle
        # Max lat - Min lat
        lat_difference = (extent[3] - extent[1])
        plt.text(extent[0], extent[1] + lat_difference * 0.005, titulo_negro,
                 horizontalalignment='left', color='white', size=7)
        plt.text(extent[2], extent[1] + lat_difference * 0.005, institucion,
                 horizontalalignment='right', color='white', size=7)

        ax = plt.gca()
        img = image.imread(LOGO)
        plt.figimage(
            img,
            25,
            100,
            # ax.figure.bbox.xmax - 160,
            # ax.figure.bbox.ymax - 70,
            zorder=1
        )

        # grabar a png
        path_imagen = path_imagenes + '/FireTemperature_ARG' + metadatos_1['fecha_img'].strftime('%Y-%m-%d_%H_%M') \
                      + '_WGS84.png'
        plt.savefig(path_imagen)
        plt.clf()
        plt.close()
        print('- finished! Time:', t.time() - start, 'seconds')
