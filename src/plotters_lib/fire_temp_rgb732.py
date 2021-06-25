import logging
import os
from datetime import datetime, timedelta
from pathlib import Path

import cartopy.crs as ccrs
import dateutil.parser
import fiona
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.image as mplimg
import numpy as np
import shapely.geometry as sgeom
import xarray
from netCDF4 import Dataset

from config.constants import SHAPEFILES, CMI, LOGO_FIRE
from config.logging_conf import GOES_LOGGER_NAME
from plotters_lib.remap_fire import remap
from wrf_api.goes_api_ingest import post_img_to_api

logger = logging.getLogger(GOES_LOGGER_NAME)


def get_geot(_extent, nlines, ncols):
    """
    Devuelve la matriz de geotransformación pasando como argumento
    la extensión geografica, la cantidad de filas y de columnas.

    extent: extensión geográfica en grados [lowerleftx, lowerlefty, upperrightx, upperrighty].
    nlines, ncols: número de filas (lines) y columnas (rows) (floats).
    return: matriz de geotransformación.
    """
    resx = (_extent[2] - _extent[0]) / ncols
    resy = (_extent[3] - _extent[1]) / nlines
    return [_extent[0], resx, 0, _extent[3], 0, -resy]


def get_files_dict(files_path: Path) -> dict:
    """Retorna un diccionario con la forma {fecha: *.nc}"""
    f_dict = {}
    nc_path_list = sorted(files_path.glob("*.nc"))
    for nc_file in nc_path_list:
        nc_dataset = xarray.open_dataset(nc_file)
        nc_timestamp = dateutil.parser.isoparse(nc_dataset.time_coverage_start)
        nc_dataset.close()
        f_dict[nc_timestamp] = nc_file
    return f_dict


def graficar_fire_temp(fecha, _canales, datos, mapas_out):  # , extent=EXTENT
    canal_07 = _canales['C07']
    canal_02 = _canales['C02']
    canal_03 = _canales['C03']
    path_fecha = fecha.strftime('%Y_%m/%d')
    files_c02: dict = get_files_dict(Path(f"{datos}/{path_fecha}/{canal_02.codigo}/"))
    files_c03: dict = get_files_dict(Path(f"{datos}/{path_fecha}/{canal_03.codigo}/"))
    files_c07: dict = get_files_dict(Path(f"{datos}/{path_fecha}/{canal_07.codigo}/"))
    path_imagenes = f"{mapas_out}/{path_fecha}/fire_rgb732"
    os.system(f"mkdir -p {path_imagenes}")
    for fecha_img in files_c02.keys():
        fire_img_name = f"FireRGB_ARG{fecha_img.strftime('%Y-%m-%d_%H_%M')}_WGS84.png"
        fire_img_path = Path(f"{path_imagenes}/{fire_img_name}")
        if fire_img_path.exists():
            continue
        try:
            blue_path: Path = files_c02[fecha_img]
            green_path: Path = files_c03[fecha_img]
            red_path: Path = files_c07[fecha_img]
        except KeyError:
            continue
        logger.info(f"Generando archivo {fire_img_name}")
        dataset_c03 = Dataset(green_path)
        dataset_c02 = Dataset(blue_path)
        dataset_c07 = Dataset(red_path)
        # Lectura de parámetros
        scale_1 = dataset_c07.variables[CMI].scale_factor
        offset_1 = dataset_c07.variables[CMI].add_offset

        scale_2 = dataset_c03.variables[CMI].scale_factor
        offset_2 = dataset_c03.variables[CMI].add_offset

        scale_3 = dataset_c02.variables[CMI].scale_factor
        offset_3 = dataset_c02.variables[CMI].add_offset
        # Obtención de fecha del archivo
        add_seconds = int(dataset_c07.variables['time_bounds'][0])
        date = datetime(2000, 1, 1, 12) + timedelta(seconds=add_seconds)
        fecha_str = date.strftime('%Y-%m-%d_%H_%M_UTC')
        date = date - timedelta(seconds=3 * 3600)
        date = date.strftime('%d-%m-%Y %H:%M')
        # Definición de extensión de mapas de la República Argentina. SETEO LIMITES dominio OHMC
        w = -75
        e = -53
        n = -25
        s = -42
        # extent = [-75, -42, -53, -25]
        extent = [w, s, e, n]  # Extensión para función remap
        #  extent_map = [-75, -53, -42, -25]
        extent_map = [w, e, s, n]  # Dominio LST

        resolution = 1.
        # Llamo a la función de reproyección con la ruta de imagen, la extensión a la que quiero el TIFF, la resolución,
        # el formato del archivo de entrada y la variable a obtener
        # RED BAND
        grid_1, array_1 = remap(str(red_path), extent, resolution, 'HDF5', CMI)
        # Genero el array de datos final
        array_1 = array_1 * scale_1 + offset_1
        # GREEN BAND
        grid_2, array_2 = remap(str(green_path), extent, resolution, 'HDF5', CMI)
        # Genero el array de datos final
        array_2 = array_2 * scale_2 + offset_2
        # BLUE BAND
        grid_3, array_3 = remap(str(blue_path), extent, resolution, 'HDF5', CMI)
        # Genero el array de datos final
        array_3 = array_3 * scale_3 + offset_3
        # PARA FIRE
        # Apply range limits for each channel (mostly important for Red channel)
        array_1 = np.maximum(array_1, 273)  # 273
        array_1 = np.minimum(array_1, 333)

        array_2 = np.maximum(array_2, 0)
        array_2 = np.minimum(array_2, 0.8)

        array_3 = np.maximum(array_3, 0)
        array_3 = np.minimum(array_3, 0.75)
        # Normalize each channel by the appropriate range of values (again, mostly important for Red channel)
        array_1 = (array_1 - 273) / (333 - 273)
        array_2 = (array_2 - 0) / (0.8 - 0)
        array_3 = (array_3 - 0) / (0.75 - 0)
        # Apply the gamma correction to Red channel.
        #   I was told gamma=0.4, but I get the right answer with gamma=2.5 (the reciprocal of 0.4)
        # R = np.power(R, 9)
        array_1 = np.power(array_1, 5)
        # The final RGB array
        # ACA invento
        geocolor = np.stack([array_1, array_2, array_3], axis=2)
        # Definición del gráfico
        fig = plt.figure(figsize=(7, 7))
        # Definición de ejes
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(extent_map, crs=ccrs.PlateCarree())
        # Agrego líneas de costa
        ax.coastlines(resolution='10m', color='black', linewidth=0.8)
        ax.stock_img()
        # Agrego grilla en la figura
        gl = ax.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--', linewidth=0.1)
        # Abro archivo shapefile con departamentos
        shpfile = f"{SHAPEFILES}/departamentos.shp"
        with fiona.open(shpfile) as records:
            geometries = [sgeom.shape(shp['geometry']) for shp in records]
            # Agrego a la figura cada uno de los departamentos
        ax.add_geometries(geometries, ccrs.PlateCarree(), edgecolor='lightgray', facecolor='none', linewidth=0.1,
                          alpha=0.5)

        shpfile = f"{SHAPEFILES}/008_límites_provinciales.shp"
        with fiona.open(shpfile) as records:
            geometries = [sgeom.shape(shp['geometry']) for shp in records]
            # Agrego a la figura cada uno de los departamentos
        ax.add_geometries(geometries, ccrs.PlateCarree(), edgecolor='lightgray', facecolor='none', linewidth=0.3,
                          alpha=0.5)
        crs = ccrs.PlateCarree()
        # Ploteo la imagen
        ax.imshow(geocolor, origin='upper', extent=extent_map, transform=crs)
        gl.top_labels = False  # Opciones de formato
        plt.title('GOES-16 RGB 732 \nCONAE-Argentina', fontweight='bold', fontsize=9, loc='left')
        plt.title(f"{fecha_str}\n{date} Hora local Argentina", fontsize=9, loc='right')
        nubes = mpatches.Patch(color='#fd6400', label='Incendios activos')
        plt.legend(handles=[nubes], loc='center', bbox_to_anchor=(0.80, 0.05), shadow='True', facecolor='#FFFFFF')
        # logo
        img_logo = mplimg.imread(LOGO_FIRE)
        ax.figure.figimage(img_logo, xo=1212.5, yo=200.0, origin='upper')
        # ax.figure.figimage(img_logo, xo=625.0, yo=1375.0, origin='upper', zorder=1)
        logger.info(f"Guardando fire_rgb732 en {fire_img_path}")
        plt.savefig(fire_img_path, bbox_inches='tight', pad_inches=0.1, dpi=300)
        plt.clf()
        plt.close()
        api_base_path = f"GOES/{path_fecha}/fire_rgb732/{fire_img_name}"
        post_img_to_api(api_base_path, fecha_img, producto="fire_rgb732", campo_prod='short_name')
