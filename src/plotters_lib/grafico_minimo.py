import datetime
import logging
import re
from pathlib import Path
from typing import List

import cartopy.crs as ccrs
import fiona
import matplotlib.image as mplimg
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pytz
import rasterio
import shapely.geometry as sgeom
from matplotlib import cm
from matplotlib.colors import ListedColormap

from config.constants import SHAPEFILES_ADEC_DEP, SHAPEFILES_ADEC_PROV, SHAPEFILES_ADEC_URUGUAY, LOGO_ADEC, \
    EXTENT_ADEC, TEMP_MIN_HOUR_UTC, MONTH_DICT
from config.logging_conf import GOES_LOGGER_NAME
from wrf_api.goes_api_ingest import post_img_to_api

logger = logging.getLogger(GOES_LOGGER_NAME)
UMBRAL_NUBOSIDAD = 50


def get_lst_geotiff_dict(lst_geotiff_list: List[Path]) -> dict:
    """Retorna un diccionario con clave=geotiff_path, valor=timestamp. Filtra los archivos que sean anteriores a la hora
     del grafico.
     """
    geotif_lst_re = re.compile(r"LST_GOES_4326_Argentina_enmascarada_(?P<timestamp>\d{4}-\d{2}-\d{2}_\d{2}_\d{2})\.tif")
    lst_geotiff_dict = {}
    for lst_geotiff_path in lst_geotiff_list:
        m = re.match(geotif_lst_re, lst_geotiff_path.name)
        if m:
            timestamp = datetime.datetime.strptime(m.group('timestamp'), "%Y-%m-%d_%H_%M")
            timestamp: datetime.datetime = pytz.utc.localize(timestamp)
            if 3 <= timestamp.hour <= TEMP_MIN_HOUR_UTC:
                lst_geotiff_dict[lst_geotiff_path] = timestamp
    return lst_geotiff_dict


def calcular_porcentaje_nubosidad(lst_filelist_sorted, nubosidad_porc, nodatavalue):
    """Retorna un numpy array con los porcentajes de pixeles ocupado por nubes en la lista de entrada
    """
    for lst_tiff in lst_filelist_sorted:
        with rasterio.open(lst_tiff, 'r') as src:
            im = src.read(1)
            nubosidad_porc[np.where(im == nodatavalue)] += 1
    nubosidad_porc = nubosidad_porc / float(len(lst_filelist_sorted)) * 100
    return nubosidad_porc


def graficar_temperaturas_minimas(lst_filelist: List[Path], path_temp_min_img: Path,
                                  ts_temp_min: datetime.datetime) -> Path:
    """Grafico de temperaturas minimas a partir de una lista de geotiff(LST)"""
    lst_geotiff_dict = get_lst_geotiff_dict(lst_filelist)
    last_timestamp: datetime.datetime = sorted(lst_geotiff_dict.values(), reverse=True)[0]
    lst_filelist_sorted = sorted(lst_geotiff_dict.keys())
    # logger.info(f"lst_filelist_sorted -> {len(lst_filelist_sorted)} elementos")
    with rasterio.open(lst_filelist_sorted[0], 'r') as src:
        min_lst = src.read(1)
        metadatos = src.meta.copy()
        nodatavalue = metadatos['nodata']
        min_lst[np.where(min_lst == nodatavalue)] = np.nan
        nubosidad_porc = np.zeros_like(min_lst)

    nubosidad_porc = calcular_porcentaje_nubosidad(lst_filelist_sorted, nubosidad_porc, nodatavalue)

    for lst_tiff in lst_filelist_sorted:
        with rasterio.open(lst_tiff, 'r') as src:
            im = src.read(1)
            im[np.where(im == nodatavalue)] = np.nan
        min_lst = np.fmin(min_lst, im)

    min_lst[np.where(min_lst == np.nan)] = nodatavalue
    min_lst[np.where(nubosidad_porc < UMBRAL_NUBOSIDAD)] = nodatavalue
    # Enmascara los valores en porcentaje (0-100) de pixeles ocupados por nubes de la lista de entrada
    path_temp_min_tif = Path(f"{path_temp_min_img.parent}/temperatura_minima.tif")

    with rasterio.open(path_temp_min_tif, 'w', **metadatos) as outf:
        outf.write(min_lst, 1)

    # Definición de extensión
    # Definición del gráfico
    plt.figure(figsize=(7, 6.66), dpi=160)
    # Definición de ejes
    ax = plt.axes(projection=ccrs.PlateCarree())  # , globe=globe))
    ax.set_extent(EXTENT_ADEC, crs=ccrs.PlateCarree())
    # Agrego líneas de costa e imagen de fondo
    ax.coastlines(resolution='10m', color='slategrey', linewidth=1)

    # Agrego grilla en la figura
    gl = ax.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--', linewidth=0.5)
    # Abro archivo shapefile con departamentos
    shpfile = Path(f"{SHAPEFILES_ADEC_DEP}/Departamentos_filtrados.shp")
    with fiona.open(shpfile) as records:
        geometries = [sgeom.shape(shp['geometry'])
                      for shp in records]
    # Agrego a la figura cada uno de los departamentos
    ax.add_geometries(geometries, ccrs.PlateCarree(),
                      edgecolor='slategrey', facecolor='none', linewidth=0.45)

    # Abro archivo shapefile con departamentos
    shpfile = Path(f"{SHAPEFILES_ADEC_PROV}/Provincia.shp")
    with fiona.open(shpfile) as records:
        geometries = [sgeom.shape(shp['geometry'])
                      for shp in records]
    # Agrego a la figura cada uno de los departamentos
    ax.add_geometries(geometries, ccrs.PlateCarree(), edgecolor='black', facecolor='none', linewidth=0.5)

    # Abro archivo shapefile con departamentos
    shpfile = Path(f"{SHAPEFILES_ADEC_URUGUAY}/Uruguay_departamentos.shp")
    with fiona.open(shpfile) as records:
        geometries = [sgeom.shape(shp['geometry'])
                      for shp in records]
    # Agrego a la figura cada uno de los departamentos
    ax.add_geometries(geometries, ccrs.PlateCarree(), edgecolor='slategrey', facecolor='none', linewidth=0.45)

    # Abro archivo shapefile con departamentos
    shpfile = Path(f"{SHAPEFILES_ADEC_URUGUAY}/Uruguay_limite_internacional.shp")
    with fiona.open(shpfile) as records:
        geometries = [sgeom.shape(shp['geometry'])
                      for shp in records]
    # Agrego a la figura cada uno de los departamentos
    ax.add_geometries(geometries, ccrs.PlateCarree(),
                      edgecolor='black', facecolor='none', linewidth=1.1)

    with rasterio.open(path_temp_min_tif, 'r') as src:
        im = src.read(1, masked=True)
    # Ploteo la imagen
    crs = ccrs.PlateCarree()

    viridis = cm.get_cmap('viridis', 5)  # Defino paleta de colores: viridis, plasma, inferno, magma
    cmap=ListedColormap(viridis(range(5)))
    cmap.set_over('white')
    cmap.set_under('#AD1457')
    cmap.set_bad('#D6D5C5')

    img_plot = ax.imshow(im, vmin=-20, vmax=20, origin='upper', extent=[-80, -52, -58, -20], cmap=cmap, transform=crs)
    # Opciones de formato
    gl.top_labels = False
    gl.xlabel_style = {'size': 6, 'color': 'black'}
    gl.ylabel_style = {'size': 6, 'color': 'black'}
    # Agrego una barra de color
    
    bounds = [-10.0, -8.0, -6.0, -4.0, -2.0, 0.0]
    plt.colorbar(img_plot,
                 extend='both',
                 ticks=bounds,
                 spacing='uniform',
                 orientation='horizontal',
                 label='Temperatura superficial (C)',
                 pad=0.05, fraction=0.05
                 )
    # plt.colorbar(img_plot, label='Temperatura superficial (C)', extend='both', orientation='horizontal',
    # pad=0.05, fraction=0.05)
    plt.title('GOES-16 (LST)\nTemperaturas Mínimas', fontweight='bold', fontsize=10, loc='left')  # Agrego un título
    local_time = last_timestamp + datetime.timedelta(hours=-3)
    local_time_str = local_time.strftime('%d ' + MONTH_DICT[local_time.month] + ' %Y %H:%M UTC-3')
    plt.title(f'{local_time_str}\nHora local Argentina', fontsize=10, loc='right')
    nubes = mpatches.Patch(color='#D6D5C5', label='Sin datos')
    plt.legend(handles=[nubes], loc='center left', bbox_to_anchor=(1.07, 0.5), shadow='True', facecolor='#ffffb6')
    img_logo = mplimg.imread(LOGO_ADEC)
    # ax.figure.figimage(img, ax.figure.bbox.xmax + 3400, ax.figure.bbox.ymax +150 , zorder=1)
    ax.figure.figimage(img_logo, xo=ax.figure.bbox.xmax - 130, yo=ax.figure.bbox.ymax - 380, origin='upper', zorder=1)
    logger.info(f"Guardando grafico de temperaturas minimas en {path_temp_min_img}")
    plt.savefig(path_temp_min_img, bbox_inches='tight', pad_inches=0.1, dpi=160)  # Guardo la imagen en archivo png
    plt.close()
    temp_min_prod_name = 'temp_min'
    api_base_path = f"GOES/{ts_temp_min.strftime('%Y_%m/%d')}/{temp_min_prod_name}/{path_temp_min_img.name}"
    post_img_to_api(api_base_path, ts_temp_min, producto=temp_min_prod_name, campo_prod='short_name')
    return path_temp_min_tif
