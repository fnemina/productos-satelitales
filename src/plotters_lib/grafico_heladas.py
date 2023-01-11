import datetime
import logging
from pathlib import Path

import cartopy.crs as ccrs
import fiona
import matplotlib.pyplot as plt
import rasterio
import shapely.geometry as sgeom
import matplotlib.patches as mpatches
import matplotlib.image as mplimg
from matplotlib import cm
from matplotlib.colors import ListedColormap
from netCDF4 import Dataset

from config.constants import SHAPEFILES_ADEC_PROV, SHAPEFILES_ADEC_DEP, SHAPEFILES_ADEC_URUGUAY, LOGO_ADEC, \
    EXTENT_ADEC, MONTH_DICT
from config.logging_conf import GOES_LOGGER_NAME
from wrf_api.goes_api_ingest import post_img_to_api

logger = logging.getLogger(GOES_LOGGER_NAME)


def graficar_heladas(lst_nc_path: str, img_path_salida: str, path_heladas_img: Path, lst_geotif_mask_path: Path,
                     img_api_path: str, timestamp: datetime.datetime):
    """
    Grafica severidad de heladas.
    
    imagen: ruta completa del archivo LST (str)
    img_path_salida: ruta donde se escriben los archivos de salida
    fecha_str: string con formato predeterminado con la fecha de la imagen. (str)

    return: fecha y hora de los archivos analizados (str)
    """
    fecha_str = timestamp.strftime('%Y-%m-%d_%H_%M')
    # Apertura del archivo LST
    file1 = Dataset(lst_nc_path)
    # Obtención de los píxeles LST y variables auxiliares
    # data1 = file1.variables['LST'][:, :]
    # data1 = data1 - 273.15

    plt.figure(figsize=(7, 6.66), dpi=160)  # Definición del gráfico

    ax = plt.axes(projection=ccrs.PlateCarree())  # Definición de ejes
    ax.set_extent(EXTENT_ADEC, crs=ccrs.PlateCarree())  # TODO: Revisar
    # Agrego líneas de costa e imagen de fondo 
    # ax.coastlines(resolution='10m', color='black', linewidth=0.8)
    ax.stock_img()
    # Agrego grilla en la figura
    gl = ax.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--', linewidth=0.5)
    # Agrego a la figura el vector de nubes
    shpfile = img_path_salida + "Mascara_nubes_vectorizada_" + fecha_str + ".shp"
    with fiona.open(shpfile) as records:
        geometries = [sgeom.shape(shp['geometry']) for shp in records]
    ax.add_geometries(geometries, ccrs.PlateCarree(), edgecolor='none', facecolor='#D6D5C5', linewidth=0.5)
    # Abro archivo shapefile con departamentos
    shpfile = f"{SHAPEFILES_ADEC_DEP}/Departamentos_filtrados.shp"
    # Agrego a la figura cada uno de los departamentos
    with fiona.open(shpfile) as records:
        geometries = [sgeom.shape(shp['geometry'])
                      for shp in records]
    ax.add_geometries(geometries, ccrs.PlateCarree(), edgecolor='slategrey', facecolor='none', linewidth=0.35)
    # Abro archivo shapefile con departamentos
    shpfile = f"{SHAPEFILES_ADEC_PROV}/Provincia.shp"
    with fiona.open(shpfile) as records:
        geometries = [sgeom.shape(shp['geometry']) for shp in records]
    # Agrego a la figura cada uno de los departamentos
    ax.add_geometries(geometries, ccrs.PlateCarree(), edgecolor='black', facecolor='none', linewidth=0.5)

    # Abro archivo shapefile con departamentos
    shpfile = Path(f"{SHAPEFILES_ADEC_URUGUAY}/Uruguay_departamentos.shp")
    with fiona.open(shpfile) as records:
        geometries = [sgeom.shape(shp['geometry'])
                      for shp in records]
    # Agrego a la figura cada uno de los departamentos
    ax.add_geometries(geometries, ccrs.PlateCarree(),
                      edgecolor='slategrey', facecolor='none', linewidth=0.45)
    # Abro archivo shapefile con departamentos
    shpfile = Path(f"{SHAPEFILES_ADEC_URUGUAY}/Uruguay_limite_internacional.shp")
    with fiona.open(shpfile) as records:
        geometries = [sgeom.shape(shp['geometry'])
                      for shp in records]
    # Agrego a la figura cada uno de los departamentos
    ax.add_geometries(geometries, ccrs.PlateCarree(),
                      edgecolor='black', facecolor='none', linewidth=1.1)

    viridis = cm.get_cmap('viridis', 5)  # Defino paleta de colores: viridis, plasma, inferno, magma
    cmap = ListedColormap(viridis(range(5)))
    cmap.set_over('white')
    cmap.set_under('#AD1457')
    # a56d5d

    # Leo y abro el archivo LST
    with rasterio.open(lst_geotif_mask_path, 'r') as src:
        im = src.read(1, masked=True)

    crs = ccrs.PlateCarree()
    # Ploteo la imagen
    img_plot = ax.imshow(im, vmin=-10.0, vmax=0.0, origin='upper', extent=[-80, -52, -58, -20], cmap=cmap,
                         transform=crs)
    # Opciones de formato
    gl.top_labels = False
    gl.xlabel_style = {'size': 6, 'color': 'black'}
    gl.ylabel_style = {'size': 6, 'color': 'black'}

    # Grafico la barra de colores
    # bounds = [-10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0]
    bounds = [-10.0, -8.0, -6.0, -4.0, -2.0, 0.0]
    plt.colorbar(img_plot,
                 extend='both',
                 ticks=bounds,
                 spacing='uniform',
                 orientation='horizontal',
                 label='Temperatura superficial (C)',
                 pad=0.05, fraction=0.05
                 )

    helada_0 = mpatches.Patch(color='#FFFFFF', label='Sin helada')
    helada_1 = mpatches.Patch(color='#fde725', label='Helada suave')
    helada_2 = mpatches.Patch(color='#5ec962', label='Helada moderada')
    helada_3 = mpatches.Patch(color='#21918c', label='Helada intensa')
    helada_4 = mpatches.Patch(color='#3b528b', label='Helada muy intensa')
    helada_5 = mpatches.Patch(color='#440154', label='Helada severa')
    helada_6 = mpatches.Patch(color='#AD1457', label='Helada extrema')
    # 281139
    nubes = mpatches.Patch(color='#D6D5C5', label='Cobertura nubosa')

    plt.legend(handles=[helada_0, helada_1, helada_2, helada_3, helada_4, helada_5, helada_6, nubes], loc='center left',
               bbox_to_anchor=(1.07, 0.5), shadow='True', facecolor='#ffffb6')

    # Agrego la fecha de la imagen
    # Obtención de fecha del archivo
    add_seconds = int(file1.variables['time_bounds'][0])
    date = datetime.datetime(2000, 1, 1, 12) + datetime.timedelta(seconds=add_seconds)
    date = date - datetime.timedelta(hours=3)
    date_str = date.strftime('%d ' + MONTH_DICT[date.month] + ' %Y %H:%M UTC-3')
    # Agrego títulos
    plt.title('GOES-16 (LST)\nMonitoreo Heladas', fontweight='bold', fontsize=10, loc='left')
    plt.title(date_str + '\nHora local Argentina', fontsize=10, loc='right')

    img_logo = mplimg.imread(LOGO_ADEC)
    ax.figure.figimage(img_logo, xo=ax.figure.bbox.xmax - 100, yo=ax.figure.bbox.ymax - 400, origin='upper', zorder=1)
    #ax.figure.figimage(img_logo, xo=ax.figure.bbox.xmax - 130, yo=ax.figure.bbox.ymax - 380, origin='upper', zorder=1)
    # Guardo la imagen en archivo png
    logger.info(f"Guardando grafico de heladas en la ruta: {path_heladas_img}")
    plt.savefig(path_heladas_img, bbox_inches='tight', pad_inches=0.1, dpi=160)
    plt.close()
    file1.close()
    post_img_to_api(img_api_path, timestamp, producto='heladas')
