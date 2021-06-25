import datetime
import logging
from pathlib import Path

import cartopy.crs as ccrs
import fiona
import geopandas as gpd
import matplotlib.image as mplimg
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import rasterio
import rasterio as rio
import shapely.geometry as sgeom
from matplotlib.colors import ListedColormap
from netCDF4 import Dataset
from osgeo import osr, gdal
from rasterio.features import shapes
from rasterio.mask import mask as mascara

from config.constants import SHAPEFILES_ADEC_DEP, SHAPEFILES_ADEC_PROV, SHAPEFILES_ADEC_URUGUAY, LOGO_ADEC, \
    EXTENT_ADEC, MONTH_DICT
from config.logging_conf import GOES_LOGGER_NAME
from plotters_lib.remap import remap
from wrf_api.goes_api_ingest import post_img_to_api

logger = logging.getLogger(GOES_LOGGER_NAME)


def get_geo_t(extent, nlines, ncols):
    """
    Devuelve la matriz de geotransformación pasando como argumento
    la extensión geografica, la cantidad de filas y de columnas.

    extent: extensión geográfica en grados [lowerleftx, lowerlefty, upperrightx, upperrighty].
    nlines, ncols: número de filas (lines) y columnas (rows) (floats).
    return: matriz de geotransformación.
    """
    resx = (extent[2] - extent[0]) / ncols
    resy = (extent[3] - extent[1]) / nlines
    return [extent[0], resx, 0, extent[3], 0, -resy]


def graficar_lst(lst_nc_path: str, bcm_nc_path: str, img_path_salida: str, lst_geotif_mask_path: Path,
                 timestamp: datetime.datetime, path_lst_img: Path, img_api_path: str):
    """
    Función principal del programa para temperatura superficial de suelo (LST) y heladas
    usando imágenes GOES-16. Esta función recibe como argumentos la dirección del archivo LST
    a analizar junto con su correspondiente archivo de BCM (Binary Cloud Mask) para enmascarar
    los datos correspondientes.


    Para el cálculo de estadísticas como así la generación de imágenes, se debe especificar
    el path de los shapefiles correspondientes a:

        Departamentos Argentina (EPSG: 4326) (.shp)
        Provincias Argentina (EPSG: 4326) (.shp)

    Genera en la ruta especificada por el argumento 'img_path_salida' los siguientes productos:
        LST Argentina geostacionario (.tif)
        BCM Argentina geostacionario (.tif)
        LST Argentina 4326 (.tif)
        BCM Argentina 4326 (.tif)
        LST Argentina con máscara de nubes 4326 (.tif)

        Shapefile BCM Argentina 4326 (.shp)

        Png LST Cordoba (.png)
        Png Heladas Cordoba (.png)


    lst_nc_path: ruta completa del archivo LST (str)
    bcm_nc_path: ruta del archivo BCM (str)
    img_path_salida: ruta donde se escriben los archivos de salida
    return: fecha y hora de los archivos analizados (str)

    """
    # Lectura del archivo y obtención del array de LST
    f_str = timestamp.strftime('%Y-%m-%d_%H_%M')
    var = "LST"
    file1 = Dataset(lst_nc_path)
    data1 = file1.variables[var][:, :]

    # Obtención de las coordenadas (x,y) en proyección geoestacionaria igual
    # a los valores de ángulo de escaneo en radianes multiplicados por la altura del satélite.
    sat_h = file1.variables['goes_imager_projection'].perspective_point_height
    x1 = file1.variables['x_image_bounds'][0] * sat_h
    x2 = file1.variables['x_image_bounds'][1] * sat_h
    y1 = file1.variables['y_image_bounds'][1] * sat_h
    y2 = file1.variables['y_image_bounds'][0] * sat_h
    goes_extent = [x1, y1, x2, y2]

    # Lectura de parámetros (escala, offset y longitud central)
    scale = file1.variables[var].scale_factor
    offset = file1.variables[var].add_offset

    # Obtención de fecha del archivo
    add_seconds = int(file1.variables['time_bounds'][0])
    date = datetime.datetime(2000, 1, 1, 12) + datetime.timedelta(seconds=add_seconds)

    date = date - datetime.timedelta(seconds=3 * 3600)
    date_str = date.strftime('%d ' + MONTH_DICT[date.month] + ' %Y %H:%M UTC-3')
    # Definición de extensión de mapas de la República Argentina
    extent = [-80, -58, -52, -20]  # Extensión para función remap
    # extent_map = [-80, -52, -58, -20]  # Extensión para Cartopy
    # % Producción de imagen TIFF a partir del dato NETCDF4 en proyección geostacionaria
    # Obtengo las dimensiones del dato
    nx = data1.data.shape[0]
    ny = data1.data.shape[1]
    # Genero una imagen GeoTiff con los datos de filas y columnas
    dst_ds = gdal.GetDriverByName('GTiff').Create(
        img_path_salida + 'LST_GOES_proy_geoestacionaria_' + f_str + '.tif', ny, nx, 1, gdal.GDT_Float32)
    # Establezco la referencia Espacial a la proyección geoestacionaria
    # con los parámetros de proyección correspondientes
    srs = osr.SpatialReference()
    srs.ImportFromProj4('+proj=geos +h=35786023.0 +a=6378137.0 +b=6356752.31414 +f=0.00335281068119356027 '
                        '+lat_0=0.0 +lon_0=-75.0 +sweep=x +no_defs')
    # Seteo la referencia producida al archivo GeoTiff
    dst_ds.SetProjection(srs.ExportToWkt())  # export coords to file
    # Obtengo la matriz de geotransformación y la seteo a la imagen tiff generada
    dst_ds.SetGeoTransform(get_geo_t(goes_extent, ny, nx))
    # Guardo en el raster creado la información de LST y escribo al disco
    dst_ds.GetRasterBand(1).WriteArray(data1)
    dst_ds.FlushCache()
    dst_ds = None  # Cierro el archivo Geotiff
    file1.close()  # Cierro el dataset
    #  Producción de imagen TIFF a partir del dato NETCDF4 a la extensión de Argentina en proyección geostacionaria
    # Elijo la resolución de la imagen de salida (en kilómetros)
    # Opcion None para misma dimensión del array que la imagen original (posibilidad de pixeles no cuadrados)
    resolution = 10.0
    # Llamo a la función de reproyección con la ruta de imagen, la extensión a la que quiero el TIFF, la resolución,
    #  el formato del archivo de entrada y la variable a obtener
    grid, array = remap(lst_nc_path, extent, resolution, 'HDF5', var)
    # Genero el array de datos final
    array = array.astype(np.uint16)
    array = array * scale + offset - 273.15
    # Seteo a 255 todos los pixeles que compartan el mismo valor que el píxel [0,0]
    # Este píxel por la extensión elegida siempre será agua y 255 es un valor proxy
    # utilizado después para definir los valores 'nodata'
    array[array == array[0, 0]] = 255
    # Correción empírica de posición para la extensión elegida (una fila para abajo, dos columnas hacia la derecha)
    array = np.roll(array, 1)
    array = np.roll(array, 2, axis=0)
    array[:, 0] = 255

    # %  Producción de imagen LST en formato TIFF a la extensión de Argentina en proyección 4326

    # Obtengo las dimensiones del dato
    nx = np.size(array, 0)
    ny = np.size(array, 1)
    # Genero una imagen GeoTiff con los datos de filas y columnas
    dst_ds = gdal.GetDriverByName('GTiff').Create(img_path_salida + 'LST_GOES_4326_Argentina_' + f_str + '.tif', ny,
                                                  nx, 1, gdal.GDT_Float32)
    # Obtengo la matriz de geotransformación y la seteo a la imagen tiff generada
    dst_ds.SetGeoTransform(get_geo_t(extent, nx, ny))
    # Establezco la referencia Espacial a la proyección WGS-84
    srs = osr.SpatialReference()
    srs.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    # Seteo la referencia producida al archivo GeoTiff
    dst_ds.SetProjection(srs.ExportToWkt())
    # Guardo en el raster creado la información de LST y escribo al disco
    dst_ds.GetRasterBand(1).WriteArray(array)
    dst_ds.FlushCache()
    dst_ds = None  # Cierro el archivo Geotiff

    # % Establezco los valores de nodata en el tiff final

    # Abro la imagen tiff producida con la librería Rasterio
    with rasterio.open(img_path_salida + 'LST_GOES_4326_Argentina_' + f_str + '.tif') as src:
        array_lst_argentina = src.read(1)
        # print(src.meta)
        lin_meta = src.meta.copy()

    nodatavalue = -999.0  # Establezco el valor de nodata como -999.0
    array_lst_argentina[array_lst_argentina > 254] = nodatavalue  # Seteo a nodata todos los valores superiores a 254
    lin_meta.update({'nodata': nodatavalue})  # Actualizo el metadato con el valor de nodata establecido

    # Defino ruta de salida y escribo el nuevo metadato con la nueva matriz con valores de nodata actualizados
    out_path_lin = img_path_salida + 'LST_GOES_4326_Argentina_' + f_str + '.tif'
    with rio.open(out_path_lin, 'w', **lin_meta) as outf:
        outf.write(array_lst_argentina, 1)

    # % Cargo el producto nubes para obtener el producto vectorizado para aplicar la máscara
    # Mismo procedimiento que LST

    imagen_nubes = bcm_nc_path
    var = 'BCM'

    file1 = Dataset(imagen_nubes)
    data1 = file1.variables['BCM'][:, :]
    sat_h = file1.variables['goes_imager_projection'].perspective_point_height
    x1 = np.float64(file1.variables['x'][:]) * np.float64(sat_h)
    y1 = np.float64(file1.variables['x'][:]) * np.float64(sat_h)

    # Exporto la imagen entera de BCM (producto binario de nubes) a TIFF
    nx = data1.data.shape[0]
    ny = data1.data.shape[1]

    dst_ds = gdal.GetDriverByName('GTiff').Create(img_path_salida + 'BCM_GOES_geoestacionaria_' + f_str + '.tif',
                                                  ny, nx, 1, gdal.GDT_Byte)
    dst_ds.SetGeoTransform(get_geo_t(goes_extent, ny, nx))
    srs = osr.SpatialReference()
    srs.ImportFromProj4('+proj=geos +h=35786023.0 +a=6378137.0 +b=6356752.31414 +f=0.00335281068119356027 '
                        '+lat_0=0.0 +lon_0=-75.0 +sweep=x +no_defs')
    dst_ds.SetProjection(srs.ExportToWkt())
    dst_ds.GetRasterBand(1).WriteArray(data1)
    dst_ds.FlushCache()
    dst_ds = None
    file1.close()

    # Reproyecto la imagen a la extensión definida a Argentina de la misma forma que LST

    resolution = 2.0
    grid, array = remap(imagen_nubes, extent, resolution, 'HDF5', var)
    array = array.astype(np.uint8)

    nx = np.size(array, 0)
    ny = np.size(array, 1)

    dst_ds = gdal.GetDriverByName('GTiff').Create(img_path_salida + 'BCM_GOES_4326_Argentina_' + f_str + '.tif', ny,
                                                  nx, 1, gdal.GDT_Byte)
    dst_ds.SetGeoTransform(get_geo_t(extent, nx, ny))
    srs = osr.SpatialReference()
    srs.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    dst_ds.SetProjection(srs.ExportToWkt())
    dst_ds.GetRasterBand(1).WriteArray(array)
    dst_ds.FlushCache()
    dst_ds = None

    # Abro la imagen y defino los valores No Data
    with rasterio.open(img_path_salida + 'BCM_GOES_4326_Argentina_' + f_str + '.tif') as src:
        array = src.read(1)
        # print(src.meta)
        lin_meta = src.meta.copy()

    nodatavalue = 254
    array[array > 250] = nodatavalue
    lin_meta.update({'nodata': nodatavalue})

    out_path_lin = img_path_salida + 'BCM_GOES_4326_Argentina_' + f_str + '.tif'
    with rasterio.open(out_path_lin, 'w', **lin_meta) as outf:
        outf.write(array, 1)
    # Abro la imagen de nubes
    with rasterio.open(img_path_salida + 'BCM_GOES_4326_Argentina_' + f_str + '.tif') as src:
        image = src.read(1)  # Leo la primera banda
        mask_value = 1  # Defino la mascara de acuerdo al valor en el raster
        if mask_value is not None:
            mask = image == mask_value
        else:  # Si no se define un valor, no se define la mascara
            mask = None
        # Se enumeran los resultados de la vectorización del raster según valores
        results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for i, (s, v) in enumerate(shapes(image, mask=mask, transform=src.transform)))
        # Se crea una lista de las geometrías obtenidas de la vectorización realizada en el paso anterior
    geoms = list(results)
    # Se crea un DataFrame a partir de la lista de geometrías
    gpd_polygonized_raster = gpd.GeoDataFrame.from_features(geoms)
    # Se exporta la vectorización de las nubes a un shapefile
    gpd_polygonized_raster.to_file(img_path_salida + "Mascara_nubes_vectorizada_" + f_str + ".shp")
    # Cargo todos las instancias de shape de la vectorización de las nubes
    with fiona.open(img_path_salida + "Mascara_nubes_vectorizada_" + f_str + ".shp", "r") as shapefile:
        shapes_nubes = [feature["geometry"] for feature in shapefile]
    # Abro la imagen original y genero una imagen con la mascara con los shapes de nubes
    with rasterio.open(img_path_salida + 'LST_GOES_4326_Argentina_' + f_str + '.tif') as src:
        out_image, out_transform = mascara(src, shapes_nubes, crop=False, invert=True)
        out_meta = src.meta
        mask_value = None
    # Exporto la imagen con la máscara aplicada
    with rasterio.open(lst_geotif_mask_path, "w", **out_meta) as \
            dest:
        dest.write(out_image)
    # Definición de extensión de mapas de la provincia de Córdoba
    # Gráfico Córdoba
    file1 = Dataset(lst_nc_path)
    # Obtención de los píxeles
    data1 = file1.variables['LST'][:, :]
    data1 = data1 - 273.15
    # Definición del gráfico
    plt.figure(figsize=(7, 6.66), dpi=160)
    # Definición de ejes
    ax = plt.axes(projection=ccrs.PlateCarree())  # , globe=globe))
    ax.set_extent(EXTENT_ADEC, crs=ccrs.PlateCarree())
    # Agrego líneas de costa e imagen de fondo
    ax.stock_img()
    # Agrego grilla en la figura
    gl = ax.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--', linewidth=0.5)
    # Abro archivo shapefile con nubes
    shpfile = img_path_salida + "Mascara_nubes_vectorizada_" + f_str + ".shp"
    with fiona.open(shpfile) as records:
        geometries = [sgeom.shape(shp['geometry'])
                      for shp in records]
    # Agrego a la figura cada uno de los vectores nubes
    ax.add_geometries(geometries, ccrs.PlateCarree(), edgecolor='none', facecolor='#D6D5C5', linewidth=0.5)
    # Abro archivo shapefile con departamentos
    shapefile_dep_path = Path(f"{SHAPEFILES_ADEC_DEP}/Departamentos_filtrados.shp")
    if not shapefile_dep_path.exists():
        raise FileNotFoundError
    shpfile = f"{SHAPEFILES_ADEC_DEP}/Departamentos_filtrados.shp"  # ToDo: Revisar
    with fiona.open(shpfile) as records:
        geometries = [sgeom.shape(shp['geometry']) for shp in records]
    # Agrego a la figura cada uno de los departamentos
    ax.add_geometries(geometries, ccrs.PlateCarree(), edgecolor='slategrey', facecolor='none', linewidth=0.45)
    # Abro archivo shapefile con departamentos
    shapefile_prov_path = Path(f"{SHAPEFILES_ADEC_PROV}/Provincia.shp")
    if not shapefile_prov_path.exists():
        raise FileNotFoundError
    shpfile = f"{SHAPEFILES_ADEC_PROV}/Provincia.shp"
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

    # Abro el archivo raster de temperatura
    with rasterio.open(lst_geotif_mask_path, 'r') as src:
        im = src.read(1, masked=True)
    # Ploteo la imagen
    crs = ccrs.PlateCarree()

    escala_color = ['#313695', '#649ac7', '#bde2ee', '#ffffbf',
                    '#fff1aa', '#fee89b', '#fedd8d', '#fece7f', '#fdbf71',
                    '#fdb062', '#fb9c59', '#f88950', '#f57547', '#ef623e', '#e64f35',
                    '#dd3d2d', '#d22c27', '#c31d27', '#b40f26', '#a50026']
    cmap = ListedColormap(escala_color)
    cmap.set_over('#800000')
    cmap.set_under('#281139')

    img_plot = ax.imshow(im, vmin=-9, vmax=51, origin='upper', extent=[-80, -52, -58, -20], cmap=cmap, transform=crs)
    # Opciones de formato
    gl.top_labels = False
    gl.xlabel_style = {'size': 6, 'color': 'black'}
    gl.ylabel_style = {'size': 6, 'color': 'black'}
    # Agrego una barra de color

    bounds = np.linspace(-9, 51, 21).tolist()
    plt.colorbar(img_plot,
                 extend='both',
                 ticks=bounds,
                 spacing='uniform',
                 orientation='horizontal',
                 label='Temperatura superficial (C)',
                 pad=0.05, fraction=0.05
                 )
    plt.title('GOES-16 (LST)', fontweight='bold', fontsize=10, loc='left')  # Agrego un título
    plt.title(date_str + '\nHora local Argentina', fontsize=10, loc='right')

    nubes = mpatches.Patch(color='#D6D5C5', label='Cobertura nubosa')

    plt.legend(handles=[nubes], loc='center left',
               bbox_to_anchor=(1.07, 0.5), shadow='True', facecolor='#ffffb6')

    img_logo = mplimg.imread(LOGO_ADEC)
    ax.figure.figimage(img_logo, xo=ax.figure.bbox.xmax - 130, yo=ax.figure.bbox.ymax - 380, origin='upper', zorder=1)
    plt.savefig(path_lst_img, bbox_inches='tight', pad_inches=0.1, dpi=160)  # Guardo la imagen en archivo png
    plt.close()
    file1.close()
    post_img_to_api(img_api_path, timestamp, producto='LST')
