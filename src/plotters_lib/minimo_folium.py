import json
import logging
import datetime
from pathlib import Path

import branca.colormap as cm
import fiona
import folium
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from folium import plugins
from rasterstats import zonal_stats
from shapely import geometry as sgeom
from shapely.geometry import shape

from config.constants import SHAPEFILES_ADEC_DEP_POL
from config.logging_conf import GOES_LOGGER_NAME
from plotters_lib.isobands_gdal import isobands
from wrf_api.goes_api_ingest import post_img_to_api

logger = logging.getLogger(GOES_LOGGER_NAME)

ESCALA_HELADAS = cm.linear.viridis.scale(-10.0, 0.0).to_step(5)


def isvalid(geom):
    """ 
    Verfica la validez de la geometría. Devuelve 1 si es válido, 0 si no.
    
    geom: geometría a verificar. 
    """
    try:
        shape(geom)
        return 1
    except Exception:
        return 0


def color_pedanias(x):
    """
    Define los colores para los distintos niveles de heladas según la pedanía.
    
    Si el valor de temperatura es nulo, o nan, devuelve transparente.
    Si el valor es mayor a 0 grados, devuelve blanco. 
    Si pasa las dos condiciones, devuelve el color correspondiente a la escala 'escala_heladas'.
   
    """

    if x['properties']['LST_mean'] is None:
        color = 'transparent'
    elif x['properties']['LST_mean'] > 0.0:
        color = 'white'
    else:
        color = ESCALA_HELADAS(x['properties']['LST_mean'])
    return color


def highlight_function(feature):
    """
    Define los parámetros visuales cuando el cursor se encuentra sobre la geometría (función de destaque visual).
    
    """
    return {
        'fillColor': color_isotermas(feature),
        'color': 'white',
        'weight': 3,
        'dashArray': '5, 5'
    }


def color_isotermas(x):
    """ 
    Define los colores para los distintos niveles de temperaturas según el vector de isotermas generado.
    
    Si el valor de isoterma es nulo, o nan, devuelve transparente.
    Si el valor es menor a 6 grados, devuelve el color del valor mínimo (-10 grados). 
    Si pasa las dos condiciones, devuelve el color correspondiente a la escala 'escala_isotermas'.
   
    """
    if x['properties']['t'] is None:
        color = 'transparent'
    elif x['properties']['t'] <= -10.0:
        color = '#AD1457'
    else:
        color = ESCALA_HELADAS(x['properties']['t'])
    return color


def obtener_estadisticas(shpfile, raster_entrada, archivo_salida, tipo_salida='shp'):
    """
    Calcula estadisticas locales (mínima, máxima y promedio) del raster ingresado en el area del shapefile indicado.
    Las estadísticas pueden ser expresadas en formato shapefile o csv.

    shpfile: ruta del shapefile. (str)
    raster_entrada: ruta del raster a analizar. (str)
    archivo_salida: archivo donde se guardan los resultados (str)
    tipo_salida: tipo del archivo de salida. Shapefile ('shp') o CSV ('csv'). Por defecto 'shp'. (str)
    """

    # Defino path del shapefile, lo abro con librería fiona y guardo en geometries los shapes correspondientes
    with fiona.open(shpfile) as records:
        geometries = [sgeom.shape(shp['geometry'])
                      for shp in records]
    # Calculo las estadísticas de valores mínimos, máximos y promedios
    zs = zonal_stats(geometries, raster_entrada)
    # Abro el shapefile usando la librería GeoPandas
    tabla_shape = gpd.read_file(shpfile)
    # Creo un Dataframe con los valores estadísticos obtenidos
    goesstats_df = pd.DataFrame(zs)
    # Renombro las columnas
    goesstats_df.rename(columns={'min': 'LST_min', 'mean': 'LST_mean', 'max': 'LST_max'}, inplace=True)
    # Concateno las dos tablas
    tabla_shape = pd.concat([tabla_shape, goesstats_df], axis=1)
    # Escribo los resultados al disco
    if tipo_salida == 'csv':
        tabla_shape.drop('geometry', axis=1).to_csv(archivo_salida)
    else:
        tabla_shape.to_file(archivo_salida)
    return


def minimo_folium(path_minimo_shp: str, path_folium_temp_min: Path) -> bool:
    """
    Genera el HTML con las estadísticas del promedio de la mínima registrada por departamento/partido, junto con las 
    isotermas correspondientes.
    
    """

    m = folium.Map(location=[-32.1, -64], zoom_start=6.5, control_scale=True, tiles=None)

    tile = 'https://ide.ign.gob.ar/geoservicios/rest/services/Mapas_IGN/mapa_topografico/MapServer/tile/{z}/{y}/{x}'
    attribute = 'Mapa del <a href="http://www.ign.gob.ar">Instituto Geográfico Nacional</a>, ' + \
                'capa de calles por colaboradores de &copy; <a href="http://openstreetmap.org">OpenStreetMap</a>'
    folium.TileLayer(tiles=tile, attr=attribute, name='IGN ').add_to(m)

    tile = 'http://wms.ign.gob.ar/geoserver/gwc/service/tms/1.0.0/capabaseargenmap@EPSG:3857@png/{z}/{x}/{-y}.png'
    attribute = "Argenmap v2 - Instituto Geográfico Nacional"
    folium.TileLayer(tiles=tile, attr=attribute, name='IGN Argenmap v2').add_to(m)

    folium.TileLayer(tiles='openstreetmap', name='OSM').add_to(m)

    collection = list(fiona.open(path_minimo_shp, 'r'))
    df1 = pd.DataFrame(collection)
    df1['isvalid'] = df1['geometry'].apply(lambda x: isvalid(x))
    df1 = df1[df1['isvalid'] == 1]
    collection = json.loads(df1.to_json(orient='records'))

    # Convert to geodataframe
    geodata_isoterma = gpd.GeoDataFrame.from_features(collection)
    geodata_isoterma.crs = {'init': 'epsg:4326'}
    geodata_isoterma['geoid'] = geodata_isoterma.index.astype(str)

    folium.features.GeoJson(
        geodata_isoterma, name='Isotermas', show=False,
        style_function=lambda x: {
            'fillColor': color_isotermas(x),
            'fillOpacity': 0.65,
            'weight': 0.5,
            'color': color_isotermas(x),
            'opacity': 0
        },
        tooltip=folium.features.GeoJsonTooltip(fields=['t'],
                                               aliases=['Isotermas Mínima [C]'],
                                               labels=True,
                                               sticky=True
                                               ),

        highlight_function=lambda x: {
            'fillColor': color_isotermas(x),
            'fillOpacity': 0.65,
            'color': 'black',
            'weight': 3,
            'dashArray': '5, 5'
        },
        embed=False,
    ).add_to(m)

    geodata = gpd.read_file(f"{path_folium_temp_min.parent}/Estadisticas.shp")
    geodata.crs = {'init': 'epsg:4326'}
    geodata['geoid'] = geodata.index.astype(str)
    geodata['LST_mean'] = geodata['LST_mean'].round(decimals=1)
    geodata['LST_mean'] = geodata['LST_mean'].apply(lambda x: 0 if x == -0 else x)

    heladas = folium.features.GeoJson(
        geodata, name='Estadísticas por departamentos/partidos', show=True,
        style_function=lambda x: {
            'fillColor': color_pedanias(x),
            'fillOpacity': 0.5,
            'weight': 1.75,
            'color': 'darkslategrey',
            'opacity': 0.5
        },
        tooltip=folium.features.GeoJsonTooltip(fields=['NAM', 'LST_mean'],
                                               aliases=['Departamento', 'Temperatura mínima promedio'],
                                               labels=True,
                                               localize=True,
                                               sticky=True
                                               ),
        highlight_function=lambda x: {
            'fillColor': color_pedanias(x),
            'fillOpacity': 0.85,
            'color': 'black',
            'weight': 3,
            'dashArray': '5, 5'
        },
        embed=False,
    ).add_to(m)

    plugins.Fullscreen().add_to(m)
    plugins.LocateControl(auto_start=True).add_to(m)
    ESCALA_HELADAS.add_to(m)
    ESCALA_HELADAS.caption = 'Escala Heladas'
    fmtr = "function(num) {return L.Util.formatNum(num, 2) + ' º ';};"
    plugins.MousePosition(position='topright', separator=' | ', prefix="Mouse:", lat_formatter=fmtr,
                          lng_formatter=fmtr).add_to(m)
    folium.LayerControl().add_to(m)
    plugins.MeasureControl(position='topright', primary_length_unit='meters', secondary_length_unit='miles',
                           primary_area_unit='sqmeters', secondary_area_unit='hectares').add_to(m)
    logger.info("Guardando HTML Folium")
    m.save(str(path_folium_temp_min))
    return True


def graficar_folium_temp_min(path_folium_temp_min: Path, path_temp_min_tif: Path, ts_temp_min: datetime.datetime):
    """

    :param path_folium_temp_min: Path absoluto del archivo HTML Folium
    :param path_temp_min_tif: Path absoluto del archivo GeoTIFF de minimas
    :return:
    """
    bounds = np.linspace(-10, 0, 6).tolist()
    # Abro el tiff de minima generado (temperatura_minima.tif) y lo guardo como numpy en la variable im
    with rasterio.open(path_temp_min_tif, 'r') as src:
        im = src.read(1)
        metadatos = src.meta.copy()
        nodatavalue = metadatos['nodata']
        im[np.where(im == nodatavalue)] = np.nan

    # Defino el shapefile de donde calcular las estadísticas, el tiff de entrada y el shp de salida
    shpfile = Path(f"{SHAPEFILES_ADEC_DEP_POL}/SHAPEFILES_ADEC_DEPARTAMENTOS_POLIGONOS_SIMPLIFICADOS.shp")
    raster_entrada = path_temp_min_tif
    archivo_salida = f"{path_folium_temp_min.parent}/Estadisticas.shp"  # definir la ruta del archivo shapefile
    # donde irian las estadisticas
    obtener_estadisticas(shpfile, raster_entrada, archivo_salida, 'shp')

    # Discretizo el tiff de minima
    minimo_discretizado_indices = np.digitize(im, np.array(bounds))
    min_discretizado = np.zeros_like(minimo_discretizado_indices)
    keys = list(range(0, len(bounds) + 1))
    values = bounds
    values.append(-999)
    new_dict = {k: v for k, v in zip(keys, values)}
    for keys in new_dict:
        min_discretizado[np.where(minimo_discretizado_indices == keys)] = new_dict[keys]
    out_path_discretizada = Path(f"{path_folium_temp_min.parent}/temperatura_minima_discretizada.tif")
    with rasterio.open(out_path_discretizada, 'w', **metadatos) as outf:
        outf.write(min_discretizado.astype('float32'), 1)

    # Genero las isotermas con la funcion isobandas
    in_file = out_path_discretizada
    band = 1
    out_folder_isotermas = str(path_folium_temp_min.parent)  # TODO: Revisar path
    out_format = 'ESRI Shapefile'
    out_shapefile_name = 'Temp_minima_discretizada_vectorizada'  # Nombre del archivo shapefile
    attr_name = 't'
    offset = 0.0
    interval = 2.0
    isobands(str(in_file), band, f"{out_folder_isotermas}/", out_format, out_shapefile_name, attr_name, offset,
             interval)

    # Genero el html con folium especificando el archivo de isotermas, y la ruta de salida
    path_minimo_shp = f"{out_folder_isotermas}/{out_shapefile_name}.shp"
    fol_status = minimo_folium(path_minimo_shp, path_folium_temp_min)
    if fol_status:
        fol_temp_min_prod_name = 'fol_temp_min'
        api_base_path = f"GOES/{ts_temp_min.strftime('%Y_%m/%d')}/{fol_temp_min_prod_name}/{path_folium_temp_min.name}"
        post_img_to_api(api_base_path, ts_temp_min, producto=fol_temp_min_prod_name, campo_prod='short_name')

    path_temp_min_tif.unlink(missing_ok=True)
