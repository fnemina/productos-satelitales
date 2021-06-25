import datetime
import glob
import logging
import os
from pathlib import Path
from typing import List
from datetime import timezone

import dateutil.parser
import xarray

from config.constants import TEMP_MIN_HOUR_UTC, NOAA_GOES16_PATH
from config.logging_conf import GOES_LOGGER_NAME
from plotters_lib.grafico_heladas import graficar_heladas
from plotters_lib.grafico_lst import graficar_lst
from plotters_lib.grafico_minimo import graficar_temperaturas_minimas
from plotters_lib.minimo_folium import graficar_folium_temp_min
from wrf_api.goes_api_ingest import post_lst_geotif_to_api

logger = logging.getLogger(GOES_LOGGER_NAME)


def delete_files(_path: Path):
    files_to_delete = ["*.tif", "*.cpg", "*.shx", "*.dbf", "*.shp", "*.prj", "*.xml"]
    for ext in files_to_delete:
        files_to_delete = glob.glob(f"{_path}/{ext}")
        for f in files_to_delete:
            os.remove(f)


def get_lst_metadata(lst_path_list: List[Path], path_fecha, _datos: str):
    lst_metadata = {}
    for lst_path in lst_path_list:
        lst_dataset = xarray.open_dataset(lst_path)
        lst_timestamp = dateutil.parser.isoparse(lst_dataset.time_coverage_start)
        acm_path_list = sorted(Path(f"{_datos}/{path_fecha}/L2-ACMF/{lst_timestamp.hour:02d}/").glob('**/*.nc'))
        try:
            lst_metadata[str(lst_path)] = {'timestamp': lst_timestamp, 'acm_path': acm_path_list[1]}
        except IndexError:
            pass
        lst_dataset.close()
    return lst_metadata


def plot_lst_heladas(fecha: datetime.datetime, _datos: str, mapas_out: str):
    path_fecha = fecha.strftime('%Y_%m/%d')
    path_geotiff_fecha = fecha.strftime('%Y-%m-%d')
    lst_nc_files = sorted(Path(f"{_datos}/{path_fecha}/L2-LSTF/").glob('**/*.nc'))
    lst_metadata = get_lst_metadata(lst_nc_files, path_fecha, NOAA_GOES16_PATH)
    logger.info(f"Generando grafico de LST y heladas. Cantidad de archivos LST = {len(lst_metadata.keys())}")
    heladas = 'heladas'
    lst = 'LST'
    lst_geotif = 'lst_geotif'
    temperaturas_minimas = 'temp_min'
    folium_temperaturas_minimas = 'fol_temp_min'
    # LST
    lst_img_path = Path(f"{mapas_out}/{path_fecha}/{lst}")
    lst_img_path.mkdir(mode=0o775, parents=True, exist_ok=True)
    # heladas
    heladas_img_path = Path(f"{mapas_out}/{path_fecha}/{heladas}")
    heladas_img_path.mkdir(mode=0o775, parents=True, exist_ok=True)
    # GeoTIFF LST
    lst_geotif_path = Path(f"{mapas_out}/{path_geotiff_fecha}/{lst_geotif}")
    lst_geotif_path.mkdir(mode=0o775, parents=True, exist_ok=True)
    lst_nc_path: str
    lst_meta: dict
    ploted_flag = False
    for lst_nc_path, lst_meta in lst_metadata.items():
        # timestamp
        timestamp: datetime.datetime = lst_meta['timestamp']
        f_str: str = timestamp.strftime('%Y-%m-%d_%H_%M')
        timestamp_path_str: str = timestamp.strftime('%Y_%m/%d')
        # LST
        lst_img_name: str = f"{lst}_ARG{f_str}_WGS84.png"
        path_lst_img = Path(f"{lst_img_path}/{lst_img_name}")
        # heladas
        heladas_img_name: str = f"{heladas}_ARG{f_str}_WGS84.png"
        path_heladas_img = Path(f"{heladas_img_path}/{heladas_img_name}")
        # lst_geotif
        lst_geotif_mask_name: str = f"LST_GOES_4326_Argentina_enmascarada_{f_str}.tif"
        lst_geotif_mask_path = Path(f"{lst_geotif_path}/{lst_geotif_mask_name}")
        if not path_lst_img.exists():
            ploted_flag = True
            img_api_path = f"GOES/{timestamp_path_str}/{lst}/{lst_img_name}"
            graficar_lst(lst_nc_path, str(lst_meta['acm_path']), f"{lst_img_path}/", lst_geotif_mask_path, timestamp,
                         path_lst_img, img_api_path)
            # publicar geotif en la API
            if lst_geotif_mask_path.exists() and path_lst_img.exists():
                lst_geotif_mask_api_path = f"GOES/{path_geotiff_fecha}/{lst_geotif}/{lst_geotif_mask_name}"
                post_lst_geotif_to_api(lst_geotif_mask_api_path, timestamp)
        if not path_heladas_img.exists():
            ploted_flag = True
            img_api_path = f"GOES/{timestamp_path_str}/{heladas}/{heladas_img_name}"
            graficar_heladas(lst_nc_path, f"{lst_img_path}/", path_heladas_img, lst_geotif_mask_path, img_api_path,
                             timestamp)
    if ploted_flag:
        delete_files(lst_img_path)
    # Temperaturas minimas
    temp_min_img_path = Path(f"{mapas_out}/{path_fecha}/{temperaturas_minimas}")
    temp_min_img_path.mkdir(mode=0o775, parents=True, exist_ok=True)
    # Temperaturas minimas Folium
    folium_temp_min_path = Path(f"{mapas_out}/{path_fecha}/{folium_temperaturas_minimas}")
    folium_temp_min_path.mkdir(mode=0o775, parents=True, exist_ok=True)
    ts_now = datetime.datetime.now(tz=timezone.utc)
    ts_temp_min = datetime.datetime(year=fecha.year, month=fecha.month, day=fecha.day, hour=TEMP_MIN_HOUR_UTC,
                                    minute=0, tzinfo=timezone.utc)
    # Temperaturas minimas
    temp_min_img_name = f"{temperaturas_minimas}_ARG{ts_temp_min.strftime('%Y-%m-%d_%H_%M')}.png"
    path_temp_min_img = Path(f"{temp_min_img_path}/{temp_min_img_name}")
    # Temperaturas minimas Folium
    folium_temp_min_name = f"{folium_temperaturas_minimas}_ARG{ts_temp_min.strftime('%Y-%m-%d_%H_%M')}.html"
    path_folium_temp_min = Path(f"{folium_temp_min_path}/{folium_temp_min_name}")
    if ts_now.hour >= TEMP_MIN_HOUR_UTC and (not path_temp_min_img.exists()):
        lst_geotiff_list = sorted(lst_geotif_path.glob("*.tif"))
        logger.info(f"Se genera grafico de minimas -> {len(lst_geotiff_list)} elementos en la lista geotiff")
        path_temp_min_tif = graficar_temperaturas_minimas(lst_geotiff_list, path_temp_min_img, ts_temp_min)
        graficar_folium_temp_min(path_folium_temp_min, path_temp_min_tif, ts_temp_min)
        delete_files(folium_temp_min_path)
    return
