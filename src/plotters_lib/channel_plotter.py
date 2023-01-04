import os
import logging

from models.canales import Canal
from config.logging_conf import GOES_LOGGER_NAME
from plotters_lib.plot import get_metadatos
from plotters_lib.graficar_canal import generar_imagen

logger = logging.getLogger(GOES_LOGGER_NAME)


def plot_channel(canal: Canal, fecha, _datos, mapas_out):
    path_fecha = fecha.strftime('%Y_%m/%d/')
    files = os.listdir(f"{_datos}/{path_fecha}{canal.codigo}")
    path_imagenes = f"{mapas_out}/{path_fecha}{canal.codigo}"
    os.system('mkdir -p ' + path_imagenes)
    for file in files:
        if not file.endswith(".nc"):
            continue
        metadatos = get_metadatos(f"{_datos}/{path_fecha}{canal.codigo}/{file}")
        if os.path.isfile(path_imagenes + '/C' + str(metadatos['icanal']) + '_ARG' +
                          metadatos['fecha_img'].strftime('%Y-%m-%d_%H_%M') + '_WGS84.png'):
            continue
        logger.info(f'Generando imagen para webmet de {file}')
        img_base_path = f"{path_fecha}{canal.codigo}/{file}"
        generar_imagen(f"{_datos}/{img_base_path}", metadatos, path_imagenes, canal)
