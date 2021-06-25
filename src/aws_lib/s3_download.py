import logging
import datetime
from pathlib import Path

import boto3
from botocore import UNSIGNED
from botocore.config import Config

from config.constants import NOAA_GOES16, NOAA_GOES16_PATH
from config.logging_conf import GOES_LOGGER_NAME


logger = logging.getLogger(GOES_LOGGER_NAME)


def descargar_canal_boto(dict_disponibles: dict, fecha: datetime.datetime, codigo, datos: str, h_sub_dir=False):
    logger.info(f'{len(dict_disponibles.keys())} links para el canal {codigo}')
    s3_resource = boto3.resource('s3', config=Config(signature_version=UNSIGNED))
    path_descarga = Path(f"{datos}/{fecha.strftime('%Y_%m/%d')}/{codigo}/")
    path_descarga.mkdir(parents=True, exist_ok=True)
    for aws_file, meta in dict_disponibles.items():
        aws_file_path = Path(aws_file)
        if not h_sub_dir:
            file_path = Path(f"{path_descarga}/{aws_file_path.name}")
        else:
            ts_hour = f"{meta['timestamp'].hour:02d}"
            sub_path = Path(f"{path_descarga}/{ts_hour}")
            sub_path.mkdir(parents=True, exist_ok=True)
            file_path = Path(f"{sub_path}/{aws_file_path.name}")
        if file_path.is_file():
            if file_path.stat().st_size == meta['size']:
                continue
        else:
            logger.info(f'Descargando {aws_file}')
            try:
                s3_resource.meta.client.download_file(NOAA_GOES16, aws_file, str(file_path))
            except Exception:
                logger.exception("Error al descargar archivo")
                continue
    return
