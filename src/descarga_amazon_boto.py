"""

Descarga de datos de Goes desde AWS.

Paquete necesario:
    - awscli (``sudo apt  install awscli``)

Tener cargado en el entorno `CHANNELS` y `PRODUCTS`.

"""
import datetime
import multiprocessing as mp
import argparse

import boto3
from botocore import UNSIGNED
from botocore.config import Config

from models.canales import canales
from config.constants import PRODUCTS, CHANNELS, NOAA_GOES16, CHANNEL_RE_DICT, FINAL_PRODUCTS, NOAA_GOES16_PATH
from config.logging_conf import GOES_LOGGER_NAME, get_logger_from_config_file
from aws_lib.s3_download import descargar_canal_boto

logger = get_logger_from_config_file(GOES_LOGGER_NAME)


def descargar(fecha_cod):
    """
    Descarga los archivos de GOES que faltan para la fecha dada
    """
    s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    fecha = datetime.datetime.strptime(fecha_cod, '%Y%m%d')
    year = fecha.strftime('%Y')
    day_jul = fecha.strftime('%j')
    pool = mp.Pool()
    for producto in PRODUCTS:
        prefix_res = s3.list_objects_v2(Bucket=NOAA_GOES16, Delimiter='/', Prefix=f'ABI-{producto}/{year}/{day_jul}/')
        prefix_list = []
        for res in prefix_res['CommonPrefixes']:
            prefix_list.append(res['Prefix'])
        dict_disponibles = {canal: {} for canal in CHANNELS}
        file_count = 0
        for prefix in prefix_list:
            objects_list = s3.list_objects_v2(Bucket=NOAA_GOES16, Delimiter='/', Prefix=prefix)
            for content in objects_list['Contents']:
                for canal, p in CHANNEL_RE_DICT.items():
                    if p.match(content['Key']):
                        file_count += 1
                        dict_disponibles[canal][content['Key']] = {'size': content['Size'],
                                                                   'timestamp': content['LastModified']}
        logger.info(f"{file_count} links para el producto {producto}")
        for canal in CHANNELS:
            logger.info(f"Descargando {canales[canal].nombre}")
            pool.apply_async(descargar_canal_boto, (dict_disponibles[canal], fecha, canal, NOAA_GOES16_PATH))

    for final_prd in FINAL_PRODUCTS:
        prefix_res = s3.list_objects_v2(Bucket=NOAA_GOES16, Delimiter='/', Prefix=f'ABI-{final_prd}/{year}/{day_jul}/')
        prefix_list = []
        for res in prefix_res['CommonPrefixes']:
            prefix_list.append(res['Prefix'])
        dict_disponibles = {}
        file_count = 0
        for prefix in prefix_list:
            objects_list = s3.list_objects_v2(Bucket=NOAA_GOES16, Delimiter='/', Prefix=prefix)
            for content in objects_list['Contents']:
                file_count += 1
                dict_disponibles[content['Key']] = {'size': content['Size'], 'timestamp': content['LastModified']}
        logger.info(f"{file_count} links para el producto {final_prd}")
        logger.info(f"Descargando {final_prd}")
        pool.apply_async(descargar_canal_boto, (dict_disponibles, fecha, final_prd, NOAA_GOES16_PATH,
                                                {'h_sub_dir': True}))

    pool.close()
    pool.join()


def main():
    parser = argparse.ArgumentParser(prog="Descargar Amazon")
    parser.add_argument("fecha", help="Día del que se descargan los datos (se descargan en horario UTC)")
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v0.7.0')
    args = parser.parse_args()
    descargar(args.fecha)


if __name__ == "__main__":
    main()
