#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue May 15 11:59:10 2018

@author: Sergio CONAE
@author: Andrés CONAE

modificado para OHMCBA
Lee los DATOS a partir de el archivo de metaDATOS generados previamente con
GenMetadato_Yaku.py
Reproyecta a WGS84, mediante Gdal dentro de Spyder
Grafica en los limites definidos
Elije la paleta en funcion del canal

"""
import argparse
import datetime
import glob
import multiprocessing as mp
import os
import time

import plotters_lib.plot as plot
from models.canales import canales
from plotters_lib.channel_plotter import plot_channel
from config.constants import PRODUCTS, CHANNELS, MAPAS, NOAA_GOES16_PATH
from config.logging_conf import GOES_LOGGER_NAME, get_logger_from_config_file
from plotters_lib.adec_plotter import plot_lst_heladas
from plotters_lib.fire_temp_rgb732 import graficar_fire_temp

logger = get_logger_from_config_file(GOES_LOGGER_NAME)


def reducir_calidad(ruta):
    """
    Reduce la calidad de las imágenes dentro de la carpeta ``ruta``
    """
    try:
        file_names = sorted(glob.glob(ruta + '/*/*.*'))

        pool = mp.Pool()
        for filename in file_names:
            pool.apply_async(os.system, ('pngquant --quality=65-80 ' + filename + ' --ext .png --force'))
        pool.close()
        pool.join()
    except:
        logger.exception("Error inesperado.")
        raise


def generar_imagenes(fecha_cod):
    """
    Genera las imágenes para la fecha dada
    """
    fecha = datetime.datetime.strptime(fecha_cod, '%Y%m%d')
    pool = mp.Pool()
    for _ in PRODUCTS:
        for canal in CHANNELS:
            if canal == "C03":
                logger.info("Generando imagenes para GeoColor")
                pool.apply_async(plot.geocolor, (fecha, canales, NOAA_GOES16_PATH, MAPAS))
                pool.apply_async(graficar_fire_temp, (fecha, canales, NOAA_GOES16_PATH, MAPAS))
            elif canal == "C01" or canal == "C02" or canal == "C07":
                continue
            else:
                logger.info(f"Generando imagenes para {canales[canal].nombre}")
                pool.apply_async(plot_channel, (canales[canal], fecha, NOAA_GOES16_PATH, MAPAS))
    logger.info(f"Generando imagenes para webmet de {canales['C13'].nombre}")
    pool.apply_async(canales['C13'].plot_webmet, (fecha, NOAA_GOES16_PATH, MAPAS))
    logger.info("Generando grafico LST y heladas")
    pool.apply_async(plot_lst_heladas, (fecha, NOAA_GOES16_PATH, MAPAS))
    logger.info("Cerrando pool.")
    pool.close()
    pool.join()
    start = time.time()
    # se reduce la calidad de los productor para la web
    reducir_calidad(MAPAS + '/' + fecha.strftime('%Y_%m/%d/'))
    logger.info(f"Tiempo de reducir calidad: {((time.time() - start) // 60)} minutos")


def main():
    parser = argparse.ArgumentParser(prog="Generar Imágenes")
    parser.add_argument("fecha", help="Día del que se descargan los DATOS (se descargan en horario UTC)")
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s v0.7.0')
    args = parser.parse_args()
    generar_imagenes(args.fecha)


if __name__ == "__main__":
    main()
