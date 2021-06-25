#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Descarga de datos de Goes desde AWS.

Paquete necesario:
    - awscli (``sudo apt  install awscli``)

Tener cargado en el entorno `CHANNELS` y `PRODUCTS`.

"""
from optparse import OptionParser
from models.canales import canales
import os
import subprocess
import datetime
import multiprocessing as mp


def descargar(fecha_cod):
    """
    Descarga los archivos de GOES que faltan para la fecha dada
    """
    fecha = datetime.datetime.strptime(fecha_cod, '%Y%m%d')
    year = fecha.strftime('%Y')
    day_jul = fecha.strftime('%j')
    DATOS = os.environ['DATOS']
    CHANNELS = os.environ['CHANNELS']
    channels = CHANNELS.split(' ')
    PRODUCTS = os.environ['PRODUCTS']
    products = PRODUCTS.split(' ')
    pool = mp.Pool()
    for producto in products:
        status, output = subprocess.getstatusoutput("""aws s3 --no-sign-request \
                             ls --recursive noaa-goes16/ABI-%s/%s/%s/ | awk \
                             '{print $3";"$4}'""" % (producto, year, day_jul))
        lista = output.split('\n')
        lista = [tuple(x.split(';')) for x in lista]
        lista_disponible = [(int(x), y.rsplit('/', 1)[0],
                            y.rsplit('/', 1)[1]) for (x, y) in lista]
        print('%s links para el producto %s' % (len(lista_disponible),
                                                producto))
        for canal in channels:
            print("Descargando %s" % canales[canal].nombre)
            pool.apply_async(canales[canal].descargar_canales,
                             (lista_disponible,
                              fecha,
                              DATOS))
    pool.close()
    pool.join()


def main():
    usage = """descarga-amazon.py [--fecha=AAAAMMDD]"""
    parser = OptionParser(usage)
    parser.add_option(
                        "--fecha",
                        dest="fecha",
                        help="Día del que se descargan los datos \
                                (se descargan en horario UTC)"
                    )
    (opts, args) = parser.parse_args()
    if not opts.fecha:
        print("Faltan parametros!")
        print(usage)
    else:
        fecha = opts.fecha
        descargar(fecha)


if __name__ == "__main__":
    main()
