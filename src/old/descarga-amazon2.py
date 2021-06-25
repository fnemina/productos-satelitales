#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Descarga de datos de Goes desde AWS.

Paquete necesario:
    - awscli (``sudo apt  install awscli``)

Tener cargado en el entorno `CHANNELS` y `PRODUCTS`.

"""
from optparse import OptionParser
import os
import subprocess
import datetime
import time


def checkAWS(fecha_cod):
    """
    Selecciona los archivos disponibles en AWS
    """
    fecha = datetime.datetime.strptime(fecha_cod, '%Y%m%d')
    year = fecha.strftime('%Y')
    day_jul = fecha.strftime('%j')
    DATOS = os.environ['GOES_DATA']
    PRODUCTS = os.environ['PRODUCTS']
    products = PRODUCTS.split(' ')

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
        descargar(lista_disponible, fecha, DATOS)


def descargar(lista_disponible, fecha, DATOS, i=0):
    """
    Consulta si el archivo ya se ha descargado, si posee el tamaño
    correspondiente y lo descarga
    """
    print('%s links para el descargar' % (len(lista_disponible)))
    path_fecha = fecha.strftime('%Y_%m/%d/')
    path_descarga = DATOS + '/GOES/' + path_fecha
    os.system('mkdir -p ' + path_descarga)
    descargo = False
    for (size, folder, file) in lista_disponible:
        if os.path.isfile(path_descarga + file):
            local_size = os.path.getsize(path_descarga + file)
            if local_size == size:
                continue
        print('Descargando %s' % file)
        os.system('aws s3 --no-sign-request cp s3://noaa-goes16/%s/%s \
                        %s/' % (folder, file, path_descarga))
        descargo = True
    if (not descargo) and i < 5:
        i = i + 1
        print('No se descargó ningún archivo, se reintentará en 1 minuto\
                (Intento %i)' % (i,))
        time.sleep(60)
        descargar(lista_disponible, fecha, DATOS, i)


def main():
    usage = """python descarga-amazon.py [--fecha=AAAAMMDD]"""
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
        checkAWS(fecha)


if __name__ == "__main__":
    main()
