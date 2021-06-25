import datetime
import logging
import os
import time
import plotters_lib.plot as plot
import errno

from config.constants import DOWNLOAD_MAX_TRY, DOWNLOAD_SLEEP_TRY, GOES_STATICS
from config.logging_conf import GOES_LOGGER_NAME

logger = logging.getLogger(GOES_LOGGER_NAME)


class Canal:
    """docstring for Canal"""

    def __init__(self, codigo, nombre, colormap, unidad, variable='CMI'):
        super(Canal, self).__init__()
        self.codigo = codigo
        self.numero = int(codigo[-2:])
        self.nombre = nombre
        self.colormap = colormap
        self.cptfile = f"{GOES_STATICS}/ctp/{colormap}"
        self.cptfile2 = f"{GOES_STATICS}/ctp/IR4AVHRR6-webmet.cpt"
        self.unidad = unidad
        self.variable = variable
        if self.numero >= 7:
            self.visible = False
        else:
            self.visible = True

    def descargar(self, lista_disponible, fecha, _datos, i=0):
        lista_interes = [x for x in lista_disponible if self.codigo in x[2]]
        print('%s links para el canal %s' % (len(lista_interes), self.codigo))
        path_fecha = fecha.strftime('%Y_%m/%d/')
        path_descarga = _datos + '/' + path_fecha + self.codigo
        os.system('mkdir -p ' + path_descarga)
        descargo = False
        for (size, folder, file) in lista_interes:
            if os.path.isfile(path_descarga + '/' + file):
                local_size = os.path.getsize(path_descarga + '/' + file)
                if local_size == size:
                    continue
            print('Descargando %s' % file)
            os.system('aws s3 --no-sign-request cp s3://noaa-goes16/%s/%s %s/' % (folder, file, path_descarga))
            descargo = True
        if (not descargo) and i < DOWNLOAD_MAX_TRY:
            i = i + 1
            print('No se descargó ningún archivo, se reintentará en 1 minuto\
                  (Intento %i para el canal %s)' % (i, self.codigo))
            time.sleep(DOWNLOAD_SLEEP_TRY)
            self.descargar(lista_disponible, fecha, _datos, i)

    def plot_webmet(self, fecha: datetime.datetime, _datos, _salidas: str):
        path_fecha = fecha.strftime('%Y_%m/%d/')
        path_fecha_salidas = fecha.strftime('%Y/%m/%d')
        path_imagenes = _salidas + '/webmet/' + path_fecha_salidas  # + self.codigo
        dir_datos = _datos + '/' + path_fecha + self.codigo
        files = os.listdir(dir_datos)
        # print(files)
        if not os.path.exists(path_imagenes):
            try:
                os.makedirs(path_imagenes, 0o755)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
        #        os.system('mkdir -p ' + path_imagenes)
        for file in files:
            if not file.endswith(".nc"):
                continue
            metadatos = plot.get_metadatos(_datos + '/' + path_fecha + self.codigo + '/' + file)
            if os.path.isfile(path_imagenes + '/GOES16' + '_' + metadatos['fecha_img'].strftime('%Y%m%dT%H%M%S') +
                              'Z_C' + str(metadatos['icanal']) + '.png'):
                continue
            logger.info(f'Generando imagen de webmet para {file}')
            plot.generar_imagen_webmet(
                _datos + '/' + path_fecha + self.codigo + '/' + file,
                metadatos,
                path_imagenes,
                self
            )


canales = {
    'C01': Canal('C01', 'Azul Visible (0.47 µm)', 'SVGAIR2_TEMP.cpt', 'Reflectancia'),
    'C02': Canal('C02', 'Rojo Visible (0.64 µm)', 'SVGAIR2_TEMP.cpt', 'Reflectancia'),
    'C03': Canal('C03', 'Vegetación (0.86 µm)', 'SVGAIR2_TEMP.cpt', 'Reflectancia'),
    'C05': Canal('C05', 'Snow/Ice (1.6 µm)', 'SVGAIR2_TEMP.cpt', 'Reflectancia'),
    'C06': Canal('C06', 'Cloud Particle Size (2.2 µm)', 'SVGAIR2_TEMP.cpt', 'Reflectancia'),
    'C07': Canal('C07', 'Shortwave Window (3.9 µm)', 'SVGAIR2_TEMP.cpt', 'Reflectancia'),
    'C08': Canal('C08', 'Vapor de Agua Nivel Alto (6.2 µm)', 'SVGAWVX_TEMP.cpt', 'Temperatura de Brillo (°C)'),
    'C09': Canal('C09', 'Vapor de Agua Nivel Medio (6.9 µm)', 'SVGAWVX_TEMP.cpt', 'Temperatura de Brillo (°C)'),
    'C10': Canal('C10', 'Vapor de Agua Nivel Bajo (7.3 µm)', 'SVGAWVX_TEMP.cpt', 'Temperatura de Brillo (°C)'),
    'C12': Canal('C12', 'Ozone Band (9.6 µm)', 'SVGAWVX_TEMP.cpt', 'Temperatura de Brillo (°C)'),
    'C13': Canal('C13', 'IR - Ventana de Onda Larga Limpia (10.3 µm)', 'IR4AVHRR6.cpt', 'Temperatura de Brillo (°C)')
}
