import os
import re

NOAA_GOES16 = 'noaa-goes16'

PICKLE_PATH_BASE = "pickles"
GOES_API_PICKLE_PATH = f"{PICKLE_PATH_BASE}/goes.pickle"
# env
GOES_DEBUG = bool(os.getenv('GOES_DEBUG', False))
MAPAS = os.getenv('MAPAS')
NOAA_GOES16_PATH = os.getenv('NOAA_GOES16_PATH')
GOES_STATICS = os.getenv('GOES_STATICS')
CHANNELS = os.getenv('CHANNELS').split(' ')
PRODUCTS = os.getenv('PRODUCTS').split(' ')
FINAL_PRODUCTS = os.getenv('FINAL_PRODUCTS').split(' ')
# EXTENT = (-75, -40, -58, -25)
EXTENT = (int(os.environ['LON_MIN']),
          int(os.environ['LAT_MIN']),
          int(os.environ['LON_MAX']),
          int(os.environ['LAT_MAX']))

EXTENT_ADEC = [-75, -53, -25, -42]
TEMP_MIN_HOUR_UTC = int(os.getenv('TEMP_MINIMAS_HOUR', 12))  # hora para grafico temperaturas minimas
# LST -> [-80, -52, -58, -20]
EXTENT_HELADAS = (int(os.environ['LON_MIN']),  # -75
                  int(os.environ['LON_MAX']),  # -58
                  int(os.environ['LAT_MAX']),  # -25
                  int(os.environ['LAT_MIN']))  # -40

EXTENT_WEBMET = (int(os.environ['LON_MIN_WEBMET']),
                 int(os.environ['LAT_MIN_WEBMET']),
                 int(os.environ['LON_MAX_WEBMET']),
                 int(os.environ['LAT_MAX_WEBMET']))

SHAPEFILES = f"{GOES_STATICS}/shapefiles"
SHAPEFILES_ADEC = f"{GOES_STATICS}/shapefiles_adec"
SHAPEFILES_ADEC_PROV = f"{SHAPEFILES_ADEC}/provincias"
SHAPEFILES_ADEC_DEP = f"{SHAPEFILES_ADEC}/departamentos"
SHAPEFILES_ADEC_DEP_POL = f"{SHAPEFILES_ADEC}/departamentos_poligonos"
SHAPEFILES_ADEC_URUGUAY = f"{SHAPEFILES_ADEC}/uruguay"
# /media/disk/goes/goes-operativo/datos/shapefiles_adec/uruguay

LOGO = f"{GOES_STATICS}/img/wug_txt2.png"
LOGO_ADEC = f"{GOES_STATICS}/img/logo_inta_gulich_ohmc_conae.png"
LOGO_FIRE = f"{GOES_STATICS}/img/logo_inta_gulich_ohmc_conae_fire_02.png"

DOWNLOAD_MAX_TRY = 2
DOWNLOAD_SLEEP_TRY = 30

CHANNEL_REGEX = r"(?P<sensor>[A-Z]{2,7})-(?P<prod>[A-Z0-9]{2,4}-[A-Z0-9]{2,8})/(?P<year>\d{4})/" \
                r"(?P<day_of_year>\d{1,3})/(?P<hour>\d{1,2})/\w{4,6}-\w{2,3}-\w{4,6}-[A-Z0-9]{2}" \
                r"(?P<channel>[A-Z][0-9]{2})_G16_s\d+_e\d+_c\d+.nc"

GENERAL_GOES_REGEX = r"(?P<sensor>[A-Z]{2,7})-(?P<prod>[A-Z0-9]{2,4}-[A-Z0-9]{2,8})/(?P<year>\d{4})/" \
                     r"(?P<day_of_year>\d{1,3})/(?P<hour>\d{1,2})/\w{4,6}-\w{2,3}-\w{4,6}-[A-Z0-9]{2}" \
                     r"(?P<channel>[A-Z][0-9]{2})?_G16_s\d+_e\d+_c\d+.nc"


def get_channel_regex(_ch):
    r_pre = "".join([
        r"[A-Z]{2,3}-(?P<prod>[A-Z0-9]{2,4}-[A-Z0-9]{2,8})/\d{4}/\d{1,3}/\d{1,2}/\w{4,6}-\w{2,3}-\w{4,6}-[A-Z0-9]{2}",
        _ch,
        r"_\w*?.nc"
    ])
    return re.compile(r_pre)


# dict_variable = {key:value for (key,value) in dictonary.items()}
CHANNEL_RE_DICT = {channel: get_channel_regex(channel) for channel in CHANNELS}

MONTH_DICT = {
    1: 'Enero',
    2: 'Febrero',
    3: 'Marzo',
    4: 'Abril',
    5: 'Mayo',
    6: 'Junio',
    7: 'Julio',
    8: 'Agosto',
    9: 'Septiembre',
    10: 'Octubre',
    11: 'Noviembre',
    12: 'Diciembre',
}

CMI = "CMI"
