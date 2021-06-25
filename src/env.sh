#!/bin/bash

#export GOES_DEBUG="True"
# WRF API
export API_BASE_URL_DICT='{"https://wrf.ohmc.com.ar/api": {"token": "53f990f391eda44c0d61ed3f013a8de0bf33e616"}, "https://wrf-beta.ohmc.com.ar/api": {"token": "53f990f391eda44c0d61ed3f013a8de0bf33e616"}}'
# DIRECTORIOS
export GOES_STATICS="/mnt/yaku/goes16-statics"
export NOAA_GOES16_PATH="/mnt/yaku/noaa-goes16"
export MAPAS="/mnt/yaku/img/GOES"
export LOGS="/mnt/yaku/goes16-logs"
# DATOS DE DESCARGA
export CHANNELS='C01 C02 C03 C07 C09 C10 C12 C13' # Ejemplo: CHANNELS='C01 C02 C03 C04 C05 C06 C07 C08 C09 C10 C11 C12 C13 C14 C15 C16'
export PRODUCTS='L2-CMIPF' # Ejemplo: PRODUCTS='L1b-RadC L1b-RadF L1b-RadM L2-CMIPC L2-CMIPF L2-CMIPM L2-MCMIPC L2-MCMIPF L2-MCMIPM'
export FINAL_PRODUCTS='L2-LSTF L2-ACMF'
# DATOS DE GRAFICO
export LON_MIN="-75"
export LAT_MIN="-40"
export LON_MAX="-58" #-50 #58
export LAT_MAX="-25"
# webmet
export LON_MIN_WEBMET="-77"
export LON_MAX_WEBMET="-51"
export LAT_MIN_WEBMET="-56"
export LAT_MAX_WEBMET="-20"
# conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate goes-ingestor-env
#conda activate goes-ingestor
