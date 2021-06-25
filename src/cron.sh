#!/bin/bash

if [ $# -lt 1 ]
then
    fecha=`date -u +"%Y%m%d"`
else
    fecha=$1
fi

# Crea las variables de entorno
source env.sh

fecha_folder=`date -d $fecha +%Y_%m/%d`

# Directorio de logs
mkdir -p $LOGS/$fecha_folder

echo "
        ################ " `date +"%Y-%m-%d %H:%M"` "################
        ################ Se descargan las imágenes de GOES ################
    "

time python descarga_amazon_boto.py $fecha > $LOGS/$fecha_folder/descarga.log 2>&1


echo "
        ################ " `date +"%Y-%m-%d %H:%M"` "################
        ################ Se crean las imágenes ################
    "

time python gen_img.py $fecha > $LOGS/$fecha_folder/imagenes.log 2>&1
