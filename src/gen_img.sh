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

#time python descarga-amazon.py --fecha $fecha > $LOGS/$fecha_folder/descarga.log


echo "
        ################ " `date +"%Y-%m-%d %H:%M"` "################
        ################ Se crean las imágenes ################
    "

time python gen_img.py $fecha > $LOGS/$fecha_folder/imagenes.log

echo "
        ################ " `date +"%Y-%m-%d %H:%M"` "################
        ################ Se copian los datos al webserver ################
    "
#rsync -r -P -e "ssh -i $HOME/.ssh/goes_eureka" $SALIDAS/mapas/  wrf@200.16.30.250:/home/wrf/datos-webwrf/img/GOES/
#rsync -r -P -e "ssh -i $HOME/.ssh/goes_eureka" $SALIDAS/webmet/  wrf@200.16.30.250:/home/wrf/datos-webwrf/img/GOES/webmet/
#rsync -r -P -e "ssh -i $HOME/.ssh/goes_eureka" $DATOS/GOES/  wrf@200.16.30.250:/home/datos/goes_16/GOES/

