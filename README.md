# GOES-operativo

Scripts para descargar y generar imágenes del satélite GOES 16.

### Descripción

Para un día, descarga todos los archivos disponibles del producto/s y banda/s solicitada/s (en caso de estar descargados ya, no lo vuelve a descargar), y genera los gráficos para todos esos archivos (en caso de ya estar generados, no los vuelve a generar).


## Pre-requicitos

1. [Anaconda](https://www.anaconda.com/download/) o [Miniconda](https://conda.io/docs/user-guide/install/linux.html)
1. [aws-cli](https://github.com/aws/aws-cli)

## Instalación

1. Clonar el repositorio en la carpeta local:

    ```
    git clone https://gitlab.unc.edu.ar/wrf/goes-operativo.git
    ```

1. Crear un entorno de conda a partir del archivo *enviroment.yml* que tiene todas las librerías de python necesarias. Puede ser a través de un gestor conda como Anaconda o Miniconda.
    
    ```
    cd goes-operativo/ejecutables/
    conda env create -f environment.yml
    ```

## Ejemplo de uso

- Para operatividad, descarga y genera lo disponible del día:

    ```
    cd ejecutables/
    ./cron.sh 
    ```

- Para descargar y generar lo de un día en particular (desde las 00UTC hasta las 24UTC de ese día):

    ```
    cd ejecutables/
    ./cron.sh AAAAMMDD
    ```
