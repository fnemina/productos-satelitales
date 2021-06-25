import datetime
import logging
from json.decoder import JSONDecodeError

import requests
from requests.exceptions import RequestException

from config.constants import GOES_DEBUG
from config.logging_conf import INGESTOR_LOGGER_NAME
from config.wrf_api_constants import API_BASE_URL_DICT, API_ROOT, TIMESTAMP_API_FORMAT

logger = logging.getLogger(INGESTOR_LOGGER_NAME)


def get_wrf_api_object_id(api_base_url, nombre, valor, campo='search'):
    try:
        r = requests.get(f"{api_base_url}/{nombre}/?{campo}={valor}")
        r.raise_for_status()
        return r.json().get('results')[0].get('id')
    except IndexError:
        logger.exception(f"No se encontro el ID de {valor}")
        raise
    except RequestException:
        logger.exception("Error al obtener objeto: ")
        raise


def create_wrf_object(api_base_url, token, nombre, payload: dict):
    headers_base = {'Authorization': f'Token {token}'}
    try:
        r = requests.post(f"{api_base_url}/{nombre}/", json=payload, headers=headers_base)
        r.raise_for_status()
    except RequestException:
        logger.exception(f"Error al crear objeto en {api_base_url} - nombre={nombre},payload={payload}")
        raise
    return r.json().get('id') if r.ok else -1


def post_img_to_api(img_base_path: str, timestamp: datetime.datetime, producto: str, campo_prod: str = 'short_name',
                    region: str = 'CBA'):
    if GOES_DEBUG:
        logger.debug("Modo DEBUG activado, skiping API ingest...")
        return
    api_file_path = f"{API_ROOT}/{img_base_path}"
    for api_base_url, meta in API_BASE_URL_DICT.items():
        try:
            product_sat_id = get_wrf_api_object_id(api_base_url, nombre='producto-satelital', valor=producto,
                                                   campo=campo_prod)
            region_id = get_wrf_api_object_id(api_base_url, nombre='region', valor=region, campo='short_name')
        except Exception:
            logger.exception(f"Error al obtener datos de la API {api_base_url}")
            continue
        timestamp_str = datetime.datetime.strftime(timestamp, TIMESTAMP_API_FORMAT)
        api_payload = {
            "path": api_file_path,
            "timestamp": timestamp_str,
            "producto_satelital": product_sat_id,
            "region": region_id
        }
        try:
            create_wrf_object(api_base_url, meta['token'], 'goes-product', api_payload)
        except Exception:
            logger.exception(f"Error al crear el objeto en la API {api_base_url}")


def post_lst_geotif_to_api(geotif_base_path: str, timestamp: datetime.datetime):
    if GOES_DEBUG:
        logger.debug("Modo DEBUG activado, skiping API ingest...")
        return
    api_file_path = f"{API_ROOT}/{geotif_base_path}"
    productos_satelitales = ["heladas", "LST"]
    for api_base_url, meta in API_BASE_URL_DICT.items():
        productos_satelitales_ids = []
        for producto_satelital in productos_satelitales:
            try:
                ps_id = get_wrf_api_object_id(api_base_url, nombre='producto-satelital', valor=producto_satelital,
                                              campo='short_name')
                productos_satelitales_ids.append(ps_id)
            except Exception:
                logger.exception(f"Error al obtener datos de la API {api_base_url}")
                continue
        if not productos_satelitales_ids:
            logger.warning(f"No se obtuvieron los id de los productos satelitales para API {api_base_url}")
            continue
        timestamp_str = datetime.datetime.strftime(timestamp, TIMESTAMP_API_FORMAT)
        api_payload = {
            "path": api_file_path,
            "timestamp": timestamp_str,
            "sat_products": productos_satelitales_ids
        }
        try:
            create_wrf_object(api_base_url, meta['token'], 'geotif', api_payload)
        except Exception:
            logger.exception(f"Error al crear el objeto en la API {api_base_url}")
