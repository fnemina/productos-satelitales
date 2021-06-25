import os
import logging
import logging.config as lc
import pathlib

import yaml

BASE_LOG_PATH = 'logs'
LOGGING_CONFIG_FILE = os.getenv('LOG_CFG', 'config/goes_logging_conf.yml')
GOES_LOGGER_NAME = 'goes'
INGESTOR_LOGGER_NAME = 'ingestor'


def get_logger_from_config_file(logger_name: str) -> logging.Logger:
    """
    Retorna logger a partir del archivo YAML

    Parameters
    ----------
    logger_name: str
        Nombre del logger

    Returns
    -------
    logging.Logger

    """
    logger_config_file = pathlib.Path(LOGGING_CONFIG_FILE)
    if not logger_config_file.exists():
        raise FileNotFoundError
    with open(logger_config_file, mode='r') as f:
        log_cfg = yaml.safe_load(f)
    lc.dictConfig(log_cfg)
    return logging.getLogger(name=logger_name)
