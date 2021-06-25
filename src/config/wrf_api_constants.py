import os
import json

API_BASE_URL_DICT: dict = json.loads(os.getenv('API_BASE_URL_DICT'))

API_ROOT = "/img"

API_RESPONSES = {
    'post_ok': 201,
    'get_ok': 200
}

TIMESTAMP_API_FORMAT = '%Y-%m-%d %H:%M'
