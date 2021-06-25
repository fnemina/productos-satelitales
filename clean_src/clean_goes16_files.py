import datetime
from pathlib import Path

CHANNELS_TO_REMOVE = ['C01', 'C02', 'C03', 'C07', 'C09', 'C10', 'C12']
PRODUCTS_TO_REMOVE = ['L2-ACMF', 'L2-LSTF']

NOAA_GOES16_PATH = '/mnt/yaku/noaa-goes16'

ts_delta = datetime.timedelta(days=-31)
timestamp = datetime.datetime.today() + ts_delta
YEAR_MONTH = f"{timestamp.strftime('%Y_%m')}"
BASE_PATH = Path(f"{NOAA_GOES16_PATH}/{YEAR_MONTH}/")


if __name__ == '__main__':

    for channel in CHANNELS_TO_REMOVE:
        print("###################################")
        print(f"Channel {channel}:")
        channel_glob = sorted(BASE_PATH.glob(f"*/{channel}/*.nc*"))
        print(f"Se eliminaran {len(channel_glob)} del producto {channel}")
        for channel_file in channel_glob:
            channel_file.unlink()
        print(f"Fin eliminación de archivos para {channel}\n")

    for product in PRODUCTS_TO_REMOVE:
        print("###################################")
        print(f"Producto {product}")
        product_glob = sorted(BASE_PATH.glob(f"*/{product}/*/*.nc*"))
        print(f"Se eliminaran {len(product_glob)} del producto {product}")
        for product_file in product_glob:
            product_file.unlink()
        print(f"Fin eliminación de archivos para {product}\n")
