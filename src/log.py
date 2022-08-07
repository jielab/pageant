import logging
from time import strftime


logger = logging.getLogger()
logger.setLevel(logging.INFO)


def stream_log_start():
    ch = logging.StreamHandler()
    logger.addHandler(ch)


def log_start(directory: str) -> str:
    log_name = f'{directory}/log_files/{strftime("%Y%m%d%H%M%S")}.log'
    fh = logging.FileHandler(log_name)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logging.info('Logging start.')
    return log_name


def gui_log_setup(hander: logging.Handler) -> None:
    logger.addHandler(hander)


def gui_log_remove(hander: logging.Handler) -> None:
    logger.removeHandler(hander)
