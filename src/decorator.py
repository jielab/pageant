from time import time
from chardet.universaldetector import UniversalDetector
from detect_delimiter import detect
from functools import wraps
import logging


progress_bar = 0


# Decorator
def use_time(process_name=None):
    """
    Record the running time of the function
    :param process_name: the name of this process
    :return: A informational log which records the time
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            start = time()
            logging.info(f'{process_name if process_name else func.__name__.replace("_", " ").capitalize()}'
                         f' start')
            fun_res = func(*args, **kwargs)
            logging.info(f'{process_name if process_name else func.__name__.replace("_", " ").capitalize()}'
                         f' used time: {time() - start:.2f}s')
            return fun_res

        return wrapper

    return decorator


def progress_value(process_value: float):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            fun_res = func(*args, **kwargs)
            global progress_bar
            progress_bar += process_value
            logging.info(f"Progress of the analysis: {progress_bar:.1f}%")
            return fun_res

        return wrapper

    return decorator


def restart_progress_bar():
    global progress_bar
    progress_bar = 0


def auto_encoding(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            fun_res = func(*args, **kwargs)
        except UnicodeDecodeError as e:
            detector = UniversalDetector()
            detector.feed(e.object)
            detector.close()
            fun_res = func(*args, **kwargs, encoding=detector.result['encoding'])
        except Exception as e:
            raise e
        return fun_res

    return wrapper


def auto_sep(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            fun_res = func(*args, **kwargs)
        except AssertionError as e:
            fun_res = func(*args, **kwargs, sep=detect(e.args[0]))
        except Exception as e:
            raise e
        return fun_res

    return wrapper
