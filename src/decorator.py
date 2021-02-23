from time import time
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
        from functools import wraps

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


def progress_value(process_value: int):
    def decorator(func):
        from functools import wraps

        @wraps(func)
        def wrapper(*args, **kwargs):
            fun_res = func(*args, **kwargs)
            global progress_bar
            progress_bar += process_value
            logging.info(f"Progress of the analysis: {progress_bar}%")
            return fun_res

        return wrapper

    return decorator


def restart_progress_bar():
    global progress_bar
    progress_bar = 0
