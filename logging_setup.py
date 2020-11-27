import logging
import tqdm
import logging.config


class Log():
    def __init__(self):
        configFile

def setup_logging(log_prefix, debug_level=logging.INFO):

    root = logging.getLogger()
    root.setLevel(debug_level)

    tqdm_handler = TqdmLoggingHandler()
    tqdm_handler.setLevel(logging.DEBUG)
    root.addHandler(tqdm_handler)


class TqdmLoggingHandler(logging.StreamHandler):
    def __init__(self):
        logging.StreamHandler.__init__(self)

    def emit(self, record):
        msg = self.format(record)
        tqdm.tqdm.write(msg)
        self.flush()

