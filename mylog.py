import logging
import os

def get_logger(f, n):

    logger = logging.getLogger(os.path.basename(f) + " - " +  n)
    logger.setLevel(logging.DEBUG)

    handler = logging.FileHandler(os.path.expanduser('~/logs/common.log'))
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)

    # add the handlers to the logger
    logger.addHandler(handler)

    return logger


