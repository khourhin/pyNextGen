import logging
import os

from colorlog import ColoredFormatter


def get_logger(f, n):

    logger = logging.getLogger(os.path.basename(f) + " - " + n)
    logger.setLevel(logging.DEBUG)

    handler = logging.FileHandler(os.path.expanduser('~/logs/common.log'))
    formatter = ColoredFormatter(
        "%(log_color)s%(asctime)s - %(name)s - %(levelname)s - %(message)s'",
        datefmt=None,
        reset=True,
        log_colors={
            'DEBUG': 'cyan',
            'INFO': 'green',
            'WARNING': 'yellow',
            'ERROR': 'red',
            'CRITICAL': 'red,bg_white',
        },
        secondary_log_colors={},
        style='%'
    )

    handler.setFormatter(formatter)

    # add the handlers to the logger
    logger.addHandler(handler)

    return logger


if __name__ == '__main__':
    logger = get_logger(__file__, __name__)
    logger.debug("This is debug")
    logger.info("This is info")
    logger.warning("This is warning")
    logger.error("This is error")
    logger.critical("This is critical")
