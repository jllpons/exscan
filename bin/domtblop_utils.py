"""
Utilities for the domtblop.py scripts.
"""

import logging
import os
import sys

from typing import (
        TextIO,
        )


logger = logging.getLogger(__name__)


def setup_logger(level: str) -> logging.Logger:
    """
    Set up a custom logger.

    Args:
        level: The logging level.

    Returns:
        The logger.
    """

    # Example: `[script.py::foo] 2021-01-01 00:00:00 - INFO - Hello, world!`
    fmt = "[%(filename)s::%(funcName)s] %(asctime)s - %(levelname)s - %(message)s"
    level = level.upper()

    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(level)
    stderr_handler.setFormatter(logging.Formatter(fmt))

    logging.basicConfig(
        level=level,
        format=fmt,
        handlers=[stderr_handler]
        )

    logger = logging.getLogger(__name__)

    # If a library is too verbose, we can use this example:
    #logging.getLogger("requests").propagate = False

    return logger


def read_input(file_path: str) -> TextIO:
    """
    Opens the input file or stdin if the file_path is "-".

    Parameters:
        file_path (str): The path to the input file.

    Returns:
        file_handle: A file handle to the opened file or stdin.

    Raises:
        FileNotFoundError: If the input file is not found.
    """

    if file_path == "-":
        return sys.stdin
    elif file_path.startswith("/dev/fd/"):
        fd = int(os.path.basename(file_path))
        return os.fdopen(fd, "r")
    else:
        if not os.path.isfile(file_path):
            logger.error(f"File not found: {file_path}")
            raise FileNotFoundError

        return open(file_path, "r")
