import gzip
import logging
import os
import re
import gffutils
import pandas as pd
from validators import url, ValidationFailure
from pydantic import BaseModel, validator
from typing import List

import configuration_files.const as const

log = logging.getLogger(__name__)

COL_ORDER = ["chrom", "start", "end", "name", "score", "strand", "frame", "attribute"]


class ConfigurationFile(BaseModel):
    run_cas_offinder: str = None
    cas_offinder_path: str = None
    run_cas_offinder_api: str = None
    run_flashfry: List[str] = None
    flashfry_path: str = None
    input_file_path: str = None
    db_list: List[str] = None
    debug: bool = None
    update_db: bool = None
    update_db_base_dir_path: str = None

    @validator("run_cas_offinder")
    def val_run_cas_offinder(cls, v):
        if v not in ["C", "G"]:
            raise ValueError("Processor type can only be C or G")
        return v

    @validator("run_cas_offinder_api")
    def val_run_cas_offinder_api(cls, v):
        test = url(v)
        if type(test) == ValidationFailure:
            raise ValueError("Not a valid address")
        return v

    @validator("run_flashfry")
    def val_run_flashfry(cls, v):
        valid_input = ["build_index", "discover", "score"]
        test = all(item in valid_input for item in v)
        if not test:
            raise ValueError("{} is not a valid input for FlashFry. Valid input are: {}".format(v, valid_input))
        return v

    @validator("db_list")
    def val_db_list(cls, v):
        valid_input = const.DB_NAME_LIST + ["all"]
        test = all(item in valid_input for item in v)
        if not test:
            raise ValueError("{} is not a valid input for FlashFry. Valid input are: {}".format(v, valid_input))
        return v

    @validator("input_file_path", "update_db_base_dir_path", "cas_offinder_path", "flashfry_path")
    def val_input_file_path(cls, v):
        if v == "":
            raise ValueError("Path is empty. Please specify a valid path")
        if not re.fullmatch(
                r"^(?!.*[\\\/]\s+)(?!(?:.*\s|.*\.|\W+)$)(?:[a-zA-Z]:)?(?:(?:[^<>:\"\|\?\*\n])+(?:\/\/|\/|\\\\|\\)?)+$",
                v):
            raise ValueError("Invalid path")
        return v


def update_cas_offinder_path(new_path):
    """
        Update the path for the Cas-Offinder executable when running locally
        Args:
            new_path: the new  path

    """
    if os.path.exists(new_path):
        log.info("Change Cas-Offinder path to: {}".format(new_path))
        const.CAS_OFFINDER_BULGE_PATH = new_path
    else:
        log.error("Given path: {} , is not valid".format(new_path))


def update_flashfry_path(new_path):
    """
        Update the path for the FlashFry executable when running locally
        Args:
            new_path: the new  path

    """
    if os.path.exists(new_path):
        log.info("Change FlashFry path to: {}".format(new_path))
        const.FLASHFRY_PATH = new_path
    else:
        log.error("Given path: {} , is not valid".format(new_path))


def init_logger(log_file=const.LOG_PATH, debug_level=logging.DEBUG, logger_name=None):
    """
    Init the logger for the program
    :param log_file: the location of the log file
    :param debug_level: the minimum log level
    :param logger_name: the name of the logger in the module
    """

    logger = logging.getLogger(logger_name)
    logger.setLevel(debug_level)

    # create file handler which logs even debug messages
    fh = logging.FileHandler(log_file)
    fh.setLevel(debug_level)
    # create formatter and add it to the handlers
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    fh.setFormatter(formatter)
    # add the handlers to logger
    logger.addHandler(fh)

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logger.addHandler(console)


def get_logger(logger_name, debug_level=logging.DEBUG):
    logging.basicConfig(
        level=debug_level,
        format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger(logger_name)

    return logger


def extract_gz_file(file_to_extract_path, file_extraction_path):
    """
    Extract gz file
    Args:
        file_to_extract_path: the path of the gz file
        file_extraction_path: the path to extract the gz file
    """
    function_name = "extract_gz_file"
    log.debug("Entering {}".format(function_name))

    extract_file = open(file_extraction_path, "wb")
    file_gz = gzip.open(file_to_extract_path, "rb")
    extract_file.writelines(file_gz)
    file_gz.close()
    extract_file.close()

    log.debug("{}: Finish extracting the file to {}".format(function_name, file_extraction_path))
