import yaml
from flask import jsonify
import sys
from datetime import timedelta
from time import perf_counter
import pandas as pd
import logging
from os import path
import warnings
from pybedtools import BedTool

from db_queries import BedQuery
from genome_db_factories import CorporateGenomeDatabaseFactory
from results_processors import DBDataFrameProcessedResult
from obj_def import OffTarget
from configuration_files.const import CAS_OFFINDER_OUTPUT_PATH, FLASHFRY_OUTPUT_PATH, YAML_CONFIG_FILE
from off_target import load_off_target_from_file, \
    load_cas_offinder_off_target, load_flashfry_off_target
from helper import get_logger
import db as db

CURRENT_DIR = path.dirname(path.abspath(__file__))
sys.path.append(path.dirname(CURRENT_DIR))
BASE_DIR = path.dirname(path.dirname(path.abspath(__file__)))

warnings.filterwarnings("ignore", category=RuntimeWarning)
pd.options.mode.chained_assignment = None

log = get_logger(logger_name=__name__, debug_level=logging.DEBUG)


class IDataExtractor(object):
    # todo: flash fry score should not be here, he should be before if no input_file
    def extract(self, db_name_list, input_file=None):
        pass


class CrisprCorporateExtractor(IDataExtractor):

    def __init__(self):
        # get from os
        self._db_factory = CorporateGenomeDatabaseFactory(YAML_CONFIG_FILE)
        pass

    def extract(self, db_name_list, genome, input_file=None):
        """
        """
        function_name = "extract_data"
        log.debug("Entering {}".format(function_name))

        off_target_df = None
        flashfry_socre = pd.DataFrame()
        time_start = perf_counter()

        if input_file:
            log.info("{}: Loading off-target from a file in: {}".format(function_name, input_file))
            if path.exists(input_file):
                off_target_df = load_off_target_from_file(input_file)
            else:
                raise Exception("There is no file in the path: {}. Please verify the file was created."
                                .format(CAS_OFFINDER_OUTPUT_PATH))
        else:
            log.info("{}: Loading local off-target".format(function_name))
            off_target_co_df, off_target_ff_df = None, None
            if path.exists(CAS_OFFINDER_OUTPUT_PATH):
                off_target_co_df = load_cas_offinder_off_target()
            if path.exists(FLASHFRY_OUTPUT_PATH):
                off_target_ff_df, flashfry_socre = load_flashfry_off_target()

            # If both
            # files was loaded, merge them
            if (off_target_co_df is not None) & (off_target_ff_df is not None):
                off_target_df = pd.concat([off_target_co_df, off_target_ff_df])
                off_target_df = off_target_df.astype({OffTarget.get_field_title("chromosome"): "category",
                                                      OffTarget.get_field_title("start"): int,
                                                      OffTarget.get_field_title("end"): int,
                                                      OffTarget.get_field_title("strand"): "category"})
                off_target_df.drop_duplicates(subset=[OffTarget.get_field_title("chromosome"),
                                                      OffTarget.get_field_title("start"),
                                                      OffTarget.get_field_title("end"),
                                                      OffTarget.get_field_title("strand")], inplace=True)

            elif off_target_co_df is not None:
                off_target_df = off_target_co_df
            elif off_target_ff_df is not None:
                off_target_df = off_target_ff_df

        # Create the bed file
        if off_target_df is not None:
            off_target_df.reset_index(drop=True, inplace=True)
            off_target_df["name"] = off_target_df.index
            off_target_bed = BedTool.from_dataframe(off_target_df).sort()
        else:
            raise Exception("Off target dataframe is empty. Please verify there is files in the output folder")

        # Initialize the result
        db_list = []
        db_results = dict()
        off_target_df = db.initialize_off_target_df(off_target_df)
        time_end = perf_counter()
        log.info("Total run for off-target initialization: {}".format(timedelta(seconds=(time_end - time_start))))

        # Start extracting information from the databases
        log.info("{}: Begin to run intersection between off-target and data".format(function_name))

        for db_name in db_name_list:
            time_start = perf_counter()
            current_db = self._db_factory.create_database("crispr", genome, db_name)
            if current_db:
                query = BedQuery(off_target_bed)
                results = current_db.select(query)
                processed_results = None
                if results is not None:
                    selected_attributes = []
                    with open(YAML_CONFIG_FILE) as f:
                        conf_yaml = yaml.load(f, Loader=yaml.FullLoader)
                        selected_attributes = conf_yaml["databases"][db_name][genome]["select_attributes"]

                    processed_results = DBDataFrameProcessedResult().process(results, selected_attributes)
                    off_target_df = off_target_df.merge(processed_results, how="left", on="off_target_id")
                    db_results["{}_result_list".format("S3")] = [{"name": db_name, "description": "Complete result",
                                                              "data": results.get_results().to_json(orient="records")}]
                # off_target_df = processed_results
                db_list.append(current_db)
                time_end = perf_counter()
                log.info("Total run for {} analyze: {}".format(current_db.get_db_name(),
                                                               timedelta(seconds=(time_end - time_start))))
            else:
                log.info("No DB was created. current DB name: {}".format("crispr"))

        log.info("Saving the results")
        time_start = perf_counter()
        ot_results = None
        with open(YAML_CONFIG_FILE) as f:
            conf_yaml = yaml.load(f, Loader=yaml.FullLoader)
            ot_results = db.save_global_off_target_results(off_target_df, flashfry_socre)
        time_end = perf_counter()
        log.info("Total run for saving: {}".format(timedelta(seconds=(time_end - time_start))))
        log.info("Clearing the result")

        return ot_results, jsonify(db_results)
