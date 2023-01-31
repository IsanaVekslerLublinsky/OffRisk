import logging
import subprocess
from datetime import timedelta
from os import path
from time import perf_counter
from bs4 import BeautifulSoup
import pandas as pd
import requests

import configuration_files.const as const
from configuration_files.const import FLASHFRY_TMP_LOCATION_PATH, FLASHFRY_DATABASE_BASE_PATH, \
    COMPLETE_GENOME_PATH, FLASHFRY_INPUT_PATH, FLASHFRY_OUTPUT_PATH, FLASHFRY_SCORE_OUTPUT_PATH, \
    CAS_OFFINDER_INPUT_FILE_PATH, CAS_OFFINDER_OUTPUT_PATH
from obj_def import OffTarget

log = logging.getLogger("Base_log")


def load_cas_offinder_off_target():
    """
    Take the off-target file from cas-offinder and convert it to be a bedtools file
    :return: pybedools object of the off-target
    """

    function_name = "convert_off_target_to_bed"
    log.debug("Entering {}".format(function_name))

    # todo: we will need to verify the off-target file format
    # Read the off-target file to a dataframe
    off_target_df = pd.read_csv(CAS_OFFINDER_OUTPUT_PATH, sep="\t")

    off_target_df = off_target_df[["Chromosome", "Position", "DNA", "Direction"]]
    off_target_df.loc[:, "end"] = off_target_df.apply(lambda row: (row.loc["Position"] + len(row.loc["DNA"])), axis=1)
    off_target_df.loc[:, "name"] = off_target_df.index
    off_target_df.loc[:, "score"] = 0
    off_target_df = off_target_df.rename(
        columns={"Chromosome": "chromosome", "Position": "start", "Direction": "strand"})
    off_target_df = off_target_df[["chromosome", "start", "end", "name", "score", "strand"]]
    return off_target_df


def load_flashfry_off_target():
    """
    Take the off-target file from flashfry and convert it to be a BedTools file

    each off-targets has are three fields, separated with underscores. The first is the actual sequence, the second
    in the number of times this off-target exists in the genome, and the third is the number of differences between
    this off-target and the target of interest. for positional information, the output will look as follows:
    individual off-target entries will have an associated list of positions within the genome, separated by the pipe
    "|" . The off-target position list will be separated from the sequences by the open and closed bracket characters
    ("<" and ">"). The position entry contains the chromosome, the start location (on the Watson strand) separated
    with a colon ":", and the orientation of the off-target separated by a carrot "^".

    Returns: T
    :return: DataFrame and pybedtools object of the off-target

    """

    function_name = "convert_off_target_flashfry_to_bed"
    log.debug("Entering {}".format(function_name))

    # todo: we will need to verify the off-target file format

    # Read the off-target file to a dataframe
    flashfry_output = pd.read_csv(FLASHFRY_OUTPUT_PATH, sep="\t")
    flashfry_score = pd.read_csv(FLASHFRY_SCORE_OUTPUT_PATH, sep="\t")

    # Keep only the row that start in the guide start position
    flashfry_output = flashfry_output[flashfry_output["start"] == 0]
    flashfry_score = flashfry_score[flashfry_score["start"] == 0]
    if (not flashfry_score.empty) | (not flashfry_output.empty):
        merge_flashfry = pd.merge(flashfry_output, flashfry_score, on=["contig", "start", "stop", "target", "context"])
        merge_flashfry.drop(["overflow_y", "orientation_y", "otCount_y"], axis=1, inplace=True)
        flashfry_score = merge_flashfry.copy()
        flashfry_score.drop(["offTargets"], axis=1, inplace=True)
    else:
        return None, None
    off_target_df_list = list()
    for i, row in merge_flashfry.iterrows():
        tmp_ot_df = pd.DataFrame(row.offTargets.split(","))
        tmp_ot_df["tmp"] = tmp_ot_df.apply(lambda s: s[0].split("<")[1].split("|"), axis=1)
        tmp_ot_df = tmp_ot_df.explode("tmp")
        tmp_ot_df["chromosome"] = tmp_ot_df["tmp"].apply(lambda s: s.split(":")[0])
        tmp_ot_df["start"] = tmp_ot_df["tmp"].apply(lambda s: s.split(":")[1].split("^")[0])
        tmp_ot_df["end"] = tmp_ot_df.apply(lambda s: int(s["start"]) + len(s[0].split("_")[0]), axis=1)
        tmp_ot_df["strand"] = tmp_ot_df["tmp"].apply(lambda s: s.split("^")[1].strip(">"))
        tmp_ot_df["strand"] = tmp_ot_df["strand"].apply(lambda s: "+" if s == "F" else "-")
        tmp_ot_df.loc[:, "name"] = ""
        tmp_ot_df.loc[:, "score"] = 0
        tmp_ot_df.loc[:, "attributes"] = "contig={};target={}".format(row["contig"], row["target"])
        # Order the DataFrame
        tmp_ot_df.drop(["tmp"], axis=1, inplace=True)
        tmp_ot_df = tmp_ot_df[["chromosome", "start", "end", "name", "score", "strand", "attributes"]]
        off_target_df_list.append(tmp_ot_df)

    off_target_to_bed = pd.concat(off_target_df_list)
    return off_target_to_bed, flashfry_score


def load_off_target_from_file(input_file):
    """
    load the off-target input file given
    Args:
        input_file: the path for the input file. need to be an absolute path

    Returns: off-target in dataframe and in off-target bed file

    """
    function_name = "load_off_target_from_file"
    log.debug("Entering {}".format(function_name))

    fields_name = OffTarget.get_fields_title()
    off_target_df = pd.read_csv(input_file, sep="\t", header=0, names=fields_name)
    off_target_df.loc[:, "name"] = off_target_df.index
    off_target_df.loc[:, "score"] = 0
    off_target_df.loc[:, "attributes"] = ""

    off_target_df = off_target_df[[OffTarget.get_field_title("chromosome"),
                                   OffTarget.get_field_title("start"), OffTarget.get_field_title("end"),
                                   "name", "score", OffTarget.get_field_title("strand"), "attributes",
                                   OffTarget.get_field_title("id"), OffTarget.get_field_title("sequence")]]
    return off_target_df


def run_external_proc(args):
    """
    Run an external program with Popen.
    Args:
        args: The argument to run, a list of string
    """

    function_name = "run_external_proc"
    log.debug("Entering {}".format(function_name))
    log.debug("{}: The following command will be run: {}".format(function_name, args))
    with subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as proc:
        while True:
            next_line = proc.stdout.readline().decode("utf-8").strip()
            if next_line == "" and proc.poll() is not None:
                break
            if next_line != "":
                log.info(next_line)

        err, output = proc.communicate()
        exit_code = proc.returncode

        if exit_code == 0:
            log.info(output.decode())
        else:
            raise Exception("Error in {}.\nCommand is : {}\nError is: {} {}".format(
                function_name, args, output.decode(), err.decode()))


def run_cas_offinder_api(server_address="http://apps-dev.crispr-il.local/off-tov/"):
    """
    Run cas offinder from a server with api.
    Important! the server must be off-tov server otherwise the request will not work

    server_address: need to be the full address with port
    """
    function_name = "run_cas_offinder_api"
    log.debug("Entering {}".format(function_name))

    with open(CAS_OFFINDER_INPUT_FILE_PATH, "r") as file:
        lines = file.readlines()
        genome_type = "human"
        pattern = lines[1].strip("\n")
        sequences = ",".join(line.strip("\n") for line in lines[2:])

        request_to_send = "{}/v1/cas-offinder-bulge/?genome={}&pattern={}&sequences={}".format(
            server_address, genome_type, pattern, sequences)
    log.info("{}: Sending the following request: {}\nThis might take a while".format(
        function_name, request_to_send))
    response = requests.get(request_to_send)
    if response.status_code == 200:
        log.info("Server returned code: {}".format(response.status_code))
        result_df = response.json()["cas-offinder-dataframe"]
        result_output = response.json()["message"]
        log.info(result_output)
        off_target_result = pd.read_json(result_df)
        off_target_result.to_csv(CAS_OFFINDER_OUTPUT_PATH, sep="\t")
    else:
        raise Exception("Server returned error code: {}: {}".format(
            response.status_code, BeautifulSoup(response.text, "html.parser").getText()))


def run_cas_offinder_locally(c_g_option):
    """
    Run Cas-OFFinder-bulge
    :param c_g_option: run on CPU or GPU
    :return:
    """
    time_start = perf_counter()
    # Verify files exist
    if not path.exists(const.CAS_OFFINDER_BULGE_PATH):
        raise Exception("No path {} exist. Please change Cas-OFFinder location".format(const.CAS_OFFINDER_BULGE_PATH))
    if not path.exists(CAS_OFFINDER_INPUT_FILE_PATH):
        raise Exception("No path {} exist. Please change Cas-OFFinder input file location.".
                        format(CAS_OFFINDER_INPUT_FILE_PATH))
    args = [const.CAS_OFFINDER_BULGE_PATH, CAS_OFFINDER_INPUT_FILE_PATH, c_g_option,
            CAS_OFFINDER_OUTPUT_PATH]

    log.info("Starting to run cas-offinder")
    run_external_proc(args)
    time_end = perf_counter()
    log.info("Finish running cas-offinder in {}".format(timedelta(seconds=(time_end - time_start))))


def run_flashfry(database_path, commands=None):
    """
    Run FlashFry with the commands received
    Args:
        database_path: the base path of the database
        commands: commands to run. can be build_index, discover or score
    """

    flashfry_header_path = "{}{}".format(database_path, FLASHFRY_DATABASE_BASE_PATH)
    log.info("FlashFry header location: {}".format(flashfry_header_path))
    if commands is None:
        commands = ["discover", "score"]
    time_start = perf_counter()
    # Check if index file already exist
    if path.exists(flashfry_header_path):
        log.info("Database already exist, it will not be build again.")
    else:
        if "build_index" in commands:
            build_index_args = ["java", "-Xmx4g", "-jar",
                                const.FLASHFRY_PATH, "index",
                                "--tmpLocation", FLASHFRY_TMP_LOCATION_PATH,
                                "--database", flashfry_header_path,
                                "--reference", COMPLETE_GENOME_PATH, "--enzyme", "spcas9ngg"]

            log.info("Starting to run FlashFry index building")
            run_external_proc(build_index_args)
            log.info("Finish running FlashFry index building")

    if "discover" in commands:
        discover_args = ["java", "-Xmx4g", "-jar", const.FLASHFRY_PATH, "discover",
                         "--database", flashfry_header_path,
                         "--fasta", FLASHFRY_INPUT_PATH,
                         "--positionOutput",
                         "--output", FLASHFRY_OUTPUT_PATH]

        log.info("Starting to run FlashFry discover")
        run_external_proc(discover_args)
        log.info("Finish running FlashFry discover")

    if "score" in commands:
        score_args = ["java", "-Xmx4g", "-jar", const.FLASHFRY_PATH, "score",
                      "--input", FLASHFRY_OUTPUT_PATH,
                      "--output", FLASHFRY_SCORE_OUTPUT_PATH,
                      "--scoringMetrics", "doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot",
                      "--database", flashfry_header_path]

        log.info("Starting to run FlashFry score")
        run_external_proc(score_args)
        log.info("Finish running FlashFry score")
    time_end = perf_counter()
    log.info("Finish running FlashFry in {}".format(timedelta(seconds=(time_end - time_start))))
