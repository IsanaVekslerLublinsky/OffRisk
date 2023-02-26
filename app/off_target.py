import logging
import subprocess
from datetime import timedelta
from time import perf_counter
from bs4 import BeautifulSoup
import pandas as pd
import requests
import tempfile
import os
import configuration_files.const as const
from configuration_files.const import FLASHFRY_TMP_LOCATION_PATH, FLASHFRY_DATABASE_BASE_PATH, \
    COMPLETE_GENOME_PATH, FLASHFRY_INPUT_PATH, FLASHFRY_OUTPUT_PATH, FLASHFRY_SCORE_OUTPUT_PATH, \
    CAS_OFFINDER_INPUT_FILE_PATH, CAS_OFFINDER_OUTPUT_PATH, BASE_DIR
from obj_def import OffTarget

log = logging.getLogger("Base_log")


def load_off_target_from_cas_offinder_and_flashfry(flashfry_output = None, flashfry_score = None,
                                                   flashfry_output_file = None, flashfry_score_output_file = None,
                                                   cas_offinder_output = None, cas_offinder_output_file = None):

    log.info("Loading local off-target")
    off_target_df = None


    off_target_co_df = load_cas_offinder_off_target(cas_offinder_output, cas_offinder_output_file)

    off_target_ff_df, flashfry_score = load_flashfry_off_target(flashfry_output, flashfry_score, flashfry_output_file, flashfry_score_output_file)

    # If both files was loaded, merge them
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
    if off_target_df is not None:
        off_target_df[OffTarget.get_field_title("chromosome")] = \
            off_target_df[OffTarget.get_field_title("chromosome")].apply(lambda s: s.strip("chr"))
        if OffTarget.get_field_title("id") not in off_target_df.columns:
            off_target_df[OffTarget.get_field_title("id")] = "."
        if OffTarget.get_field_title("sequence") not in off_target_df.columns:
            off_target_df[OffTarget.get_field_title("sequence")] = "."

    return off_target_df

def load_cas_offinder_off_target(cas_offinder_output = None, cas_offinder_output_file = None):
    """
    Take the off-target file from cas-offinder and convert it to be a bedtools file
    :return: pybedools object of the off-target
    """

    function_name = "load_cas_offinder_off_target"
    log.debug("Entering {}".format(function_name))

    log.debug(cas_offinder_output)
    if cas_offinder_output is None:
        if cas_offinder_output_file is None:
            cas_offinder_output = pd.DataFrame()
        else:
            cas_offinder_output = pd.read_csv(cas_offinder_output_file, sep="\t")

    if cas_offinder_output.empty:
        return None

    # todo: we will need to verify the off-target file format
    # Read the off-target file to a dataframe
    off_target_df = cas_offinder_output

    off_target_df = off_target_df[["Chromosome", "Position", "DNA", "Direction"]]
    off_target_df.loc[:, "end"] = off_target_df.apply(lambda row: (row.loc["Position"] + len(row.loc["DNA"])), axis=1)
    off_target_df.loc[:, "name"] = off_target_df.index
    off_target_df.loc[:, "score"] = 0
    off_target_df = off_target_df.rename(
        columns={"Chromosome": "chromosome", "Position": "start", "Direction": "strand"})
    off_target_df = off_target_df[["chromosome", "start", "end", "name", "score", "strand"]]
    return off_target_df


def load_flashfry_off_target(flashfry_output = None, flashfry_score = None, flashfry_output_file = None, flashfry_score_output_file = None):
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

    function_name = "load_flashfry_off_target"
    log.debug("Entering {}".format(function_name))

    # todo: we will need to verify the off-target file format

    # Read the off-target file to a dataframe
    if flashfry_output is None:
        if flashfry_output_file is None:
            flashfry_output = pd.DataFrame()
        else:
            flashfry_output = pd.read_csv(flashfry_output_file, sep="\t")

    if flashfry_score is None:
        if flashfry_score_output_file is None:
            flashfry_score = pd.DataFrame()
        else:
            flashfry_score = pd.read_csv(flashfry_score_output_file, sep="\t")

    if flashfry_score.empty and flashfry_output.empty:
        return None, None

    # Keep only the row that start in the guide start position
    if not flashfry_output.empty:
        flashfry_output = flashfry_output[flashfry_output["start"] == 0]

    if not flashfry_score.empty:
        flashfry_score = flashfry_score[flashfry_score["start"] == 0]

    merge_flashfry = pd.merge(flashfry_output, flashfry_score, on=["contig", "start", "stop", "target", "context"])
    merge_flashfry.drop(["overflow_y", "orientation_y", "otCount_y"], axis=1, inplace=True)
    flashfry_score = merge_flashfry.copy()
    flashfry_score.drop(["offTargets"], axis=1, inplace=True)

    off_target_df_list = list()
    for i, row in merge_flashfry.iterrows():
        tmp_ot_df = pd.DataFrame(row.offTargets.split(","))
        tmp_ot_df["tmp"] = tmp_ot_df.apply(lambda s: s[0].split("<")[1].split("|"), axis=1)
        tmp_ot_df = tmp_ot_df.explode("tmp")
        tmp_ot_df["chromosome"] = tmp_ot_df["tmp"].apply(lambda s: s.split(":")[0])
        tmp_ot_df["start"] = tmp_ot_df["tmp"].apply(lambda s: s.split(":")[1].split("^")[0])
        tmp_ot_df["end"] = tmp_ot_df.apply(lambda s: int(s["start"]) + len(s[0].split("_")[0]), axis=1)
        tmp_ot_df["strand"] = tmp_ot_df["tmp"].apply(lambda s: s.split("^")[1].strip(">"))
        tmp_ot_df["strand"] = tmp_ot_df["strand"].apply(lambda s: "+" if s == "F" else "se-")
        tmp_ot_df["mismatch"] = tmp_ot_df.apply(lambda s: "mismatch={}".format(s[0].split("_")[2].split("<")[0]), axis=1)
        tmp_ot_df["occurrence"] = tmp_ot_df.apply(lambda s: "occurrence={}".format(s[0].split("_")[1]), axis=1)
        tmp_ot_df.loc[:, "name"] = ""
        tmp_ot_df.loc[:, "score"] = 0
        tmp_ot_df.loc[:, "attributes"] = "contig={};target={}".format(row["contig"], row["target"])
        tmp_ot_df.loc[:, "attributes"] = tmp_ot_df.apply(
            lambda s: ";".join([str(s["attributes"]), str(s["mismatch"]), str(s["occurrence"])]), axis=1)
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


def run_cas_offinder_api(cas_offinder_input_file, server_address=""):
    """
    Run cas offinder from a server with api.
    Important! the server must be off-risk server otherwise the request will not work

    server_address: need to be the full address with port
    """
    function_name = "run_cas_offinder_api"
    log.debug("Entering {}".format(function_name))

    with open(cas_offinder_input_file, "r") as file:
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
    else:
        raise Exception("Server returned error code: {}: {}".format(
            response.status_code, BeautifulSoup(response.text, "html.parser").getText()))

    return off_target_result

def run_cas_offinder_locally(c_g_option, pattern, seqs, docker_path_to_genome):
    """
    Run Cas-OFFinder-bulge
    :param c_g_option: run on CPU or GPU
    :return:
    """

    time_start = perf_counter()
    cas_offinder_output = None

    # Verify files exist
    if not os.path.exists(const.CAS_OFFINDER_BULGE_PATH):
        raise Exception("No path {} exist. Please change Cas-OFFinder location".format(const.CAS_OFFINDER_BULGE_PATH))

    try:
        fd_in, path_in = write_cas_offinder_input(pattern, seqs, docker_path_to_genome)
        temp_cas_offinder_output_file = tempfile.NamedTemporaryFile(mode='w', delete=False,
                                                     dir='{}/app/configuration_files/'.format(BASE_DIR))

        log.debug("cas_offinder_input_file_path: {}".format(path_in))

        log.debug("cas_offinder_output_file_path: {}".format(temp_cas_offinder_output_file.name))
        temp_cas_offinder_output_file.close()

        args = [const.CAS_OFFINDER_BULGE_PATH, path_in, c_g_option,
                temp_cas_offinder_output_file.name]

        log.info("Starting to run cas-offinder")
        run_external_proc(args)

        cas_offinder_output = pd.read_csv(temp_cas_offinder_output_file.name, sep="\t")
        time_end = perf_counter()
        log.info("Finish running cas-offinder in {}".format(timedelta(seconds=(time_end - time_start))))


    finally:
        os.remove(temp_cas_offinder_output_file.name)

    return cas_offinder_output

def write_cas_offinder_input(pattern, seqs, docker_path_to_genome, input_file=None):
    """
    Create input file for cas-offinder
    Args:
        pattern: desired pattern including PAM site and optional DNA or RNA bulge sizes, separated by spaces.
        seqs: query sequences and maximum mismatch numbers, separated by spaces, each sequence is separate with ","
        genome_type: which genome type to run. default is human GRCH38
        input_file: path to input file. if None then create the path here.

    Returns: io file and path
    """

    fd_in = None
    if input_file is None:
        fd_in, input_file = tempfile.mkstemp()
    with open(input_file, "w") as tmp_in:
        tmp_in.write("{}\n".format(docker_path_to_genome))
        tmp_in.write("{}".format(pattern))

        for line in seqs.split(","):
            tmp_in.write("\n")
            tmp_in.write(line)
    return fd_in, input_file

def run_flashfry(database_path, command=None, **kwargs):
    """
    Run FlashFry with the commands received
    Args:
        database_path: the base path of the database
        commands: commands to run. can be build_index, discover or score
    """

    if command is None:
        raise Exception('A Command should be one of the following: [discover, score, index]')

    flashfry_header_path = "{}{}".format(database_path, FLASHFRY_DATABASE_BASE_PATH)

    if os.path.exists(flashfry_header_path):
        log.info("Database already exist, it will not be build again.")
    else:
        _run_flashfry_index(flashfry_header_path)

    time_start = perf_counter()

    if command == 'index':
        _run_flashfry_index(flashfry_header_path)

    elif command == "discover" :
        return _run_flashfry_discover(flashfry_header_path, kwargs.get('sites',[]))

    elif command == "score":
        return _run_flashfry_score(flashfry_header_path, **kwargs)


    time_end = perf_counter()



    log.info("Finish running FlashFry in {}".format(timedelta(seconds=(time_end - time_start))))


def _run_flashfry_index(flashfry_header_path):
    build_index_args = ["java", "-Xmx4g", "-jar",
                        const.FLASHFRY_PATH, "index",
                        "--tmpLocation", FLASHFRY_TMP_LOCATION_PATH,
                        "--database", flashfry_header_path,
                        "--reference", COMPLETE_GENOME_PATH, "--enzyme", "spcas9ngg"]

    log.info("Starting to run FlashFry index building")
    run_external_proc(build_index_args)
    log.info("Finish running FlashFry index building")
def _run_flashfry_discover(flashfry_header_path, sites):

    flashfry_output = None

    log.info("FlashFry header location: {}".format(flashfry_header_path))

    try:
        temp_flashfry_input_file = tempfile.NamedTemporaryFile(mode='w', delete=False,
                                                     dir='{}/app/configuration_files/'.format(BASE_DIR))
        log.info("flashfry_input_file_path: {}".format(temp_flashfry_input_file.name))
        temp_flashfry_output_file = tempfile.NamedTemporaryFile(mode='w', delete=False,
                                                     dir='{}/app/configuration_files/'.format(BASE_DIR))
        log.info("flashfry_output_file_path: {}".format(temp_flashfry_output_file.name))

        temp_flashfry_input_file.writelines([">sequence_{}\n{}".format(i, site["sequence"]) for i, site in enumerate(sites)])

        temp_flashfry_input_file.close()
        temp_flashfry_output_file.close()

        discover_args = ["java", "-Xmx4g", "-jar", const.FLASHFRY_PATH, "discover",
                         "--database", flashfry_header_path,
                         "--fasta", temp_flashfry_input_file.name,
                         "--positionOutput",
                         "--output", temp_flashfry_output_file.name]

        log.info("Starting to run FlashFry discover")
        run_external_proc(discover_args)
        log.info("Finish running FlashFry discover")

        flashfry_output = pd.read_csv(temp_flashfry_output_file.name, sep="\t")

    finally:
        os.remove(temp_flashfry_input_file.name)
        os.remove(temp_flashfry_output_file.name)

    return flashfry_output

def _run_flashfry_score(flashfry_header_path, flashfry_output_file = None, flashfry_output_df = None):

    flashfry_score = None
    temp_flashfry_output_file = None

    if flashfry_output_df is None and flashfry_output_file is None:
        raise Exception('No FlashFry output to score')

    log.info("FlashFry header location: {}".format(flashfry_header_path))

    try:
        if flashfry_output_file is None:
            temp_flashfry_output_file = tempfile.NamedTemporaryFile(mode='w', delete=False,
                                                                    dir='{}/app/configuration_files/'.format(BASE_DIR))
            log.info("temp_flashfry_output_file: {}".format(temp_flashfry_output_file.name))

            flashfry_output_df.to_csv(temp_flashfry_output_file, sep='\t', index=False)
            flashfry_output_file = temp_flashfry_output_file.name
            temp_flashfry_output_file.close()

        temp_flashfry_score_output_path = tempfile.NamedTemporaryFile(mode='w', delete=False,
                                                         dir='{}/app/configuration_files/'.format(BASE_DIR))

        log.info("temp_flashfry_score_output_path: {}".format(temp_flashfry_score_output_path.name))

        score_args = ["java", "-Xmx4g", "-jar", const.FLASHFRY_PATH, "score",
                      "--input", flashfry_output_file,
                      "--output", temp_flashfry_score_output_path.name,
                      "--scoringMetrics", "doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot",
                      "--database", flashfry_header_path]

        log.info("Starting to run FlashFry score")
        run_external_proc(score_args)
        log.info("Finish running FlashFry score")


        flashfry_score = pd.read_csv(temp_flashfry_score_output_path.name, sep="\t")

    finally:
        os.remove(temp_flashfry_score_output_path.name)
        if temp_flashfry_output_file is not None:
            os.remove(temp_flashfry_output_file.name)


    return flashfry_score

