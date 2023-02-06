import logging
import os
import re
import subprocess
import sys
import tempfile
from datetime import timedelta
from time import perf_counter
import pandas as pd
import yaml
from flask import Flask, request, make_response, jsonify
# import uvicorn
# from fastapi import FastAPI, Request, Response
# from fastapi.responses import JSONResponse
from pydantic_webargs import webargs
import warnings

from configuration_files.const import FLASHFRY_INPUT_PATH, CAS_OFFINDER_OUTPUT_PATH, CAS_OFFINDER_INPUT_FILE_PATH, \
    FLASHFRY_OUTPUT_PATH, FLASHFRY_SCORE_OUTPUT_PATH, YAML_CONFIG_FILE
from off_target import run_flashfry, run_cas_offinder_locally
from db import update_database_base_path, get_database_path
from helper import get_logger
from off_tov import extract_data
from obj_def import OffTargetList, AllDbResult, OtResponse, SitesList, DB_NAME_LIST, FlashFrySite
warnings.filterwarnings("ignore", category=RuntimeWarning)
pd.options.mode.chained_assignment = None

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(CURRENT_DIR))


BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

app = Flask(__name__)
app.config["JSON_SORT_KEYS"] = False
# app = FastAPI()

log = get_logger(logger_name=__name__, debug_level=logging.DEBUG)

conf_yaml = None
with open(YAML_CONFIG_FILE, "r") as f:
    conf_yaml = yaml.load(f, Loader=yaml.FullLoader)

update_database_base_path(conf_yaml["databases"]["base_path"])


@app.route("/")
# @app.get("/")
def hello():
    """
    Home page
    Returns: a message

    """
    return "Welcome to OffRisk server"


@app.route("/v1/off-target-analyze/", methods=["POST"])
@webargs(body=OffTargetList)
def off_target_analyze(**kwargs):
    # @app.post("/v1/off-target-analyze/")
    # def off_target_analyze(body: OffTargetList):
    """
    Analyze off-target request. Receive off-target sites.
    Returns: OtResponse object - off-target analysis

    """
    time_start = perf_counter()
    body = OffTargetList(**kwargs["payload"])
    body = body.dict()
    dbs = body["db_list"]
    log.info("Got new request: {}".format(body))

    if body["organism"] not in conf_yaml["genomes"]:
        log.info("request_id: {} - the organism {} is not supported".format(body["request_id"], body["organism"]))
        return OtResponse(request_id=body["request_id"],
                          flashfry_score="[]",
                          off_target_df="[]",
                          all_result=AllDbResult(),
                          time=0).dict()

    try:
        tempinput_file = tempfile.NamedTemporaryFile(mode='w', delete=False,
                                                     dir='{}/app/configuration_files/'.format(BASE_DIR))
        log.info("input_file_path: {}".format(tempinput_file.name))
        input_file_df = pd.DataFrame([dict(s) for s in body["off_targets"]])
        input_file_df.to_csv(tempinput_file, sep="\t", index=False)
        tempinput_file.close()
        response = analyze(dbs, body["organism"], body["request_id"], tempinput_file.name, time_start)

    finally:
        os.remove(tempinput_file.name)

    time_end = perf_counter()
    log.info("Total run: {}".format(timedelta(seconds=(time_end - time_start))))
    return response


@app.route("/v1/on-target-analyze/", methods=["POST"])
@webargs(body=SitesList)
def analyze_ot(**kwargs):
    # @app.post("/v1/on-target-analyze/")
    # def analyze_ot(body: SitesList):
    """
    Analyze on-target request. Receive sequences to search for off-target with FlashFry and Cas-Offinder,
    and then analyze them

    Returns: OtResponse object - off-target analysis
    """
    # Get the dbs
    time_start = perf_counter()
    body = SitesList(**kwargs["payload"])
    body = body.dict()
    dbs = body["db_list"]
    log.info("Got new request: {}".format(body))
    tools_list = body["search_tools"]

    if "cas_offinder" in tools_list:
        # Create the input file for cas-offinder
        seqs = ",".join(["{} {}".format(s["sequence"], s["mismatch"]) for s in body["sites"]])
        pattern = "{} {} {}".format(body["pattern"], body["pattern_dna_bulge"], body["pattern_rna_bulge"])
        genome_type = conf_yaml["cas_offinder"]["default_genome"]

        path_out = CAS_OFFINDER_OUTPUT_PATH
        if os.path.exists(path_out):
            log.info("Removing already existed output file")
            os.remove(path_out)

        path_in = CAS_OFFINDER_INPUT_FILE_PATH
        # Run Cas-Offinder
        try:
            fd_in, path_in = write_cas_offinder_input(pattern, seqs, genome_type, input_file=path_in)
            run_cas_offinder_locally("C")
            log.info("Finish running cas-offinder")
            os.remove(path_in) if os.path.exists(path_in) else None
        except Exception as e:
            log.error("An error has occurred while running cas-offinder {}".format(e))
            log.debug("input file is: {}, output file is: {}".format(path_in, path_out))

    if "flashfry" in tools_list:
        try:
            run_flashfry_from_server(body=body)
            # Remove input file
            log.info("Finish running FlashFry")
            os.remove(FLASHFRY_INPUT_PATH) if os.path.exists(FLASHFRY_INPUT_PATH) else None
        except Exception as e:
            log.error("An error has occurred while running cas-offinder {}".format(e))

    # analyze
    response = analyze(dbs, "human", body["request_id"], time_start=time_start)
    time_end = perf_counter()
    log.info("Total run: {}".format(timedelta(seconds=(time_end - time_start))))

    # Remove output files
    log.info("Cleaning the environment")
    os.remove(FLASHFRY_OUTPUT_PATH) if os.path.exists(FLASHFRY_OUTPUT_PATH) else None
    os.remove(FLASHFRY_SCORE_OUTPUT_PATH) if os.path.exists(FLASHFRY_SCORE_OUTPUT_PATH) else None
    os.remove(CAS_OFFINDER_OUTPUT_PATH) if os.path.exists(CAS_OFFINDER_OUTPUT_PATH) else None
    return response


def analyze(dbs, genome, request_id, input_file=None, time_start=perf_counter()):
    """
    Analyze the off-target with extract_data function
    Args:
        genome: The genome type (eg human)
        time_start: The time analyze started to run
        dbs: which db to analyze
        request_id: the request ID for this analyzing.
        input_file: path to input file with off-targets

    Returns:

    """
    if "all" in dbs:
        db_name_list = DB_NAME_LIST
    else:
        db_name_list = dbs

    if genome == 'human':
        off_t_result, all_result = extract_data(db_name_list=db_name_list, input_file=input_file)
        time_end = perf_counter()
        all_db_result = AllDbResult(**all_result.json)
        total_time = timedelta(seconds=(time_end - time_start))
        response = OtResponse(request_id=request_id,
                              flashfry_score=off_t_result["flashfry_score"],
                              off_targets=off_t_result["off_targets"],
                              all_result=all_db_result,
                              time=total_time.total_seconds())


    return make_response(jsonify(response.dict()), 200)
    # return response


@app.route("/v1/flashfry/", methods=["POST"])
@webargs(body=FlashFrySite)
def flashfry(**kwargs):
    # @app.post("/v1/flashfry/")
    # def flashfry(body: FlashFrySite):
    """
    Run FlashFry for discover and score option.
    Returns: Dataframe with the result from discover and score
    """
    time_start = perf_counter()
    body = FlashFrySite(**kwargs["payload"])
    body = body.dict()
    log.info("Got new request: {}".format(body))
    try:
        run_flashfry_from_server(body)
        flashfry_output = pd.read_csv(FLASHFRY_OUTPUT_PATH, sep="\t")
        flashfry_score = pd.read_csv(FLASHFRY_SCORE_OUTPUT_PATH, sep="\t")
        time_end = perf_counter()
        log.info("Total run: {}".format(timedelta(seconds=(time_end - time_start))))
        response = {"flashfry-discover": flashfry_output.to_json(), "flashfry-score": flashfry_score.to_json(),
                    "message": "Finish running FlashFry"}
        return make_response(jsonify(response), 200)
        # return response
    except Exception as e:
        message = "An error has occurred while running FlashFry {}".format(e)
        log.error(message)
        response = {"flashfry-dataframe": pd.DataFrame(), "message": message}
        return make_response(jsonify(response), 400)
        # response = {"flashfry-dataframe": pd.DataFrame().to_json(), "message": message}
        # return JSONResponse(content=response, status_code=400)

    finally:
        os.remove(FLASHFRY_OUTPUT_PATH) if os.path.exists(FLASHFRY_OUTPUT_PATH) else None
        os.remove(FLASHFRY_SCORE_OUTPUT_PATH) if os.path.exists(FLASHFRY_SCORE_OUTPUT_PATH) else None
        os.remove(FLASHFRY_INPUT_PATH) if os.path.exists(FLASHFRY_INPUT_PATH) else None


def run_flashfry_from_server(body):
    """
    Run FlashFry
    Args:
        body: FlashFrySite object
    """

    # Create the input file for flashfry
    flashfry_input_list = list()
    for i, site in enumerate(body["sites"]):
        next_line = "" if i == 0 else "\n"
        flashfry_input_list.append("{}>sequence_{}\n{}".format(next_line, i, site["sequence"]))

    flashfry_input_file_path = FLASHFRY_INPUT_PATH
    with open(flashfry_input_file_path, "w") as tmp_in:
        for line in flashfry_input_list:
            tmp_in.write(line)

    # Run FlashFry
    run_flashfry(database_path=get_database_path(), commands=["discover", "score"])


@app.route("/v1/cas-offinder-bulge/", methods=["GET"])
def cas_offinder_bulge():
    # @app.get("/v1/cas-offinder-bulge/")
    # def cas_offinder_bulge(request: Request):
    """
    Run cas-offinder-bulge
    Returns: on success A json format file of the result dataframe and the output from the program. on failure the error
    """

    log.info("Got new request for cas-offinder bulge. Starting to work")
    return run_cas_offinder_server(request.args, "bulge")
    # todo: fix this issue - convert request ars from Flask to Fastapi
    # return run_cas_offinder_server(request.query_params, "bulge")


def run_cas_offinder_server(received_request, cas_offinder_type="default"):
    """
    run cas-offinder or cas-offinder-bulge
    Args:
        received_request: the request body that was received
        cas_offinder_type: default will run cas-offinder and bulge will run cas-offinder-bulge

    Returns: on success A json format file of the result dataframe and the output from the program. on failure the error
    """

    cas_offinder_proc = ["cas-offinder"]
    if cas_offinder_type == "bulge":
        cas_offinder_proc = ["python", "/app/tools/cas-offinder-bulge"]

    log.debug("Running {}".format(cas_offinder_proc))
    validate_received_request(received_request)
    genome_type = received_request.get("genome")
    pattern = received_request.get("pattern")
    seqs = received_request.get("sequences")
    off_target_df = pd.DataFrame()
    path_in, path_out = "", ""
    try:
        fd_in, path_in = write_cas_offinder_input(pattern, seqs, genome_type)
        path_out = "{}.out".format(path_in)
        run_external_proc(cas_offinder_proc + [path_in, "C", path_out])
        if os.path.exists(path_out):
            off_target_df = pd.read_csv(path_out, sep="\t")
            log.info("Finish running cas-offinder")
    except Exception as e:
        message = "An error has occurred while running cas-offinder {}".format(e)
        log.error(message)
        response = {"cas-offinder-dataframe": pd.DataFrame().to_json(), "message": message}
        return make_response(jsonify(response), 400)
        # response = {"cas-offinder-dataframe": pd.DataFrame().to_json(), "message": message}
        # return JSONResponse(content=response, status_code=400)
    finally:
        os.remove(path_in) if os.path.exists(path_in) else None
        os.remove(path_out) if os.path.exists(path_out) else None
    response = {"cas-offinder-dataframe": off_target_df.to_json(), "message": "Finish running cas-offinder-buge"}
    return make_response(jsonify(response), 200)
    # return response


def write_cas_offinder_input(pattern, seqs, genome_type, input_file=None):
    """
    Create input file for cas-offinder
    Args:
        pattern: desired pattern including PAM site and optional DNA or RNA bulge sizes, separated by spaces.
        seqs: query sequences and maximum mismatch numbers, separated by spaces, each sequence is separate with ","
        genome_type: which genome type to run. default is human GRCH38
        input_file: path to input file. if None then create the path here.

    Returns: io file and path
    """
    if genome_type not in conf_yaml["genomes"]:
        log.info("genome_type {} is not supported".format(genome_type))
        raise Exception("No such genome {}. Allowed genome: {}".format(genome_type, list(conf_yaml["genomes"].keys())))

    docker_path_to_genome = conf_yaml["genomes"][genome_type]
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


def validate_received_request(received_request):
    """
    Validate the received request for Cas-Offinder
    Args:
        received_request: has the following fields
        genome - which genome to work on. for example "human"
        pattern -  desired pattern including PAM site and optional DNA or RNA bulge sizes, separated by spaces.
        For example "NNNNNNNNNNNNNNNNNNNNNGG 0 0"
        sequences -  query sequences and maximum mismatch numbers, separated by spaces, each sequence is
        separate with ",". For example: "CTTAAGAATACGCGTAGTCGAGG 4"
    Returns: Raise an error when validation fails

    """
    # Validate genome
    genome_list = ["human"]
    if not received_request.get("genome") in genome_list:
        raise Exception("No such genome {}. Allowed genome: {}".format(received_request.get("genome"), genome_list))
    # Validate pattern
    pattern = received_request.get("pattern")
    if not re.fullmatch(r"^[AGTCRYSWKMBDHVN]*\s[0-9]+\s[0-9]+", pattern):
        raise Exception("Not a valid pattern. Please refer to cas-offinder documentation.")
    # Validate sequences
    seqs = received_request.get("sequences")
    if seqs:
        seqs = seqs.split(",")
        for seq in seqs:
            if not re.fullmatch(r"^[AGTCRYSWKMBDHVN]*\s[0-9]+(\s)*[A-Za-z0-9]*", seq):
                raise Exception("Not a valid sequence. Please refer to cas-offinder documentation.")


def run_external_proc(args):
    """
    Run an external program with Popen.
    Args:
        args: The argument to run, a list of string

    Returns: return the output of the program

    """
    log.info("The following command will be run: {}\nThis process will take some time.".format(args))
    result = subprocess.run(args, capture_output=True, text=True)
    result_message = "STDOUT: {}\nSTDERR: {}".format(result.stdout, result.stderr)
    log.info(result_message)
    if result.returncode != 0:
        raise Exception("Error while running Cas-offinder")
    return result_message


#  this doesn't get executed inside the flask docker
if __name__ == "__main__":
    log.info("Starting to run off-risk server")

    # Only for debugging while developing
    app.run(host="0.0.0.0", debug=True, port=8002)
    # uvicorn.run("app.main:app", host="0.0.0.0", port=8002)
