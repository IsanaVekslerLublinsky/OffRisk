import datetime
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
APP_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOG_PATH = "{}/log/run.log".format(BASE_DIR)
CURRENT_DATA_TIME = str(datetime.datetime.now().strftime("%Y-%m-%d-%H:%M"))
CONF_FILE = "{}/configuration_files/conf_param.json".format(APP_DIR)
DB_NAME_LIST = ["gencode", "mirgene", "remapepd", "enhanceratlas", "pfam", "targetscan", "omim", "humantf",
                "protein_atlas", "rbp", "cosmic"]

YAML_CONFIG_FILE = "{}/configuration_files/off-risk-config.yaml".format(APP_DIR)
"""
This file store default values and path for different variables used by this program.
You can change the default values from here.
"""
# Cas-OFFinder - default path in docker
CAS_OFFINDER_BULGE_PATH = "/app/tools/cas-offinder-bulge"

# Cas-OFFinder was installed from 'conda install -c bioconda cas-offinder'.
# For other uses please specify the path for cas-offinder bellow
CAS_OFFINDER_PATH = "/app/tools/cas-offinder"

CAS_OFFINDER_INPUT_FILE_PATH = "{}/configuration_files/cas_offinder_input_bulge.txt".format(APP_DIR)
CAS_OFFINDER_OUTPUT_PATH = "{}/off_target_output/cas_offinder_output_bulge.txt".format(APP_DIR)


# FlashFry - - default path in docker
FLASHFRY_PATH = "/app/tools/FlashFry-assembly-1.12.jar"
FLASHFRY_INPUT_PATH = "{}/configuration_files/flashfry_input.fa".format(APP_DIR)
FLASHFRY_TMP_LOCATION_PATH = "/app/tmp/"
FLASHFRY_DATABASE_BASE_PATH = "/FlashFry/cas9ngg_database"
FLASHFRY_OUTPUT_PATH = "{}/off_target_output/flashfry_output.output".format(APP_DIR)
FLASHFRY_SCORE_OUTPUT_PATH = "{}/off_target_output/flashfry_output.output.scored".format(APP_DIR)

# Database location
COMPLETE_GENOME_PATH = "/databases/genome/hg38.fa"
COMPLETE_GENOME_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
