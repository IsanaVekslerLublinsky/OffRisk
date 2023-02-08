import subprocess
import sys
import argparse
import logging

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

log = logging.getLogger(__name__)


DOCKER_COMPOSE_PATH = ""
SERVER_DOCKER_PATH = ""
UI_DOCKER_PATH = ""
DB_PATH = ""


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


def main(argv):

    parser = argparse.ArgumentParser(description='Docker action: export, run')
    parser.add_argument('--export', action="store_true",
                        help='Export current docker UI, docker server and db. Save to zip file')
    args = parser.parse_args(argv)
    if args.export:
        log.info("Starting to export")
        run_external_proc(["docker-compose", "build"])
        run_external_proc(["docker", "save", "-o", "off-risk-server.tar", "off-risk-server:latest"])
        run_external_proc(["docker", "save", "-o", "off-risk-ui.tar", "off-risk-ui:latest"])
        run_external_proc(["tar", "cvzf", "off-risk-db.tar.gz", "<path>"])
        run_external_proc(["tar", "cvzf", "off-risk-files.tar.gz", "off-risk-server.tar", "off-risk-ui.tar",
                           "off-risk-ui.tar"])


if __name__ == '__main__':
    main(sys.argv[1:])
