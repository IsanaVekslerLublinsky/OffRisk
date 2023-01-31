import logging
import yaml
from db import GffDB
from helper import get_logger

log = get_logger(logger_name=__name__, debug_level=logging.DEBUG)


class AGenomeDatabaseFactory:

    def create_database(self, db_type, genome, db_name):
        pass


class CorporateGenomeDatabaseFactory(AGenomeDatabaseFactory):

    def __init__(self, conf_file):
        self._conf_file = conf_file
        pass

    def _get_db_attributes(self, genome, db_name):
        location = None
        columns = None
        with open(self._conf_file) as f:
            conf_yaml = yaml.load(f, Loader=yaml.FullLoader)
            if db_name in conf_yaml["databases"] and genome in conf_yaml["databases"][db_name]:
                location = conf_yaml["databases"][db_name][genome]["path"]
                columns = conf_yaml["databases"][db_name][genome]["columns"]

        return location, columns

    def create_database(self, db_type, genome, db_name):
        db = None
        if db_type == "crispr":
            (location, columns) = self._get_db_attributes(genome, db_name)
            if location is not None and columns is not None:
                db = GffDB(db_name, location, columns)
        elif db_type == "gil-ad":
            # db = DB()
            pass
        return db
