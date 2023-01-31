from __future__ import annotations
import logging
import pandas as pd
from abc import abstractmethod
from db import GffDB, DataFrameDB
from helper import get_logger
log = get_logger(logger_name=__name__, debug_level=logging.DEBUG)


class ADBQuery(object):

    @abstractmethod
    def execute_gffdb(self, db: GffDB):
        pass

    @abstractmethod
    def execute_dataframedb(self, db: DataFrameDB):
        pass


class BedQuery(ADBQuery):

    def __init__(self, query):
        self._query = query

    @property
    def query(self):
        return self._query

    def execute_gffdb(self, db: GffDB):
        results = None
        intersection_bed = db.db_bed.intersect(self.query, wb=True, wa=True, sorted=True)
        # Convert intersection result to dataframe
        if intersection_bed.count() != 0:
            intersection_group = intersection_bed.to_dataframe(header=None, names=db.columns)
            results = intersection_group

        else:
            # if count is 0, it means that there is not result for the intersection
            log.info("There are no result from {} intersection".format(db.name))

        return results

    def execute_dataframedb(self, db: DataFrameDB):
        return None


class DataFrameQuery(ADBQuery):
    def __init__(self, query):
        self._query = query

    def _separate_attributes(self, df_to_separate):
        """
        Separate the attributes to different columns
        :param df_to_separate: A dictionary of the result to separate. key is index, value is dataframe
        :return: A dictionary the result intersection with attribute separate. key is index, value is dataframe
        """

        function_name = "separate_attributes"
        log.debug("Entering {}".format(function_name))

        # for result_key in off_target_result_dict:
        # Separate the attributes column
        if not df_to_separate.empty:
            df_to_separate.loc[:, "separate"] = df_to_separate.loc[:, "attributes"].apply(
                lambda s: self._separate_attributes(s))
            df_attributes = pd.DataFrame(df_to_separate["separate"].values.tolist(), index=df_to_separate.index)
            result = pd.concat([df_to_separate, df_attributes], axis=1).drop("separate", axis=1)
            return result
        else:
            return df_to_separate

    def execute_gffdb(self, db: GffDB):
        results = None
        intersection_bed = self.db_bed.intersect(db.db_bed, wb=True, wa=True, sorted=True)
        # Convert intersection result to dataframe
        if intersection_bed.count() != 0:
            intersection_group = intersection_bed.to_dataframe(header=None, names=self.columns_name)
            if self.separate_attributes:
                intersection_group = self._separate_attributes(intersection_group)
            results = intersection_group
            self.pr_df.append({"name": self.db_name, "description": "Complete result",
                               "data": intersection_group.to_json(orient="columns")})
        else:
            # if count is 0, it means that there is not result for the intersection
            log.info("There is no result from {} intersection".format(self.db_name))

        return results

    def execute_dataframedb(self, db: DataFrameDB):
        return None
