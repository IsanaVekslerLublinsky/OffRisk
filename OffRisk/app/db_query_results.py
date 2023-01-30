import pandas as pd


class IDBResult(object):
    # return pandas dataframe
    def get_results(self):
        pass


class DBDataFrameResult(IDBResult):

    def __init__(self, results):
        assert isinstance(results, pd.DataFrame)
        self.result_df = results

    def get_results(self):
        return self.result_df
