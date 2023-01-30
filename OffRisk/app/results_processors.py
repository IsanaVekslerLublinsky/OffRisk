import logging
import pandas as pd
from db import separate_attributes
from helper import get_logger

log = get_logger(logger_name=__name__, debug_level=logging.DEBUG)


class IDBProcessResult:
    def __init__(self):
        pass

    def process(self, db_result):
        pass


class DBDataFrameProcessedResult(IDBProcessResult):

    def __init__(self):
        pass

    def process(self, db_result, selected_attributes):
        tmp_df = separate_attributes(db_result.get_results())

        tmp_df = self._segment_filtered_result(tmp_df)

        dict_column_rename = {item['original_name']: item['return_name'] for item in selected_attributes}

        agg_dict = dict()
        for attribute in dict_column_rename.keys():
            agg_dict[attribute] = lambda x: list(x)

        tmp_df = tmp_df.groupby("off_target_id", as_index=False).agg(agg_dict)
        tmp_df = tmp_df.rename(columns=dict_column_rename)

        return tmp_df

    def _segment_filtered_result(self, df):
        transcripts_types = ['ncRNA', 'nmRNA', 'sRNA', 'smnRNA', 'tRNA', 'sRNA', 'mRNA', 'pcRNA', 'rRNA', '5S rRNA',
                             '5.8S rRNA',
                             'SSU rRNA', 'LSU rRNA', 'NoRC RNA', 'pRNA', '6S RNA', 'SsrS RNA', 'aRNA', 'asRNA',
                             'asmiRNA',
                             'cis-NAT', 'crRNA', 'tracrRNA', 'CRISPR RNA', 'DD RNA', 'diRNA', 'dsRNA', 'endo-siRNA',
                             'exRNA',
                             'gRNA', 'hc-siRNA', 'hcsiRNA', 'hnRNA', 'RNAi', 'lincRNA', 'lncRNA', 'miRNA', 'mrpRNA',
                             'nat-siRNA',
                             'natsiRNA', 'OxyS RNA', 'piRNA', 'qiRNA', 'rasiRNA', 'RNase MRP', 'RNase P', 'scaRNA',
                             'scnRNA',
                             'scRNA', 'scRNA', 'SgrS RNA', 'shRNA', 'siRNA', 'SL RNA', 'SmY RNA', 'snoRNA', 'snRNA',
                             'snRNP',
                             'SPA lncRNA', 'SRP RNA', 'ssRNA', 'stRNA', 'tasiRNA', 'tmRNA', 'uRNA', 'vRNA', 'vtRNA',
                             'Xist RNA',
                             'Y RNA', 'NATs', 'pre-mRNA', 'circRNA', 'msRNA', 'cfRNA', 'transcript']

        group_gencode = df.groupby(["off_target_id", "gene"])
        index_list_to_remove = list()
        for group_name, group_df in group_gencode:
            df_tmp = pd.DataFrame()
            available_segment = group_df["segment"].unique().tolist()
            final_level = ["CDS", "five_prime_UTR", "three_prime_UTR"]
            if final_level in available_segment:
                df_tmp = group_df[group_df["segment"] in final_level]
            elif "exon" in available_segment:
                df_tmp = group_df[group_df["segment"] == "exon"]
            elif "transcript" in available_segment:
                df_tmp = group_df[group_df["segment"].isin(transcripts_types)]
            elif "gene" in available_segment:
                df_tmp = group_df[group_df["segment"] == "gene"]

            if group_name[1] == 'None':
                df_index_to_remove = group_df.index.tolist()
                index_list_to_remove.extend(df_index_to_remove)

            elif not df_tmp.empty:
                df_index_to_remove = group_df.loc[~group_df.index.isin([df_tmp.index[0]])].index.tolist()
                index_list_to_remove.extend(df_index_to_remove)


        return df.drop(index=index_list_to_remove)
