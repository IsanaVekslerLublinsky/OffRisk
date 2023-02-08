from __future__ import annotations
import numpy as np
import pandas as pd
import logging
import os
import re
import wget
from flask import jsonify, Flask
from pybedtools import BedTool
from configuration_files.const import BASE_DIR, COMPLETE_GENOME_URL, COMPLETE_GENOME_PATH
from abc import abstractmethod
from helper import get_logger, extract_gz_file

log = get_logger(logger_name=__name__, debug_level=logging.DEBUG)
database_base_path = "{}/databases".format(BASE_DIR)


class Db(object):

    def __init__(self, name, file_path=None):
        """
        Initialize default db
        Args:
            name: The name of the DB
            file_path: a path for all DB files
        """
        self.db_name = name
        self.file_path = file_path
        self.complete_result = pd.DataFrame()

        self.separate_attributes = False
        self.columns_name = None
        self.dtype_to_intersect = {}
        self.column_name_for_intersection = "gene_ensembl_id"
        self.db_bed = None
        self.db_df = pd.DataFrame()
        self.pr_df = list()

        if not os.path.exists(self.file_path):
            log.error("File for db {} in {} does not exist".format(self.db_name, self.file_path))
            self.is_data_loaded = 0
            return
        self.load_data()

    def load_data(self):
        """
        Load the data
        """
        pass

    def analyze(self, off_target_bed):
        """
        Analyze the off-target with the db
        Args:
            off_target_bed: bed file of the off-target

        Returns: save the result in "complete_result"

        """
        if type(self.db_bed) is not BedTool:
            log.error("The data for {} was not loaded".format(self.db_name))
            return
        log.info("Starting to analyze {}".format(self.db_name))
        # Run BedTools intersection
        intersection_bed = self.db_bed.intersect(off_target_bed, wb=True, wa=True, sorted=True)
        # Convert intersection result to dataframe
        if intersection_bed.count() != 0:
            intersection_group = intersection_bed.to_dataframe(header=None, names=self.columns_name,
                                                               dtype=self.dtype_to_intersect)
            if self.separate_attributes:
                intersection_group = separate_attributes(intersection_group)
            self.complete_result = intersection_group
            self.pr_df.append({"name": self.db_name, "description": "Complete result",
                               "data": intersection_group.to_json(orient="columns")})
        else:
            # if count is 0, it means that there is not result for the intersection
            log.info("There is no result from {} intersection".format(self.db_name))

    def get_db_result(self):
        """
        get analyze result
        Returns: dataframe of the complete result

        """
        return self.complete_result

    def get_db_name(self):
        """

        Returns: db name

        """
        return self.db_name

    def get_pr_df(self):
        """

        Returns: pr_df which is the process result

        """
        return self.pr_df

    def update_db(self):
        """
        Update the database if possible
        Returns:

        """
        log.info("Update function is not implemented for {}".format(self.db_name))

    def process_result(self, off_target_df):
        return off_target_df


class GencodeDb(Db):

    def __init__(self, file_path, final_columns):
        file_path = "{}/{}".format(database_base_path, file_path)
        super().__init__("GENCODE", file_path)
        self.url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gff3.gz"
        self.separate_attributes = True
        self.columns_name = ["chromosome", "source", "segment", "start", "end", "score", "strand", "frame",
                             "attributes", "ot_chromosome", "ot_start", "ot_end", "off_target_id", "risk_score",
                             "ot_strand", "ot_attributes", "ot_id", "ot_seq"]
        self.dtype_to_intersect = {"chromosome": str, "ot_chromosome": str}
        self.final_columns = final_columns

    def load_data(self):
        """
        Load the GENCODE file v42 from GRCh38. Downloaded on 6.11.2022.
        :return: pybedtools object with GENCODE information
        """
        self.db_bed = BedTool(self.file_path)

    def update_db(self):
        gencode_file = wget.download(url=self.url, out="{}/GENCODE".format(database_base_path))
        if os.path.exists(self.file_path["GENCODE"]):
            os.rename(self.file_path["GENCODE"], "{}.old".format(self.file_path["GENCODE"]))
        extract_gz_file(gencode_file, self.file_path["GENCODE"])
        log.info("{} was downloaded to: {}".format(self.db_name, gencode_file))

    def process_result(self, off_target_df):
        """
        Returns:

        """
        if len(self.complete_result) != 0:
            # Convert gene id to be without the decimal
            self.complete_result["gene_id"] = self.complete_result["gene_id"].apply(
                lambda x: re.sub("(\:.*$)", "", re.sub("^[^:]*:", "", re.sub("\.\d+", "", x))))

            self.complete_result.rename(columns={"gene_name": "gene_symbol", "gene_id": "gene_ensembl_id",
                                                 "transcript_id": "transcript_ensembl_id",
                                                 "transcript_name": "transcript_symbol"}, inplace=True)

            self.segment_filtered_result()

            # Update final result
            column_to_keep = sorted(set(self.final_columns).intersection(self.complete_result.columns),
                                    key=lambda x: self.final_columns.index(x))
            self.complete_result = self.complete_result[column_to_keep]
            self.pr_df[0] = {"name": self.db_name, "description": "Complete result",
                             "data": self.complete_result.to_json(orient="records")}

            # Update relevant fields in global off_target dataframe
            intersection_group = self.complete_result.groupby("off_target_id", as_index=False).agg(
                {"gene_ensembl_id": lambda x: list(set(x)), "gene_symbol": lambda x: list(set(x)),
                 "segment": lambda x: list(set(x))})
            off_target_df = merge_off_target_information(off_target_df, intersection_group,
                                                         ["gene_ensembl_id", "gene_symbol", "segment"])
        return off_target_df

    def segment_filtered_result(self):
        group_gencode = self.complete_result.groupby(["off_target_id", "gene_ensembl_id"])
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
                df_tmp = group_df[group_df["segment"] == "transcript"]
            elif "gene" in available_segment:
                df_tmp = group_df[group_df["segment"] == "gene"]
            if len(df_tmp.index) != 0:
                df_index_to_remove = group_df.loc[~group_df.index.isin([df_tmp.index[0]])].index.tolist()
                index_list_to_remove.extend(df_index_to_remove)
        self.complete_result.drop(index=index_list_to_remove, inplace=True)


class MirGeneDB(Db):

    def __init__(self, file_path, final_columns):
        file_path = "{}/{}".format(database_base_path, file_path)
        super().__init__("MirGene", file_path)
        self.url = "https://mirgenedb.org/static/data/hsa/hsa-all.bed"

        self.columns_name = ["chromosome", "start", "end", "mir_symbol", "score", "strand", "ot_chromosome",
                             "ot_start", "ot_end", "off_target_id", "risk_score", "ot_strand",
                             "ot_attributes", "ot_id", "ot_seq"]
        self.dtype_to_intersect = {"chromosome": str, "ot_chromosome": str}
        self.final_columns = final_columns

    def load_data(self):
        """
        Load the MirGene db. Version 2.1, downloaded on 6/11/2022
        :return: pybedtools object with MirGene information
        """
        log.debug("Loading MirGene data")
        self.db_bed = BedTool(self.file_path)

    def update_db(self):
        """
        update MirGene DB
        Returns:

        """
        if os.path.exists(self.file_path["MirGene"]):
            os.rename(self.file_path["MirGene"], "{}.old".format(self.file_path["MirGene"]))
        mirgene_file = wget.download(url=self.url, out=self.file_path["MirGene"])
        log.info("{} was downloaded to: {}".format(self.db_name, mirgene_file))

    def process_result(self, off_target_df):
        """
        Process the result of MirGeneDB.
        Add information to global off-target dataframe
        """
        if len(self.complete_result.index) != 0:

            # from each group take one group type
            final_df_index_list = list()
            self.complete_result[["group_name", "group_type"]] = self.complete_result["mir_symbol"].str.rsplit(
                "-", n=1, expand=True)
            self.complete_result["group_range"] = self.complete_result[["start", "end"]].apply(
                lambda s: s.end - s.start + 1, axis=1)
            df_group = self.complete_result.groupby("group_name")
            for name, group in df_group:
                group = group.sort_values("group_range")
                index = list(group.index)[0]
                for row in group.itertuples():
                    if (row.ot_start >= row.start) & (row.ot_end <= row.end):
                        index = row.Index
                        break
                final_df_index_list.append(index)
            df_index_to_remove = self.complete_result.loc[~self.complete_result.index.isin(final_df_index_list)]. \
                index.tolist()
            self.complete_result.drop(index=df_index_to_remove, inplace=True)

            # Update final result
            self.complete_result = self.complete_result[self.final_columns]
            self.pr_df[0] = {"name": self.db_name, "description": "Complete result",
                             "data": self.complete_result.to_json(orient="records")}

            # Update relevant fields in global off_target dataframe
            off_target_mirgene = self.complete_result.groupby("off_target_id").agg(
                {"mir_symbol": lambda x: list(set(x))})
            off_target_mirgene.rename(columns={"mir_symbol": "mir_gene"}, inplace=True)
            off_target_df = merge_off_target_information(off_target_df, off_target_mirgene, ["mir_gene"])

        return off_target_df


class ReMapEPD(Db):

    def __init__(self, file_path, final_columns):
        file_path = "{}/{}".format(database_base_path, file_path)
        super().__init__("ReMapEPD", file_path)
        self.separate_attributes = True
        self.columns_name = ["chromosome", "start", "end", "attributes", "score", "strand", "ot_chromosome",
                             "ot_start", "ot_end", "off_target_id", "risk_score", "ot_strand", "ot_attributes",
                             "ot_id", "ot_seq"]
        self.dtype_to_intersect = {"chromosome": str, "ot_chromosome": str}
        self.final_columns = final_columns

    def load_data(self):
        """
        Load the EPD file.
        Downloaded on the 7/11/2022
        :return: pybedtools object with epd information
        """
        log.debug("Loading ReMap and EPD data")
        self.db_bed = BedTool(self.file_path)

    def process_result(self, off_target_df):
        """
        Process the result of ReMap and EPD.
        Add information to global off-target dataframe
        """
        if len(self.complete_result.index) != 0:
            # Update final result
            self.complete_result = self.complete_result[self.final_columns]
            self.pr_df[0] = {"name": self.db_name, "description": "Complete result",
                             "data": self.complete_result.to_json(orient="records")}

            # Update relevant fields in global off_target dataframe
            off_target_remap_epd = self.complete_result.groupby("off_target_id").agg(
                {"gene_ensembl_id": lambda x: list(set(x))})
            off_target_remap_epd.rename(columns={"gene_ensembl_id": "remap_epd_gene_ensembl_id"},
                                        inplace=True)
            off_target_df = merge_off_target_information(off_target_df, off_target_remap_epd,
                                                         ["remap_epd_gene_ensembl_id"])
        return off_target_df


class EnhancerAtlas(Db):

    def __init__(self, file_path, final_columns):
        file_path = "{}/{}".format(database_base_path, file_path)
        super().__init__("EnhancerAtlas", file_path)
        self.columns_name = ["chr_enhancer", "enh_start", "enh_stop", "name", "ot_chromosome",
                             "ot_start", "ot_end", "off_target_id", "risk_score", "ot_strand",
                             "ot_attributes", "ot_id", "ot_seq"]
        self.dtype_to_intersect = {"chr_enhancer": str, "ot_chromosome": str}
        self.final_columns = final_columns

    def load_data(self):
        """
        Load the EnhancerAtlas db.
        Downloaded on 7/11/2022
        :return: pybedtools object with EnhancerAtlas information
        """
        log.debug("Loading Enhancer Atlas data")
        self.db_bed = BedTool(self.file_path)

    def process_result(self, off_target_df):
        """
        Process the result of MirGeneDB.
        Add information to global off-target dataframe
        """
        if len(self.complete_result.index) != 0:
            self.complete_result["name"] = self.complete_result["name"].apply(lambda s: s.split(","))
            self.complete_result = self.complete_result.explode("name")
            self.complete_result.reset_index(drop=True, inplace=True)
            self.complete_result["name"] = self.complete_result["name"].apply(lambda s: s.strip("[\'\']").strip(" \'"))
            column_name = ["gene_ensembl_id", "gene_symbol", "gene_chrom", "gene_start", "gene_strand",
                           "enhancer_gene_score", "tissue"]
            new_information = self.complete_result["name"].str.split("$", expand=True)
            new_information.columns = column_name
            self.complete_result = pd.concat([self.complete_result[:], new_information[:]], axis=1)

            # Update final result
            self.complete_result = self.complete_result[self.final_columns]
            self.pr_df[0] = {"name": self.db_name, "description": "Complete result",
                             "data": self.complete_result.to_json(orient="records")}

            # Update relevant fields in global off_target dataframe
            off_target_enhancer_atlas = self.complete_result.groupby("off_target_id").agg(
                {"gene_ensembl_id": lambda x: list(set(x))})
            off_target_enhancer_atlas.rename(columns={"gene_ensembl_id": "enhancer_atlas_gene_ensembl_id"},
                                             inplace=True)
            off_target_df = merge_off_target_information(off_target_df, off_target_enhancer_atlas,
                                                         ["enhancer_atlas_gene_ensembl_id"])

        return off_target_df


class Pfam(Db):

    def __init__(self, file_path, final_columns):
        file_path = "{}/{}".format(database_base_path, file_path)
        super().__init__("Pfam", file_path)
        self.separate_attributes = True
        self.columns_name = ["chromosome", "start", "end", "attributes", "ot_chromosome",
                             "ot_start", "ot_end", "off_target_id", "risk_score", "ot_strand",
                             "ot_attributes", "ot_id", "ot_seq"]
        self.dtype_to_intersect = {"chromosome": str, "ot_chromosome": str}
        self.final_columns = final_columns

    def load_data(self):
        """
        Load the Pfam db. Downloaded on 7/11/2022.
        :return: pybedtools object with Pfam protein domains information
        """
        log.debug("Loading Pfam protein domain data")
        self.db_bed = BedTool(self.file_path)

    def process_result(self, off_target_df):
        if len(self.complete_result.index) != 0:
            # Update final result
            self.complete_result = self.complete_result[self.final_columns]
            self.pr_df[0] = {"name": self.db_name, "description": "Complete result",
                             "data": self.complete_result.to_json(orient="records")}

            # Update relevant fields in global off_target dataframe
            off_target_pfam = self.complete_result.groupby("off_target_id").agg(
                {"pfam_domain_name": lambda x: list(set(x))})
            off_target_pfam.rename(columns={"pfam_domain_name": "pfam_protein_domains"}, inplace=True)
            off_target_df = merge_off_target_information(off_target_df, off_target_pfam, ["pfam_protein_domains"])
        return off_target_df


class TargetScan(Db):

    def __init__(self, file_path, final_columns):
        file_path = "{}/{}".format(database_base_path, file_path)
        super().__init__("TargetScan", file_path)
        self.columns_name = ["chromosome", "start", "end", "name", "score", "strand", "ot_chromosome",
                             "ot_start", "ot_end", "off_target_id", "risk_score", "ot_strand",
                             "ot_attributes", "ot_id", "ot_seq"]
        self.dtype_to_intersect = {"chromosome": str, "ot_chromosome": str}
        self.final_columns = final_columns

    def load_data(self):
        """
        Load the TargetScan db.
        :return: pybedtools object with TargetScan information
        """
        log.debug("Loading TargetScan protein domain data")
        self.db_bed = BedTool(self.file_path)

    def process_result(self, off_target_df):
        if len(self.complete_result.index) != 0:
            self.complete_result[["gene_symbol", "mir_symbol", "gene_ensembl_id"]] = \
                self.complete_result["name"].str.rsplit(":", expand=True)

            # Update final result
            self.complete_result = self.complete_result[self.final_columns]
            self.pr_df[0] = {"name": self.db_name, "description": "Complete result",
                             "data": self.complete_result.to_json(orient="records")}

            # Update relevant fields in global off_target dataframe
            off_target_targetscan = self.complete_result.groupby("off_target_id").agg(
                {"mir_symbol": lambda x: list(set(x))})
            off_target_targetscan.rename(columns={"mir_symbol": "targetscan"}, inplace=True)
            off_target_df = merge_off_target_information(off_target_df, off_target_targetscan, ["targetscan"])
        return off_target_df


class OmimDb(Db):

    def __init__(self, file_path, final_columns):
        file_path = "{}/{}".format(database_base_path, file_path)
        super().__init__("OMIM", file_path)
        self.enhancer_atlas = pd.DataFrame()
        self.remap_epd = pd.DataFrame()
        self.final_columns = final_columns

    def load_data(self):
        """
        Load the OMIM db. downloaded on 24/12/2020.
        :return: dataframe with OMIM information
        """
        self.db_df = pd.read_csv(self.file_path, sep="\t")

    def process_result(self, off_target_df):
        if len(self.complete_result.index) != 0:
            self.complete_result = self.complete_result.astype({"disease_related": "string",
                                                                "inheritance_model": "string"})
            self.complete_result[["disease_related", "inheritance_model"]].replace(np.nan, "", inplace=True)

            # Update final result
            self.complete_result = self.complete_result[self.final_columns]
            self.pr_df[0] = {"name": self.db_name, "description": "Complete result",
                             "data": self.complete_result.to_json(orient="records")}

            # Update relevant fields in global off_target dataframe
            intersection_group = self.complete_result.groupby("off_target_id", as_index=False).agg(
                {"disease_related": lambda x: list(set(x)),
                 "inheritance_model": lambda x: list(set(x))})
            off_target_df = merge_off_target_information(off_target_df,
                                                         intersection_group.astype({"off_target_id": int}),
                                                         ["disease_related", "inheritance_model"])

            # Update fields for Enhancer Atlas
            if len(self.enhancer_atlas.index) != 0:
                self.enhancer_atlas = self.enhancer_atlas.astype({"disease_related": "string",
                                                                  "inheritance_model": "string"})
                self.enhancer_atlas[["disease_related", "inheritance_model"]].replace(np.nan, "", inplace=True)
                intersection_group = self.enhancer_atlas.groupby("off_target_id", as_index=False).agg(
                    {"disease_related": lambda x: list(set(x)),
                     "inheritance_model": lambda x: list(set(x))})
                intersection_group = intersection_group.rename(
                    columns={"disease_related": "enhancer_atlas_disease_related",
                             "inheritance_model": "enhancer_atlas_inheritance_model"})
                off_target_df = merge_off_target_information(off_target_df,
                                                             intersection_group.astype({"off_target_id": int}),
                                                             ["enhancer_atlas_disease_related",
                                                              "enhancer_atlas_inheritance_model"])

            # Update fields for ReMap EPD
            if len(self.remap_epd.index) != 0:
                self.remap_epd = self.remap_epd.astype({"disease_related": "string",
                                                        "inheritance_model": "string"})
                self.remap_epd[["disease_related", "inheritance_model"]].replace(np.nan, "", inplace=True)
                intersection_group = self.remap_epd.groupby("off_target_id", as_index=False).agg(
                    {"disease_related": lambda x: list(set(x)),
                     "inheritance_model": lambda x: list(set(x))})
                intersection_group = intersection_group.rename(
                    columns={"disease_related": "remap_epd_disease_related",
                             "inheritance_model": "remap_epd_inheritance_model"})
                off_target_df = merge_off_target_information(off_target_df,
                                                             intersection_group.astype({"off_target_id": int}),
                                                             ["remap_epd_disease_related",
                                                              "remap_epd_inheritance_model"])
        return off_target_df


class HumanTFDb(Db):

    def __init__(self, file_path, final_columns):
        file_path = "{}/{}".format(database_base_path, file_path)
        super().__init__("HumanTF", file_path)
        human_tf_url = "http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Homo_sapiens_TF"
        human_tf_cofactor_url = \
            "http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Homo_sapiens_TF_cofactors"
        self.url = {"TF_PATH": human_tf_url, "TF_COFACTORS_PATH": human_tf_cofactor_url}
        self.final_columns = final_columns

    def load_data(self):
        """
        Load the HumanTFDb db. downloaded on 7/11/2022.
        :return: dataframe with HumanTFDb information
        """
        log.debug("Loading Human TF data")
        self.db_df = pd.read_csv(self.file_path, sep="\t")

    def process_result(self, off_target_df):
        if len(self.complete_result.index) != 0:
            # Update final result
            self.complete_result = self.complete_result[self.final_columns]
            self.pr_df[0] = {"name": self.db_name, "description": "Complete result",
                             "data": self.complete_result.to_json(orient="records")}

            # Update relevant fields in global off_target dataframe
            off_target_df = merge_off_target_information(off_target_df,
                                                         self.complete_result[["off_target_id", "HumanTF_source"]].astype(
                                                             {"off_target_id": int, "HumanTF_source": str}),
                                                         ["HumanTF_source"],
                                                         "string")
        return off_target_df

    def update_db(self):
        if os.path.exists(self.file_path["TF_PATH"]):
            os.rename(self.file_path["TF_PATH"], "{}.old".format(self.file_path["TF_PATH"]))
        human_tf_file = wget.download(url=self.url["TF_PATH"], out=self.file_path["TF_PATH"])
        log.info("{} cofactors was downloaded to: {}".format(self.db_name, human_tf_file))

        if os.path.exists(self.file_path["TF_COFACTORS_PATH"]):
            os.rename(self.file_path["TF_COFACTORS_PATH"], "{}.old".format(
                self.file_path["TF_COFACTORS_PATH"]))
        human_tf_co_file = wget.download(url=self.url["TF_COFACTORS_PATH"],
                                         out=self.file_path["TF_COFACTORS_PATH"])
        log.info("{} was downloaded to: {}".format(self.db_name, human_tf_co_file))


class ProteinAtlas(Db):

    def __init__(self, file_path):
        file_path = "{}/{}".format(database_base_path, file_path)
        super().__init__("Protein_Atlas", file_path)

    def load_data(self):
        """
        Load the ProteinA tlas file.
        Downloaded on the 7/11/2022
        :return: dataframe with Protein Atlas information
        """
        log.debug("Loading Protein Atlas data")
        self.db_df = pd.read_csv(self.file_path)

    def process_result(self, off_target_df):
        """
        Process the result of Protein Atlas.
        Add information to global off-target dataframe
        """
        if len(self.complete_result.index) != 0:
            df_to_vector = self.complete_result.fillna("None")

            # create dictionary with value to integer mappings
            value_to_int = {value: i for i, value in
                            enumerate(["Not representative", "None", "Not detected", "Low", "Medium", "High"])}

            df_to_vector = df_to_vector.replace(value_to_int).drop(
                ["gene_ensembl_id", "gene_symbol", "off_target_id"], axis=1)
            df_to_vector["expression_information"] = df_to_vector.apply(
                lambda row: ",".join(row.values.astype(str)), axis=1)
            self.complete_result = pd.concat([self.complete_result, df_to_vector], axis=1)
            off_target_protein_atlas = self.complete_result.astype(
                {"off_target_id": int, "expression_information": "string", "gene_symbol": "string"})
            off_target_protein_atlas["expression_information"] = \
                off_target_protein_atlas[["gene_symbol", "expression_information"]].apply(lambda row: "{}:({})".format(
                    row["gene_symbol"], row["expression_information"]), axis=1)

            off_target_protein_atlas = off_target_protein_atlas.groupby("off_target_id", as_index=False).agg(
                {"expression_information": lambda x: list(set(x))})
            off_target_df = merge_off_target_information(off_target_df, off_target_protein_atlas,
                                                         ["expression_information"])

        return off_target_df


class RBP(Db):

    def __init__(self, file_path, final_columns):
        file_path = "{}/{}".format(database_base_path, file_path)
        super().__init__("RBP", file_path)
        self.final_columns = final_columns

    def load_data(self):
        """
        Load the RBP file.
        Downloaded on the 7/11/2022
        :return: dataframe with RBP information
        """
        log.debug("Loading RBP data")
        self.db_df = pd.read_csv(self.file_path, sep="\t")

    def process_result(self, off_target_df):
        if len(self.complete_result.index) != 0:
            # Update final result
            self.complete_result = self.complete_result[self.final_columns]
            self.pr_df[0] = {"name": self.db_name, "description": "Complete result",
                             "data": self.complete_result.to_json(orient="records")}

            # Update relevant fields in global off_target dataframe
            off_target_cosmic = self.complete_result.astype(
                {"off_target_id": int}).groupby("off_target_id").agg(
                {"gene_ensembl_id": lambda x: list(set(x))})
            off_target_cosmic.rename(columns={"gene_ensembl_id": "rbp_gene_ensembl_id"}, inplace=True)
            off_target_df = merge_off_target_information(off_target_df, off_target_cosmic, ["rbp_gene_ensembl_id"])
        return off_target_df


class COSMIC(Db):

    def __init__(self, file_path, final_columns):
        file_path = "{}/{}".format(database_base_path, file_path)
        super().__init__("COSMIC", file_path)
        self.enhancer_atlas = pd.DataFrame()
        self.remap_epd = pd.DataFrame()
        self.final_columns = final_columns

    def load_data(self):
        """
        Load the COSMIC file.
        Downloaded on the 7/11/2022
        :return: dataframe with COSMIC information
        """
        log.debug("Loading COSMIC data")
        self.db_df = pd.read_csv(self.file_path, sep="\t")

    def process_result(self, off_target_df):
        """
        Process the result of COSMIC.
        Add information to global off-target dataframe
        """
        if len(self.complete_result.index) != 0:
            # Update final result
            self.complete_result = self.complete_result[self.final_columns]
            self.pr_df[0] = {"name": self.db_name, "description": "Complete result",
                             "data": self.complete_result.to_json(orient="records")}

            # Update relevant fields in global off_target dataframe
            off_target_cosmic = self.complete_result.astype(
                {"off_target_id": int}).groupby("off_target_id").agg(
                {"Role in Cancer": lambda x: list(set(x))})
            off_target_cosmic.rename(columns={"Role in Cancer": "cancer_related"}, inplace=True)
            off_target_df = merge_off_target_information(off_target_df, off_target_cosmic, ["cancer_related"])

        # Update fields for Enhancer Atlas
        if len(self.enhancer_atlas.index) != 0:
            off_target_cosmic = self.enhancer_atlas.astype(
                {"off_target_id": int}).groupby("off_target_id").agg(
                {"Role in Cancer": lambda x: list(set(x))})
            off_target_cosmic.rename(columns={"Role in Cancer": "enhancer_atlas_cancer_related"}, inplace=True)
            off_target_df = merge_off_target_information(off_target_df, off_target_cosmic,
                                                         ["enhancer_atlas_cancer_related"])

        # Update fields for ReMap EPD
        if len(self.remap_epd.index) != 0:
            off_target_cosmic = self.remap_epd.astype(
                {"off_target_id": int}).groupby("off_target_id").agg(
                {"Role in Cancer": lambda x: list(set(x))})
            off_target_cosmic.rename(columns={"Role in Cancer": "remap_epd_cancer_related"}, inplace=True)
            off_target_df = merge_off_target_information(off_target_df, off_target_cosmic,
                                                         ["remap_epd_cancer_related"])
        return off_target_df


class ANonHumanDB(object):

    def __init__(self, name, location, columns):
        self._name = name
        self._location = location
        self._columns = columns
        self._db_bed = BedTool(self.location)

    @abstractmethod
    def select(self, query):
        pass

    @property
    def name(self):
        return self._name

    @property
    def location(self):
        return self._location

    @property
    def columns(self):
        return self._columns

    @property
    def db_bed(self):
        return self._db_bed

    def get_db_name(self):
        return self.name


def separate_attribute(line):
    """
    for each attribute separate by ';' the attribute itself separate by '=' in the format 'key=value' separate it to
    different columns
    :param line: attribute line
    :return: A dictionary of the attributes
    """

    function_name = "separate_attribute"
    dict_result = dict()
    if line:
        first_split = line.split(";")

        for s in first_split:
            second_split = s.split("=")
            if len(second_split) > 1:
                if second_split[0] in dict_result.keys():
                    raise Exception("ERROR in {}: The same key exist twice - {}".format(function_name, second_split[0]))
                else:
                    dict_result[second_split[0]] = second_split[1]
    return dict_result


def separate_attributes(df_to_separate):
    """
    Separate the attributes to different columns
    :param df_to_separate: A dictionary of the result to separate. key is index, value is dataframe
    :return: A dictionary the result intersection with attribute separate. key is index, value is dataframe
    """

    function_name = "separate_attributes"
    log.debug("Entering {}".format(function_name))

    # for result_key in off_target_result_dict:
    # Separate the attributes column
    if len(df_to_separate.index) != 0:
        df_to_separate.loc[:, "separate"] = df_to_separate.loc[:, "attributes"].apply(lambda s: separate_attribute(s))
        df_attributes = pd.DataFrame(df_to_separate["separate"].values.tolist(), index=df_to_separate.index)
        result = pd.concat([df_to_separate, df_attributes], axis=1).drop("separate", axis=1)
        return result
    else:
        return df_to_separate


def analyze_with_id_list(current_db, off_target_df, ensembl_id_list, db_name, column_to_search):
    log.info("Starting to analyze {}".format(current_db.db_name))
    if len(current_db.db_df) == 0:
        log.error("The data for {} was not loaded".format(current_db.db_name))
        return
    db_result = current_db.db_df[current_db.db_df[current_db.column_name_for_intersection].isin(ensembl_id_list)]
    if current_db.separate_attributes:
        db_result = separate_attributes(db_result)

    # Add the off_target_id to the complete result
    if len(db_result.index) != 0:
        db_result["off_target_id"] = ""
        if column_to_search in off_target_df.columns:
            # For each off_target add the id to the complete result if it is the ensemble_id
            for result in off_target_df.itertuples():
                if type(getattr(result, column_to_search)) is list:
                    result_id_list = list(getattr(result, column_to_search))
                    db_result.loc[:, "off_target_id"] = db_result.apply(
                        lambda x: "{},".format(result.off_target_id)
                        if x[current_db.column_name_for_intersection] in result_id_list
                        else x["off_target_id"], axis=1)
            db_result["off_target_id"] = db_result["off_target_id"].str.strip(",")
            db_result.reset_index(inplace=True, drop=True)
            setattr(current_db, db_name, db_result)
            if db_name == "complete_result":
                current_db.pr_df.append({"name": current_db.db_name, "description": "Complete result",
                                         "data": current_db.complete_result.to_json(orient="records")})
        else:
            log.info("No Gene ensembl ID therefore there are no result for {}".format(current_db.db_name))


def initialize_off_target_df(off_target_df):
    """
    Args:
        off_target_df: the global_off_target_df
    """
    off_target_df = off_target_df.rename(columns={"name": "off_target_id"})
    return off_target_df


def merge_off_target_information(left_df, right_df, columns_name, column_type="list"):
    """
    Args:
        right_df: right dataframe to merge
        left_df: usually the off_target_df
        columns_name: a list of column type to add to the left dataframe
        column_type: data type of the columns
    """
    df_result = left_df.merge(right_df, how="left", on="off_target_id")
    if column_type == "list":
        for column in columns_name:
            df_result[column] = df_result[column].apply(lambda x: [] if isinstance(x, float) else x)
            df_result[column] = df_result[column].apply(lambda x: [i for i in x if i not in [np.nan, ""]])

    return df_result


def save_global_off_target_results(off_target_df, flashfry_score, columns_order=None):
    """
    Args
    save the result as a json file that can be loaded later.
    :return:
    """
    off_target_df = off_target_df.astype({"chromosome": "string"})
    if columns_order:
        for value in columns_order:
            if value not in off_target_df.columns:
                off_target_df[value] = ""
        off_target_df = off_target_df[columns_order]
    save_result = {"off_targets": off_target_df.to_json(orient="records"),
                   "flashfry_score": flashfry_score.to_json(orient="records")}
    return save_result


def calculate_score(off_target_df):
    off_target_df["risk_score"] = ""
    off_target_complete_col = ["chromosome", "start", "end", "off_target_id", "score", "strand",
                               "attributes", "id", "sequence", "gene_ensembl_id", "gene_symbol",
                               "segment", "mir_gene", "remap_epd_gene_ensembl_id",
                               "enhancer_atlas_gene_ensembl_id", "disease_related",
                               "inheritance_model", "enhancer_atlas_disease_related",
                               "enhancer_atlas_inheritance_model", "remap_epd_disease_related",
                               "remap_epd_inheritance_model", "HumanTF_source",
                               "expression_information", "cancer_related",
                               "enhancer_atlas_cancer_related", "remap_epd_cancer_related", "risk_score"]
    off_target_current_col = off_target_df.columns.tolist()
    missing_col = [x for x in off_target_complete_col if x not in off_target_current_col]

    off_target_df[missing_col] = ""
    for row in off_target_df.itertuples():
        if "exon" in row.segment:
            if (row.disease_related != "") or (row.inheritance_model != "") or (row.cancer_related != ""):
                off_target_df.loc[row.Index, "risk_score"] = "High_coding"
            else:
                off_target_df.loc[row.Index, "risk_score"] = "Medium_coding"
        elif "transcript" in row.segment:
            off_target_df.loc[row.Index, "risk_score"] = "Low_coding"
        # todo : There is a bug here!
        elif row.remap_epd_gene_ensembl_id or row.enhancer_atlas_gene_ensembl_id:
            if row.enhancer_atlas_cancer_related or row.enhancer_atlas_inheritance_model or \
                    row.enhancer_atlas_disease_related or row.remap_epd_cancer_related or \
                    row.remap_epd_inheritance_model or row.remap_epd_disease_related:
                off_target_df.loc[row.Index, "risk_score"] = "Medium_regulatory"
            else:
                off_target_df.loc[row.Index, "risk_score"] = "Low_regulatory"
    return off_target_df


def add_db(db_list, db):
    """
    Add new DB to db_list
    Args:
        db_list: list of dbs
        db: db to add
    """

    db_list.append(db)


def update_database(db_list):
    """
    Update the existing databases
    Database that will be update: genome, GENCODE, HumanTF, mirGene.
    OMIM need to be brought manually
    """
    function_name = "update_database"
    log.debug("Entering {}".format(function_name))

    complete_genome = wget.download(url=COMPLETE_GENOME_URL, out="{}/genome".format(database_base_path))
    if os.path.exists(COMPLETE_GENOME_PATH):
        os.rename(COMPLETE_GENOME_PATH, "{}.old".format(COMPLETE_GENOME_PATH))
    extract_gz_file(complete_genome, COMPLETE_GENOME_PATH)
    log.info("{}: genome was downloaded to: {}".format(function_name, complete_genome))

    for db in db_list:
        db.update_db()


def save_db_result(db_list):
    """
    Save all the dataframe complete_result for all db
    Returns: jsonify object

    """
    app = Flask(__name__)
    with app.app_context():
        all_result = dict()
        for db in db_list:
            result_db_list = db.get_pr_df()
            db_name = db.get_db_name()
            if result_db_list:
                log.info("Saving {}".format(db_name))
                all_result.update({"{}_result_list".format(db.get_db_name()): result_db_list})
            else:
                log.info("No result for {}".format(db_name))
        return jsonify(all_result)


def update_database_base_path(new_base_path):
    """
    Update the base path for the databases
    Args:
        new_base_path: the new base path

    """
    global database_base_path
    if os.path.isdir(new_base_path):
        log.info("Change database base path to: {}".format(new_base_path))
        database_base_path = new_base_path
    else:
        log.error("Given path: {} , is not a valid directory".format(new_base_path))


def get_database_path():
    """

    Returns: the base path for the databases

    """
    global database_base_path
    return database_base_path
