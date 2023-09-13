import os
import pprint
import re
import subprocess
import sys
from datetime import timedelta
from time import perf_counter

import requests
from pybedtools import BedTool
import pandas as pd
import numpy as np
from BCBio.GFF import GFFExaminer
import mygene
import ssl
from tqdm import tqdm
import argparse

ssl._create_default_https_context = ssl._create_unverified_context

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HG19TOHG38_CONVERTOR = ""
LIFTOVER_PATH = ""
pd.set_option('display.max_columns', 7)

# Global function
def run_external_proc(args):
    """
    Run an external program with Popen.
    Args:
        args: The argument to run, a list of string

    Returns: return the output of the program

    """
    print("The following command will be run: {}\nThis process will take some time.".format(args))
    result = subprocess.run(args, capture_output=True, text=True, shell=True)
    result_message = "STDOUT: {}\nSTDERR: {}".format(result.stdout, result.stderr)
    print(result_message)
    return result_message


ENSG_SYMBOL_MAP = ""


def create_ensg_symbol_mapping():
    gencode_bed = BedTool(GENCODE_OUTPUT_PATH_FINAL)
    gencode_df = gencode_bed.to_dataframe(comment="#")
    gencode_df = separate_attributes(gencode_df)
    final_df = gencode_df[["gene_id", "gene_name"]]
    final_df.loc[:,"gene_id"] = final_df.loc[:,"gene_id"].apply(
        lambda x: re.sub("(\:.*$)", "", re.sub("^[^:]*:", "", re.sub("\.\d+", "", x))))
    final_df.drop_duplicates(inplace=True)
    final_df.to_csv(ENSG_SYMBOL_MAP, sep="\t", index=False)


# GENCODE
GENCODE_INPUT_PATH = ""
GENCODE_OUTPUT_PATH_FINAL = ""

GENCODE_MAPPING = ""
GENCODE_FILTER_ROWS = f"{BASE_DIR}/scripts/Remove.csv"


def separate_attribute(line):
    """
    for each attribute separate by ';' the attribute itself separate by '=' in the format 'key=value' separate it to
    different columns
    :param line: attribute line
    :return: A dictionary of the attributes
    """
    dict_result = dict()
    if line:
        first_split = line.split(";")

        for s in first_split:
            second_split = s.split("=")
            if len(second_split) > 1:
                if second_split[0] in dict_result.keys():
                    raise Exception("ERROR - The same key exist twice - {}".format(second_split[0]))
                else:
                    dict_result[second_split[0]] = second_split[1]
    return dict_result


def separate_attributes(df_to_separate):
    """
    Separate the attributes to different columns
    :param df_to_separate: A dictionary of the result to separate. key is index, value is dataframe
    :return: A dictionary the result intersection with attribute separate. key is index, value is dataframe
    """
    # for result_key in off_target_result_dict:
    # Separate the attributes column
    if not df_to_separate.empty:
        df_to_separate.loc[:, "separate"] = df_to_separate.loc[:, "attributes"].apply(lambda s: separate_attribute(s))
        df_attributes = pd.DataFrame(df_to_separate["separate"].values.tolist(), index=df_to_separate.index)
        result = pd.concat([df_to_separate, df_attributes], axis=1).drop("separate", axis=1)
        return result
    else:
        return df_to_separate


def load_assembly_report_to_csv():
    go_genome_mapping = pd.read_csv(GENCODE_MAPPING, comment="#", sep="\t",
                                    names=["Sequence-Name", "Sequence-Role", "Assigned-Molecule",
                                           "Assigned-Molecule-Location/Type", "GenBank-Accn", "Relationship",
                                           "RefSeq-Accn", "Assembly-Unit", "Sequence-Length", "UCSC-style-name"])
    go_genome_mapping.rename(columns={}, inplace=True)
    return go_genome_mapping


def create_mapping_from_gencode(chrom_names_orig):
    """
    Create mapping name for GoGenome according to GENCODE naming
    Args:
        chrom_names_orig: unique list of old naming from gencode

    Returns: Dataframe with mapping from old name from GENCODE to new one

    """
    go_genome_mapping = load_assembly_report_to_csv()

    # Create mapping dictionary
    name_mapping_dict = {}
    for chrom_name in chrom_names_orig.itertuples():
        name_mapping_dict[chrom_name[1]] = {"candidates_lst": [], "sequence_name": []}

    # Add mapping names for each of the original chrom names
    for chrom_name_tup in chrom_names_orig.itertuples():
        chrom_name = chrom_name_tup[1]
        for mapping_row in go_genome_mapping.itertuples():
            sequence_name = mapping_row[1]
            gen_bank_curr = mapping_row[5]
            ucsc_curr = mapping_row[10]
            if chrom_name == gen_bank_curr:
                name_mapping_dict[chrom_name]["candidates_lst"] += [gen_bank_curr]
                name_mapping_dict[chrom_name]["sequence_name"] += [sequence_name]
            elif chrom_name == ucsc_curr:
                name_mapping_dict[chrom_name]["candidates_lst"] += [ucsc_curr]
                name_mapping_dict[chrom_name]["sequence_name"] += [sequence_name]

    chrom_names = chrom_names_orig.copy(deep=True)
    chrom_names["correct_names"] = chrom_names[0].apply(
        lambda x: name_mapping_dict[x]["sequence_name"][0] if name_mapping_dict[x]["sequence_name"] else x)
    chrom_names.rename(columns={0: "prev_names"}, inplace=True)
    return chrom_names


def replace_chrom_names():
    """
    Replace all the chromosome name according to GoGenome mapping in GENCODE
    Returns: Dataframe for GENCODE with new chromosome names
    """

    gencode_bed = BedTool(GENCODE_INPUT_PATH)
    gencode_df = gencode_bed.to_dataframe(comment="#")

    # Retrieve the mapping for old naming to new
    chrom_names_orig = pd.DataFrame(pd.unique(gencode_df["seqname"]))
    chrom_names = create_mapping_from_gencode(chrom_names_orig)

    # Replace chromosome names
    new_gencode_df = gencode_df.merge(chrom_names, left_on="seqname", right_on="prev_names")
    new_gencode_df.drop(columns=["seqname", "prev_names"], inplace=True)
    new_gencode_df.rename(columns={"correct_names": "seqname"}, inplace=True)
    return new_gencode_df


def filter_rows(new_gencode_df):
    """
    Remove all records that are not relevant for analysis for better performance.
    Records are determine base on the Remove.csv file
    Args:
        new_gencode_df: Tha dataframe from which to remove the records
    Returns: given dataframe without the records
    """
    print("separate attributes")
    # Dataframe might be to big to handle in amazon worksapce - separate it.
    a6 = len(new_gencode_df)
    a5 = int(a6 * 5 / 6)
    a4 = int(a6 * 4 / 6)
    a3 = int(a6 * 3 / 6)
    a2 = int(a6 * 2 / 6)
    a1 = int(a6 / 6)
    a0 = 0

    # Extract attributes
    df1 = separate_attributes(new_gencode_df[a0:a1])
    df2 = separate_attributes(new_gencode_df[a1:a2])
    df3 = separate_attributes(new_gencode_df[a2:a3])
    df4 = separate_attributes(new_gencode_df[a3:a4])
    df5 = separate_attributes(new_gencode_df[a4:a5])
    df6 = separate_attributes(new_gencode_df[a5:a6])
    data = [df1, df2, df3, df4, df5, df6]
    new_gencode_df = pd.concat(data)
    print("Filter rows")
    rows_to_filter_df = pd.read_csv(GENCODE_FILTER_ROWS)
    print("Filter rows 1")
    # Remove all records with the type         ~(new_gencode_df.gene_type.str.contains("pseudogene")) | (
    new_gencode_df = new_gencode_df[
        ~(new_gencode_df.gene_type.str.contains("pseudogene")) | (
            new_gencode_df.transcript_type.str.contains("pseudogene"))]
    new_gencode_df.reset_index(drop=True, inplace=True)
    print("Filter rows 2")
    # Remove all records that match in the file Remove.csv
    indexes_to_remove = []
    for row_gencode in new_gencode_df.itertuples():
        for row_remove in rows_to_filter_df.itertuples():
            if row_gencode.gene_type == row_remove.gene_type and \
                    row_gencode.transcript_type == row_remove.transcript_type and \
                    row_gencode.feature == row_remove.feature:
                indexes_to_remove.append(row_gencode.Index)
                break
    new_gencode_df = new_gencode_df[~(new_gencode_df.index.isin(indexes_to_remove))]
    new_gencode_df.loc[:, "gene_id"] = new_gencode_df.loc[:, "gene_id"].apply(
        lambda x: re.sub("(\:.*$)", "", re.sub("^[^:]*:", "", re.sub("\.\d+", "", x))))

    new_gencode_df.fillna("", inplace=True)

    new_gencode_df["new_attributes"] = ""
    new_gencode_df["new_attributes"] = new_gencode_df.apply(
        lambda s: f"gene_id={s.gene_id}" if s.gene_id != "" else f"{s.new_attributes}", axis=1)

    new_gencode_df["new_attributes"] = new_gencode_df.apply(
        lambda s: f"{s.new_attributes};gene_name={s.gene_name}" if s.gene_name != "" else f"{s.new_attributes}", axis=1)

    new_gencode_df["new_attributes"] = new_gencode_df.apply(
        lambda s: f"{s.new_attributes};gene_type={s.gene_type}" if s.gene_type != "" else f"{s.new_attributes}", axis=1)

    new_gencode_df["new_attributes"] = new_gencode_df.apply(
        lambda s: f"{s.new_attributes};transcript_id={s.transcript_id}"
        if s.transcript_id != "" else f"{s.new_attributes}", axis=1)

    new_gencode_df["new_attributes"] = new_gencode_df.apply(
        lambda  s: f"{s.new_attributes};transcript_name={s.transcript_name}"
        if s.transcript_name != "" else f"{s.new_attributes}", axis=1)

    new_gencode_df["new_attributes"] = new_gencode_df.apply(
        lambda s: f"{s.new_attributes};transcript_type={s.transcript_type}"
        if s.transcript_type != "" else f"{s.new_attributes}",axis=1)

    new_gencode_df["new_attributes"] = new_gencode_df.apply(
        lambda s: f"{s.new_attributes};protein_id={s.protein_id}" if s.protein_id != "" else f"{s.new_attributes}",
        axis=1)

    new_gencode_df["new_attributes"] = new_gencode_df.apply(
        lambda s: f"{s.new_attributes};exon_number={s.exon_number}" if s.exon_number != "" else f"{s.new_attributes}",
        axis=1)

    new_gencode_df["new_attributes"] = new_gencode_df.apply(
        lambda s: f"{s.new_attributes};exon_id={s.exon_id}" if s.exon_id != "" else f"{s.new_attributes}", axis=1)

    new_gencode_df.drop(columns=['ID', 'gene_id', 'gene_type', 'gene_name',
       'level', 'tag', 'Parent', 'transcript_id', 'transcript_type',
       'transcript_name', 'transcript_support_level', 'havana_transcript',
       'exon_number', 'exon_id', 'hgnc_id', 'havana_gene', 'ont', 'protein_id',
       'ccdsid', 'artif_dupl'], inplace=True)

    new_gencode_df = new_gencode_df[
        ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "new_attributes"]]
    return new_gencode_df


def preprocess_gencode():
    """
    Version 42, Downloaded on 6/11/2022
    """
    new_gencode_df = replace_chrom_names()
    new_gencode_df = filter_rows(new_gencode_df)
    gencode = BedTool.from_dataframe(new_gencode_df).sort()
    gencode.saveas(GENCODE_OUTPUT_PATH_FINAL)


# MirGeneDB
MIRGENEDB_HG37_INPUT_PATH = ""
MIRGENEDB_OUTPUT_PATH = ""


def preprocess_mirgenedb():
    """
    Version 2.1, Downloaded on 6/11/2022
    Load the file to Bedtools and use the built in sort. Save in the output path
    """
    mirgene_bed = BedTool(MIRGENEDB_HG37_INPUT_PATH).sort()
    mirgene_df = mirgene_bed.to_dataframe(comment="#")
    mirgene_df["chrom"] = mirgene_df["chrom"].apply(lambda x: x.strip("chr"))
    mirgene = BedTool.from_dataframe(mirgene_df).sort()
    mirgene.saveas(MIRGENEDB_OUTPUT_PATH)


# Remap and EPD
REMAP_INPUT_PATH = ""
EPD_H_INPUT_PATH = ""
EPD_H_NC_INPUT_PATH = ""
EPD_H_NEW_PATH = ""
EPD_H_NC_NEW_PATH = ""
MERGE_EPD_PATH = ""
REMAP_EPD_OUTPUT_PATH = ""


def process_epd(epd_path, coding):
    if coding:
        epd_search_link = "https://epd.epfl.ch/cgi-bin/get_doc?db=hgEpdNew&format=genome&entry={}"
    else:
        epd_search_link = "https://epd.epfl.ch/cgi-bin/get_doc?db=hsNCEpdNew&format=genome&entry={}"
    epd_df = pd.read_csv(epd_path, sep="\t", header=None)
    epd_df[8] = ""
    mg = mygene.MyGeneInfo()
    for row in tqdm(epd_df.itertuples(), total=epd_df.shape[0]):
        df_pandas = pd.read_html(epd_search_link.format(row._4), flavor="bs4", thousands=".")[1]
        gene_symbol = df_pandas[df_pandas[0] == "Gene Symbol:"][1].to_string(index=False)
        gene_ensembl_id = df_pandas[df_pandas[0] == "Ensembl:"][1].to_string(index=False)
        if gene_ensembl_id == "Series([], )":
            out = mg.query(gene_symbol, scopes="symbol", fields="ensembl.gene", species="human")
            if out.get("total") != 0:
                for result_item in out.get("hits"):
                    if "ensembl" in result_item:
                        gene_ensembl_id = result_item.get("ensembl")
                        if type(gene_ensembl_id) is list:
                            gene_ensembl_id = gene_ensembl_id[0].get("gene")
                        else:
                            gene_ensembl_id = gene_ensembl_id.get("gene")
                        break
        epd_df.loc[row.Index, 8] = gene_ensembl_id
    epd_df[8] = epd_df[[3, 8]].apply(lambda s: "{};epd_gene_symbol={};epd_coding={}".format(s[8], s[3],coding), axis=1)
    return epd_df


def preprocess_remap_epd(file_location_epd, output_file_path):
    remap = BedTool(REMAP_INPUT_PATH)

    epd_df = pd.read_csv(file_location_epd, sep="\t", header=None)
    epd = BedTool.from_dataframe(epd_df)
    # Intersect ReMap and EPD databases
    intersection_result = remap.intersect(epd, wa=True, wb=True, sorted=True)
    print("Number of records: {}".format(intersection_result.count()))
    intersection_result_df = intersection_result.to_dataframe(header=None)
    intersection_result_df.drop([4, 5, 6, 7, 8, 9, 10, 11, 13], axis=1, inplace=True)
    intersection_result_df = process_result_intersection(intersection_result_df)
    intersection_result_df[0] = intersection_result_df[0].apply(lambda x: x.strip("chr"))
    final_bed = BedTool.from_dataframe(intersection_result_df).sort()
    print("Number of records: {}".format(final_bed.count()))
    final_bed.saveas(output_file_path)


def process_result_intersection(result_to_process):
    final_df = result_to_process.fillna(value={3: ".", 12: "."})
    final_df.loc[:, 3] = final_df.loc[:, 3].apply(lambda s: ("remap={}".format(s)) if s != "." else s)
    final_df.loc[:, 12] = final_df.loc[:, 12].apply(lambda s: ("gene_ensembl_id={}".format(s)) if s != "." else s)
    final_df.loc[:, "attributes"] = final_df.apply(lambda s: ";".join([str(s[3]), str(s[12])]), axis=1)
    final_df.drop([3, 12], axis=1, inplace=True)

    final_df.loc[:, "score"] = 0
    final_df.loc[:, "frame"] = ""

    final_df = final_df[[0, 1, 2, "attributes", "score", 14, "frame"]]
    return final_df


def preprocess_remapepd():
    epd_h_df = process_epd(EPD_H_INPUT_PATH, 1)
    epd_h_df_bed = BedTool.from_dataframe(epd_h_df[[0, 1, 2, 8, 4, 5]]).sort()
    epd_h_df_bed.saveas(EPD_H_NEW_PATH)
    epd_h_nc_df = process_epd(EPD_H_NC_INPUT_PATH, 0)
    epd_h_nc_bed = BedTool.from_dataframe(epd_h_nc_df[[0, 1, 2, 8, 4, 5]]).sort()
    epd_h_nc_bed.saveas(EPD_H_NC_NEW_PATH)

    concatenated = pd.concat([epd_h_df, epd_h_nc_df])
    final_bed = BedTool.from_dataframe(concatenated[[0, 1, 2, 8, 4, 5]]).sort()
    final_bed.saveas(MERGE_EPD_PATH)
    preprocess_remap_epd(MERGE_EPD_PATH, REMAP_EPD_OUTPUT_PATH)


# EnhancerAtlas
TISSUE_OR_CELL_STR_NAMES = "A375 	A549 	AML_blast 	Astrocyte 	BJ 	Bronchia_epithelial Caco-2 	Calu-3 	CD14+ 	" \
                           "CD19+ 	CD20+ 	CD34+ CD36+ 	CD4+ 	CD8+ 	Cerebellum 	CUTLL1 	DOHH2 ECC-1 	" \
                           "ESC_neuron 	Esophagus 	Fetal_heart 	Fetal_kidney 	Fetal_muscle_leg Fetal_placenta 	" \
                           "Fetal_small_intestine 	Fetal_spinal_cord Fetal_stomach 	Fetal_thymus 	FT246 FT33 	" \
                           "GM10847 	GM12878 	GM12891 	GM12892 	GM18505 GM18526 	GM18951 	" \
                           "GM19099 	GM19193 	GM19238 	GM19239 GM19240 	H1 	H9 	HCC1954 	HCT116 " \
                           "	HEK293T HEK293 	Hela-S3 	Hela 	HepG2 	HFF 	HL-60 hMADS-3 	HMEC 	hNCC 	" \
                           "HSMM 	HT1080 	HT29 HUVEC 	IMR90 	Jurkat 	K562 	Kasumi-1 	KB Keratinocyte 	" \
                           "Left_ventricle 	LHCN-M2 	Liver 	LNCaP-abl 	LNCaP Lung 	MCF-7 	MCF10A 	ME-1 	" \
                           "Melanocyte 	melanoma Mesendoderm 	MS1 	Myotube 	Namalwa 	NB4 	NHDF NHEK 	" \
                           "NHLF 	NKC 	OCI-Ly7 	Osteoblast 	Ovary PANC-1 	Pancreas 	Pancreatic_islet 	" \
                           "PBMC 	PC3 	PrEC SGBS_adipocyte 	SK-N-SH 	SK-N-SH_RA 	Skeletal_muscle 	" \
                           "Small_intestine 	Sperm Spleen 	T47D 	T98G 	th1 	Thymus 	U2OS VCaP 	ZR75-30"
BASE_URL = "http://www.enhanceratlas.org/data/AllEPs/hs/{}_EP.txt"
ENHANCER_ATLAS_HG19 = ""
ENHANCER_ATLAS_HG38 = ""
ENHANCER_ATLAS_OUTPUT = ""


def create_tissue_or_cell_list(string_to_convert):
    """
    Convert from tab delimiter string to list of tissue or cell name
    Args:
        string_to_convert: const string .

    Returns: list of tissue or cell names.

    """
    string_to_convert = string_to_convert.replace("\t", "")
    names_str_list = string_to_convert.split(" ")
    return names_str_list


def create_enhancer_atlas_bed_file_hg19():
    """
        Preprocess preprocess_enhanceratlas. Get from a URL all the files relevant to Homo sapience.
        Merge all files to one file to be loaded ad BED.
        Downloaded pn 7/11/2022
        The original format:
        For human, the data format in the file listed as "A:B-C_D$E$F$G$H	I":
        A - Chromosome of enhancer.
        B - The starting position of enhancer.
        C - The ending position of enhancer.
        D - Gene ensembl ID.
        E - Gene ensembl Name.
        F - chromosome of gene.
        G - the position of the gene transcription start site.
        H - The strand of DNA the gene located.
        I - The predition score of the enhancer-gene interaction.

        Returns: save the file to enhanceratlas_db.csv, tab delimiter, with no header and index.

        """

    tissue_or_cell_name_list = create_tissue_or_cell_list(TISSUE_OR_CELL_STR_NAMES)
    final_db = pd.DataFrame()

    for tissue_or_cell in tissue_or_cell_name_list:
        response = requests.get(BASE_URL.format(tissue_or_cell))
        data = response.text

        data = data.replace(":", ";")
        data = data.replace("\t", "$")
        data_list = data.split("\r\n")
        current_db = pd.DataFrame(data_list)
        current_db[0] = current_db[0].apply(lambda s: s.replace("-", ";", 1))
        current_db[0] = current_db[0].apply(lambda s: s.replace("_", ";", 1))
        current_db = current_db[0].str.split(";", expand=True)

        current_db.rename(inplace=True, columns={0: "chr_enhancer", 1: "enh_start", 2: "enh_stop", 3: "name"})
        current_db["tissue_or_cell_name"] = tissue_or_cell
        current_db["name"] = current_db[["name", "tissue_or_cell_name"]].apply(
            lambda row: "$".join(row.values.astype(str)), axis=1)
        current_db.dropna(inplace=True, thresh=4)
        final_db = pd.concat([final_db, current_db])

    # final_db.to_csv("enhancerAtlas_hg19.csv", sep="\t")
    db_bed = BedTool.from_dataframe(final_db[["chr_enhancer", "enh_start", "enh_stop", "name"]], na_rep=".")
    db_bed = db_bed.sort()
    db_bed.saveas(ENHANCER_ATLAS_HG19)


def preprocess_enhanceratlas():
    create_enhancer_atlas_bed_file_hg19()
    # Convert with LiftOver the cordinate from hg19 to hg38
    run_external_proc([f"{LIFTOVER_PATH} "
                       f"{ENHANCER_ATLAS_HG19} "
                       f"{HG19TOHG38_CONVERTOR} "
                       f"{ENHANCER_ATLAS_HG38} "
                       f"enhancer_atlas_unMapped"])

    enhancer_atlas_hg38_bed = BedTool(ENHANCER_ATLAS_HG38).sort()
    enhancer_atlas_hg38_df = enhancer_atlas_hg38_bed.to_dataframe(comment="#")

    # df1 = enhancer_atlas_hg38_df["name"].str.split("$", expand=True)
    # new_db = pd.concat([enhancer_atlas_hg38_df[:], df1[:]], axis=1)
    #
    # new_db["chrom"] = new_db["chrom"].apply(lambda x: x.strip("chr"))
    # final_db_group = new_db.groupby(
    #     ["chrom", "start", "end", 0, 1, 2, 3, 4, 5]) \
    #     .agg({"name": lambda x: x.tolist(), 6: lambda x: x.tolist()}) \
    #     .reset_index()

    # Dataframe might is to big to handle in amazon worksapce - separate it.
    a4 = len(enhancer_atlas_hg38_df)
    a3 = int(a4 * 3 / 4)
    a2 = int(a4 * 2 / 4)
    a1 = int(a4 / 4)
    a0 = 0

    # Extract attributes
    df1 = enhancer_atlas_hg38_df[a0:a1]["name"].str.split("$", expand=True)

    df2 = enhancer_atlas_hg38_df[a1:a2]["name"].str.split("$", expand=True)
    df3 = enhancer_atlas_hg38_df[a2:a3]["name"].str.split("$", expand=True)
    df4 = enhancer_atlas_hg38_df[a3:a4]["name"].str.split("$", expand=True)

    data = [df1, df2, df3, df4]
    new_enhancer_atlas_hg38_df = pd.concat(data)
    new_db = pd.concat([enhancer_atlas_hg38_df[:], new_enhancer_atlas_hg38_df[:]], axis=1)

    final_db_group = new_db.groupby(
        ["chrom", "start", "end", "name", 0, 1, 2, 3, 4, 5]) \
        .agg({6: lambda x: x.tolist()}) \
        .reset_index()
    final_db_group["chrom"] = final_db_group["chrom"].apply(lambda x: str(x).strip("chr"))

    db_bed = BedTool.from_dataframe(final_db_group[["chrom", "start", "end", "name"]], na_rep=".")
    db_bed = db_bed.sort()
    db_bed.saveas(ENHANCER_ATLAS_OUTPUT)
    return final_db_group


# OMIM
GENEMAP2_INPUT_PATH = ""
MIM2GENE_INPUT_PATH = ""
OMIM_CSV_OUTPUT_PATH = ""
OMIM_BED_OUTPUT_PATH = ""


def preprocess_omim():
    """
    Load genemap2.txt and mim2gene.txt. Downloaded on 24/12/2020
    Returns:

    """
    print("Loading OMIM data")
    genemap2_df = pd.read_csv(GENEMAP2_INPUT_PATH, sep="\t", skiprows=3,
                              usecols=(
                                  "# Chromosome", "Genomic Position Start", "Genomic Position End", "MIM Number",
                                  "Phenotypes"))
    # Remove all rows that do not contains any data
    genemap2_df = genemap2_df[~genemap2_df["# Chromosome"].astype(str).str.startswith("#")]

    genemap2_df_bed = genemap2_df.astype(
        {"Genomic Position Start": int, "Genomic Position End": int, "Phenotypes": "string", "MIM Number": int})

    genemap2_df_bed.rename(columns={"# Chromosome": "chromosome", "Genomic Position Start": "start",
                                    "Genomic Position End": "end", "MIM Number": "omim_id",
                                    "Phenotypes": "disease_related"}, inplace=True)
    genemap2_df_bed.dropna(subset=["disease_related"], inplace=True)
    # Merger the Ensembl ID and add it to attributes as the same format in gff
    mim2gene_pd = pd.read_csv(MIM2GENE_INPUT_PATH, sep="\t", skiprows=4)
    mim2gene_pd.rename(columns={"# MIM Number": "omim_id",
                                "Ensembl Gene ID (Ensembl)": "gene_ensembl_id"}, inplace=True)

    genemap2_df_bed = genemap2_df_bed.merge(mim2gene_pd, on="omim_id")
    genemap2_df_bed.rename(columns={"Approved Gene Symbol (HGNC)": "gene_symbol"}, inplace=True)

    genemap2_df_bed["inheritance_model"] = ""
    for row in genemap2_df_bed.itertuples():
        current_phenotypes = row.disease_related
        current_inheritance_model = ""
        if current_phenotypes:
            if ("X-linked dominant" in current_phenotypes) | ("X-linked recessive" in current_phenotypes):
                if "X-linked dominant" in current_phenotypes:
                    current_inheritance_model = current_inheritance_model + "X-linked dominant "
                    current_phenotypes = current_phenotypes.replace("X-linked dominant", "")
                if "X-linked recessive" in current_phenotypes:
                    current_inheritance_model = current_inheritance_model + "X-linked recessive "
                    current_phenotypes = current_phenotypes.replace("X-linked recessive", "")
            elif "X-linked" in current_phenotypes:
                current_inheritance_model = current_inheritance_model + "X-linked "
                current_phenotypes = current_phenotypes.replace("X-linked", "")
            if "Autosomal dominant" in current_phenotypes:
                current_inheritance_model = current_inheritance_model + "Autosomal dominant "
                current_phenotypes = current_phenotypes.replace("Autosomal dominant", "")
            if "Autosomal recessive" in current_phenotypes:
                current_inheritance_model = current_inheritance_model + "Autosomal recessive"
                current_phenotypes = current_phenotypes.replace("Autosomal recessive", "")
            genemap2_df_bed.loc[row.Index, "inheritance_model"] = current_inheritance_model.strip()
            genemap2_df_bed.loc[row.Index, "disease_related"] = current_phenotypes.rstrip(" ,{}")

    genemap2_df_bed.to_csv(OMIM_CSV_OUTPUT_PATH, sep="\t", index=False,
                           columns=["gene_ensembl_id", "gene_symbol", "omim_id", "disease_related",
                                    "inheritance_model"])

    genemap2_df_bed.rename(columns={"omim_id": "name", "disease_related": "Phenotypes"}, inplace=True)

    genemap2_df_bed.loc[:, "score"] = 0
    genemap2_df_bed.loc[:, "strand"] = "+"
    genemap2_df_bed.loc[:, "frame"] = ""

    genemap2_df_bed.loc[:, "Phenotypes"] = genemap2_df_bed.loc[:, "Phenotypes"].apply(
        lambda s: ("Phenotypes={}".format(s.replace(";", ":"))) if s != "nan" else None)
    genemap2_df_bed = genemap2_df_bed.fillna(value={"attributes": ".", "gene_ensembl_id": ".", "gene_symbol": "."})
    genemap2_df_bed.loc[:, "gene_symbol"] = genemap2_df_bed.loc[:, "gene_symbol"].apply(
        lambda s: ("gene_symbol={}".format(s)) if s != "." else s)
    genemap2_df_bed.loc[:, "gene_ensembl_id"] = genemap2_df_bed.loc[:, "gene_ensembl_id"].apply(
        lambda s: ("gene_ensembl_id={}".format(s)) if s != "." else s)
    genemap2_df_bed.loc[:, "attributes"] = genemap2_df_bed.apply(
        lambda s: ";".join([str(s["Phenotypes"]), str(s["gene_ensembl_id"]), str(s["gene_symbol"]),
                            str(s["inheritance_model"])]), axis=1)
    genemap2_df_bed.loc[:, "attributes"] = genemap2_df_bed.loc[:, "attributes"].apply(
        lambda s: s.rstrip(".;").lstrip(";."))
    genemap2_df_bed = genemap2_df_bed.rename(columns={"test": "attributes"})
    genemap2_df_bed = genemap2_df_bed[["chromosome", "start", "end", "name", "score", "strand", "frame", "attributes"]]
    genemap2_df_bed.loc[:, "attributes"] = genemap2_df_bed.loc[:, "attributes"].replace("", ".")

    BedTool.from_dataframe(genemap2_df_bed, na_rep=".").saveas(OMIM_BED_OUTPUT_PATH)


# Human TF
TF_INPUT_PATH = ""
TF_COFACTOR_INPUT_PATH = ""
HUMANTF_OUTPUT_PATH = ""


def preprocess_humantf():
    """
    Downloade om 7/11/2022
    Returns:

    """
    print("Loading Human TF data")
    tf_list_df = pd.read_csv(TF_INPUT_PATH, sep="\t")
    tf_list_df["source"] = "TF"
    tf_cofactors_list_df = pd.read_csv(TF_COFACTOR_INPUT_PATH, sep="\t")
    tf_cofactors_list_df["source"] = "TF cofactors"
    db_df = pd.concat([tf_list_df, tf_cofactors_list_df])
    db_df.rename(columns={"Symbol": "gene_symbol", "Ensembl": "gene_ensembl_id", }, inplace=True)
    db_df.drop(["Species", "Protein", "Entrez ID"], axis=1, inplace=True)
    db_df.to_csv(HUMANTF_OUTPUT_PATH, sep="\t", index=False)


# Protein Atlas
PA_INPUT_PATH = ""
PA_OUTPUT_PATH = ""
def get_info_from_tsv(tsv_file, k):
    """
    get_info_from_tsv function:
    receives tsv file as input, AND an integer K which indicates how to filter the tissue-cell types.
    for example if tissue A and cell type B appears together 40 times in the input file, and K = 50
    so A,B will be dropped from the result file.
    """

    print("-----Begin filtering & creating the result matrix-----")
    df = pd.read_csv(tsv_file, sep="\t", header=0)
    df.rename(columns={"Gene": "gene_ensembl_id", "Gene name": "gene_symbol"}, inplace=True)
    df = df.astype({"gene_symbol": "string", "Tissue": "string", "Cell type": "string"})
    # get amount of unique genes in file.
    genes = df["gene_ensembl_id"].nunique()
    # get amount of unique tissue - cell type combinations in file.
    combinations = df.groupby(["Tissue", "Cell type"]).size().reset_index().rename(columns={0: "count"})
    s = combinations.loc[combinations["count"] >= k]
    lst = ["gene_ensembl_id", "gene_symbol"]
    genes_dictionary = {}
    # build a list contains all "chosen" tissue-cell type combinations.
    for index, row in s.iterrows():
        lst.append(row["Tissue"] + " " + row["Cell type"])
    # build gene dictionary, key: gene, value: list of all combinations appeared with the gene.
    for index, row in df.iterrows():
        try:
            if (row["gene_ensembl_id"], row["gene_symbol"]) not in genes_dictionary:
                genes_dictionary[(row["gene_ensembl_id"], row["gene_symbol"])] = [(row[2] + " " + row[3], row[4])]
            else:
                genes_dictionary[(row["gene_ensembl_id"], row["gene_symbol"])].append((row[2] + " " + row[3], row[4]))
        except Exception as e:
            print("Error: {}".format(e))
            continue
    new_df = pd.DataFrame(columns=lst)
    # run over each gene in the gene dictionary
    for gene in genes_dictionary.keys():
        curr_row = [gene[0], gene[1]]
        # check for common tissue-cell type combo`s
        for key in lst[2:]:
            flag = False
            for values_in_gene in genes_dictionary[gene]:
                if key == str(values_in_gene[0]):
                    flag = True
                    curr_row.append(values_in_gene[1])
                    break
            if not flag:
                curr_row.append("NA")
        new_df.loc[len(new_df)] = curr_row
    new_df.to_csv(PA_OUTPUT_PATH, index=False)
    print("-----DONE-----")


FILTER_INTEGER = 50


def preprocess_proteinatlas():
    """
    Downloaded on 7/11/2022
    Returns:

    """
    get_info_from_tsv(PA_INPUT_PATH, FILTER_INTEGER)


# RBP
RBP_INPUT_PATH = ""
RBP_OUTPUT_PATH = ""


def preprocess_rbp():
    """"
    Downloaded on 7/11/2022
    """
    db_df = pd.read_excel(RBP_INPUT_PATH, index_col=0, header=None, skiprows=2)
    db_df.reset_index(inplace=True)
    col_name = ["gene_ensembl_id", "Essential Genes", "Splicing regulation", "Spliceosome", "RNA modification",
                "3' end processing", "rRNA processing", "Ribosome & basic translation", "RNA stability & decay",
                "microRNA processing", "RNA localization", "RNA export", "Translation regulation",
                "tRNA regulation",
                "mitochondrial RNA regulation", "Viral RNA regulation", "snoRNA / snRNA / telomerase",
                "P-body / stress granules", "Exon Junction Complex"]
    column_map = dict(zip(range(1, len(col_name) + 1), col_name))
    db_df = db_df.rename(columns=column_map)
    db_df = db_df.rename(columns={0: "gene_symbol"})
    col_list = db_df.columns.to_list()
    col_list.insert(0, col_list.pop(1))
    db_df = db_df[col_list]
    db_df.drop(db_df.iloc[:, 20:], inplace=True, axis=1)
    db_df.to_csv(RBP_OUTPUT_PATH, sep="\t", index=False)


# COSMIC
COSMIC_INPUT_PATH = ""
COSMIC_OUTPUT_PATH = ""


def separate_id(line):
    if line:
        id_list = line.split(",")
        ensembl_id = [item for item in id_list if "ENSG" in item]
        if len(ensembl_id) > 0:
            if len(ensembl_id) > 1:
                print(
                    "Size of the ensembl ID list is {} when it should be 1. returning only the first one".format(
                        len(ensembl_id)))
            return ensembl_id[0].split(".")[0]


def preprocess_cosmic():
    """
    Downloaded on 7/11/2022
    """
    db_df = pd.read_csv(COSMIC_INPUT_PATH)
    db_df.drop(["Entrez GeneId", "Genome Location", "Tier", "Hallmark", "Chr Band", "Cancer Syndrome", "Tissue Type",
                "Mutation Types", "Translocation Partner", "Other Germline Mut", "Other Syndrome"],
               axis=1, inplace=True)
    db_df.rename(columns={"Gene Symbol": "gene_symbol"}, inplace=True)
    db_df["gene_ensembl_id"] = db_df["Synonyms"].apply(lambda s: separate_id(str(s)))
    db_df.drop(["Synonyms"], axis=1, inplace=True)
    print("Cosmic unique values in column Role in Cancer before {}".format(db_df["Role in Cancer"].unique().tolist()))
    db_df.dropna(subset=["Role in Cancer"], inplace=True)
    print("Cosmic unique values in column Role in Cancer after {}".format(db_df["Role in Cancer"].unique().tolist()))
    db_df.to_csv(COSMIC_OUTPUT_PATH, sep="\t", index=False)


PFAM_INPUT_PATH = ""
PFAM_GENCODE_INTERSECT_PATH = ""


# Pfam Protein Domain
def preprocess_pfam_protein_domains():
    """
    Downloaded on 7/11/2022
    """
    pfam_df = pd.read_csv(PFAM_INPUT_PATH, delimiter="\t")
    new_pfam = pfam_df.copy(deep=True)[0:0][["chrom", "chromStart", "chromEnd", "name", "score", "strand"]]
    print("separate all the locations into separate row")
    # separate all the locations into separate row
    for row in pfam_df.itertuples():
        block_starts = row.chromStarts.split(",")[:-1]
        block_sizes = row.blockSizes.split(",")[:-1]
        for b_index, block_size in enumerate(block_sizes):
            new_row = {"chrom": row.chrom, "chromStart": int(row.chromStart) + int(block_starts[b_index]),
                       "chromEnd": int(row.chromStart) + int(block_starts[b_index]) + int(block_size), "name": row.name,
                       "score": row.score, "strand": row.strand}
            new_pfam.loc[len(new_pfam.index)] = new_row

    # new_pfam.to_csv("PFAM_WITH_bLOCKS.csv", index=False)
    #
    # new_pfam = pd.read_csv("PFAM_WITH_bLOCKS.csv")

    print("Map to NCBI name")
    chrom_names_orig = pd.DataFrame(pd.unique(new_pfam["chrom"]))
    chrom_names = create_mapping_from_gencode(chrom_names_orig)

    # Replace chromosome names
    new_pfam = new_pfam.merge(chrom_names, left_on="chrom", right_on="prev_names")
    new_pfam.drop(columns=["chrom", "prev_names"], inplace=True)
    new_pfam.rename(columns={"correct_names": "chromosome"}, inplace=True)

    new_pfam = new_pfam[["chromosome", "chromStart", "chromEnd", "name", "score", "strand"]]
    # Remove records of "_ML" since it is not in GENCODE
    chromosome_to_remove = list(set([item for item in new_pfam["chromosome"].unique().tolist() if "_ML" in item]))
    df_index_to_remove = new_pfam[new_pfam.chromosome.isin(chromosome_to_remove)].index.tolist()
    new_pfam.drop(index=df_index_to_remove, inplace=True)
    pfam_merged_bed = BedTool.from_dataframe(new_pfam).sort()

    print("Retrive Ensembl ID from GENCODE")
    # Load GENCODE
    gencode = BedTool(GENCODE_OUTPUT_PATH_FINAL)
    gencode_df = gencode.to_dataframe()
    gencode_df = gencode_df[gencode_df.feature == "CDS"]
    gencode_separated = separate_attributes(gencode_df)
    gencode_separated = gencode_separated[
        ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "gene_id", "gene_name",
         "transcript_id"]]
    gencode_sep_bed = BedTool.from_dataframe(gencode_separated)

    intersection_result = pfam_merged_bed.intersect(gencode_sep_bed, wa=True, wb=True)
    print("Number of records: {}".format(intersection_result.count()))
    intersection_result_df = intersection_result.to_dataframe(header=None, low_memory=False)

    intersection_result_df[14] = intersection_result_df[14].apply(
        lambda x: re.sub("(\:.*$)", "", re.sub("^[^:]*:", "", re.sub("\.\d+", "", x))))
    intersection_result_df[14] = intersection_result_df[14].apply(lambda s: "gene_ensembl_id={}".format(s))
    intersection_result_df[15] = intersection_result_df[15].apply(lambda s: "gene_symbol={}".format(s))
    intersection_result_df[3] = intersection_result_df[3].apply(lambda s: "pfam_domain_name={}".format(s))

    intersection_result_df[17] = intersection_result_df.apply(lambda s: ";".join([str(s[14]), str(s[15]), str(s[3])]),
                                                              axis=1)
    intersection_result_df = intersection_result_df.iloc[:, [0, 1, 2, 17]].drop_duplicates()
    db_bed = BedTool.from_dataframe(intersection_result_df, na_rep=".").sort()
    db_bed.saveas(PFAM_GENCODE_INTERSECT_PATH)


# TargetScan
TARGET_SCAN_INPUT_PATH = ""
TARGET_SCAN_OUTPUT_PATH = ""


def preprocess_target_scan():
    """
    Pre process targetScan
    """
    run_external_proc([f"{LIFTOVER_PATH} "
                       f"{TARGET_SCAN_INPUT_PATH} "
                       f"{HG19TOHG38_CONVERTOR} "
                       f"{TARGET_SCAN_OUTPUT_PATH} "
                       f"target_scan_unMapped"])
    targetscan_hg38_bed = BedTool(TARGET_SCAN_OUTPUT_PATH).sort()
    targetscan_hg38_df = targetscan_hg38_bed.to_dataframe(comment="#")
    # Update chrom name to NCBI format
    targetscan_hg38_df["chrom"] = targetscan_hg38_df["chrom"].apply(lambda x: x.strip("chr"))

    # Add gene symbol
    ensg_symbol_df = pd.read_csv(ENSG_SYMBOL_MAP, sep="\t")
    targetscan_hg38_df["gene_symbol"] = targetscan_hg38_df["name"].apply(lambda s: s.split(":")[0])
    targetscan_hg38_df = targetscan_hg38_df.merge(ensg_symbol_df, left_on="gene_symbol", right_on="gene_name")
    targetscan_hg38_df["name"] = targetscan_hg38_df.apply(lambda row: f"{row['name']}:{row['gene_id']}", axis=1)
    targetscan_hg38_df.drop(columns=['thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts',
                                     "gene_symbol", "gene_name", "gene_id"], inplace=True)
    db_bed = BedTool.from_dataframe(targetscan_hg38_df, na_rep=".")
    db_bed = db_bed.sort()
    db_bed.saveas(TARGET_SCAN_OUTPUT_PATH)


def present_gencode_file_stat():
    # examine the data we have
    examiner = GFFExaminer()
    in_handle = open(GENCODE_INPUT_PATH)
    available_limits = examiner.available_limits(in_handle)
    parent_child_map = examiner.parent_child_map(in_handle)

    pprint.pprint(available_limits)
    pprint.pprint(parent_child_map)
    in_handle.close()


def parse_arg(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--list", nargs="+", help="DB that need to preprocess", required=True)
    args = parser.parse_args(argv)
    return args


def main(argv):
    args = parse_arg(argv[1:])

    if "gencode" in args.list:
        print("Starting preprocess for gencode")
        time_start = perf_counter()
        preprocess_gencode()
        time_end = perf_counter()
        print("Total run for GENCODE: {}".format(timedelta(seconds=(time_end - time_start))))

    if "mirgenedb" in args.list:
        print("Starting preprocess for mirgenedb")
        time_start = perf_counter()
        preprocess_mirgenedb()
        time_end = perf_counter()
        print("Total run for MirGeneDB: {}".format(timedelta(seconds=(time_end - time_start))))

    if "remapepd" in args.list:
        print("Starting preprocess for remapepd")
        time_start = perf_counter()
        preprocess_remapepd()
        time_end = perf_counter()
        print("Total run for ReMap_EPD: {}".format(timedelta(seconds=(time_end - time_start))))

    if "enhanceratlas" in args.list:
        print("Starting preprocess for enhanceratlas")
        time_start = perf_counter()
        preprocess_enhanceratlas()
        time_end = perf_counter()
        print("Total run ofr Enhancer Atlas: {}".format(timedelta(seconds=(time_end - time_start))))

    if "omim" in args.list:
        print("Starting preprocess for omim")
        time_start = perf_counter()
        preprocess_omim()
        time_end = perf_counter()
        print("Total run for OMIM: {}".format(timedelta(seconds=(time_end - time_start))))

    if "humantf" in args.list:
        print("Starting preprocess for humantf")
        time_start = perf_counter()
        preprocess_humantf()
        time_end = perf_counter()
        print("Total run for HumanTF: {}".format(timedelta(seconds=(time_end - time_start))))

    if "proteinatlas" in args.list:
        print("Starting preprocess for proteinatlas")
        time_start = perf_counter()
        preprocess_proteinatlas()
        time_end = perf_counter()
        print("Total run for Protein Atlas: {}".format(timedelta(seconds=(time_end - time_start))))

    if "rbp" in args.list:
        print("Starting preprocess for rbp")
        time_start = perf_counter()
        preprocess_rbp()
        time_end = perf_counter()
        print("Total run for RBP: {}".format(timedelta(seconds=(time_end - time_start))))

    if "cosmic" in args.list:
        print("Starting preprocess for cosmic")
        time_start = perf_counter()
        preprocess_cosmic()
        time_end = perf_counter()
        print("Total run for COSMIC: {}".format(timedelta(seconds=(time_end - time_start))))

    if "pfam" in args.list:
        print("Starting preprocess for Pfam")
        time_start = perf_counter()
        preprocess_pfam_protein_domains()
        time_end = perf_counter()
        print("Total run for Pfam: {}".format(timedelta(seconds=(time_end - time_start))))

    if "targetscan" in args.list:
        print("Starting preprocess for TargetScan")
        time_start = perf_counter()
        preprocess_target_scan()
        time_end = perf_counter()
        print("Total run for TargetScan: {}".format(timedelta(seconds=(time_end - time_start))))
    # create_ensg_symbol_mapping()


if __name__ == "__main__":
    main(sys.argv)
