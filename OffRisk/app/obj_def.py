import re
from pydantic import BaseModel, Json, validator, Field
from typing import List, Union

from configuration_files.const import DB_NAME_LIST


class ProcessedResult(BaseModel):
    """
    Object definition for each saved DB
    """
    name: str
    description: str = None
    type: str = "DataFrame"
    data: Union[Json, List, str] = None


class AllNoneHumanDbResult(BaseModel):
    """
    Object definition for structure of a database object
    """
    S3_result_list: List[ProcessedResult] = None


class AllDbResult(BaseModel):
    """
    Object definition for structure of all databases object
    """
    GENCODE_result_list: List[ProcessedResult] = None
    MirGene_result_list: List[ProcessedResult] = None
    ReMapEPD_result_list: List[ProcessedResult] = None
    EnhancerAtlas_result_list: List[ProcessedResult] = None
    Pfam_result_list: List[ProcessedResult] = None
    TargetScan_result_list: List[ProcessedResult] = None
    OMIM_result_list: List[ProcessedResult] = None
    HumanTF_result_list: List[ProcessedResult] = None
    Protein_Atlas_result_list: List[ProcessedResult] = None
    RBP_result_list: List[ProcessedResult] = None
    COSMIC_result_list: List[ProcessedResult] = None


class OtResponse(BaseModel):
    """
    Object definition for response from the server with the result for each off-target
    """
    request_id: int
    off_targets: Json  # Dataframe as json object. This is the complete summary
    flashfry_score: Json  # If FlashFry was used, Dataframe with the result of the score
    all_result: AllDbResult  # All the databases inforamtion
    time: float  # Totoal time for running the analysis


class OtNoneHumanResponse(BaseModel):
    """
    Object definition for response from the server with the result for each off-target
    """
    request_id: int
    off_targets: Json
    flashfry_score: Json
    all_result: AllNoneHumanDbResult
    time: float

    @validator("off_targets")
    def val_off_targets(cls, v):
        for item in v:
            if item['chromosome'].startswith('chr'):
                item['chromosome'] = item['chromosome'][3:]
        return v


class OffTarget(BaseModel):
    """
    Object definition for Off-target location for off-target analysis
    """
    chromosome: Union[str, int] = Field(..., title="chromosome")
    start: int = Field(..., title="start")
    end: int = Field(..., title="end")
    strand: str = Field(None, title="strand")
    id: int = Field(None, title="id") # Field for CRISPRIL collaboration
    sequence: str = Field(None, title="sequence") # Field for CRISPRIL collaboration

    @classmethod
    def get_fields_name(cls, alias=False):
        return cls.schema(alias).get("properties")

    @classmethod
    def get_fields_title(cls):
        return [field['title'] for field in cls.schema(False).get("properties").values()]

    @classmethod
    def get_field_title(cls, field_name):
        fields = cls.schema(False).get("properties")
        return fields[field_name]['title'] if field_name in fields else None

    @validator("end")
    def val_end(cls, v, values):
        assert isinstance(v, int)
        if "start" in values and v < values["start"]:
            raise ValueError("End need to be larger then start")
        return v

    @validator("strand")
    def val_strand(cls, v):
        assert v in ["+", "-"]
        return v

    @validator("sequence")
    def val_seq(cls, v):
        if v:
            if not re.fullmatch(r"^[AaGcTtCcUu]*$", v):
                raise ValueError("Got invalid string for site. Only A, T, U, C and G are allowed")
            if not len(v) < 24:
                raise ValueError("seq can not be longer then 24")
        return v


class OffTargetList(BaseModel):
    """
    Object definition for List of OffTarget. This is the server input.
    """
    request_id: int
    organism: str = 'human'
    off_targets: List[OffTarget]
    on_target: OffTarget = None
    db_list: List[str] = ["all"]

    @validator("db_list")
    def val_db_list(cls, v):
        db = DB_NAME_LIST + ["all"]
        test = all(item in db for item in v)
        if not test:
            raise ValueError("Not a valid DB")
        return v


class Site(BaseModel):
    """
    Object definition for on-target search
    """

    sequence: str  # include sequences and pam
    mismatch: int = 4

    @validator("sequence")
    def val_site(cls, v):
        if not re.fullmatch(r"^[AaGcTtCcUu]*$", v):
            raise ValueError("Got invalid string for site. Only A, T, U, C and G are allowed")
        if len(v) == 0:
            raise ValueError("Got empty sequence. Please specify A sequence")
        return v

    @validator("mismatch")
    def val_mismatch(cls, v):
        if v > 20:
            raise ValueError("Number is greater then 20. Only number Between 0 to 20 are allowed")
        if v < 0:
            raise ValueError("Number is less then 0. Only number Between 0 to 20 are allowed")
        return v


class SitesList(BaseModel):
    """
    Object definition for List of Sites. This is The input for on-target search.
    """

    request_id: int
    pattern: str = "NNNNNNNNNNNNNNNNNNNNNGG"
    pattern_dna_bulge: int = 0
    pattern_rna_bulge: int = 0
    sites: List[Site]
    db_list: List[str] = ["all"]
    search_tools: List[str] = ["flashfry"]

    @validator("db_list")
    def val_db_list(cls, v):
        db = DB_NAME_LIST + ["all"]
        test = all(item in db for item in v)
        if not test:
            raise ValueError("Not a valid DB")
        return v

    @validator("pattern")
    def val_pattern(cls, v):
        if not re.fullmatch(r"[AGTCRYSWKMBDHVN]*$", v):
            raise ValueError("Not a valid pattern. Please refer for Cas-Offinder documentation")
        return v

    @validator("pattern_dna_bulge", "pattern_rna_bulge")
    def val_bulge(cls, v):
        if v > 20:
            raise ValueError("Number is greater then 20. Only number Between 0 to 20 are allowed")
        if v < 0:
            raise ValueError("Number is less then 0. Only number Between 0 to 20 are allowed")
        return v

    @validator("search_tools")
    def val_search_tool(cls, v):
        tools = ["flashfry", "cas_offinder"]
        test = all(item in tools for item in v)
        if not test:
            raise ValueError("Not a valid tools. Only {} are supported".format(tools))
        return v


class FlashFrySite(BaseModel):
    request_id: int
    sites: List[Site]
