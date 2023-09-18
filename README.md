
# OffRisk User guide

OffRisk is a pipeline that comprise of two dockers – server docker and web UI docker.

The main function of the server is to analyze off-target and return the result of the analysis from
different databases.

The web UI enables users to easily run analysis with offRisk server for two use cases:

1. Off-target analysis – input from the user for off-target location and returning the results of
    the analysis.
2. On-target analysis – user input gRNA sequences, and then using <cite> [Cas-OFFinder][1] </cite>, <cite> [FlashFry][2] </cite> or <cite> [CRISPRitz][3] </cite>
    the pipeline will search for the off-target locations and then analysis them and return the
    results.

for usage on how is best work with Off-Risk, please refer the [Off-Risk UI repository](https://github.com/gili311/OffRisk-ui)

## Getting Started

<details>

<summary>Installation</summary><blockquote>

<details>

<summary>Docker</summary>

### Docker Installation

**For your convenience, we have also prepared an [docker installation video](https://www.dropbox.com/scl/fi/mie1rnb8b8r6oumjcs6yo/Docker-Installation-Made-with-Clipchamp.mp4?rlkey=6fcsm29adhfhws9vjq570tck1&dl=0) that you can refer to for guidance.**

To work with the docker files, docker and docker-compose must be installed.
https://docs.docker.com/get-docker/
#### If downloading for windows use the “WSL 2 backend” option

- Enable the WSL 2 feature on Windows. For detailed instructions, refer to the [Microsoft documentation](https://learn.microsoft.com/en-us/windows/wsl/install).
- Download and install the [Linux kernel update package](https://learn.microsoft.com/en-us/windows/wsl/install-manual#step-4---download-the-linux-kernel-update-package).
- Set WSL 2 as your default version, as it is needed by the docker.

Restart computer and run “Docker desktop” from the start menu / desktop shortcut (if you
created one)

This will require 8 GB. It contains the database folder, this user guide and docker-compose.yml file.


OffRisk dockers are on DockerHub: 

- [Off-Risk-Server](https://hub.docker.com/r/talmalulbgu/off-risk-server)
- [Off-Risk-UI](https://hub.docker.com/r/talmalulbgu/off-risk-ui)

</details>

<details>
<summary>Databases</summary>

### Databases Installation

**For your convenience, we have also prepared an [full Off-Risk installation video](https://www.dropbox.com/scl/fi/p8aqz75opxv73ro0ktoi7/Off-Risk-Installation-Made-with-Clipchamp.mp4?rlkey=qgomzpn9rkrh55bjc1auaardh&dl=0) that you can refer to for guidance.**

All necessary database files, including in the *Supported Databases* list bellow, except for OMIM and COSMIC. OMIM and COSMIC require licenses, can be found in the [OffRisk-db.zip archive](https://doi.org/10.5281/zenodo.8289271) on <cite> [Zenodo][4] </cite>.</br>
The total disk space required for this archive is approximately 10.5 GB.

#### Supported Databases:
 
- GENCODE
- MirGeneDB
- ReMapEPD
- EnhancerAtlas 2.0
- Pfam
- TargetScan 8.0
- OMIM
- HumanTF 3.0
- Protein Atlas
- RBP
- COSMIC

#### Database Files:

All necessary files for these databases are included in the [OffRisk-db.zip archive](https://doi.org/10.5281/zenodo.8289271), except for OMIM and COSMIC. OMIM and COSMIC require licenses. Therefore, users interested in using them should follow these steps:

#### COSMIC Database:

1. Download the COSMIC file from [https://cancer.sanger.ac.uk/cosmic/download](https://cancer.sanger.ac.uk/cosmic/download). Look for the "Cancer Gene Census" file.

2. Under "Tools", navigate to the database conversion page in OffRisk.

3. Upload the downloaded COSMIC file.

4. Click the "Pre-process COSMIC" button.

5. Download the result file with the given name.

6. Place the downloaded file in the "COSMIC" folder in your database directory. Name it "cosmic.csv".

#### OMIM Database:

1. Download two OMIM files from [https://www.omim.org/contact](https://www.omim.org/contact). You'll need both "mim2gene.txt" and "genemap2.txt".

2. Under "Tools", navigate to the database conversion page in OffRisk.

3. Upload both OMIM files.

4. Click the "Pre-process OMIM" button.

5. Download the result file with the given name.

6. Place the downloaded file in the "OMIM" folder in your database directory. Name it "omim.csv".

</details>
<details>

<summary>Off-Risk</summary>

### Download off-risk-sever and off-risk-ui dockers

**For your convenience, we have also prepared an [full Off-Risk installation video](https://www.dropbox.com/scl/fi/p8aqz75opxv73ro0ktoi7/Off-Risk-Installation-Made-with-Clipchamp.mp4?rlkey=qgomzpn9rkrh55bjc1auaardh&dl=0) that you can refer to for guidance.**

**Installing Off-Risk with Docker:**

1. Start by installing Docker on your computer if you haven't already. You can download it from [Docker's official website](https://www.docker.com/get-started).

2. Download the necessary databases for Off-Risk.

3. Next, you'll want to download the Off-Risk server and UI Docker containers. To do this, follow these steps:

   - Download the `docker-compose.yml` file from the [OffRisk repository](https://github.com/gili311/OffRisk/blob/main/docker/docker-compose.yml).
   
   - Your `docker-compose.yml` file should resemble the following configuration:

     ```yaml
     version: "3.9"  # Optional, required for Docker version v1.27.0 and later
     services:
       off-risk-server:
         ports:
           - "8123:80"
         volumes:
           - <path_to_database_folder>/databases:/databases
         image: almaliahbgu/off_risk:off-risk-server

       off-risk-ui:
         ports:
           - "8501:8501"
         image: almaliahbgu/off_risk:off-risk-ui
     networks:
       default:
         external:
           name: OffRisk-net
     ```
   
   - Replace `<path_to_database_folder>` under the `volumes` section with the actual path to your database folder.

4. Once you have the `docker-compose.yml` file, open Docker Desktop on your computer.

5. Open your terminal (e.g., Command Prompt or Terminal).

6. Navigate to the location where you downloaded the `docker-compose.yml` file using the `cd` command, e.g., `cd <path_to_docker-compose.yml_folder>`.

7. In the terminal, enter the following command: `docker-compose up`.</br>
This command will start the download and setup of the Off-Risk server and UI Docker containers on your computer.

By following these steps, you'll have Off-Risk installed and ready to use on your machine.

The off-risk-server and off-risk-ui dockers are also available to download without the docker-compose.yml file from:

</details></blockquote>

</details>




## Usage

&ensp; The web UI docker must communicate with the server docker.

&ensp; The server docker can work as standalone and integrate with other tools by using API.

&ensp; Validation on the input is done using python module pydantic and have the following structure:




## API Documentation

### `POST /v1/flashfry/`

&ensp; This endpoint is used to perform a FlashFry analysis for a list of target sites.

&ensp; **Request Body:**
 ```json
 {
   "request_id": 123,
   "sites": [
     {
       "sequence": "ATGCATGCATGC",
       "mismatch": 3
     },
     {
       "sequence": "CGTACGTACGTA",
       "mismatch": 2
     }
   ]
 }
 ```
 - `request_id`: An integer representing the request identifier for this on-target search.
 - `sites`: A list of objects representing target sites and their mismatch limits.
   - `sequence`: A sequence

### `POST /v1/on-target-analyze/`

&ensp; This endpoint is used to perform an on-target analysis for a list of target sites.

&ensp; **Request Body:**
```json
{
  "request_id": 456,
  "pam": "NGG",
  "downstream": true,
  "pattern_dna_bulge": 0,
  "pattern_rna_bulge": 0,
  "sites": [
    {
      "sequence": "ATGCATGCATGC",
      "mismatch": 4
    }
  ],
  "db_list": ["all"],
  "search_tools": ["flashfry"]
}
```
- `request_id`: An integer identifier for the on-target analysis request.
- `pam`: The Protospacer Adjacent Motif (PAM) pattern used for identifying target sites. (`default: NGG`)
- `downstream`: A boolean indicating whether to include downstream sites.
- `pattern_dna_bulge` and `pattern_rna_bulge`: Integers representing allowable DNA and RNA bulge patterns. (`default: 0 and 0`)
- `sites`: A list of objects representing target sites and their mismatch limits. Each site object includes:
  - `sequence`: The DNA sequence of the target site, **excluding the PAM**.
  - `mismatch`: The number of allowed mismatches for the target site. (`default: 4`)
- `db_list`: A list of supported databases to search in. Supported databases include: ["gencode", "mirgene", "remapepd", "enhanceratlas", "pfam", "targetscan", "omim", "humantf", "protein_atlas", "rbp", "cosmic"]. (`default: ["all"]`)
- `search_tools`: A list of search tools to be used. available search tools are ["cas_offinder", "flashfry", "crispritz"]. (`default: ["flashfry"]`)




### `POST /v1/off-target-analyze/`

&ensp;This endpoint is used to perform an off-target analysis for a list of target sites.


&ensp;**Request Body:**
```json
{
  "request_id": 789,
  "organism": "human",
  "off_targets": [
    {
      "chromosome": "chr1",
      "start": 100,
      "end": 150,
      "strand": "+",
      "sequence": "ATGCATGCATGC"
    },
    {
      "chromosome": "chr2",
      "start": 200,
      "end": 250,
      "strand": "-",
      "sequence": "CGTACGTACGTA"
    }
  ],
  "on_target": {
    "chromosome": "chr3",
    "start": 300,
    "end": 350,
    "strand": "+",
    "sequence": "TGCATGCATGCA"
  },
  "db_list": ["all"]
}

```
- `request_id`: An integer identifier for the off-target analysis request.
- `organism`: The organism for which the off-target analysis is performed.
- `off_targets`: A list of objects representing off-target locations. Each off-target object includes:
  - `chromosome`: The chromosome identifier of the off-target site.
  - `start and end`: The start and end positions of the off-target site.
  - `strand`: The DNA strand on which the off-target site is located (+ or -).
  - `sequence`: The DNA sequence of the off-target site.
- `on_target`: An object representing the on-target location. The on-target object includes:
  - `chromosome`: The chromosome identifier of the off-target site.
  - `start and end`: The start and end positions of the off-target site.
  - `strand`: The DNA strand on which the off-target site is located (+ or -).
  - `sequence`: The DNA sequence of the off-target site.
- `db_list`: A list of supported databases to search in. Supported databases include: ["gencode", "mirgene", "remapepd", "enhanceratlas", "pfam", "targetscan", "omim", "humantf", "protein_atlas", "rbp", "cosmic"].

## Contact Us

Your feedback is incredibly valuable to us. If you have any specific suggestions or encounter any issues, please don't hesitate to share them with us via the GitHub issue tracker.</br>
We are genuinely eager to hear your thoughts and ideas to continually improve our tool.


[1]: https://doi.org/10.1093/bioinformatics/btu048 

[2]: https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-018-0545-0

[3]: https://doi.org/10.1093/bioinformatics/btz867

[4]: https://www.re3data.org/repository/r3d100010468

