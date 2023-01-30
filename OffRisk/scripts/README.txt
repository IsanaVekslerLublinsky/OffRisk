Preprocess db

For each DB there is relevant path (input, output, etc) - search for the DB and update the paths.
Supported DB:
gencode, mirgenedb, remapepd, enhanceratlas, pfam, targetscan, omim, humantf, proteinatlas, rbp, cosmic,

Example for preprocessing gencode:
"python preprocess.py -l gencode"


The following tools were used:
liftover:
1. Download the binaries. can be found here: https://genome.sph.umich.edu/wiki/LiftOver#Binary_liftOver_tool
2. change the path to it in the file (LIFTOVER_PATH)


