[![Build](https://github.com/crispr-il/off-tov/actions/workflows/ci.yml/badge.svg)](https://github.com/crispr-il/off-tov/actions/workflows/ci.yml)
[![Deploy to EKS](https://github.com/crispr-il/off-tov/actions/workflows/cd.yml/badge.svg)](https://github.com/crispr-il/off-tov/actions/workflows/cd.yml)


# off-tov

<img src="https://user-images.githubusercontent.com/553010/109801117-b576df00-7c26-11eb-99c1-fdd210ddafdf.png" width="300" height="300" />



### Build the docker locally
~~~
docker build --no-cache -t off-tov-server ./ -f xeon.flask.Dockerfile
~~~

### Run the docker
~~~
docker run -p8123:80 -v/mnt/gogenome/off-tov/databases:/databases off-tov-server
~~~
If you are running with off-tov UI as well you need to connect them to the same network.
To create the network and then connect the docker to it:

~~~
docker network create off-tov-net
docker run -p8123:80 --net off-tov-net -v/mnt/gogenome/off-tov/databases:/databases off-tov-server
~~~

For debug:
~~~
docker run -ti --entrypoint bash -v/mnt/gogenome/off-tov/databases:/databases off-tov-server
~~~


### content of input.txt (see where the genome is)
~~~
/host/home/zozo123/genomes/GCF_000004515.5_Glycine_max_v2.1_genomic_for_CRISPRIL.fna
NNNNNNNNNNNNNNNNNNNNNRG
GGCCGACCTGTCGCTGACGCNNN 5
CGCCAGCGTCAGCGACAGGTNNN 5
ACGGCGCCAGCGTCAGCGACNNN 5
GTCGCTGACGCTGGCGCCGTNNN 5
~~~

### inside the docker run:
~~~
./cas-offinder /host/path/to/input.txt C ./output.txt
~~~

## Build the docker:
~~~
docker build -t cas-offinder-local-server ./ -f xeon.flask.Dockerfile
~~~

## Run the docker server:
if gogenome is mounted run:
~~~
docker run -p8123:80 -v/mnt/gogenome/:/gogenome off-tov-server
~~~
For debug:
~~~
docker run -ti --entrypoint bash -v/mnt/gogenome/:/gogenome off-tov-server
~~~


else run:
~~~
docker run -p8123:80 -v/:/host off-tov-server
~~~
for debug: 
~~~
docker run -ti --entrypoint bash -v/:/host off-tov-server
~~~
