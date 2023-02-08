# offRisk

A new pipeline that aims to help design gRNA by analyzing the off-targets location on the genome by using biological
information in different levels of DNA, RNA and protein gathered from different databases, and presenting the results to
the user

OffRisk is running with the UI OffRisk-UI, and both need to be on the same network.
To create the network:

~~~
docker network create off-risk-net
~~~

## Build the docker:

In order to build the docker you need to make sure the full path for the folder to be build is 
in the context field in docker-compose.yaml
~~~
docker-compose build
~~~

## Run the docker server:
The path for the volume with the database file need to be updates in the docker-compose.yaml:

for example:

<path_to_database_folder>:/databases -> /user/gilad/off-risk/databases:/databases

~~~
docker-compose up
~~~