# Flask serving via nginx
FROM tiangolo/uwsgi-nginx-flask:python3.8
# FROM tiangolo/uvicorn-gunicorn:python3.8
#RUN echo "deb http://ftp.us.debian.org/debian testing main contrib non-free" >> /etc/apt/sources.list

RUN apt-get update && apt-get install --no-install-recommends -y \
    apt-utils \
    git \
    unzip \
    zip \
    tar \
    curl \
    wget \
    xz-utils \
    alien \
    clinfo \
    software-properties-common \
    libc6-dev libc6 libc-bin \
    build-essential \
    pkg-config \
    cmake \
    ca-certificates \
    gnupg \
    gpg-agent \
    bedtools \
    ;


ARG DEBIAN_FRONTEND=noninteractive
ARG APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1

RUN wget https://github.com/intel/compute-runtime/releases/download/21.13.19438/intel-gmmlib_20.4.1_amd64.deb && \
    wget https://github.com/intel/intel-graphics-compiler/releases/download/igc-1.0.6748/intel-igc-core_1.0.6748_amd64.deb && \
    wget https://github.com/intel/intel-graphics-compiler/releases/download/igc-1.0.6748/intel-igc-opencl_1.0.6748_amd64.deb && \
    wget https://github.com/intel/compute-runtime/releases/download/21.13.19438/intel-opencl_21.13.19438_amd64.deb && \
    wget https://github.com/intel/compute-runtime/releases/download/21.13.19438/intel-ocloc_21.13.19438_amd64.deb && \
    wget https://github.com/intel/compute-runtime/releases/download/21.13.19438/intel-level-zero-gpu_1.0.19438_amd64.deb && \
    dpkg -i *.deb && \
    wget -qO - https://repositories.intel.com/graphics/intel-graphics.key | apt-key add - && \
    apt-add-repository 'deb [arch=amd64] https://repositories.intel.com/graphics/ubuntu focal main'

WORKDIR /tmp
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB && \
    apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB && \
    rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB && \
    echo "deb https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list && \
    wget -qO - https://repositories.intel.com/graphics/intel-graphics.key | apt-key add - && \
    apt-add-repository 'deb [arch=amd64] https://repositories.intel.com/graphics/ubuntu focal main' && \
    apt-get update -y && \
    apt-get install -y --no-install-recommends \
    intel-oneapi-runtime-ccl \
    intel-oneapi-runtime-compilers \
    intel-oneapi-runtime-dal \
    intel-oneapi-runtime-dnnl \
    intel-oneapi-runtime-dpcpp-cpp \
    intel-oneapi-runtime-dpcpp-library \
    intel-oneapi-runtime-fortran \
    intel-oneapi-runtime-ipp \
    intel-oneapi-runtime-ipp-crypto \
    intel-oneapi-runtime-libs \
    intel-oneapi-runtime-mkl \
    intel-oneapi-runtime-mpi \
    intel-oneapi-runtime-opencl \
    intel-oneapi-runtime-openmp \
    intel-oneapi-runtime-tbb \
    intel-oneapi-runtime-vpl

# Add FlashFry relevant files
ENV SDKMAN_DIR /root/.sdkman
ENV JAVA_VERSION 8.0.282-open

RUN rm /bin/sh && ln -s /bin/bash /bin/sh && \
    curl -s "https://get.sdkman.io" | bash && \
    chmod a+x "$SDKMAN_DIR/bin/sdkman-init.sh" && \
    set -x \
    && echo "sdkman_auto_answer=true" > $SDKMAN_DIR/etc/config \
    && echo "sdkman_auto_selfupdate=false" >> $SDKMAN_DIR/etc/config \
    && echo "sdkman_insecure_ssl=false" >> $SDKMAN_DIR/etc/config

WORKDIR $SDKMAN_DIR
RUN [[ -s "$SDKMAN_DIR/bin/sdkman-init.sh" ]] && source "$SDKMAN_DIR/bin/sdkman-init.sh" && exec "$@" && \
    source /root/.bashrc && \
    source "$SDKMAN_DIR/bin/sdkman-init.sh" && sdk list java && sdk install java $JAVA_VERSION

# Add FlashFry and cas-offinder
WORKDIR /app/tools
RUN git clone https://github.com/hyugel/cas-offinder-bulge . && \
    chmod 0777 cas-offinder-bulge && chmod 0777 cas-offinder && \
    wget https://github.com/mckennalab/FlashFry/releases/download/1.12/FlashFry-assembly-1.12.jar && \
    chmod 0777 FlashFry-assembly-1.12.jar


ENV PATH="/app/tools:${PATH}"
ENV PATH="/root/.sdkman/candidates/java/current/bin:${PATH}"

COPY ./requirements.txt /off-tov/requirements.txt
RUN pip install -r /off-tov/requirements.txt

COPY ./app /off-tov/app
COPY ./timeout.conf  /etc/nginx/conf.d/timeout.conf
WORKDIR /app/tmp
WORKDIR /off-tov/log
WORKDIR /off-tov/app/off_target_output
WORKDIR /off-tov/app

RUN rm -rf /app/tools/cas-offinder.zip