# Builds a Docker image for VALENCE https://github.com/VALENCE-software/VALENCE
# docker build . -t <MYIMAGE>
# docker run -it <MYIMAGE> bash
# Authors:
# Murat Keceli <keceli@gmail.com>

FROM keceli/chembox
LABEL maintainer "Murat Keceli <keceli@gmail.com>"

ENV MPICH_VERSION=3.2.1

# Install system packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        gcc \
        gfortran \
        git \
        curl \
        ca-certificates \
        wget && \
    apt-get clean && \
    apt-get autoremove && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    mkdir -p /container && \
    cd /container && \
    wget http://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz && \
    tar xf mpich-${MPICH_VERSION}.tar.gz && \
    rm -f  mpich-${MPICH_VERSION}.tar.gz  && \
    cd mpich-${MPICH_VERSION} && \
    ./configure --prefix=/container/mpich-${MPICH_VERSION}/install --disable-wrapper-rpath && \
    make install 
ENV PATH=$PATH:/container/mpich-${MPICH_VERSION}/install/bin
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/container/mpich-${MPICH_VERSION}/install/lib

RUN cd /container && \
    git clone https://github.com/VALENCE-software/VALENCE.git && \
    cd VALENCE && \
    ./install-simint.sh && \
    SEQUENTIAL=false make && \
    conda env update -n base -f environment.yml
