FROM ubuntu:20.04

# ARGUMENTS
###########
ARG OPENSTRUCTURE_VERSION="2.3.1"
ARG SRC_FOLDER="/usr/local/src"
ARG CPUS_FOR_MAKE=2
ARG OPENMM_VERSION="7.1.1"
ARG OPENMM_INCLUDE_PATH="/usr/local/openmm/include/"
ARG OPENMM_LIB_PATH="/usr/local/openmm/lib/"
ARG DEBIAN_FRONTEND=noninteractive

# INSTALL SYSTEM DEPS
#####################
RUN apt-get update -y && apt-get install -y cmake \
                                            g++ \
                                            gfortran \
                                            wget \
                                            tar \
                                            libsqlite3-dev \
                                            sip-dev \
                                            libtiff-dev \
                                            libfftw3-dev \
                                            libeigen3-dev \
                                            libboost-all-dev \
                                            libpng-dev \
                                            python3-all \
                                            python3-numpy \
                                            python3-scipy \
                                            python3-pandas \
                                            doxygen \
                                            swig \
                                            clustalw \
                                            locales && \
                                            # CLEANUP
                                            rm -rf /var/lib/apt/lists/*

# INSTALL OPENMM
################
RUN cd ${SRC_FOLDER} && \
    wget -O openmm-${OPENMM_VERSION}.tar.gz -nc https://github.com/pandegroup/openmm/archive/${OPENMM_VERSION}.tar.gz && \
    mkdir ${SRC_FOLDER}/openmm-${OPENMM_VERSION} && \
    tar xf openmm-${OPENMM_VERSION}.tar.gz -C ${SRC_FOLDER}/openmm-${OPENMM_VERSION} --strip-components=1 && \
    mkdir -p ${SRC_FOLDER}/openmm-${OPENMM_VERSION}/build && \
    cd ${SRC_FOLDER}/openmm-${OPENMM_VERSION}/build && \
    cmake .. && make -j $CPUS_FOR_MAKE && make install && \
    cd ${SRC_FOLDER}/openmm-${OPENMM_VERSION}/build/python && \
    python3 setup.py build && python3 setup.py install && \
    rm ${SRC_FOLDER}/openmm-${OPENMM_VERSION}.tar.gz && \
    rm -rf ${SRC_FOLDER}/openmm-${OPENMM_VERSION}

# INSTALL OST
#############
RUN cd ${SRC_FOLDER} && \
    # copy ost release
    wget -O openstructure-${OPENSTRUCTURE_VERSION}.tar.gz -nc https://git.scicore.unibas.ch/schwede/openstructure/-/archive/${OPENSTRUCTURE_VERSION}/openstructure-${OPENSTRUCTURE_VERSION}.tar.gz && \
    mkdir openstructure-${OPENSTRUCTURE_VERSION} && \
    tar xf openstructure-${OPENSTRUCTURE_VERSION}.tar.gz -C ${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION} --strip-components=1 && \
    mkdir -p ${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION}/build && \
    cd ${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION}/build && \
    cmake .. -DOPTIMIZE=ON \
             -DENABLE_MM=ON \
             -DCOMPILE_TMTOOLS=1 \
             -DUSE_NUMPY=1 \
             -DOPEN_MM_LIBRARY=$OPENMM_LIB_PATH/libOpenMM.so \
             -DOPEN_MM_INCLUDE_DIR=$OPENMM_INCLUDE_PATH \
             -DOPEN_MM_PLUGIN_DIR=$OPENMM_LIB_PATH/plugins \
             -DENABLE_GFX=ON \
             -DENABLE_GUI=OFF \
             -DENABLE_INFO=OFF \
             -DCMAKE_C_FLAGS="-isystem /usr/include/boost/ -isystem ${OPENMM_INCLUDE_PATH}/include" \
             -DCMAKE_CXX_FLAGS="-isystem /usr/include/boost/ -isystem ${OPENMM_INCLUDE_PATH}/include" && \
    make -j ${CPUS_FOR_MAKE} && \
    wget ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz && \
    stage/bin/chemdict_tool create components.cif.gz compounds.chemlib pdb && stage/bin/chemdict_tool update ../modules/conop/data/charmm.cif compounds.chemlib charmm && \
    cmake .. -DCOMPOUND_LIB=${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION}/build/compounds.chemlib && \
             make -j ${CPUS_FOR_MAKE} && make check && make install && \
    rm ${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION}.tar.gz
    #rm -rf ${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION}

# INSTALL MUSTANG
#################

RUN wget https://lcb.infotech.monash.edu/mustang/mustang_v3.2.3.tgz &&\
    tar -zxvf mustang_v3.2.3.tgz && \
    cd MUSTANG_v3.2.3 && \
    make


# INSTALL PYTHON PAKS
#####################
COPY python_paks.txt .
RUN apt update -qq && \
    apt install -y --no-install-recommends software-properties-common dirmngr && \
    add-apt-repository universe && \
    apt install -y python3-pip && \
    pip3 install -r python_paks.txt

# INSTALL R
############
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    apt install -y --no-install-recommends r-base

COPY install_rpaks.R .
RUN Rscript install_rpaks.R

# ENVIRONMENT
#############
WORKDIR /home
ENV OST_ROOT="/usr/local"
ENV PYTHONPATH="/usr/local/lib64/python3.8/site-packages"
ENV PYTHONPATH="/usr/local/src/openstructure-2.3.1/build/stage/lib64/python3.8/site-packages:$PYTHONPATH"
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib64:/usr/local/openmm/lib/"
ENV OPENSTRUCTURE_VERSION=$OPENSTRUCTURE_VERSION
ENV PATH="$OST_ROOT:${PATH}"
ENV PATH="/MUSTANG_v3.2.3/bin/:${PATH}"

EXPOSE 8900
CMD jupyter lab --no-browser --port 8900 --allow-root --ip 0.0.0.0

