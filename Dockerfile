FROM rocker/r-base:4.4

ARG DEBIAN_FRONTEND=noninteractive
ARG GH_PAT='NOT_SET'

ADD . /SpectraWrapper

RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    uuid-dev \
    libxml2-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    libsqlite3-dev \
    pkg-config \
    git-all \
    wget \
    libbz2-dev \
    zlib1g-dev \
    python3-dev \
    libffi-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev && \
    mkdir /GW_Python && \
    cd /GW_Python && \
    wget https://github.com/python/cpython/archive/refs/tags/v3.8.10.tar.gz && \
    tar -zxvf v3.8.10.tar.gz && \
    cd cpython-3.8.10 && \
    ./configure --prefix=/GW_Python && \
    cd /GW_Python/cpython-3.8.10 && \
    make && \
    make install && \
    /GW_Python/bin/pip3 --no-cache-dir install numpy scipy scikit-learn scikit-misc matplotlib tqdm sympy setuptools pandas pyyaml  && \
    /GW_Python/bin/pip3 --no-cache-dir install scSpectra && \
    /GW_Python/bin/pip3 --no-cache-dir install dill && \
    chmod -R 777 /GW_Python


#install R and dependencies

RUN apt-get update && apt-get install -y r-base r-base-dev && \
    if [ "${GH_PAT}" != 'NOT_SET' ]; then \
        echo 'Setting GH_PAT'; \
        export GITHUB_PAT="${GITHUB_TOKEN}"; \
    fi && \
    Rscript -e "install.packages('Matrix', repos='http://cran.us.r-project.org', version='1.6.4')" && \
    Rscript -e "install.packages(c('devtools', 'remotes', 'Seurat', 'SeuratObject', 'data.table', 'jsonlite', 'readr', 'BiocManager'))" && \
    Rscript -e "BiocManager::install('DropletUtils', ask = FALSE, upgrade = 'always')"

#installing RIRA using remotes/devtools install_github() was giving me a
#"malformed DESCRIPTION" error on the encoding line of DESCRIPTION.
#makes very little sense to me as to why (the characters are fine/UTF-8).
#trying to build from source through git instead.

RUN cd / && \
    git clone https://github.com/BimberLab/RIRA.git && \
    cd /RIRA && \
    R CMD build . && \
    #TODO: remove ls once I have a handle on base file structure in /RIRA
    ls && \
    R CMD INSTALL --build *.tar.gz && \
	  rm -Rf /tmp/downloaded_packages/ /tmp/*.rds

#build SpectraWrapper
RUN cd /SpectraWrapper && \
    R CMD build . && \
	  R CMD INSTALL --build *.tar.gz && \
	  rm -Rf /tmp/downloaded_packages/ /tmp/*.rds

ENV NUMBA_CACHE_DIR=/work/numba_cache
ENV MPLCONFIGDIR=/work/mpl_cache

ENTRYPOINT ["/bin/bash"]
