FROM ubuntu:jammy-20231128

LABEL maintainer "Daniel Park <dpark@broadinstitute.org>"

# non-interactive session just for build
ARG DEBIAN_FRONTEND=noninteractive

# update apt database and install R apt repo; install all desired packages
RUN apt-get update && \
  apt-get -y -qq install software-properties-common && \
  apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
  add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/' && \
  apt-get update && \
  apt-get -y -qq install \
    less nano vim git wget curl jq zstd pigz parallel locales \
    gnupg libssl-dev libcurl4-openssl-dev \
    libgsl-dev libxml2 libxml2-dev \
    imagemagick libmagick++-dev \
    texlive-base texlive-latex-recommended texlive texlive-latex-extra texlive-extra-utils texlive-fonts-extra \
    fonts-roboto \
    r-base r-base-dev r-cran-devtools \
    r-cran-tidyverse r-cran-extradistr \
    r-cran-rcpp r-cran-rcppgsl r-cran-rcppparallel \
    r-cran-segmented r-cran-pixmap r-cran-ape r-cran-seqinr r-cran-ade4 \
  && apt-get clean

# Set default locale to en_US.UTF-8
RUN locale-gen en_US.UTF-8
ENV LANG="en_US.UTF-8" LANGUAGE="en_US:en" LC_ALL="en_US.UTF-8"

# Install necessary R dependencies for reconstructR that don't have apt packages
RUN R -e "for (lib in c( 'LaplacesDemon', 'kmer', 'phylogram', 'aphid', 'insect' )) { install.packages(lib, dependencies=TRUE); library(lib, character.only=TRUE) }"
# Rfast in CRAN is broken, install from github
RUN R -e "devtools::install_github('RfastOfficial/Rfast', dependencies=TRUE); library(Rfast)"

# Install reconstructR R package
COPY . /opt/reconstructR
RUN R -e "devtools::install_local('/opt/reconstructR', dependencies=TRUE, upgrade='never'); library(reconstructR)"

# Bash prompt
ENV PATH="/opt/reconstructR/scripts:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
CMD ["/bin/bash"]
