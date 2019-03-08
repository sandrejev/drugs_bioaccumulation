FROM rocker/r-ver:3.4.4
LABEL maintainer="Sergej Andrejev <sandrejev@gmail.com>"

COPY *.R  /drugs_bioaccumulation/
COPY data  /drugs_bioaccumulation/data

RUN apt-get -y update \
  && apt-get install -y \
    git-core \
    libcurl4-openssl-dev \
    libgmp-dev \
    libssl-dev \
    libxml2-dev \
    make \
    default-jre \
    default-jdk \
    libpcre++-dev \
    liblzma-dev \ 
    libbz2-dev \
    libz-dev \
    libnetcdf-dev \
  && R CMD javareconf

RUN ["install2.r", "versions"]
RUN ["Rscript", \
    "-e", "source('https://bioconductor.org/biocLite.R'); biocLite(c('pcaMethods', 'impute', 'mzR', 'MSnbase', 'xcms'), suppressUpdates=T)", \
    "-e", "versions::install.versions('curl', '3.1')", \
    "-e", "install.packages('devtools')", \
    "-e", "devtools::install_version('htmltools', version='0.3.6')", \
    "-e", "devtools::install_version('XML', version='3.98-1.11')", \
    "-e", "devtools::install_version('formatR', version='1.5')", \
    "-e", "devtools::install_version('caTools')"]

RUN ["Rscript", \
    "-e", "versions::install.versions('rJava', '0.9-10')", \
    "-e", "versions::install.versions('scales', '1.0.0')", \
    "-e", "versions::install.versions('imputeLCMD', '2.0')", \
    "-e", "versions::install.versions('norm', '1.0-9.5')", \
    "-e", "versions::install.versions('tmvtnorm', '1.4-10')", \
    "-e", "versions::install.versions('RSQLite', '2.1.1')"]

RUN ["Rscript", \
    "-e", "versions::install.versions('gmm', '1.6-2')", \
    "-e", "versions::install.versions('sandwich', '2.5-0')", \
    "-e", "versions::install.versions('Matrix', '1.2-12')", \
    "-e", "versions::install.versions('mvtnorm', '1.0-8')", \
    "-e", "versions::install.versions('gplots', '3.0.1')"]

RUN ["Rscript", \
    "-e", "versions::install.versions('Rcpp', '1.0.0')", \
    "-e", "versions::install.versions('reshape2', '1.4.3')", \
    "-e", "versions::install.versions('stringr', '1.3.1')", \
    "-e", "versions::install.versions('xlsx', '0.6.1')", \
    "-e", "versions::install.versions('readr', '1.1.1')"]

RUN ["Rscript", \
    "-e", "versions::install.versions('dplyr', '0.7.6')", \
    "-e", "versions::install.versions('networkD3', '0.4')", \
    "-e", "versions::install.versions('bindrcpp', '0.2.2')", \
    "-e", "versions::install.versions('vegan', '2.5-2')", \
    "-e", "versions::install.versions('lattice', '0.20-35')"]

RUN ["Rscript", \
    "-e", "versions::install.versions('permute', '0.9-4')", \
    "-e", "versions::install.versions('labdsv', '1.8-0')", \
    "-e", "versions::install.versions('cluster', '2.0.6')", \
    "-e", "versions::install.versions('MASS', '7.3-50')", \
    "-e", "versions::install.versions('mgcv', '1.8-23')"]

RUN ["Rscript", \
    "-e", "versions::install.versions('nlme', '3.1-131')", \
    "-e", "versions::install.versions('widyr', '0.1.1')", \
    "-e", "versions::install.versions('RColorBrewer', '1.1-2')", \
    "-e", "versions::install.versions('tidyr', '0.8.1')", \
    "-e", "versions::install.versions('colorspace', '1.3-2')"]

RUN ["Rscript", \
    "-e", "versions::install.versions('rjson', '0.2.20')", \
    "-e", "versions::install.versions('futile.logger', '1.4.3')", \
    "-e", "versions::install.versions('rstudioapi', '0.7')", \
    "-e", "versions::install.versions('fansi', '0.4.0')", \
    "-e", "versions::install.versions('codetools', '0.2-15')"]

RUN ["Rscript", \
    "-e", "versions::install.versions('doParallel', '1.0.11')", \
    "-e", "versions::install.versions('jsonlite', '1.5')", \
    "-e", "versions::install.versions('stevedore', '0.9.1')", \
    "-e", "versions::install.versions('broom', '0.5.0')", \
    "-e", "versions::install.versions('semver', '0.2.0')"]

RUN ["Rscript", \
    "-e", "versions::install.versions('backports', '1.1.2')", \
    "-e", "versions::install.versions('assertthat', '0.2.0')", \
    "-e", "versions::install.versions('lazyeval', '0.2.1')", \
    "-e", "versions::install.versions('cli', '1.0.1')"]

RUN ["Rscript", \
    "-e", "versions::install.versions('igraph', '1.2.2')", \
    "-e", "versions::install.versions('gtable', '0.2.0')", \
    "-e", "versions::install.versions('glue', '1.3.0')", \
    "-e", "versions::install.versions('RANN', '2.6.1')", \
    "-e", "versions::install.versions('MALDIquant', '1.18')"]

RUN ["Rscript", \
    "-e", "versions::install.versions('RJSONIO', '1.3-0')", \
    "-e", "versions::install.versions('gdata', '2.18.0')", \
    "-e", "versions::install.versions('iterators', '1.0.10')", \
    "-e", "versions::install.versions('xlsxjars', '0.6.1')", \
    "-e", "versions::install.versions('gtools', '3.8.1')"]

RUN ["Rscript", \
    "-e", "versions::install.versions('zoo', '1.8-3')", \
    "-e", "versions::install.versions('hms', '0.4.2')", \
    "-e", "versions::install.versions('lambda.r', '1.2.3')", \
    "-e", "versions::install.versions('yaml', '2.2.0')", \
    "-e", "versions::install.versions('stringi', '1.2.4')"]

RUN ["Rscript", \
    "-e", "versions::install.versions('foreach', '1.4.4')", \
    "-e", "versions::install.versions('rlang', '0.3.1')", \
    "-e", "versions::install.versions('pkgconfig', '2.0.2')", \
    "-e", "versions::install.versions('bitops', '1.0-6')"]

RUN ["Rscript", \
    "-e", "versions::install.versions('purrr', '0.2.5')", \
    "-e", "versions::install.versions('bindr', '0.1.1')", \
    "-e", "versions::install.versions('htmlwidgets', '1.2')", \
    "-e", "versions::install.versions('labeling', '0.3')", \
    "-e", "versions::install.versions('cowplot', '0.9.3')"]

RUN ["Rscript", \
    "-e", "versions::install.versions('tidyselect', '0.2.4')", \
    "-e", "versions::install.versions('magrittr', '1.5')", \
    "-e", "versions::install.versions('R6', '2.3.0')", \
    "-e", "versions::install.versions('pillar', '1.3.1')", \
    "-e", "versions::install.versions('whisker', '0.3-2')"]

RUN ["Rscript", \
    "-e", "versions::install.versions('withr', '2.1.2')", \
    "-e", "versions::install.versions('survival', '2.41-3')", \
    "-e", "versions::install.versions('tibble', '1.4.2')", \
    "-e", "versions::install.versions('crayon', '1.3.4')", \
    "-e", "versions::install.versions('futile.options', '1.0.1')"]

RUN ["Rscript", \
    "-e", "versions::install.versions('KernSmooth', '2.23-15')", \
    "-e", "versions::install.versions('utf8', '1.1.4')", \
    "-e", "versions::install.versions('digest', '0.6.18')", \
    "-e", "versions::install.versions('munsell', '0.5.0')", \
    "-e", "versions::install.versions('remotes', '1.1.1')", \
    "-e", "versions::install.versions('tidyverse', '1.2.1')"]
RUN ["installGithub.r", "ramnathv/rCharts@8d3fe35be5a4b41907d800944a14ac0d193c5b1a", "slowkow/ggrepel@3a14848bacbd5a4ce6a5a857f16af898e0eb1bfc"]

WORKDIR /drugs_bioaccumulation
RUN ["Rscript", "/drugs_bioaccumulation/run_everything.R"]
CMD ["Rscript", "/drugs_bioaccumulation/run_everything.R"]
