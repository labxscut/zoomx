#this should be working with ubuntu:latest
#NOTE: developer tested the commands by folloing actions:
#      docker pull ubuntu:latest
#      docker images
#      docker run --memory=2g -i -t ubuntu:latest /bin/bash
#      docker run --memory=2g -i -t id /bin/bash      
#      docker-machine scp Dockerfile main:~/zoomx
#      docker build --memory=2g ~/zoomx
#      docker start --memory=2g 301176b69086
#      docker exec -it 301176b69086 /bin/bash
FROM ubuntu:latest

MAINTAINER Charlie Xia <lxia.usc@gmail.com>

WORKDIR $HOME/
RUN mkdir $HOME/setup
RUN mkdir $HOME/bin
RUN export PATH=$PATH:$HOME/bin
RUN ulimit -s unlimited

### install git and dev files ###
RUN apt-get update && apt-get install -y build-essential software-properties-common
RUN apt-get install -y zlib1g-dev libncurses5-dev wget git unzip libcurl4-openssl-dev libxml2-dev
RUN apt-get install -y bzip2 libbz2-dev liblzma-dev

### python, has to precede bedtools ###
RUN apt-get install -y python libpython-dev python-pip
RUN apt-get install -y python-numpy python-scipy
### install python packages ###
RUN pip install -U numpy
RUN pip install -U scipy
RUN pip install -U tables
RUN pip install -U six
RUN pip install -U pandas
RUN pip install -U pysam
RUN pip install -U pybedtools

### install bedtools ###
RUN cd $HOME/setup && git clone https://github.com/arq5x/bedtools.git
RUN cd $HOME/setup/bedtools && make && cp bin/* $HOME/bin

### install samtools ###
RUN cd $HOME/setup && git clone https://github.com/samtools/samtools.git
RUN cd $HOME/setup && git clone https://github.com/samtools/htslib.git
RUN cd $HOME/setup/samtools && make && cp samtools $HOME/bin

### R ###
RUN apt-get install -y r-base 
### install R packages ###
RUN mkdir -p $(Rscript -e 'cat(Sys.getenv("R_LIBS_USER"))')
#CRAN: optparse, rjson, RColorBrewer
RUN Rscript -e 'install.packages("optparse",repos="http://cran.us.r-project.org")'  
RUN Rscript -e 'install.packages("rjson",repos="http://cran.us.r-project.org")'  
RUN Rscript -e 'install.packages("RColorBrewer",repos="http://cran.us.r-project.org")'  
#Bioconductor: GenomicRanges, ggbio, GenomeInfoDb
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("rhdf5",ask=F)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("GenomicRanges",ask=F)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("GenomeInfoDb",ask=F)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("BSgenome.Hsapiens.UCSC.hg19",ask=F)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("BSgenome.Hsapiens.UCSC.hg38",ask=F)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("Rsamtools",ask=F)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("Biostrings",ask=F)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("ggbio",ask=F)'

### install zoomx ###
RUN cd $HOME/setup && git clone https://bitbucket.org/charade/zoomx.git
RUN cd $HOME/setup/zoomx && python setup.py install
RUN cd $HOME/setup/zoomx/test && ./test_dep.sh

### run test scripts ### 
#RUN export PATH=$PATH:$HOME/bin
#RUN cd $HOME/setup/zoomx/test
#RUN wget https://s3-us-west-2.amazonaws.com/lixiabucket/zoomx.test.full.tgz
#RUN tar -zxvf zoomx.test.full.tgz
#RUN zoomxStage germline10.cfg
#RUN zoomxStage somatic10.cfg
