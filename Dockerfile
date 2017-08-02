#NOTE: developer tested the commands by folloing actions:
#      docker pull ubuntu:latest
#      docker images
#      docker run --memory=2g -i -t ubuntu:latest /bin/bash
#      docker run --memory=2g -i -t id /bin/bash      
#      docker-machine scp Dockerfile main:~/zoomx
#      docker build --memory=2g ~/zoomx
#      docker start --memory=2g 301176b69086
#      docker exec -it 301176b69086 /bin/bash
#NOTE: merging RUN as file become stable as every RUN creates a commit which has a limit
FROM ubuntu:xenial

MAINTAINER Charlie Xia <lxia.usc@gmail.com>

WORKDIR $HOME/
RUN ulimit -s unlimited
RUN mkdir $HOME/setup && mkdir $HOME/bin

### install git and dev files ###
#NOTE: RUN echo "PATH=$PATH:$HOME/bin" >>/etc/enviroment #not working, why?
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN gpg -a --export E084DAB9 >$HOME/rstudio.key
RUN apt-key add $HOME/rstudio.key
RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" >>/etc/apt/sources.list
RUN echo "PATH=$PATH:$HOME/bin" >>/etc/enviroment
RUN apt-get update 
RUN apt-get install -y build-essential software-properties-common
RUN apt-get install -y bzip2 libbz2-dev liblzma-dev openssl libssl-dev
RUN apt-get install -y zlib1g-dev libncurses5-dev wget git unzip libcurl4-openssl-dev libxml2-dev

### python, has to precede bedtools ###
RUN apt-get install -y python libpython-dev python-pip
RUN apt-get install -y python-numpy python-scipy
### install python packages ###
RUN pip install -U numpy scipy tables six pandas pysam pybedtools

### install bedtools ###
RUN cd $HOME/setup && git clone https://github.com/arq5x/bedtools.git
RUN cd $HOME/setup/bedtools && make && cp bin/* $HOME/bin

### install samtools ###
RUN cd $HOME/setup && git clone https://github.com/samtools/samtools.git
RUN cd $HOME/setup && git clone https://github.com/samtools/htslib.git
RUN cd $HOME/setup/samtools && make && cp samtools $HOME/bin

### R ###
RUN apt-get install -y r-base r-base-dev 
### install R packages ###
RUN mkdir -p $(Rscript -e 'cat(Sys.getenv("R_LIBS_USER"))')
#Bioconductor: GenomicRanges, ggbio, GenomeInfoDb, BSgenome.Hsapiens.UCSC.hg19, BSgenome.Hsapiens.UCSC.hg38, Biostrings
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("BSgenome.Hsapiens.UCSC.hg19",ask=F)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("BSgenome.Hsapiens.UCSC.hg38",ask=F)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("GenomicRanges",ask=F)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("GenomeInfoDb",ask=F)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("Biostrings",ask=F)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("Rsamtools",ask=F)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("biovizBase",ask=F)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("ensembldb",ask=F)'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("ggbio",ask=F)'
#CRAN: optparse, rjson, RColorBrewer, devtools, igraph, cowplot, ggplot2
RUN Rscript -e 'install.packages("optparse",repos="http://cran.us.r-project.org")'  
RUN Rscript -e 'install.packages("igraph",repos="http://cran.us.r-project.org")'  
RUN Rscript -e 'install.packages("rjson",repos="http://cran.us.r-project.org")'  
RUN Rscript -e 'install.packages("ggplot2",repos="http://cran.us.r-project.org")'  
RUN Rscript -e 'install.packages("RColorBrewer",repos="http://cran.us.r-project.org")'  
RUN Rscript -e 'install.packages("cowplot",repos="http://cran.us.r-project.org")'  
RUN Rscript -e 'install.packages("devtools",repos="http://cran.us.r-project.org")'  

### install zoomx ###
RUN cd $HOME/setup && git clone https://bitbucket.org/charade/zoomx.git
RUN cd $HOME/setup/zoomx && python setup.py install
RUN cd $HOME/setup/zoomx/test && ./test_dep.sh

### run test scripts ### 
RUN cd $HOME/setup/zoomx/test && wget https://s3-us-west-2.amazonaws.com/lixiabucket/zoomx.test.full.tgz -O zoomx.test.full.tgz && tar -zxvf zoomx.test.full.tgz
RUN export "PATH=$PATH:$HOME/bin" && cd $HOME/setup/zoomx/test && ./germline10.run
RUN export "PATH=$PATH:$HOME/bin" && cd $HOME/setup/zoomx/test && ./somatic10.run

### mount your server dir for data analysis
#docker run -d --name="foo" -v "/home/lixia:/home/lixia" ubuntu
