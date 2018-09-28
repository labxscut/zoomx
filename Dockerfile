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
FROM charade/xlibbox:zero

MAINTAINER Charlie Xia <lxia.usc@gmail.com>

WORKDIR $HOME/
ENV PATH "PATH=$PATH:$HOME/bin"
RUN ulimit -s unlimited

### install zoomx ###
RUN cd $HOME/setup && git lfs clone https://bitbucket.org/charade/zoomx.git
RUN cd $HOME/setup/zoomx && python setup.py install
RUN cd $HOME/setup/zoomx/test && ./test_dep.sh

### run test scripts ### 
RUN cd $HOME/setup/zoomx/test && wget https://s3-us-west-2.amazonaws.com/lixiabucket/zoomx.test.full.tgz -O zoomx.test.full.tgz && tar -zxvf zoomx.test.full.tgz
RUN export "PATH=$PATH:$HOME/bin" && cd $HOME/setup/zoomx/test && ./germline10.run
RUN export "PATH=$PATH:$HOME/bin" && cd $HOME/setup/zoomx/test && ./somatic10.run

### mount your server dir for data analysis
#docker run -d --name="foo" -v "/home/lixia:/home/lixia" ubuntu
