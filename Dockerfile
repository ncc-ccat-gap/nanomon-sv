FROM python:3.7.3
MAINTAINER aokada <aokada@ncc.go.jp>

# RUN: build
WORKDIR /usr/local

RUN apt-get -y update && \
    apt-get install -y wget unzip zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev

RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar jxvf samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    ./configure && make && make install

RUN pip install pysam==0.15.2

RUN wget https://github.com/ncc-ccat-gap/nanomon-sv/archive/v0.0.1.zip && \
    unzip v0.0.1.zip && \
    rm -f v0.0.1.zip

CMD ["/bin/bash"]
