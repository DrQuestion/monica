FROM ubuntu:xenial

RUN apt-get clean all && apt-get update && apt-get install -y -q git wget perl \
    software-properties-common zlibc zlib1g-dev

RUN add-apt-repository -y ppa:deadsnakes/ppa

RUN apt-get update && apt install python3.7 python3-pip

WORKDIR /opt/

RUN git clone  https://github.com/DrQuestion/monica.git

RUN pip3 install -r ./monica/requirements.txt

RUN ln -s /opt/monica/monica/monica.py /usr/local/bin/monica



