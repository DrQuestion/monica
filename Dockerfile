FROM ubuntu:bionic

RUN apt-get clean all && apt-get update && apt-get install -y -q git wget perl \
    software-properties-common zlibc zlib1g-dev python3-pip

WORKDIR /opt/

RUN git clone  https://github.com/DrQuestion/monica.git

RUN pip3 install -r ./monica/requirements.txt

RUN ln -s /opt/monica/monica/monica.py /usr/local/bin/monica



