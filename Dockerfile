FROM ubuntu:xenial

RUN apt-get clean all && apt-get update && apt-get install -y -q git wget perl \
    python3.7 software-properties-common python3-pip zlibc zlib1g-dev

WORKDIR /opt/

RUN git clone  https://github.com/DrQuestion/monica.git

RUN pip3 install -r ./monica/requirements.txt

RUN ln -s /opt/monica/monica/monica.py /usr/local/bin/monica



