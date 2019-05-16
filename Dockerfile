FROM ubuntu:bionic

ENV TZ=Europe/Minsk

RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get clean all && apt-get update && apt-get install -y -q git wget perl \
    software-properties-common zlibc zlib1g-dev python3-pip python3-matplotlib 

WORKDIR /opt/

RUN git clone  https://github.com/DrQuestion/monica.git

RUN pip3 install -r ./monica/requirements.txt

RUN ln -s /opt/monica/monica/monica.py /usr/local/bin/monica

RUN apt-get update && apt-get install -y python3-tk



