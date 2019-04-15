MAINTAINER Michael M. Weinstein, Zymo Research
LABEL version="1.0.0"

FROM python:3.6.7

WORKDIR /
RUN cd /opt &&\
    mkdir trimParameterFinder16s &&\

COPY . /opt/trimParameterFinder16s

RUN pip3 install --upgrade pip &&\
    pip3 install -r requirements.txt

CMD ["python3", "/opt/trimParameterFinder16s/analyzeFastqFolder.py"]