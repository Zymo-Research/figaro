FROM python:3.6.7

WORKDIR /
RUN pip3 install --upgrade pip &&\
    pip3 install matplotlib==3.0.2 &&\
    pip3 install numpy==1.13.1 &&\
    pip3 install scipy==1.2.1 &&\
    cd /opt &&\
    mkdir trimParameterFinder16s

COPY . /opt/trimParameterFinder16s

CMD ["python3", "/opt/trimParameterFinder16s/analyzeFastqFolder.py"]