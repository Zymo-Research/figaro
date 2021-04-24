FROM python:3.6.7

MAINTAINER Michael M. Weinstein, Zymo Research
LABEL version="1.0.0"

WORKDIR /
RUN cd /opt &&\
    mkdir figaro

COPY . /opt/figaro

RUN cd opt/figaro &&\
    pip3 install --upgrade pip &&\
    pip3 install -r requirements.txt

ENV PYTHONUNBUFFERED=1

CMD ["python3", "/opt/figaro/figaro/figaro.py"]