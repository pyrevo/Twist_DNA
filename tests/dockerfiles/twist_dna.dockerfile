FROM python:3.8

RUN mkdir /Twist_DNA

COPY * /Twist_DNA/

WORKDIR /Twist_DNA

RUN pip install -r requirements.txt
