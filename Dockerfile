FROM ubuntu:19.04

COPY . /app
WORKDIR /app

ENV LANG en_GB.UTF-8

# default bokeh server port
EXPOSE 5006/tcp

# use mirrors for apt, much faster than default
RUN sed -i -e 's/http:\/\/archive.ubuntu.com\/ubuntu/mirror:\/\/mirrors.ubuntu.com\/mirrors.txt/' /etc/apt/sources.list

# install python3 and nplinker deps
RUN apt update 
RUN apt install -y python3 python3-pip
RUN pip3 install -v biopython jupyter numpy pandas matplotlib tornado bokeh pytz pyparsing python-dateutil 
RUN pip3 install -v terminado wcwidth seaborn scipy sklearn psutil toml networkx xdg

# can remove last bit if using port 5006 externally too
CMD cd /app && bokeh serve webapp/npapp --dev --allow-websocket-origin=localhost:5010 
