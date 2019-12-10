FROM ubuntu:19.04

COPY . /app
WORKDIR /app

ENV LANG en_GB.UTF-8

# default bokeh server port
EXPOSE 5006/tcp

# use mirrors for apt, much faster than default
RUN sed -i -e 's/http:\/\/archive.ubuntu.com\/ubuntu/mirror:\/\/mirrors.ubuntu.com\/mirrors.txt/' /etc/apt/sources.list

# install python3 and nplinker deps
RUN apt-get update 
RUN apt-get install -y python3 python3-pip
RUN pip3 install biopython numpy pandas tornado pytz pyparsing python-dateutil 
RUN pip3 install scipy sklearn psutil toml networkx xdg
RUN pip3 install -v /app/bokeh-1.3.5.dev3+5.ge5c1c99e1.dirty-py3-none-any.whl

# can remove last bit if using port 5006 externally too
CMD cd /app && bokeh serve webapp/npapp --dev --allow-websocket-origin=localhost:5010 
