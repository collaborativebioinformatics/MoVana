# Python base image [pre-built ubuntu-base python env]
FROM python:3.12-slim

# reduce package overhead
ENV DEBIAN_FRONTEND=noninteractive \
    TZ=UTC \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    APP_HOME=/opt/MoVana

# install dependencies
RUN apt-get update &&  apt-get install -y --no-install-recommends build-essential \
        wget \
        git \
        bedtools \
        bcftools \
        openjdk-17-jdk \
        gfortran \
        zlib1g-dev && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# workspace
RUN mkdir -p $APP_HOME

# Setup Working directory`
WORKDIR $APP_HOME

# Copy the WDL script into the container
COPY . $APP_HOME

# Cromwell (make sure Java is installed)
RUN mkdir -p cromwell && wget https://github.com/broadinstitute/cromwell/releases/download/90/cromwell-90.jar && mv cromwell-90.jar cromwell/

# python pkgs
#RUN python3 -m venv 0env && source 0env/bin/activate && \
RUN python3 -m pip install --upgrade pip setuptools wheel && python3 -m pip install --no-cache-dir pandas numpy matplotlib

# Set login shell
CMD ["/bin/bash -c"]
