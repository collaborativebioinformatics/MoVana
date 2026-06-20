# ubuntu 24.04 base image
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=UTC \
    APP_HOME=/opt/MoVana \
    PATH=/opt/venv/bin:$PATH \
    UV_PYTHON_INSTALL_DIR=/opt/python \
    UV_LINK_MODE=copy

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        zlib1g-dev \
        wget \
        ca-certificates \
	bedtools \
        bcftools \
        openjdk-17-jre-headless&& \
    rm -rf /var/lib/apt/lists/*

# uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /usr/local/bin/

RUN uv python install 3.12
COPY requirements.txt /tmp/requirements.txt
RUN uv venv /opt/venv --python 3.12 && \
    uv pip install --python /opt/venv/bin/python -r /tmp/requirements.txt

# cromwell
ARG CROMWELL_VERSION=90
ARG CROMWELL_SHA256=d90e46f60f430ff627222c97b950c43f1ededc992619e0aeceaa334690d06073
RUN wget -q -O $APP_HOME/cromwell/cromwell.jar \
        "https://github.com/broadinstitute/cromwell/releases/download/${CROMWELL_VERSION}/cromwell-${CROMWELL_VERSION}.jar" && \
    echo "${CROMWELL_SHA256}  /tmp/cromwell.jar" | sha256sum -c -

# working dir
WORKDIR $APP_HOME

# Run as non-root - movana
RUN useradd --create-home --shell /bin/bash movana && \
    chown -R movana:movana $APP_HOME
USER movana

COPY --chown=movana:movana . $APP_HOME

# bash
CMD ["bash"]
