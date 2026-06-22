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
COPY --from=ghcr.io/astral-sh/uv:0.9.7 /uv /uvx /usr/local/bin/

RUN uv python install 3.12
COPY requirements.txt /tmp/requirements.txt
RUN uv venv /opt/venv --python 3.12 && \
    uv pip install --python /opt/venv/bin/python -r /tmp/requirements.txt

# working dir
WORKDIR $APP_HOME

# cromwell
ARG CROMWELL_VERSION=92
ARG CROMWELL_SHA256=e0e3a050d4124e81369a79059e5774142b2f06bd89df4a0b035f559db85cedf5
RUN wget -q -O $APP_HOME/cromwell.jar https://github.com/broadinstitute/cromwell/releases/download/${CROMWELL_VERSION}/cromwell-${CROMWELL_VERSION}.jar && \
    echo "${CROMWELL_SHA256}  $APP_HOME/cromwell.jar" | sha256sum -c -

# Run as non-root - movana
RUN useradd --create-home --shell /bin/bash movana && \
    chown -R movana:movana $APP_HOME
USER movana

COPY --chown=movana:movana . $APP_HOME

# bash
CMD ["bash"]
