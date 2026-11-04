ARG UBUNTU_VERSION
FROM us-central1-docker.pkg.dev/dnastack-asap-parkinsons/workflow-images/util:1.2.0 as scripts

FROM ubuntu:${UBUNTU_VERSION}
ENV UBUNTU_VERSION "${UBUNTU_VERSION}"

LABEL MAINTAINER="DNAstack <support@dnastack.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -qq update \
	&& apt-get -qq install \
		build-essential \
		wget \
		curl \
		git \
		time \
		xxd \
		zlib1g-dev \
		libffi-dev \
		libssl-dev \
		pkg-config \
		libbz2-dev \
		liblzma-dev \
		libreadline-dev \
		libsqlite3-dev \
		libncursesw5-dev \
		tk-dev \
		uuid-dev \
	&& rm -rf /var/lib/apt/lists/*

# gcloud sdk; needed to upload output files
ARG GCLOUD_CLI_VERSION
ENV GCLOUD_CLI_VERSION "${GCLOUD_CLI_VERSION}"
RUN wget "https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz" \
	&& tar -zxvf "google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz" --directory /opt \
	&& rm "google-cloud-cli-${GCLOUD_CLI_VERSION}-linux-x86_64.tar.gz"

ENV PATH "${PATH}:/opt/google-cloud-sdk/bin"

ARG PYTHON3_VERSION
ENV PYTHON3_VERSION "${PYTHON3_VERSION}"
RUN curl -O https://www.python.org/ftp/python/${PYTHON3_VERSION}/Python-${PYTHON3_VERSION}.tar.xz \
	&& tar -xvf Python-${PYTHON3_VERSION}.tar.xz --directory /opt/ \
	&& rm Python-${PYTHON3_VERSION}.tar.xz
RUN cd /opt/Python-${PYTHON3_VERSION} \
	&& ./configure \
	&& make \
	&& make altinstall

RUN update-alternatives --install /usr/bin/python3 python3 /opt/Python-${PYTHON3_VERSION}/python 1 \
	&& update-alternatives --set python3 /opt/Python-${PYTHON3_VERSION}/python

# Build context: repo root (docker build -f docker/cell2location_tools/Dockerfile .)
COPY docker/cell2location_tools/requirements.txt /opt/requirements.txt
RUN python3 -m pip cache purge
RUN python3 -m pip install -r /opt/requirements.txt

COPY --from=scripts /opt/scripts /opt/scripts
COPY docker/cell2location_tools/scripts /opt/scripts
COPY utils /opt/scripts/utils
ENV PATH "${PATH}:/opt/scripts:/opt/scripts/main"
ENV PYTHONPATH /opt/scripts
