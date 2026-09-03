FROM ubuntu:22.04

ARG MOSDEPTH_VERSION=0.3.14
ARG IGV_VERSION=2.16.0

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        samtools \
        bcftools \
        bedtools \
        tabix \
        python3 \
        python3-pip \
        python-is-python3 \
        git \
        ca-certificates \
        curl \
        unzip \
        xvfb \
        xauth \
        fontconfig \
        fonts-dejavu-core \
        libxrender1 \
        libxtst6 \
        libxi6 \
        libxext6 \
    && rm -rf /var/lib/apt/lists/*

# mosdepth: use the official static Linux binary (pinned)
RUN curl -sSL -o /usr/local/bin/mosdepth \
        "https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth" \
    && chmod +x /usr/local/bin/mosdepth

# IGV (WithJava bundle: ships its own JDK, so no separate Java needed).
# Extracts to /IGV_Linux_${IGV_VERSION}/
RUN curl -sSL -o /tmp/igv.zip \
        "https://data.broadinstitute.org/igv/projects/downloads/2.16/IGV_Linux_${IGV_VERSION}_WithJava.zip" \
    && unzip -q /tmp/igv.zip -d / \
    && rm /tmp/igv.zip

# Default IGV genome to hg38 so `new`/startup never loads the built-in hg19 default
RUN mkdir -p /root/igv && printf 'DEFAULT_GENOME_KEY=hg38\n' > /root/igv/prefs.properties

RUN pip3 install --no-cache-dir \
        "matplotlib>=3.5" \
        "pillow>=9" \
        "numpy>=1.21" \
        "pandas>=1.3"

ARG VI_COMMIT=868f90801760c17b874ada8af9961ca8afef2ad5
RUN git clone https://github.com/yuliamostovoy/variant-interpretation.git \
        /src/variant-interpretation \
    && git -C /src/variant-interpretation checkout ${VI_COMMIT}

# Point HTS_PATH to the correct location so samtools can stream gs:// BAMs
ENV HTS_PATH=/usr/lib/x86_64-linux-gnu/htslib

ENV MPLBACKEND=Agg

WORKDIR /src/variant-interpretation
