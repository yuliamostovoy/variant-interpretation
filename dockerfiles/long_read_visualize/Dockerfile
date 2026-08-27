# long_read_visualize_docker
#
# Image for the long-read visualization workflow's compute-heavy tasks:
#   - LongReadDepthPlot.depth_plot : samtools (subset remote BAM), bedtools (merge/windows),
#                                    mosdepth (per-window depth), python+matplotlib (plot)
#   - RunVisualizePlotsLongRead.concat_plots : python + Pillow (stack IGV over depth)
#
# gcloud/gsutil are present so tasks can auth to and stream gs:// BAMs
# (`export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)`),
# matching the pattern used by the existing GATK-SV rdtest image.

FROM google/cloud-sdk:slim

ARG MOSDEPTH_VERSION=0.3.14

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        bedtools \
        samtools \
        tabix \
        python3 \
        python3-pip \
        git \
        ca-certificates \
        curl \
    && rm -rf /var/lib/apt/lists/*

# mosdepth: use the official static Linux binary (pinned)
RUN curl -sSL -o /usr/local/bin/mosdepth \
        "https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth" \
    && chmod +x /usr/local/bin/mosdepth

# Python plotting deps (headless matplotlib + Pillow for image stacking)
RUN pip3 install --no-cache-dir --break-system-packages \
        "matplotlib>=3.5" \
        "pillow>=9" \
        "numpy>=1.21"

# Repo scripts at the conventional path
RUN git clone --depth 1 --branch main \
        https://github.com/yuliamostovoy/variant-interpretation.git \
        /src/variant-interpretation

ENV MPLBACKEND=Agg

WORKDIR /src/variant-interpretation
