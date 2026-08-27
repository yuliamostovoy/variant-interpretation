# long_read_visualize_docker
#
# Self-contained image for ALL script/compute tasks of the long-read visualization
# workflow (RunVisualizePlotsLongRead.wdl), so the pipeline depends only on this image
# plus a stock gatk-sv sv_base_mini image:
#   - ReformatVariants.reformat_variants        : python3 + bgzip (reformat_variants_for_visualization.py)
#   - CreateIgvBamPlots.update_sample_bam_bai    : python3 + pandas (renameCrams.py)
#   - IgvBamPlots.runIGV_*                        : IGV + xvfb + python (makeigvpesr.py)
#   - LongReadDepthPlot.depth_plot               : samtools + bedtools + mosdepth + matplotlib
#   - RunVisualizePlotsLongRead.concat_plots     : python + Pillow
#
# gcloud/gsutil are present so tasks can auth to and stream gs:// BAMs
# (`export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)`),
# matching the pattern used by the existing GATK-SV rdtest image.
#
# Scripts land at /src/variant-interpretation/scripts/ to match the WDL task paths.

FROM google/cloud-sdk:slim

ARG MOSDEPTH_VERSION=0.3.14
ARG IGV_VERSION=2.16.0

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        bedtools \
        samtools \
        tabix \
        python3 \
        python3-pip \
        python-is-python3 \
        git \
        ca-certificates \
        curl \
        unzip \
        xvfb \
        fontconfig \
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
# Extracts to /IGV_Linux_${IGV_VERSION}/, matching the igv.sh path the WDL invokes.
RUN curl -sSL -o /tmp/igv.zip \
        "https://data.broadinstitute.org/igv/projects/downloads/2.16/IGV_Linux_${IGV_VERSION}_WithJava.zip" \
    && unzip -q /tmp/igv.zip -d / \
    && rm /tmp/igv.zip

# Python deps: plotting (matplotlib + Pillow) and dataframes (pandas for renameCrams.py /
# makeigvpesr.py, which also uses numpy)
RUN pip3 install --no-cache-dir --break-system-packages \
        "matplotlib>=3.5" \
        "pillow>=9" \
        "numpy>=1.21" \
        "pandas>=1.3"

# Repo scripts at the conventional path
RUN git clone --depth 1 --branch main \
        https://github.com/yuliamostovoy/variant-interpretation.git \
        /src/variant-interpretation

ENV MPLBACKEND=Agg

WORKDIR /src/variant-interpretation
