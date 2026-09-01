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
# gs:// BAMs are streamed directly by htslib/samtools (built with libcurl) using an OAuth
# token the tasks fetch from the GCE metadata server -- no gcloud/gsutil SDK required.
#
# Scripts land at /src/variant-interpretation/scripts/ to match the WDL task paths.

FROM ubuntu:22.04

ARG MOSDEPTH_VERSION=0.3.14
ARG IGV_VERSION=2.16.0

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        samtools \
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
# Extracts to /IGV_Linux_${IGV_VERSION}/, matching the igv.sh path the WDL invokes.
RUN curl -sSL -o /tmp/igv.zip \
        "https://data.broadinstitute.org/igv/projects/downloads/2.16/IGV_Linux_${IGV_VERSION}_WithJava.zip" \
    && unzip -q /tmp/igv.zip -d / \
    && rm /tmp/igv.zip

# Default IGV genome to hg38 so `new`/startup never loads the built-in hg19 default
RUN mkdir -p /root/igv && printf 'DEFAULT_GENOME_KEY=hg38\n' > /root/igv/prefs.properties

# Python deps: plotting (matplotlib + Pillow) and dataframes (pandas for renameCrams.py /
# makeigvpesr.py, which also uses numpy)
RUN pip3 install --no-cache-dir \
        "matplotlib>=3.5" \
        "pillow>=9" \
        "numpy>=1.21" \
        "pandas>=1.3"

# Repo scripts at the conventional path. Pin a commit so bumping VI_COMMIT busts this
# layer's cache — a plain clone of main is cached by command string and would rebuild stale.
ARG VI_COMMIT=1b3305b2a0c6e7ce022149f32d9cd73113dca91e
RUN git clone https://github.com/yuliamostovoy/variant-interpretation.git \
        /src/variant-interpretation \
    && git -C /src/variant-interpretation checkout ${VI_COMMIT}

# Debian/Ubuntu ships the htslib URL plugins (hfile_gcs, hfile_libcurl, ...) in the
# multiarch dir, but the packaged libhts's compiled-in plugin path points at /usr/local.
# Point HTS_PATH at the real location so samtools can stream gs:// BAMs (otherwise every
# remote open fails with "Protocol not supported").
ENV HTS_PATH=/usr/lib/x86_64-linux-gnu/htslib

ENV MPLBACKEND=Agg

WORKDIR /src/variant-interpretation
