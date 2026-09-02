version 1.0

##########################################################################################
##
## IGV screenshots from long-read BAMs: subsets each BAM to the plotted regions with
## `samtools view -b`, then generates IGV batch scripts (long-read mode: no `viewaspairs`,
## reference genome loaded via `--genome`) and captures snapshots under xvfb.
##
##########################################################################################

import "Structs2.wdl"

workflow IGV {
    input{
        File varfile
        String family
        File ped_file
        Array[String] samples
        Boolean file_localization
        Boolean requester_pays
        Boolean long_read
        Array[File]? bams_localize
        Array[File]? bais_localize
        String igv_genome
        Int igv_max_window
        String buffer
        String igv_docker
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_igv
        RuntimeAttr? runtime_attr_localize_reads
        Array[String]? bams_parse
        Array[String]? bais_parse
        File? updated_sample_bam_bai
        File? sample_bam_bai
        Array[File] annotation_beds = []
        Array[String] annotation_names = []
        File? gene_track
        File? gene_track_index
    }

    if (file_localization) {
        Array[File] bams_localize_ = select_first([bams_localize])
        Array[File] bais_localize_ = select_first([bais_localize])
        File sample_bam_bai_ = select_first([sample_bam_bai])

        if (requester_pays){
        # move the reads nearby -- handles requester_pays and makes cross-region transfers just once
            scatter(i in range(length(bams_localize_))) {
                call LocalizeReads as LocalizeReadsLocalize{
                    input:
                        reads_path = bams_localize_[i],
                        reads_index = bais_localize_[i],
                        runtime_attr_override = runtime_attr_localize_reads
                }
            }
        }

        call runIGV_whole_genome_localize {
            input:
                varfile = varfile,
                family = family,
                ped_file = ped_file,
                samples = samples,
                bams = select_first([LocalizeReadsLocalize.output_file, bams_localize_]),
                bais = select_first([LocalizeReadsLocalize.output_index, bais_localize_]),
                sample_bam_bai = sample_bam_bai_,
                buffer = buffer,
                igv_genome = igv_genome,
                igv_max_window = igv_max_window,
                long_read = long_read,
                annotation_beds = annotation_beds,
                annotation_names = annotation_names,
                gene_track = gene_track,
                gene_track_index = gene_track_index,
                igv_docker = igv_docker,
                runtime_attr_override = runtime_attr_igv
        }
    }

    if (!(file_localization)) {
        Array[String] bams_parse_ = select_first([bams_parse])
        Array[String] bais_parse_ = select_first([bais_parse])
        File updated_sample_bam_bai_ = select_first([updated_sample_bam_bai])

        call runIGV_whole_genome_parse{
            input:
                varfile = varfile,
                family = family,
                ped_file = ped_file,
                samples = samples,
                updated_sample_bam_bai = updated_sample_bam_bai_,
                buffer = buffer,
                igv_genome = igv_genome,
                igv_max_window = igv_max_window,
                long_read = long_read,
                annotation_beds = annotation_beds,
                annotation_names = annotation_names,
                gene_track = gene_track,
                gene_track_index = gene_track_index,
                igv_docker = igv_docker,
                runtime_attr_override = runtime_attr_igv
        }
    }

    output{
        File tar_gz_pe = select_first([runIGV_whole_genome_localize.pe_plots, runIGV_whole_genome_parse.pe_plots])
    }
}

task runIGV_whole_genome_localize{
        input{
            File varfile
            String igv_genome
            Int igv_max_window
            String family
            File ped_file
            Array[String] samples
            Array[File] bams
            Array[File] bais
            File sample_bam_bai
            String buffer
            Boolean long_read
            Array[File] annotation_beds = []
            Array[String] annotation_names = []
            File? gene_track
            File? gene_track_index
            String igv_docker
            RuntimeAttr? runtime_attr_override
        }

    Float input_size = size(select_all([varfile, ped_file, gene_track, gene_track_index]), "GB") + size(bams, "GB") + size(bais, "GB") + size(annotation_beds, "GB")
    # Peak IGV RAM is driven by how many BAM tracks are held at once (load-once) and the
    # capped read depth, NOT the variant count (variants render sequentially). Scale the
    # heap by sample count, with a hard ceiling so a large family can't run up cloud cost.
    # The JVM's -Xmx is set at runtime from the VM's actual RAM (see command block).
    Int n_samples = length(samples)
    Float base_mem_gb = 6.0
    Float mem_per_sample_gb = 2.0
    Float mem_ceiling_gb = 20.0
    Float mem_uncapped_gb = base_mem_gb + mem_per_sample_gb * n_samples
    Float dynamic_mem_gb = if mem_uncapped_gb < mem_ceiling_gb then mem_uncapped_gb else mem_ceiling_gb

    RuntimeAttr default_attr = object {
                                      mem_gb: dynamic_mem_gb,
                                      disk_gb: ceil(20 + input_size * 2),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
            set -euo pipefail
            mkdir pe_igv_plots
            head -n+1 ~{ped_file} > family_ped.txt
            grep -w ~{family} ~{ped_file} >> family_ped.txt
            python3 /src/variant-interpretation/scripts/renameCramsLocalize.py --ped family_ped.txt --scc ~{sample_bam_bai}
            cut -f4 changed_sample_crai_cram.txt > bams.txt

            while read sample bai bam new_bam new_bai
            do
                mv $bam $new_bam
                mv $bai $new_bai
            done<changed_sample_crai_cram.txt

            # gene track: prefer a user-supplied bgzipped + tabix-indexed file (IGV then
            # region-queries only the visible locus); otherwise index it once here
            GENES_ARG=""
            GENE_TRACK="~{default='' gene_track}"
            GENE_TBI="~{default='' gene_track_index}"
            if [ -n "$GENE_TRACK" ]; then
                if [ -n "$GENE_TBI" ]; then
                    # place the index beside the gzip unless it's already there (Cromwell may
                    # localize both to the same dir, where $GENE_TRACK.tbi == $GENE_TBI)
                    if [ "$GENE_TBI" != "$GENE_TRACK.tbi" ]; then
                        ln -sf "$GENE_TBI" "$GENE_TRACK.tbi"
                    fi
                    GENES_ARG="--genes $GENE_TRACK"
                else
                    low=$( echo "$GENE_TRACK" | tr '[:upper:]' '[:lower:]' )
                    case "$low" in
                        *.gtf|*.gtf.gz|*.gff|*.gff.gz|*.gff3|*.gff3.gz)
                            zcat -f "$GENE_TRACK" | grep -v '^#' | sort -k1,1 -k4,4n | bgzip > genes.idx.gtf.gz
                            tabix -p gff genes.idx.gtf.gz
                            GENES_ARG="--genes genes.idx.gtf.gz" ;;
                        *.bed|*.bed.gz)
                            zcat -f "$GENE_TRACK" | grep -v '^#' | sort -k1,1 -k2,2n | bgzip > genes.idx.bed.gz
                            tabix -p bed genes.idx.bed.gz
                            GENES_ARG="--genes genes.idx.bed.gz" ;;
                        *)
                            GENES_ARG="--genes $GENE_TRACK" ;;
                    esac
                fi
            fi

            # one IGV batch (and one JVM) for all of this family's variants
            python /src/variant-interpretation/scripts/makeigvpesr.py -v ~{varfile} -fam_id ~{family} -samples ~{sep="," samples} -crams bams.txt -p ~{ped_file} -o pe_igv_plots -b ~{buffer} -i pe.all.txt -bam pe.all.sh -m ~{igv_max_window} --genome ~{igv_genome} ~{true="--long_read" false="" long_read} $GENES_ARG --annotation_beds ~{sep=" " annotation_beds} --annotation_names ~{sep=" " annotation_names}
            bash pe.all.sh

            # Size the IGV JVM heap to this VM's actual RAM (mem_gb is set per-family in the
            # runtime block). Leave ~20% headroom for xvfb/OS/samtools.
            MEM_KB=$(awk '/MemTotal/{print $2}' /proc/meminfo)
            XMX_MB=$(( MEM_KB * 80 / 100 / 1024 ))
            sed -i -E "s/-Xmx[0-9]+[mMgG]/-Xmx${XMX_MB}m/g" /IGV_Linux_2.16.0/igv.sh
            # Terminate the JVM immediately on OutOfMemoryError; otherwise IGV catches it
            # internally and the process never exits.
            export _JAVA_OPTIONS="-XX:+ExitOnOutOfMemoryError ${_JAVA_OPTIONS:-}"

            # Retry the batch: the hg38 reference is fetched from igv.org over HTTP and can
            # fail transiently. set -e does not trip on a failure in the until-condition.
            n=0
            until xvfb-run --server-args="-screen 0, 1920x1080x24" bash /IGV_Linux_2.16.0/igv.sh -g ~{igv_genome} -b pe.all.txt; do
                n=$((n+1))
                if [ $n -ge 3 ]; then echo "IGV batch failed after $n attempts" >&2; exit 1; fi
                echo "IGV batch attempt $n failed; retrying in $((n*30))s..." >&2
                sleep $((n*30))
            done
            tar -czf ~{family}_pe_igv_plots.tar.gz pe_igv_plots

        >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: igv_docker
    }
    output{
        File pe_plots="~{family}_pe_igv_plots.tar.gz"
        Array[File] pe_txt = glob("pe.*.txt")
        Array[File] pe_sh = glob("pe.*.sh")
        Array[File] new_varfiles = glob("new.varfile.*.bed")
        }
    }

task runIGV_whole_genome_parse{
    input{
        File varfile
        String igv_genome
        Int igv_max_window
        String family
        File ped_file
        Array[String] samples
        File updated_sample_bam_bai
        String buffer
        Boolean long_read
        Array[File] annotation_beds = []
        Array[String] annotation_names = []
        File? gene_track
        File? gene_track_index
        String igv_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([varfile, ped_file, gene_track, gene_track_index]), "GB") + size(annotation_beds, "GB")
    # Peak IGV RAM is driven by how many BAM tracks are held at once (load-once) and the
    # capped read depth, NOT the variant count (variants render sequentially). Scale the
    # heap by sample count, with a hard ceiling so a large family can't run up cloud cost.
    # The JVM's -Xmx is set at runtime from the VM's actual RAM (see command block).
    Int n_samples = length(samples)
    Float base_mem_gb = 6.0
    Float mem_per_sample_gb = 2.0
    Float mem_ceiling_gb = 20.0
    Float mem_uncapped_gb = base_mem_gb + mem_per_sample_gb * n_samples
    Float dynamic_mem_gb = if mem_uncapped_gb < mem_ceiling_gb then mem_uncapped_gb else mem_ceiling_gb

    RuntimeAttr default_attr = object {
                                      mem_gb: dynamic_mem_gb,
                                      disk_gb: ceil(20 + input_size * 2),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
            set -euo pipefail
            mkdir pe_igv_plots
            cat ~{varfile} | cut -f1-3 | awk '{if (($3-$2)+int(($3-$2)*1.5)>=~{igv_max_window}) print $1"\t"$2-~{buffer}"\t"$2+~{buffer} "\n" $1"\t"$3-~{buffer}"\t"$3+~{buffer};
                else print $1"\t"($2-int(($3-$2)*0.25))-~{buffer}"\t"$3+int(($3-$2)*0.25)+~{buffer}}' | sort -k1,1 -k2,2n | bgzip -c > regions.bed.gz
            tabix -p bed regions.bed.gz
            # OAuth token from the GCE metadata server (no gcloud SDK needed); htslib
            # reads gs:// BAMs via libcurl using GCS_OAUTH_TOKEN
            export GCS_OAUTH_TOKEN=$(curl -s -H "Metadata-Flavor: Google" \
                "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token" \
                | python3 -c "import sys,json;print(json.load(sys.stdin)['access_token'])")
            # Bill a project for requester-pays buckets: htslib sends it as X-Goog-User-Project,
            # and we add the same header to gcs_cp. Use this VM's compute project; ignored for
            # non-requester-pays buckets.
            export GCS_REQUESTER_PAYS_PROJECT=$(curl -s -H "Metadata-Flavor: Google" \
                "http://metadata.google.internal/computeMetadata/v1/project/project-id")
            # gs://bucket/obj -> local file, via the GCS XML API with a bearer token
            gcs_cp () {
                p="${1#gs://}"
                curl -sf -H "Authorization: Bearer $GCS_OAUTH_TOKEN" \
                    -H "X-Goog-User-Project: $GCS_REQUESTER_PAYS_PROJECT" \
                    -o "$2" "https://storage.googleapis.com/${p}"
            }
            # subset each remote BAM to the plotted regions
            while read sample bai bam new_bam new_bai
            do
                gcs_cp "$bai" "$( basename $bam | sed 's/\.bam$/.bai/g' )"
                # name the subset by sample id so makeigvpesr maps each track to its sample
                samtools view -h -b -o $sample.bam $bam -L regions.bed.gz -M
                samtools index $sample.bam
            done<~{updated_sample_bam_bai}
            ls *.bam > bams.txt

            # gene track: prefer a user-supplied bgzipped + tabix-indexed file (IGV then
            # region-queries only the visible locus); otherwise index it once here
            GENES_ARG=""
            GENE_TRACK="~{default='' gene_track}"
            GENE_TBI="~{default='' gene_track_index}"
            if [ -n "$GENE_TRACK" ]; then
                if [ -n "$GENE_TBI" ]; then
                    # place the index beside the gzip unless it's already there (Cromwell may
                    # localize both to the same dir, where $GENE_TRACK.tbi == $GENE_TBI)
                    if [ "$GENE_TBI" != "$GENE_TRACK.tbi" ]; then
                        ln -sf "$GENE_TBI" "$GENE_TRACK.tbi"
                    fi
                    GENES_ARG="--genes $GENE_TRACK"
                else
                    low=$( echo "$GENE_TRACK" | tr '[:upper:]' '[:lower:]' )
                    case "$low" in
                        *.gtf|*.gtf.gz|*.gff|*.gff.gz|*.gff3|*.gff3.gz)
                            zcat -f "$GENE_TRACK" | grep -v '^#' | sort -k1,1 -k4,4n | bgzip > genes.idx.gtf.gz
                            tabix -p gff genes.idx.gtf.gz
                            GENES_ARG="--genes genes.idx.gtf.gz" ;;
                        *.bed|*.bed.gz)
                            zcat -f "$GENE_TRACK" | grep -v '^#' | sort -k1,1 -k2,2n | bgzip > genes.idx.bed.gz
                            tabix -p bed genes.idx.bed.gz
                            GENES_ARG="--genes genes.idx.bed.gz" ;;
                        *)
                            GENES_ARG="--genes $GENE_TRACK" ;;
                    esac
                fi
            fi

            # one IGV batch (and one JVM) for all of this family's variants
            python /src/variant-interpretation/scripts/makeigvpesr.py -v ~{varfile} -fam_id ~{family} -samples ~{sep="," samples} -crams bams.txt -p ~{ped_file} -o pe_igv_plots -b ~{buffer} -i pe.all.txt -bam pe.all.sh -m ~{igv_max_window} --genome ~{igv_genome} ~{true="--long_read" false="" long_read} --status_labels $GENES_ARG --annotation_beds ~{sep=" " annotation_beds} --annotation_names ~{sep=" " annotation_names}
            bash pe.all.sh

            # Size the IGV JVM heap to this VM's actual RAM (mem_gb is set per-family in the
            # runtime block). Leave ~20% headroom for xvfb/OS/samtools.
            MEM_KB=$(awk '/MemTotal/{print $2}' /proc/meminfo)
            XMX_MB=$(( MEM_KB * 80 / 100 / 1024 ))
            sed -i -E "s/-Xmx[0-9]+[mMgG]/-Xmx${XMX_MB}m/g" /IGV_Linux_2.16.0/igv.sh
            # Terminate the JVM immediately on OutOfMemoryError; otherwise IGV catches it
            # internally and the process never exits.
            export _JAVA_OPTIONS="-XX:+ExitOnOutOfMemoryError ${_JAVA_OPTIONS:-}"

            # Retry the batch: the hg38 reference is fetched from igv.org over HTTP and can
            # fail transiently. set -e does not trip on a failure in the until-condition.
            n=0
            until xvfb-run --server-args="-screen 0, 1920x1080x24" bash /IGV_Linux_2.16.0/igv.sh -g ~{igv_genome} -b pe.all.txt; do
                n=$((n+1))
                if [ $n -ge 3 ]; then echo "IGV batch failed after $n attempts" >&2; exit 1; fi
                echo "IGV batch attempt $n failed; retrying in $((n*30))s..." >&2
                sleep $((n*30))
            done
            tar -czf ~{family}_pe_igv_plots.tar.gz pe_igv_plots

        >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: igv_docker
    }
    output{
        File pe_plots="~{family}_pe_igv_plots.tar.gz"
        Array[File] pe_txt = glob("pe.*.txt")
        Array[File] pe_sh = glob("pe.*.sh")
        Array[File] new_varfiles = glob("new.varfile.*.bed")
        Array[File] bamfiles = glob("*bam")
        }
    }

task LocalizeReads {
  input {
    File reads_path
    File reads_index
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(60.0),
                                  cpu: 2,
                                  preemptible: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }

  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu, runtime_default.cpu])
    preemptible: select_first([runtime_override.preemptible, runtime_default.preemptible])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: "ubuntu:18.04"
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command {
      set -exuo pipefail

      cp ~{reads_path} $(basename ~{reads_path})
      cp ~{reads_index} $(basename ~{reads_index})
  }
  output {
    File output_file = basename(reads_path)
    File output_index = basename(reads_index)
  }
}
