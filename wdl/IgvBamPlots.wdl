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
        File reference
        File reference_index
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
                reference = reference,
                reference_index = reference_index,
                igv_max_window = igv_max_window,
                long_read = long_read,
                annotation_beds = annotation_beds,
                annotation_names = annotation_names,
                gene_track = gene_track,
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
                reference = reference,
                reference_index = reference_index,
                igv_max_window = igv_max_window,
                long_read = long_read,
                annotation_beds = annotation_beds,
                annotation_names = annotation_names,
                gene_track = gene_track,
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
            File reference
            File reference_index
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
            String igv_docker
            RuntimeAttr? runtime_attr_override
        }

    Float input_size = size(select_all([varfile, ped_file, bams, bais]), "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size),
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

            i=0
            while read -r line
            do
                let "i=$i+1"
                echo "$line" > new.varfile.$i.bed
                python /src/variant-interpretation/scripts/makeigvpesr.py -v new.varfile.$i.bed -fam_id ~{family} -samples ~{sep="," samples} -crams bams.txt -p ~{ped_file} -o pe_igv_plots -b ~{buffer} -i pe.$i.txt -bam pe.$i.sh -m ~{igv_max_window} --genome ~{reference} ~{true="--long_read" false="" long_read} ~{"--genes " + gene_track} --annotation_beds ~{sep=" " annotation_beds} --annotation_names ~{sep=" " annotation_names}
                bash pe.$i.sh
                xvfb-run --server-args="-screen 0, 1920x1080x24" bash /IGV_Linux_2.16.0/igv.sh -b pe.$i.txt
            done < ~{varfile}
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
        File reference
        File reference_index
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
        String igv_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([varfile, ped_file]), "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size),
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
            # subset each remote BAM to the plotted regions
            while read sample bai bam new_bam new_bai
            do
                gsutil cp $bai $( basename $bam | sed 's/\.bam$/.bai/g' )
                export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
                # name the subset by sample id so makeigvpesr maps each track to its sample
                samtools view -h -b -o $sample.bam $bam -L regions.bed.gz -M
                samtools index $sample.bam
            done<~{updated_sample_bam_bai}
            ls *.bam > bams.txt

            i=0
            while read -r line
            do
                let "i=$i+1"
                echo "$line" > new.varfile.$i.bed
                python /src/variant-interpretation/scripts/makeigvpesr.py -v new.varfile.$i.bed -fam_id ~{family} -samples ~{sep="," samples} -crams bams.txt -p ~{ped_file} -o pe_igv_plots -b ~{buffer} -i pe.$i.txt -bam pe.$i.sh -m ~{igv_max_window} --genome ~{reference} ~{true="--long_read" false="" long_read} --status_labels ~{"--genes " + gene_track} --annotation_beds ~{sep=" " annotation_beds} --annotation_names ~{sep=" " annotation_names}
                bash pe.$i.sh
                xvfb-run --server-args="-screen 0, 1920x1080x24" bash /IGV_Linux_2.16.0/igv.sh -b pe.$i.txt
            done < ~{varfile}
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
