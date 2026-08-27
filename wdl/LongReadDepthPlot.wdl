version 1.0

##########################################################################################
##
## Component 4: long-read read-depth track. Replaces RdVisualization.wdl (GATK-SV RdTest).
##
## Instead of GATK-SV batch bincov matrices / medianfiles, depth is computed on demand with
## mosdepth over the plotted CNV loci from each sample's BAM. Normalization is local
## (per-sample median over flanking windows), so no batch/median/outlier inputs are needed.
##
## Output shape (per-variant PNGs, tarred) matches the RD track so the downstream
## integrate + concatenate steps are reused unchanged.
##
##########################################################################################

import "Structs2.wdl"

workflow LongReadDepthPlot {
    input {
        String prefix
        File bed                 # canonical 6-col varfile: chrom,start,end,ID,svtype,samples (bgzipped, header)
        File ped_file
        File? fam_ids
        File sample_bam_bai      # sample <tab> bai <tab> bam
        Int? flank
        Int? window
        String sv_base_mini_docker
        String long_read_visualize_docker
        RuntimeAttr? runtime_attr_depth
        RuntimeAttr? runtime_attr_create_bed
    }

    Int flank_ = select_first([flank, 5000])
    Int window_ = select_first([window, 250])

    if (defined(fam_ids)) {
        File fam_ids_ = select_first([fam_ids])
        Array[String] family_ids = transpose(read_tsv(fam_ids_))[0]
    }
    if (!(defined(fam_ids))) {
        call generate_families{
            input:
                bed = bed,
                ped_file = ped_file,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_create_bed
        }
    }

    scatter (family in select_first([family_ids, generate_families.families])){
        call generate_per_family_bed{
            input:
                bed = bed,
                family = family,
                ped_file = ped_file,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_create_bed
        }

        call depth_plot{
            input:
                family = family,
                per_family_bed = generate_per_family_bed.bed_file,
                ped_file = ped_file,
                sample_bam_bai = sample_bam_bai,
                flank = flank_,
                window = window_,
                prefix = prefix,
                long_read_visualize_docker = long_read_visualize_docker,
                runtime_attr_override = runtime_attr_depth
        }
    }

    call integrate_depth_plots{
        input:
            depth_tar = depth_plot.plots,
            prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_create_bed
    }

    output{
        File Plots = integrate_depth_plots.plot_tar
    }
}

task generate_families{
    input {
        File bed
        File ped_file
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(select_all([bed, ped_file]), "GB")
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
        cat ~{bed} | gunzip | tail -n+2 | cut -f 1-6 | grep 'DEL\|DUP' | cut -f6 | tr ',' '\n' | sort -u > samples.txt
        grep -w -f samples.txt ~{ped_file} | cut -f1 | sort -u > families.txt
        >>>

    output{
        Array[String] families = read_lines("families.txt")
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: sv_base_mini_docker
    }
}

task generate_per_family_bed {
    input{
        File bed
        String family
        File ped_file
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([bed, ped_file]), "GB")
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
        cat ~{ped_file} | grep -w ~{family} | cut -f2 | sort -u > samples_in_family.txt
        # keep DEL/DUP rows carried by a member of this family; preserve all 6 columns
        cat ~{bed} | gunzip | tail -n+2 | cut -f1-6 | grep 'DEL\|DUP' | grep -w -f samples_in_family.txt > per_family_bed.bed || true
        >>>

    output {
        File bed_file = "per_family_bed.bed"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task depth_plot {
    input{
        String family
        File per_family_bed
        File ped_file
        File sample_bam_bai
        Int flank
        Int window
        String prefix
        String long_read_visualize_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(select_all([per_family_bed, sample_bam_bai, ped_file]), "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(20 + input_size),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail
        mkdir rd_plots

        # family members and their bam/bai
        grep -w ^~{family} ~{ped_file} | cut -f2 | sort -u > fam_samples.txt
        grep -w -f fam_samples.txt ~{sample_bam_bai} > fam_scc.txt   # sample <tab> bai <tab> bam

        # regions (+/- flank) and fine windows for mosdepth
        cut -f1-3 ~{per_family_bed} \
            | awk -v F=~{flank} '{s=$2-F; if(s<0)s=0; print $1"\t"s"\t"$3+F}' \
            | sort -k1,1 -k2,2n | bedtools merge -i - > regions.bed
        bedtools makewindows -b regions.bed -w ~{window} > windows.bed

        while read sample bai bam; do
            gsutil cp $bai $( basename $bam ).bai || true
            export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
            samtools view -b -o $sample.bam $bam -L regions.bed -M
            samtools index $sample.bam
            mosdepth --by windows.bed --no-per-base -n $sample $sample.bam
        done < fam_scc.txt

        python3 /src/variant-interpretation/scripts/plot_longread_depth.py \
            --bed ~{per_family_bed} \
            --ped ~{ped_file} \
            --family ~{family} \
            --flank ~{flank} \
            --depth-dir . \
            --outdir rd_plots

        tar -czf rd_plots.tar.gz rd_plots/
    >>>

    output {
        File plots = "rd_plots.tar.gz"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: long_read_visualize_docker
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task integrate_depth_plots{
    input {
        Array[File] depth_tar
        String prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(depth_tar, "GB")
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
        mkdir ~{prefix}_rd_plots
        while read file; do
            tar -zxf ${file}
            mv rd_plots/*  ~{prefix}_rd_plots/ || true
        done < ~{write_lines(depth_tar)};

        tar -czf ~{prefix}_rd_plots.tar.gz ~{prefix}_rd_plots
    >>>

    output{
        File plot_tar = "~{prefix}_rd_plots.tar.gz"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
