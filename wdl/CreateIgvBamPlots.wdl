version 1.0

##########################################################################################
##
## Long-read / BAM adaptation of CreateIgvCramPlots.wdl.
##
## Fans out over families (or single samples), subsets the per-family aligned reads, and
## produces IGV screenshots at each variant locus using the BAM/long-read leaf workflow.
## The GATK-SV complex-SV bed splitting (updateCpxBed) is intentionally omitted: long-read
## SV callers (Sniffles2 / pbsv) do not emit that CPX format.
##
##########################################################################################

import "IgvBamPlots.wdl" as igv_plots
import "Structs2.wdl"

workflow IGV_all_samples {
    input {
        File ped_file
        File? fam_ids
        File sample_bam_bai
        File varfile
        File reference
        File reference_index
        Int igv_max_window
        Boolean file_localization
        Boolean requester_pays
        Boolean long_read
        String prefix
        String buffer
        String sv_base_mini_docker
        String igv_docker
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_update_scc
        RuntimeAttr? runtime_attr_run_igv
        RuntimeAttr? runtime_attr_igv
        RuntimeAttr? runtime_attr_cpx
    }

    if (defined(fam_ids)) {
        File fam_ids_ = select_first([fam_ids])
        Array[String] family_ids = transpose(read_tsv(fam_ids_))[0]
    }

    if (!(defined(fam_ids))) {
        call generate_families{
            input:
                varfile = varfile,
                ped_file = ped_file,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_run_igv
        }
    }
    scatter (family in select_first([family_ids, generate_families.families])){
        call generate_per_family_sample_bam_bai{
            input:
                family = family,
                ped_file = ped_file,
                sample_bam_bai = sample_bam_bai,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_run_igv
            }

        call update_sample_bam_bai{
            input:
                family = family,
                ped_file = ped_file,
                sample_bam_bai = generate_per_family_sample_bam_bai.subset_sample_bam_bai,
                bais_files = generate_per_family_sample_bam_bai.per_family_bais_strings,
                bams_files = generate_per_family_sample_bam_bai.per_family_bams_strings,
                variant_interpretation_docker = variant_interpretation_docker,
                runtime_attr_override = runtime_attr_update_scc
        }

        call generate_per_family_bed{
            input:
                varfile = varfile,
                samples = update_sample_bam_bai.per_family_samples,
                family = family,
                ped_file = ped_file,
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_run_igv
        }

        if (file_localization){
            call igv_plots.IGV as IGV_localize {
                input:
                    varfile = generate_per_family_bed.per_family_varfile,
                    family = family,
                    ped_file = ped_file,
                    samples = update_sample_bam_bai.per_family_samples,
                    file_localization = file_localization,
                    requester_pays = requester_pays,
                    long_read = long_read,
                    igv_max_window = igv_max_window,
                    bams_localize = generate_per_family_sample_bam_bai.per_family_bams_files,
                    bais_localize = generate_per_family_sample_bam_bai.per_family_bais_files,
                    sample_bam_bai = generate_per_family_sample_bam_bai.subset_sample_bam_bai,
                    buffer = buffer,
                    reference = reference,
                    reference_index = reference_index,
                    igv_docker = igv_docker,
                    variant_interpretation_docker = variant_interpretation_docker,
                    runtime_attr_igv = runtime_attr_igv
            }
        }

        if (!(file_localization)){
            call igv_plots.IGV as IGV_parse {
                input:
                    varfile = generate_per_family_bed.per_family_varfile,
                    family = family,
                    ped_file = ped_file,
                    file_localization = file_localization,
                    requester_pays = requester_pays,
                    long_read = long_read,
                    igv_max_window = igv_max_window,
                    bams_parse = generate_per_family_sample_bam_bai.per_family_bams_strings,
                    bais_parse = generate_per_family_sample_bam_bai.per_family_bais_strings,
                    samples = update_sample_bam_bai.per_family_samples,
                    updated_sample_bam_bai = update_sample_bam_bai.changed_sample_bam_bai,
                    buffer = buffer,
                    reference = reference,
                    reference_index = reference_index,
                    igv_docker = igv_docker,
                    variant_interpretation_docker = variant_interpretation_docker,
                    runtime_attr_igv = runtime_attr_igv
            }
        }
    }

    call integrate_igv_plots{
        input:
            igv_tar = select_all(flatten([IGV_localize.tar_gz_pe, IGV_parse.tar_gz_pe])),
            prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_run_igv
    }

    output{
        File tar_gz_pe = integrate_igv_plots.plot_tar
    }
}

task generate_families{
    input {
        File varfile
        File ped_file
        String sv_base_mini_docker
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
        cat ~{varfile} | gunzip | tail -n+2 | cut -f6 | tr ',' '\n' | sort -u > samples.txt #must have header line
        grep -w -f samples.txt ~{ped_file} | cut -f1 | sort -u  > families.txt
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

task generate_per_family_sample_bam_bai{
    input {
        String family
        File ped_file
        File sample_bam_bai
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(select_all([sample_bam_bai, ped_file]), "GB")
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
        grep -w ^~{family} ~{ped_file} | cut -f2 > samples_list.txt
        grep -w -f samples_list.txt ~{sample_bam_bai} > subset_sample_bam_bai.txt
        cut -f1 subset_sample_bam_bai.txt > samples.txt
        cut -f2 subset_sample_bam_bai.txt > bai.txt
        cut -f3 subset_sample_bam_bai.txt > bam.txt
        >>>

    output{
        Array[String] per_family_samples = read_lines("samples.txt")
        Array[File] per_family_bams_files = read_lines("bam.txt")
        Array[File] per_family_bais_files = read_lines("bai.txt")
        Array[String] per_family_bams_strings = read_lines("bam.txt")
        Array[String] per_family_bais_strings = read_lines("bai.txt")
        File subset_sample_bam_bai = "subset_sample_bam_bai.txt"
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

task update_sample_bam_bai{
    input {
        String family
        File ped_file
        File sample_bam_bai
        Array[String] bams_files
        Array[String] bais_files
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(select_all([sample_bam_bai, ped_file]), "GB")
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
        head -n+1 ~{ped_file} > family_ped.txt
        grep -w ^~{family} ~{ped_file} >> family_ped.txt
        python3 /src/variant-interpretation/scripts/renameCrams.py --ped family_ped.txt --scc ~{sample_bam_bai}
        mv changed_sample_crai_cram.txt changed_sample_bam_bai.txt
        cut -f1 changed_sample_bam_bai.txt > samples.txt
        cut -f5 changed_sample_bam_bai.txt > bai.txt
        cut -f4 changed_sample_bam_bai.txt > bam.txt
        >>>

    output{
        Array[String] per_family_samples = read_lines("samples.txt")
        Array[File] per_family_bams_files = read_lines("bam.txt")
        Array[File] per_family_bais_files = read_lines("bai.txt")
        Array[String] per_family_bams_strings = read_lines("bam.txt")
        Array[String] per_family_bais_strings = read_lines("bai.txt")
        File changed_sample_bam_bai = "changed_sample_bam_bai.txt"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: variant_interpretation_docker
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

task generate_per_family_bed{
    input {
        File varfile
        Array[String] samples
        String family
        File ped_file
        String sv_base_mini_docker
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
    String filename = basename(varfile, ".bed")

    command <<<
        set -euo pipefail
        cat ~{varfile} | gunzip | cut -f1-6 > updated_varfile.bed
        grep -f ~{write_lines(samples)} updated_varfile.bed | cut -f1-6 | awk '{print $1,$2,$3,$4,$5,$6}' | sed -e 's/ /\t/g' > ~{filename}.~{family}.bed
        >>>

    output{
        File per_family_varfile = "~{filename}.~{family}.bed"
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

task integrate_igv_plots{
    input {
        Array[File] igv_tar
        String prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(igv_tar, "GB")
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
        mkdir ~{prefix}_igv_plots
        while read file; do
            tar -zxf ${file}
            mv pe_igv_plots/*  ~{prefix}_igv_plots/
        done < ~{write_lines(igv_tar)};
        tar -czf ~{prefix}_igv_plots.tar.gz ~{prefix}_igv_plots
    >>>

    output{
        File plot_tar = "~{prefix}_igv_plots.tar.gz"
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
