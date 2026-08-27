version 1.0

##########################################################################################
##
## Component 2: input adapter for the long-read visualization workflow.
##
## Takes a user-curated variant list (the limited subset to be plotted) and emits the
## canonical bgzipped 6-column BED (chrom,start,end,ID,svtype,samples) with a header line
## that the IGV / depth tracks consume.
##
##########################################################################################

import "Structs2.wdl"

workflow ReformatVariants {
    input {
        File variant_list
        String prefix
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    call reformat_variants {
        input:
            variant_list = variant_list,
            prefix = prefix,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File varfile = reformat_variants.varfile
    }
}

task reformat_variants {
    input {
        File variant_list
        String prefix
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(variant_list, "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size * 1.5),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail

        python3 /src/variant-interpretation/scripts/reformat_variants_for_visualization.py \
            --input ~{variant_list} \
            --output ~{prefix}.variants_for_visualization.bed

        bgzip ~{prefix}.variants_for_visualization.bed
    >>>

    output {
        File varfile = "~{prefix}.variants_for_visualization.bed.gz"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}
