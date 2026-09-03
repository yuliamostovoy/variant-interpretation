version 1.0

##########################################################################################
##
## Input adapter: takes a user-curated variant list and emits the canonical bgzipped
## 6-column BED (chrom,start,end,ID,svtype,samples) with a header line that the IGV / depth
## tracks consume.
##
##########################################################################################

import "Structs2.wdl"

workflow ReformatVariants {
    input {
        File variant_list
        File? variant_vcf
        File? variant_vcf_index
        String prefix
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    call reformat_variants {
        input:
            variant_list = variant_list,
            variant_vcf = variant_vcf,
            variant_vcf_index = variant_vcf_index,
            prefix = prefix,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File varfile = reformat_variants.varfile
        File genotypes = reformat_variants.genotypes
    }
}

task reformat_variants {
    input {
        File variant_list
        File? variant_vcf
        File? variant_vcf_index
        String prefix
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([variant_list, variant_vcf]), "GB")
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

        # Optional per-sample genotypes for the pedigree glyph. Reference the VCF as a File
        # placeholder (not concatenated into a String) so Cromwell localizes it and bcftools
        # reads the local copy -- coercing a File to String yields the raw gs:// path, which
        # htslib then fails to open ("Permission denied").
        VCF="~{default='' variant_vcf}"
        GT_ARGS=()
        if [ -n "$VCF" ]; then
            bcftools query \
                -f '%CHROM\t%POS\t%END\t%ID\t%INFO/SVTYPE[\t%SAMPLE=%GT]\n' \
                "$VCF" > gts.raw.tsv
            GT_ARGS=(--genotypes-raw gts.raw.tsv --genotypes-out ~{prefix}.genotypes.tsv)
        fi

        python3 /src/variant-interpretation/scripts/reformat_variants_for_visualization.py \
            --input ~{variant_list} \
            --output ~{prefix}.variants_for_visualization.bed \
            "${GT_ARGS[@]}"

        bgzip ~{prefix}.variants_for_visualization.bed
        # ensure the genotypes output always exists (empty when no VCF was supplied)
        touch ~{prefix}.genotypes.tsv
    >>>

    output {
        File varfile = "~{prefix}.variants_for_visualization.bed.gz"
        File genotypes = "~{prefix}.genotypes.tsv"
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
