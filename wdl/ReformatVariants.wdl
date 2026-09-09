version 1.0

##########################################################################################
##
## Input adapter: takes a user-curated variant list and emits the canonical bgzipped BED
## (chrom,start,end,ID,svtype,samples,svlen) with a header line that the IGV / depth tracks
## consume. svlen is the allele length in bp (from REF/ALT via the VCF when available).
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

        # Reference the VCF/index as File placeholders (not concatenated into a String) so
        # Cromwell localizes them and bcftools reads the local copies -- coercing a File to
        # String yields the raw gs:// path, which htslib then fails to open ("Permission denied").
        VCF="~{default='' variant_vcf}"
        VCFIDX="~{default='' variant_vcf_index}"
        INFO_ARGS=()
        GT_ARGS=()
        if [ -n "$VCF" ]; then
            # Restrict the VCF queries to the plotted loci when an index is available (make it
            # discoverable next to the localized VCF, which Cromwell may place apart); otherwise
            # fall back to a full scan.
            REGION=()
            if [ -n "$VCFIDX" ]; then
                case "$VCFIDX" in
                    *.csi) IDXLINK="${VCF}.csi" ;;
                    *)     IDXLINK="${VCF}.tbi" ;;
                esac
                # only create the link when the index isn't already sitting next to the VCF
                # (Cromwell may localize both into the same dir, making src == dst)
                if [ ! "$VCFIDX" -ef "$IDXLINK" ]; then
                    ln -sf "$VCFIDX" "$IDXLINK"
                fi
                awk 'BEGIN{FS=OFS="\t"} $2 ~ /^[0-9]+$/ {s=$2-1; if(s<0)s=0; print $1, s, $3}' \
                    ~{variant_list} | sort -k1,1 -k2,2n > var_regions.bed
                REGION=(-R var_regions.bed)
            fi

            # REF/ALT per variant (multiallelics split) for allele-length header sizing. Uses
            # only always-defined fields, so it works on SNV/indel VCFs lacking SVTYPE/SVLEN.
            bcftools view "${REGION[@]}" "$VCF" | bcftools norm -m- \
                | bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%REF\t%ALT\n' > variant_info.tsv
            INFO_ARGS=(--variant-info variant_info.tsv)

            # Per-sample genotypes for the pedigree glyph. SVTYPE is undefined in SNV/indel VCFs
            # and querying an undefined tag aborts bcftools, so emit a literal '.' in that case.
            if bcftools view -h "$VCF" | grep -q '##INFO=<ID=SVTYPE,'; then SVT='%INFO/SVTYPE'; else SVT='.'; fi
            bcftools query "${REGION[@]}" \
                -f "%CHROM\t%POS\t%END\t%ID\t${SVT}[\t%SAMPLE=%GT]\n" \
                "$VCF" > gts.raw.tsv
            GT_ARGS=(--genotypes-raw gts.raw.tsv --genotypes-out ~{prefix}.genotypes.tsv)
        fi

        python3 /src/variant-interpretation/scripts/reformat_variants_for_visualization.py \
            --input ~{variant_list} \
            --output ~{prefix}.variants_for_visualization.bed \
            "${INFO_ARGS[@]}" \
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
