version 1.0

##########################################################################################
##
## Component 1: top-level long-read visualization workflow.
##
## Long-read (PacBio HiFi) adaptation of RunVisualizePlots.wdl for data NOT processed with
## GATK-SV. Variants come from Sniffles2/pbsv (SV), HiFiCNV (CNV) and DeepVariant->GLNexus
## (SNV/indel); the user supplies a curated subset as a bcftools-query TSV.
##
## Tracks (vs the GATK-SV original):
##   - IGV reads track over aligned BAMs (kept, generalized to long reads).
##   - mosdepth-based depth track for DEL/DUP (replaces the GATK-SV RdTest bincov track).
##   - The GATK-SV PE/SR-evidence track and complex-SV bed splitting are dropped.
##
##########################################################################################

import "Structs2.wdl"
import "ReformatVariants.wdl" as reformat
import "CreateIgvBamPlots.wdl" as igv_bam
import "LongReadDepthPlot.wdl" as depth

workflow VisualizePlotsLongRead {
    input {
        File variant_list          # curated subset (bcftools-query TSV: chrom,POS0,END,ID,SVTYPE,samples)
        File pedfile
        File sample_bam_bai        # sample <tab> bai <tab> bam
        File reference
        File reference_index
        String prefix
        File? fam_ids

        Boolean run_IGV = true
        Boolean run_depth = true
        Boolean file_localization = false
        Boolean requester_pays = false
        Boolean long_read = true

        Int? igv_max_window
        String? buffer
        Int? depth_flank
        Int? depth_window

        String sv_base_mini_docker
        String igv_docker
        String variant_interpretation_docker
        String long_read_visualize_docker

        RuntimeAttr? runtime_attr_reformat
        RuntimeAttr? runtime_attr_run_igv
        RuntimeAttr? runtime_attr_igv
        RuntimeAttr? runtime_attr_update_scc
        RuntimeAttr? runtime_attr_depth
        RuntimeAttr? runtime_attr_create_bed
        RuntimeAttr? runtime_attr_concat
    }

    String buffer_ = select_first([buffer, "500"])
    Int igv_max_window_ = select_first([igv_max_window, 150000])

    # Component 2: normalize the curated variant list into the canonical bgzipped BED
    call reformat.ReformatVariants as reformat_variants {
        input:
            variant_list = variant_list,
            prefix = prefix,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_override = runtime_attr_reformat
    }

    # Component 3: IGV reads track
    if (run_IGV) {
        call igv_bam.IGV_all_samples as igv_plots {
            input:
                ped_file = pedfile,
                fam_ids = fam_ids,
                sample_bam_bai = sample_bam_bai,
                varfile = reformat_variants.varfile,
                reference = reference,
                reference_index = reference_index,
                igv_max_window = igv_max_window_,
                file_localization = file_localization,
                requester_pays = requester_pays,
                long_read = long_read,
                prefix = prefix,
                buffer = buffer_,
                sv_base_mini_docker = sv_base_mini_docker,
                igv_docker = igv_docker,
                variant_interpretation_docker = variant_interpretation_docker,
                runtime_attr_run_igv = runtime_attr_run_igv,
                runtime_attr_igv = runtime_attr_igv,
                runtime_attr_update_scc = runtime_attr_update_scc
        }
    }

    # Component 4: mosdepth depth track for DEL/DUP
    if (run_depth) {
        call depth.LongReadDepthPlot as depth_plots {
            input:
                prefix = prefix,
                bed = reformat_variants.varfile,
                ped_file = pedfile,
                fam_ids = fam_ids,
                sample_bam_bai = sample_bam_bai,
                reference = reference,
                reference_index = reference_index,
                flank = depth_flank,
                window = depth_window,
                sv_base_mini_docker = sv_base_mini_docker,
                long_read_visualize_docker = long_read_visualize_docker,
                runtime_attr_depth = runtime_attr_depth,
                runtime_attr_create_bed = runtime_attr_create_bed
        }
    }

    # Stack IGV (top) + depth (bottom) per variant when both tracks ran
    if (run_IGV && run_depth) {
        call concat_plots {
            input:
                igv_tar = select_first([igv_plots.tar_gz_pe]),
                depth_tar = select_first([depth_plots.Plots]),
                prefix = prefix,
                long_read_visualize_docker = long_read_visualize_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File? igv_plots_tar = igv_plots.tar_gz_pe
        File? depth_plots_tar = depth_plots.Plots
        File? combined_plots_tar = concat_plots.combined_tar
    }
}

task concat_plots {
    input {
        File igv_tar
        File depth_tar
        String prefix
        String long_read_visualize_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([igv_tar, depth_tar]), "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size * 2),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail
        mkdir -p igv_in depth_in ~{prefix}_igv_depth_plots

        tar -zxf ~{igv_tar} -C igv_in
        tar -zxf ~{depth_tar} -C depth_in

        # the track tars each wrap a single top-level dir; collect the PNGs from within
        mkdir -p igv_pngs depth_pngs
        find igv_in -name '*.png' -exec cp -n {} igv_pngs/ \;
        find depth_in -name '*.png' -exec cp -n {} depth_pngs/ \;

        python3 /src/variant-interpretation/scripts/concat_igv_depth.py \
            --igv-dir igv_pngs \
            --depth-dir depth_pngs \
            --outdir ~{prefix}_igv_depth_plots

        tar -czf ~{prefix}_igv_depth_plots.tar.gz ~{prefix}_igv_depth_plots
    >>>

    output {
        File combined_tar = "~{prefix}_igv_depth_plots.tar.gz"
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
