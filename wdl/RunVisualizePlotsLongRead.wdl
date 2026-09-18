version 1.0

##########################################################################################
##
## Top-level long-read (PacBio HiFi) visualization workflow. The user supplies a curated
## subset of variants as a bcftools-query TSV; it produces:
##   - an IGV reads track over aligned BAMs, and
##   - a mosdepth-based depth track for DEL/DUP,
## combined per variant with an SV-info header.
##
##########################################################################################

import "Structs2.wdl"
import "ReformatVariants.wdl" as reformat
import "CreateIgvBamPlots.wdl" as igv_bam
import "LongReadDepthPlot.wdl" as depth

workflow VisualizePlotsLongRead {
    input {
        File variant_list          # curated subset (bcftools-query TSV: chrom,POS,END,ID,SVTYPE,samples)
        File? variant_vcf          # optional source VCF; supplies per-sample genotypes for the pedigree glyph
        File? variant_vcf_index
        File pedfile
        Array[String] sample_ids
        Array[String] bams
        Array[String] bais
        String igv_genome = "hg38"
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
        Float? depth_flank_frac
        Int? depth_window
        Int? depth_target_bins      # target bins per plotted region; window scales up for large events
        Int? depth_min_svlen        # skip the depth plot for DEL/DUP shorter than this (default 1 kb)
        # optional per-sample genome-wide median coverage, aligned by index to sample_ids/bams.
        # When supplied the depth track normalizes each sample by its own median so chrX/chrY
        # ploidy is correct; otherwise it falls back to a local per-region normalization.
        Array[Float]? sample_median_coverages
        # optional reference BEDs to highlight on the depth plots (e.g. N-gaps, segdups),
        # with matching labels in the same order
        Array[File] annotation_beds = []
        Array[String] annotation_names = []
        # optional gene annotation shown as an IGV gene track; supply a bgzipped, tabix-indexed
        # file as gene_track plus its .tbi as gene_track_index (otherwise it is indexed at runtime)
        File? gene_track
        File? gene_track_index

        String sv_base_mini_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52"
        String long_read_visualize_docker = "quay.io/ymostovoy/lr-visualize:latest"

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

    # Assemble the sample <tab> bai <tab> bam manifest the sub-workflows expect from the
    # per-individual arrays (index second, bam third).
    call make_sample_manifest {
        input:
            sample_ids = sample_ids,
            bais = bais,
            bams = bams,
            sv_base_mini_docker = sv_base_mini_docker
    }
    File sample_bam_bai = make_sample_manifest.manifest

    # Zip the per-sample median coverages (aligned to sample_ids) into a 'sample <tab> median'
    # file for the depth track; skipped entirely when no medians are provided.
    if (defined(sample_median_coverages)) {
        call make_median_manifest {
            input:
                sample_ids = sample_ids,
                median_coverages = select_first([sample_median_coverages]),
                sv_base_mini_docker = sv_base_mini_docker
        }
    }

    # normalize the curated variant list into the canonical bgzipped BED
    call reformat.ReformatVariants as reformat_variants {
        input:
            variant_list = variant_list,
            variant_vcf = variant_vcf,
            variant_vcf_index = variant_vcf_index,
            prefix = prefix,
            variant_interpretation_docker = long_read_visualize_docker,
            runtime_attr_override = runtime_attr_reformat
    }

    # IGV reads track
    if (run_IGV) {
        call igv_bam.IGV_all_samples as igv_plots {
            input:
                ped_file = pedfile,
                fam_ids = fam_ids,
                sample_bam_bai = sample_bam_bai,
                varfile = reformat_variants.varfile,
                igv_genome = igv_genome,
                igv_max_window = igv_max_window_,
                file_localization = file_localization,
                requester_pays = requester_pays,
                long_read = long_read,
                annotation_beds = annotation_beds,
                annotation_names = annotation_names,
                gene_track = gene_track,
                gene_track_index = gene_track_index,
                prefix = prefix,
                buffer = buffer_,
                sv_base_mini_docker = sv_base_mini_docker,
                igv_docker = long_read_visualize_docker,
                variant_interpretation_docker = long_read_visualize_docker,
                runtime_attr_run_igv = runtime_attr_run_igv,
                runtime_attr_igv = runtime_attr_igv,
                runtime_attr_update_scc = runtime_attr_update_scc
        }
    }

    # mosdepth depth track for DEL/DUP
    if (run_depth) {
        call depth.LongReadDepthPlot as depth_plots {
            input:
                prefix = prefix,
                bed = reformat_variants.varfile,
                ped_file = pedfile,
                fam_ids = fam_ids,
                sample_bam_bai = sample_bam_bai,
                median_coverage_file = make_median_manifest.manifest,
                annotation_beds = annotation_beds,
                annotation_names = annotation_names,
                flank = depth_flank,
                flank_frac = depth_flank_frac,
                min_svlen = depth_min_svlen,
                window = depth_window,
                target_bins = depth_target_bins,
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
                varfile = reformat_variants.varfile,
                pedfile = pedfile,
                genotypes = reformat_variants.genotypes,
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

task make_sample_manifest {
    input {
        Array[String] sample_ids
        Array[String] bais
        Array[String] bams
        String sv_base_mini_docker
    }

    command <<<
        set -euo pipefail
        paste ~{write_lines(sample_ids)} ~{write_lines(bais)} ~{write_lines(bams)} > sample_bam_bai.tsv
    >>>

    output {
        File manifest = "sample_bam_bai.tsv"
    }

    runtime {
        cpu: 1
        memory: "1 GiB"
        disks: "local-disk 10 HDD"
        bootDiskSizeGb: 8
        docker: sv_base_mini_docker
        preemptible: 2
        maxRetries: 1
    }
}

task make_median_manifest {
    input {
        Array[String] sample_ids
        Array[Float] median_coverages
        String sv_base_mini_docker
    }

    command <<<
        set -euo pipefail
        # one median per sample, in the same order as sample_ids
        paste ~{write_lines(sample_ids)} <(printf '%s\n' ~{sep=" " median_coverages}) > sample_median.tsv
        # guard against a length mismatch silently truncating the join
        if [ "$(wc -l < ~{write_lines(sample_ids)})" -ne "$(wc -l < sample_median.tsv)" ]; then
            echo "ERROR: sample_median_coverages length does not match sample_ids" >&2
            exit 1
        fi
    >>>

    output {
        File manifest = "sample_median.tsv"
    }

    runtime {
        cpu: 1
        memory: "1 GiB"
        disks: "local-disk 10 HDD"
        bootDiskSizeGb: 8
        docker: sv_base_mini_docker
        preemptible: 2
        maxRetries: 1
    }
}

task concat_plots {
    input {
        File igv_tar
        File depth_tar
        File varfile
        File pedfile
        File genotypes
        String prefix
        String long_read_visualize_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([igv_tar, depth_tar, varfile, pedfile, genotypes]), "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      # this task holds both gathered plot tars, extracts them,
                                      # copies the PNGs out (a 2nd copy), writes the combined PNGs,
                                      # then re-tars -> peak disk ~5x the tar size. The old 10+2x
                                      # under-provisioned and ran out of disk on the last large batch.
                                      disk_gb: ceil(20 + input_size * 6),
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
            --varfile ~{varfile} \
            --ped ~{pedfile} \
            --genotypes ~{genotypes} \
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
