#!/usr/bin/env nextflow

// Include workflows
include { ref_mod } from "./fasta_ref_x_mod.nf"
include { long_read as long_ref_pacbio; long_read as long_ref_ont; long_read as long_mod_pacbio; long_read as long_mod_ont} from "./long_read.nf"
include { short_read as short_ref; short_read as short_mod } from "./short_read.nf"
include { truvari_comparison } from "./vcf_comparison.nf"
include { qc } from "./subworkflows.nf"

include { compare_unmapped; compare_unmapped as compare_unmapped_ont; compare_unmapped as compare_unmapped_pacbio } from "../modules/mapping.nf"
include { nanoplot as nanoplot_pacbio; nanoplot as nanoplot_ont } from "../modules/qc.nf"
include { restructure_sv_tbl; create_empty_tbl as create_ont_tbl; create_empty_tbl as create_asm_tbl; create_empty_tbl as create_pacbio_tbl; create_empty_tbl as create_short_tbl } from "../modules/sv_calling.nf"


workflow analysis {
    
    take:
        params_json

    main:

        // Parse validated_params.json into a value channel (reusable across all branches)
        pmap = params_json.map { f -> new groovy.json.JsonSlurper().parse(f) }.first()

        pmap | map { json ->
            def lines = ["", "╔══════════════════════════════════════════════════════╗",
                             "║           validated_params.json — loaded params      ║",
                             "╚══════════════════════════════════════════════════════╝",
                             "  [switches]",
                             "  run_ref_x_mod:      ${json.run_ref_x_mod}",
                             "  run_illumina:       ${json.run_illumina}",
                             "  run_nanopore:       ${json.run_nanopore}",
                             "  run_pacbio:         ${json.run_pacbio}",
                             "  run_vcf_annotation: ${json.run_vcf_annotation}",
                             "  contig_file_size:   ${json.contig_file_size}",
                             "  validation_timestamp: ${json.validation_timestamp}",
                             "  [files]",
                             "  ref_fasta:    ${json.ref_fasta_validated   ?: 'null'}",
                             "  mod_fasta:    ${json.mod_fasta_validated   ?: 'null'}",
                             "  ref_plasmid:  ${json.ref_plasmid_fasta     ?: 'null'}",
                             "  mod_plasmid:  ${json.mod_plasmid_fasta     ?: 'null'}",
                             "  gff:          ${json.gff                   ?: 'null'}",
                             "  [read lists]",
                             "  illumina_fastqs: ${json.illumina_fastqs ?: []}",
                             "  ont_fastqs:      ${json.ont_fastqs      ?: []}",
                             "  ont_bams:        ${json.ont_bams        ?: []}",
                             "  pacbio_fastqs:   ${json.pacbio_fastqs   ?: []}",
                             "  pacbio_bams:     ${json.pacbio_bams     ?: []}",
                             "  contig_files:    ${json.contig_files    ?: []}",
                             ""]
            log.info lines.join("\n")
        }

        // Pipeline switches
        run_ref_x_mod       = pmap.map { it.run_ref_x_mod }
        run_illumina        = pmap.map { it.run_illumina }
        run_nanopore        = pmap.map { it.run_nanopore }
        run_pacbio          = pmap.map { it.run_pacbio }
        run_vcf_annotation  = pmap.map { it.run_vcf_annotation }
        contig_file_size    = pmap.map { it.contig_file_size }
        validation_timestamp = pmap.map { it.validation_timestamp }

        // Single-file paths — emits only when the path is present, empty channel otherwise
        ref_fasta   = pmap.filter { it.ref_fasta_validated   }.map { file(it.ref_fasta_validated) }
        mod_fasta   = pmap.filter { it.mod_fasta_validated   }.map { file(it.mod_fasta_validated) }
        ref_plasmid = pmap.filter { it.ref_plasmid_fasta     }.map { file(it.ref_plasmid_fasta)   }
        mod_plasmid = pmap.filter { it.mod_plasmid_fasta     }.map { file(it.mod_plasmid_fasta)   }
        gff         = pmap.filter { it.gff                   }.map { file(it.gff)                 }

        // Read file list channels (empty list when absent)
        illumina_fastqs = pmap.flatMap { it.illumina_fastqs ?: [] }.map { file(it) }
        ont_fastqs      = pmap.flatMap { it.ont_fastqs      ?: [] }.map { file(it) }
        ont_bams        = pmap.flatMap { it.ont_bams        ?: [] }.map { file(it) }
        pacbio_fastqs   = pmap.flatMap { it.pacbio_fastqs   ?: [] }.map { file(it) }
        pacbio_bams     = pmap.flatMap { it.pacbio_bams     ?: [] }.map { file(it) }
        contigs_ch    = pmap.flatMap { it.contig_files    ?: [] }.map { file(it) }

        // --- Assembly pipeline — only when run_ref_x_mod = true ---
        // Gate both inputs through the flag so processes never start when false
        active_contigs  = run_ref_x_mod.combine(contigs_ch)
            .filter { flag, f -> flag }.map { flag, f -> f }
        active_ref_asm  = run_ref_x_mod.combine(ref_fasta)
            .filter { flag, f -> flag }.map { flag, f -> f }

        ref_mod(active_ref_asm, active_contigs)
        create_asm_tbl(run_ref_x_mod.filter { !it }.map { "assembly" })

        // --- PacBio long-read pipeline — only when run_pacbio = true ---
        active_pacbio = run_pacbio.combine(pacbio_fastqs)
            .filter { flag, f -> flag }.map { flag, f -> f }

        nanoplot_pacbio(active_pacbio, "pacbio")
        long_ref_pacbio(active_pacbio, ref_fasta, "map-pb", ref_plasmid, "pacbio/long-ref")
        long_mod_pacbio(active_pacbio, mod_fasta, "map-pb", mod_plasmid, "pacbio/long-mod")
        compare_unmapped_pacbio(long_ref_pacbio.out.unmapped_fastq, long_mod_pacbio.out.unmapped_fastq, "pacbio")
        create_pacbio_tbl(run_pacbio.filter { !it }.map { "pb" })

        // --- Nanopore long-read pipeline — only when run_nanopore = true ---
        active_ont = run_nanopore.combine(ont_fastqs)
            .filter { flag, f -> flag }.map { flag, f -> f }

        nanoplot_ont(active_ont, "ont")
        long_ref_ont(active_ont, ref_fasta, "map-ont", ref_plasmid, "ont/long-ref")
        long_mod_ont(active_ont, mod_fasta, "map-ont", mod_plasmid, "ont/long-mod")
        compare_unmapped_ont(long_ref_ont.out.unmapped_fastq, long_mod_ont.out.unmapped_fastq, "ont")
        create_ont_tbl(run_nanopore.filter { !it }.map { "ont" })

        // --- Illumina short-read pipeline — only when run_illumina = true ---
        active_illumina = run_illumina.combine(illumina_fastqs)
            .filter { flag, f -> flag }.map { flag, f -> f }

        illumina_reads = active_illumina
            .map { fobj ->
                def matcher = fobj.name =~ /^(.+?)(?:[_\.](S[0-9]+_L[0-9]+_)?(R[12]|[12]))?\.f(ast)?q\.gz$/
                matcher.matches() ? [matcher[0][1], fobj] : null
            }
            .filter { it != null }
            .groupTuple(sort: true)

        qc(illumina_reads, "illumina/qc_trimming") | set { trimmed }

        short_ref(trimmed, ref_fasta, "illumina/short-ref", ref_plasmid, run_vcf_annotation, gff)
        short_mod(trimmed, mod_fasta, "illumina/short-mod", mod_plasmid, run_vcf_annotation, gff)
        compare_unmapped(short_ref.out.unmapped_fastq, short_mod.out.unmapped_fastq, "short")
        create_short_tbl(run_illumina.filter { !it }.map { "short" })

        // --- Collect SV tables from all active branches ---
        // Each branch contributes either real results OR an empty placeholder,
        // so .out is always safe to reference (one of the two always ran).
        sv_tbl = ref_mod.out.sv_tbl.map         { tbl -> tuple('assembly', tbl) }
            .mix(create_asm_tbl.out.map          { tbl -> tuple('assembly', tbl) })
            .mix(long_ref_pacbio.out.sv_tbl.map  { tbl -> tuple('pb',       tbl) })
            .mix(create_pacbio_tbl.out.map        { tbl -> tuple('pb',       tbl) })
            .mix(long_ref_ont.out.sv_tbl.map     { tbl -> tuple('ont',      tbl) })
            .mix(create_ont_tbl.out.map           { tbl -> tuple('ont',      tbl) })
            .mix(short_ref.out.sv_tbl.map        { tbl -> tuple('short',    tbl) })
            .mix(create_short_tbl.out.map         { tbl -> tuple('short',    tbl) })

        supp_reads_ch = long_ref_pacbio.out.supp_reads
            .mix(long_ref_ont.out.supp_reads)

        // --- Truvari comparison (when enabled and multiple VCF sources available) ---
        if (params.run_truvari) {
            active_vcfs = ref_mod.out.sv_vcf
                .mix(long_ref_pacbio.out.sv_vcf)
                .mix(long_ref_ont.out.sv_vcf)
                .mix(short_ref.out.sv_vcf)

            truvari_comparison(ref_fasta, active_vcfs)
        }

        // --- Final SV table aggregation ---
        script = file("${workflow.projectDir}/modules/utils/create_sv_output.py")

        def tbl_channel = sv_tbl.collect().map { list ->
            def tagged = list.collate(2)
            def asm     = tagged.find { it[0] == 'assembly' }[1]
            def long_pb = tagged.find { it[0] == 'pb'       }[1]
            def long_ont = tagged.find { it[0] == 'ont'     }[1]
            def sht     = tagged.find { it[0] == 'short'    }[1]
            tuple(asm, long_ont, long_pb, sht)
        }

        restructure_sv_tbl(script, tbl_channel, supp_reads_ch.collect().ifEmpty(file('NO_FILE')))
}