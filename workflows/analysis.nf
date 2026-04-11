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
        // Parse validated params JSON at runtime (value channel for multiple consumption)
        pmap = params_json.map { json ->
            new groovy.json.JsonSlurper().parse(json.toFile())
        }.first()

        // Core genome file channels (single-item value channels)
        ref_fasta = pmap.map { file(it.ref_fasta_validated) }
        mod_fasta = pmap.map { file(it.mod_fasta_validated) }

        ref_plasmid = pmap.map { it.ref_plasmid_fasta ? [file(it.ref_plasmid_fasta)] : [] }
        mod_plasmid = pmap.map { it.mod_plasmid_fasta ? [file(it.mod_plasmid_fasta)] : [] }

        // --- Assembly pipeline (run_ref_x_mod && contig_file_size >= 1) ---
        contigs_ch = pmap
            .filter { it.run_ref_x_mod && it.contig_file_size >= 1 }
            .flatMap { it.contig_files.collect { f -> file(f) } }

        ref_mod(ref_fasta, contigs_ch)

        // Empty assembly table when ref_x_mod is not active
        create_asm_tbl(pmap.filter { !it.run_ref_x_mod }.map { "assembly" })

        // Helper: convert list of fastq paths to named tuple channel
        def toNamedFastq = { f ->
            def fobj = file(f)
            def name = fobj.name.replaceFirst(/\.fastq\.gz$/, '')
            tuple(name, fobj)
        }

        // --- PacBio long-read pipeline ---
        pacbio_fastqs = pmap
            .filter { it.run_pacbio }
            .flatMap { it.pacbio_fastqs }
            .map(toNamedFastq)

        nanoplot_pacbio(pacbio_fastqs, "pacbio")

        long_ref_pacbio(pacbio_fastqs, ref_fasta, "map-pb", ref_plasmid, "pacbio/long-ref")
        long_mod_pacbio(pacbio_fastqs, mod_fasta, "map-pb", mod_plasmid, "pacbio/long-mod")
        compare_unmapped_pacbio(long_ref_pacbio.out.unmapped_fastq, long_mod_pacbio.out.unmapped_fastq, "pacbio")

        // Empty pacbio table when not active
        create_pacbio_tbl(pmap.filter { !it.run_pacbio }.map { "pb" })

        // --- Nanopore long-read pipeline ---
        ont_fastqs = pmap
            .filter { it.run_nanopore }
            .flatMap { it.ont_fastqs }
            .map(toNamedFastq)

        nanoplot_ont(ont_fastqs, "ont")

        long_ref_ont(ont_fastqs, ref_fasta, "map-ont", ref_plasmid, "ont/long-ref")
        long_mod_ont(ont_fastqs, mod_fasta, "map-ont", mod_plasmid, "ont/long-mod")
        compare_unmapped_ont(long_ref_ont.out.unmapped_fastq, long_mod_ont.out.unmapped_fastq, "ont")

        // Empty ont table when not active
        create_ont_tbl(pmap.filter { !it.run_nanopore }.map { "ont" })

        // --- Illumina short-read pipeline ---
        illumina_reads = pmap
            .filter { it.run_illumina }
            .flatMap { it.illumina_fastqs }
            .map { f ->
                def fobj = file(f)
                def matcher = fobj.name =~ /^(.+?)(?:[_\.](S[0-9]+_L[0-9]+_)?(R[12]|[12]))?\.f(ast)?q\.gz$/
                matcher.matches() ? [matcher[0][1], fobj] : null
            }
            .filter { it != null }
            .groupTuple(sort: true)

        qc(illumina_reads, "illumina/qc_trimming") | set { trimmed }

        // VCF annotation flag and GFF from validated params
        run_vcf_annotation = pmap.map { it.run_vcf_annotation }
        gff = pmap.map { it.gff ? file(it.gff) : file('NO_GFF') }

        short_ref(trimmed, ref_fasta, "illumina/short-ref", ref_plasmid, run_vcf_annotation, gff)
        short_mod(trimmed, mod_fasta, "illumina/short-mod", mod_plasmid, run_vcf_annotation, gff)
        compare_unmapped(short_ref.out.unmapped_fastq, short_mod.out.unmapped_fastq, "short")

        // Empty short table when not active
        create_short_tbl(pmap.filter { !it.run_illumina }.map { "short" })

        // --- Collect SV tables from all branches (tagged with explicit labels) ---
        sv_tbl = ref_mod.out.sv_tbl.map { tbl -> tuple('assembly', tbl) }
            .mix(create_asm_tbl.out.map { tbl -> tuple('assembly', tbl) })
            .mix(long_ref_pacbio.out.sv_tbl.map { tbl -> tuple('pb', tbl) })
            .mix(create_pacbio_tbl.out.map { tbl -> tuple('pb', tbl) })
            .mix(long_ref_ont.out.sv_tbl.map { tbl -> tuple('ont', tbl) })
            .mix(create_ont_tbl.out.map { tbl -> tuple('ont', tbl) })
            .mix(short_ref.out.sv_tbl.map { tbl -> tuple('short', tbl) })
            .mix(create_short_tbl.out.map { tbl -> tuple('short', tbl) })

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
            def asm = tagged.find { it[0] == 'assembly' }[1]
            def long_pb = tagged.find { it[0] == 'pb' }[1]
            def long_ont = tagged.find { it[0] == 'ont' }[1]
            def sht = tagged.find { it[0] == 'short' }[1]
            tuple(asm, long_ont, long_pb, sht)
        }

        restructure_sv_tbl(script, tbl_channel, supp_reads_ch.collect().ifEmpty(file('NO_FILE')))
}