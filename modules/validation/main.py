#!/usr/bin/env python3

import argparse
import sys
import traceback
from pathlib import Path

# Add modules/validation/ to sys.path so utils/ is importable.
sys.path.insert(0, str(Path(__file__).parent))

from validation_pkg import (
    ConfigManager,
    GenomeValidator,
    ReadValidator,
    FeatureValidator,
    ReadXReadSettings,
    GenomeXGenomeSettings,
    ValidationReport,
    validate_genome,
    validate_reads,
    validate_feature,
    setup_logging,
    get_logger,
    readxread_validation,
    genomexgenome_validation
)
from validation_pkg.exceptions import ValidationError
from validation_pkg.utils.formats import CodingType, GenomeFormat
from utils.ref_defragment import defragment_reference

import utils.nextflow_params_handler as nf_params


def main():
    parser = argparse.ArgumentParser(description="Validation pipeline for genomic input files")
    parser.add_argument("config_path", help="Path to config.json")
    parser.add_argument("--threads",          type=int,  help="Number of threads (overrides config.json)")
    parser.add_argument("--validation-level", choices=["strict", "trust", "minimal"], help="Validation depth (overrides config.json)")
    parser.add_argument("--logging-level",    choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Log verbosity (overrides config.json)")
    parser.add_argument("--type",             dest="organism_type", choices=["prokaryote", "eukaryote"], help="Organism type (overrides config.json)")
    parser.add_argument("--force-defragment-ref", action="store_true", default=False, help="Merge fragmented reference contigs (unsupported workaround)")
    parsed = parser.parse_args()

    config_path = Path(parsed.config_path).resolve()

    cli_options = {}
    if parsed.threads            is not None: cli_options["threads"]             = parsed.threads
    if parsed.validation_level   is not None: cli_options["validation_level"]    = parsed.validation_level
    if parsed.logging_level      is not None: cli_options["logging_level"]       = parsed.logging_level
    if parsed.organism_type      is not None: cli_options["type"]                = parsed.organism_type
    if parsed.force_defragment_ref:           cli_options["force_defragment_ref"] = True

    output_dir = Path.cwd()
    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = config_path.parent.parent / "outputs" / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    logger = None
    log_file = None

    # Setup logging
    log_filename = "validation.log"
    try:
        logger = setup_logging(console_level='DEBUG', log_file=logs_dir / log_filename)
    except (PermissionError, OSError) as e:
        logger = setup_logging(console_level='DEBUG')
        logger.warning(f"Could not write log file ({e}); logging to console only")
    log_file = getattr(logger, 'log_file', None)

    # ========================================================================
    # Step 1: Read and validate config
    # ========================================================================
    config = None
    try:
        config = ConfigManager.load(config_path, cli_options=cli_options or None)
    except Exception as e:
        logger.error(f"Loading a config file failed: {e}")
        return 1

    # Override sub-config output dirs to CWD. ConfigManager sets output_dir from
    # the config file location, which is wrong inside a Nextflow work directory.
    _all_sub_configs = [
        config.ref_genome, config.mod_genome,
        config.ref_plasmid, config.mod_plasmid,
        config.ref_feature,
    ] + list(config.reads or [])
    for sub_cfg in _all_sub_configs:
        if sub_cfg is not None:
            sub_cfg.output_dir = output_dir

    # ========================================================================
    # Step 1.5 (optional): Defragment reference if force_defragment_ref is set
    # Priority: config.json > CLI arg > default (false)
    # ========================================================================
    force_defragment = config.force_defragment_ref
    if force_defragment:
        logger.warning("=" * 70)
        logger.warning("UNSUPPORTED WORKAROUND: --force-defragment-ref is ACTIVE")
        logger.warning("=" * 70)
        logger.warning(
            "The reference genome is highly fragmented and is NOT suitable "
            "for this pipeline. Merging contigs is a workaround only."
        )
        logger.warning(
            "All downstream results (inter-genome alignment, feature "
            "coordinate mapping, variant calling) may be INCORRECT or "
            "MEANINGLESS when run on an artificially merged reference."
        )
        logger.warning(
            "EFSA pipeline does NOT support fragmented references. "
            "Do NOT use these results for regulatory submissions or "
            "biological conclusions without expert review."
        )
        logger.warning(
            f"Original reference: {config.ref_genome.filepath} "
            "— consider obtaining a properly assembled genome instead."
        )
        logger.warning("=" * 70)
        try:
            merged_fasta, join_tsv = defragment_reference(config.ref_genome)
        except Exception as e:
            logger.error(f"Defragmentation failed: {e}")
            return 1
        config.ref_genome.filepath = merged_fasta
        config.ref_genome.filename = merged_fasta.name
        config.ref_genome.coding_type = CodingType.NONE
        config.ref_genome.detected_format = GenomeFormat.FASTA
        config.ref_genome._extract_basename()
        logger.warning(
            f"Reference replaced for this run. "
            f"Merged file: {merged_fasta} | Join order: {join_tsv}"
        )
        logger.warning(
            "Proceeding with validation on merged reference. "
            "Results should be interpreted with extreme caution."
        )

    # ========================================================================
    # Step 2: Edit settings for each validator
    # ========================================================================

    # Settings for reference genome
    try:
        ref_genome_settings = GenomeValidator.Settings(
            plasmids_to_one=True,
            main_longest=True,
            coding_type=None,
            output_filename_suffix='ref',
            replace_id_with_incremental='chr',
            min_sequence_length=100
        )

        # Settings for modified genome
        mod_genome_settings = GenomeValidator.Settings(
            plasmids_to_one=False,
            coding_type=None,
            output_filename_suffix='mod',
            replace_id_with_incremental='chr',
            min_sequence_length=100
        )

        # Settings for ref plasmid genomes (if you have them)
        ref_plasmid_settings = GenomeValidator.Settings(
            is_plasmid=True,
            plasmids_to_one=True,
            coding_type=None,
            output_filename_suffix='ref_plasmid'
        )

        # Settings for mod plasmid genomes (if you have them)
        mod_plasmid_settings = GenomeValidator.Settings(
            is_plasmid=True,
            plasmids_to_one=True,
            coding_type=None,
            output_filename_suffix='mod_plasmid'
        )

        # Settings for reads
        reads_settings = ReadValidator.Settings(
            coding_type='gz',
            outdir_by_ngs_type=True
        )

        # Settings for reference features
        ref_feature_settings = FeatureValidator.Settings(
            sort_by_position=False,
            check_coordinates=False,
            replace_id_with='chr',
            coding_type=None,
            output_filename_suffix='ref'
        )

        # Inter genome validation settings (using defaults)
        genomexgenome_settings = GenomeXGenomeSettings(
            characterize=True,
            same_sequence_ids=False,
            same_number_of_sequences=False
        )

        # Inter read validation settings (using defaults)
        readxread_settings = ReadXReadSettings()

    except Exception as e:
        logger.error(f"Setting up validators failed: {e}")
        return 1


    # ========================================================================
    # Step 3: Run validation using functional API
    # ========================================================================
    report_filename = "report.txt"
    report = ValidationReport(logs_dir / report_filename)
    fatal_errors: list[str] = []

    def register_required_failure(label: str, exc: Exception) -> None:
        message = f"{label} validation failed: {exc}"
        logger.error(message)
        fatal_errors.append(message)

    def register_missing_output(label: str, result) -> None:
        output_file = getattr(result, "output_file", None) if result is not None else None
        if not output_file:
            message = f"{label} validation produced no usable output file"
            logger.error(message)
            fatal_errors.append(message)

    # Validate reference genome (required)
    ref_genome_res = None
    if hasattr(config, 'ref_genome') and config.ref_genome:
        try:
            ref_genome_res = validate_genome(config.ref_genome, ref_genome_settings)
            report.write(ref_genome_res, file_type="genome")
            register_missing_output("ref_genome", ref_genome_res)
        except ValidationError as e:
            register_required_failure("ref_genome", e)

    # Validate modified genome (optional)
    mod_genome_res = None
    if hasattr(config, 'mod_genome') and config.mod_genome:
        try:
            mod_genome_res = validate_genome(config.mod_genome, mod_genome_settings)
            report.write(mod_genome_res, file_type="genome")
        except ValidationError as e:
            logger.error(f"Optional mod_genome validation failed: {e}")

    # Validate plasmid genomes (optional)
    ref_plasmid_res = None
    if hasattr(config, 'ref_plasmid') and config.ref_plasmid:
        try:
            ref_plasmid_res = validate_genome(config.ref_plasmid, ref_plasmid_settings)
            report.write(ref_plasmid_res, file_type="genome")
        except ValidationError as e:
            logger.error(f"Optional ref_plasmid validation failed: {e}")

    mod_plasmid_res = None
    if hasattr(config, 'mod_plasmid') and config.mod_plasmid:
        try:
            mod_plasmid_res = validate_genome(config.mod_plasmid, mod_plasmid_settings)
            report.write(mod_plasmid_res, file_type="genome")
        except ValidationError as e:
            logger.error(f"Optional mod_plasmid validation failed: {e}")

    # Inter-genome validation — only if both genomes validated successfully and mod is not fragmented
    genomexgenome_res = None
    if mod_genome_res is not None and ref_genome_res is not None and not getattr(mod_genome_res, 'fragmented', False):
        try:
            genomexgenome_res = genomexgenome_validation(ref_genome_res, mod_genome_res, genomexgenome_settings, mod_plasmid_res)
            report.write(genomexgenome_res, file_type="genomexgenome")
        except ValidationError as e:
            logger.error(f"Inter-genome validation failed: {e}")
    else:
        logger.info("Inter-genome validation skipped")

    # Validate reads (required)
    reads_res = None
    if hasattr(config, 'reads') and config.reads:
        try:
            reads_res = validate_reads(config.reads, reads_settings)
            report.write(reads_res, file_type="read")
            for read_result in reads_res:
                register_missing_output("reads", read_result)
        except ValidationError as e:
            register_required_failure("reads", e)

    # Add interread validation — skip when all reads are BAM (pairing check is meaningless)
    readxread_res = None
    fastq_reads = [r for r in (reads_res or []) if getattr(r, "input_format", None) != "bam"]
    if reads_res is not None and fastq_reads:
        try:
            readxread_res = readxread_validation(fastq_reads, readxread_settings)
            report.write(readxread_res, file_type="readxread")
        except ValidationError as e:
            logger.error(f"Inter-read validation failed: {e}")
    else:
        logger.info("Inter-read validation skipped")

    # Validate features (optional — non-fatal)
    ref_feature_res = None
    if hasattr(config, 'ref_feature') and config.ref_feature and not force_defragment:
        try:
            ref_feature_res = validate_feature(config.ref_feature, ref_feature_settings)
            report.write(ref_feature_res, file_type="feature")
        except ValidationError as e:
            logger.error(f"Optional ref_feature validation failed: {e}")

    if fatal_errors:
        report.add_fatal_errors(fatal_errors)

    report.flush(format='text')
    print(f"Log file: {log_file}")

    # ========================================================================
    # Step 4: Write validated_params.json for Nextflow (-params-file)
    # ========================================================================
    validation_results = {
        "ref_genome":    ref_genome_res,
        "mod_genome":    mod_genome_res,
        "ref_plasmid":   ref_plasmid_res,
        "mod_plasmid":   mod_plasmid_res,
        "genomexgenome": genomexgenome_res,
        "reads":         reads_res,
        "ref_feature":   ref_feature_res,
    }
    if force_defragment:
        logger.warning(
            "force_defragment_ref is active: GFF validation for the reference is "
            "skipped. Feature coordinates are not meaningful on a defragmented reference."
        )
    repo_root = config_path.parent.parent.parent
    params = nf_params.build_params(validation_results, base_dir=repo_root, organism_type=config.type.value)
    nf_params.write_params(params, output_dir / "validated_params.json")

    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as e:
        logger = get_logger()
        logger.error(f"✗ Fatal error: {e}")
        logger.debug(traceback.format_exc())
        if len(sys.argv) >= 2:
            actual_log_file = getattr(get_logger(), 'log_file', None) or (Path(sys.argv[1]).resolve().parent.parent / "valid" / "validation.log")
            print(f"Log file: {actual_log_file}")
        sys.exit(1)
