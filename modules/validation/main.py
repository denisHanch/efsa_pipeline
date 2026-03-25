#!/usr/bin/env python3

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
    readxread_validation,
    genomexgenome_validation
)
from validation_pkg.exceptions import ValidationError
from validation_pkg.utils.formats import CodingType, GenomeFormat
from utils.ref_defragment import defragment_reference
from validation_utils.logger import setup_logging, ValidationLogger

import nextflow_params_handler as nf_params

def _parse_cli_args(argv):
    """Parse CLI arguments and return (config_path, cli_options)."""
    args = list(argv)
    cli_options = {}

    # Extract boolean flag
    if '--force-defragment-ref' in args:
        cli_options['force_defragment_ref'] = True
        args.remove('--force-defragment-ref')

    # Extract key-value options
    i = 0
    remaining = []
    while i < len(args):
        if args[i] == '--threads' and i + 1 < len(args):
            val = args[i + 1]
            cli_options['threads'] = None if val == 'auto' else int(val)
            i += 2
        elif args[i] == '--validation-level' and i + 1 < len(args):
            cli_options['validation_level'] = args[i + 1]
            i += 2
        elif args[i] == '--logging-level' and i + 1 < len(args):
            cli_options['logging_level'] = args[i + 1].upper()
            i += 2
        elif args[i] == '--type' and i + 1 < len(args):
            cli_options['type'] = args[i + 1]
            i += 2
        else:
            remaining.append(args[i])
            i += 1

    return remaining, cli_options


def main():
    # Check command line arguments
    args, cli_options = _parse_cli_args(sys.argv[1:])

    if not args:
        print("Usage: python main.py <config_path> [--force-defragment-ref]")
        print("       [--threads N|auto] [--validation-level LEVEL]")
        print("       [--logging-level LEVEL] [--type TYPE]")
        print("\nExample:")
        print("  python main.py config.json")
        return 1

    # Derive output directory from config path (mirrors ConfigManager._setup_output_directory)
    config_path = Path(args[0]).resolve()
    output_dir = config_path.parent.parent / "valid"
    logger = None
    log_file = None

    # Setup logging
    try:
        logger = setup_logging(console_level='DEBUG', log_file=output_dir / "validation.log")
    except (PermissionError, OSError) as e:
        logger = setup_logging(console_level='DEBUG')
        logger.warning(f"Could not write log file ({e}); logging to console only")
    log_file = getattr(logger, 'log_file', None)

    # ========================================================================
    # Step 1: Read and validate config
    # ========================================================================
    config = None
    try:
        config = ConfigManager.load(config_path, cli_options=cli_options)
    except Exception as e:
        logger.error(f"Loading a config file failed: {e}")
        return 1

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

        # Settings for plasmid genomes (if you have them)
        plasmid_settings = GenomeValidator.Settings(
            is_plasmid=True,
            plasmids_to_one=True,
            coding_type=None,
            output_filename_suffix='plasmid'
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

        # Settings for modified features
        mod_feature_settings = FeatureValidator.Settings(
            sort_by_position=False,
            check_coordinates=False,
            replace_id_with='chr',
            coding_type=None,
            output_filename_suffix='mod'
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
    report = ValidationReport(output_dir / "report.txt")
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

    # Inter-genome validation — only if both genomes validated successfully and mod is not fragmented
    genomexgenome_res = None
    if mod_genome_res is not None and ref_genome_res is not None and not getattr(mod_genome_res, 'fragmented', False):
        try:
            genomexgenome_res = genomexgenome_validation(ref_genome_res, mod_genome_res, genomexgenome_settings)
            report.write(genomexgenome_res, file_type="genomexgenome")
        except ValidationError as e:
            logger.error(f"Inter-genome validation failed: {e}")
    else:
        logger.info("Inter-genome validation skipped")

    # Validate plasmid genomes (optional)
    if hasattr(config, 'ref_plasmid') and config.ref_plasmid:
        try:
            res = validate_genome(config.ref_plasmid, plasmid_settings)
            report.write(res, file_type="genome")
        except ValidationError as e:
            logger.error(f"Optional ref_plasmid validation failed: {e}")

    if hasattr(config, 'mod_plasmid') and config.mod_plasmid:
        try:
            res = validate_genome(config.mod_plasmid, plasmid_settings)
            report.write(res, file_type="genome")
        except ValidationError as e:
            logger.error(f"Optional mod_plasmid validation failed: {e}")

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

    # Add interread validation
    readxread_res = None
    if reads_res is not None:
        try:
            readxread_res = readxread_validation(reads_res, readxread_settings)
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

    if hasattr(config, 'mod_feature') and config.mod_feature:
        try:
            res = validate_feature(config.mod_feature, mod_feature_settings)
            report.write(res, file_type="feature")
        except ValidationError as e:
            logger.error(f"Optional mod_feature validation failed: {e}")

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
        "genomexgenome": genomexgenome_res,
        "reads":         reads_res,
        "ref_feature":   ref_feature_res,
    }
    if force_defragment:
        logger.warning(
            "force_defragment_ref is active: GFF validation for the reference is "
            "skipped. Feature coordinates are not meaningful on a defragmented "
            "reference — run_vcf_annotation will be disabled."
        )
    params = nf_params.build_params(validation_results, force_defragment_ref=force_defragment)
    nf_params.write_params(params, output_dir / "validated_params.json")

    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as e:
        _fallback_logger = ValidationLogger()
        _fallback_logger.error(f"✗ Fatal error: {e}")
        _fallback_logger.debug(traceback.format_exc())
        if len(sys.argv) >= 2:
            actual_log_file = (Path(sys.argv[1]).resolve().parent.parent / "valid" / "validation.log")
            print(f"Log file: {actual_log_file}")
        sys.exit(1)
