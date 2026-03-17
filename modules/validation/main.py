#!/usr/bin/env python3

import sys
import traceback
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

from pathlib import Path

import nextflow_params as nf_params

def main():
    # Check command line arguments
    if len(sys.argv) < 2:
        print("Usage: python main.py <config_path>")
        print("\nExample:")
        print("  python main.py config.json")
        return 1

    # Derive output directory from config path (mirrors ConfigManager._setup_output_directory)
    config_path = Path(sys.argv[1]).resolve()
    output_dir = config_path.parent.parent / "valid"
    logger = None
    log_file = None

    # Setup logging,
    try:
        logger = setup_logging(console_level='DEBUG', log_file=output_dir / "validation.log")
        log_file = logger.log_file
    except Exception as e:
        print(f"Failed to setup logging: {e}")
        return 1


    # ========================================================================
    # Step 1: Read and validate config
    # ========================================================================
    config = None
    try:
        config = ConfigManager.load(config_path)
    except Exception as e:
        logger.error(f"Loading a config file failed: {e}")
        return 1


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


    # ========================================================================
    # Step 3: Run validation using functional API
    # ========================================================================
    report = ValidationReport(output_dir / "report.txt")

    # Validate reference genome (required)
    ref_genome_res = None
    if hasattr(config, 'ref_genome') and config.ref_genome:
        try:
            ref_genome_res = validate_genome(config.ref_genome, ref_genome_settings)
            report.write(ref_genome_res, file_type="genome")
        except ValidationError as e:
            logger.error(f"Required ref_genome validation failed: {e}")

    # Validate modified genome (optional)
    mod_genome_res = None
    if hasattr(config, 'mod_genome') and config.mod_genome:
        try:
            mod_genome_res = validate_genome(config.mod_genome, mod_genome_settings)
            report.write(mod_genome_res, file_type="genome")
        except ValidationError as e:
            logger.warning(f"Optional mod_genome validation failed: {e}")

    # Inter-genome validation — only if both genomes validated successfully
    genomexgenome_res = None
    if mod_genome_res is not None:
        try:
            genomexgenome_res = genomexgenome_validation(ref_genome_res, mod_genome_res, genomexgenome_settings)
            report.write(genomexgenome_res, file_type="genomexgenome")
        except ValidationError as e:
            logger.warning(f"Inter-genome validation failed: {e}")

    # Validate plasmid genomes (optional)
    if hasattr(config, 'ref_plasmid') and config.ref_plasmid:
        try:
            res = validate_genome(config.ref_plasmid, plasmid_settings)
            report.write(res, file_type="genome")
        except ValidationError as e:
            logger.warning(f"Optional ref_plasmid validation failed: {e}")

    if hasattr(config, 'mod_plasmid') and config.mod_plasmid:
        try:
            res = validate_genome(config.mod_plasmid, plasmid_settings)
            report.write(res, file_type="genome")
        except ValidationError as e:
            logger.warning(f"Optional mod_plasmid validation failed: {e}")

    # Validate reads (required)
    try:
        reads_res = validate_reads(config.reads, reads_settings)
        report.write(reads_res, file_type="read")
    except ValidationError as e:
        logger.warning(f"Read validation failed: {e}")

    # Add interread validation
    try:
        readxread_res = readxread_validation(reads_res, readxread_settings)
        report.write(readxread_res, file_type="readxread")
    except ValidationError as e:
        logger.warning(f"Inter-read validation failed: {e}")

    # Validate features (optional — non-fatal)
    ref_feature_res = None
    if hasattr(config, 'ref_feature') and config.ref_feature:
        try:
            ref_feature_res = validate_feature(config.ref_feature, ref_feature_settings)
            report.write(ref_feature_res, file_type="feature")
        except ValidationError as e:
            logger.warning(f"Optional ref_feature validation failed: {e}")

    if hasattr(config, 'mod_feature') and config.mod_feature:
        try:
            res = validate_feature(config.mod_feature, mod_feature_settings)
            report.write(res, file_type="feature")
        except ValidationError as e:
            logger.warning(f"Optional mod_feature validation failed: {e}")

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
    params = nf_params.build_params(validation_results)
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
            actual_log_file = get_logger().log_file or (Path(sys.argv[1]).resolve().parent.parent / "valid" / "validation.log")
            print(f"Log file: {actual_log_file}")
        sys.exit(1)
