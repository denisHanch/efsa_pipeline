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

sys.path.insert(0, str(Path(__file__).parent))
import nextflow_schema as nf_schema

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

    # Setup logging,
    logger = setup_logging(console_level='DEBUG', log_file=output_dir / "validation.log")
    log_file = logger.log_file

    # ========================================================================
    # Step 1: Read and validate config
    # ========================================================================
    config = ConfigManager.load(config_path)

    # ========================================================================
    # Step 2: Edit settings for each validator
    # ========================================================================

    # Settings for reference genome
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
        error_n_sequences=5,
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

    # ========================================================================
    # Step 3: Run validation using functional API
    # ========================================================================
    report = ValidationReport(output_dir / "report.txt")

    # Validate reference genome (required — fatal on failure)
    ref_genome_res = validate_genome(config.ref_genome, ref_genome_settings)
    report.write(ref_genome_res, file_type="genome")

    # Validate modified genome (optional — non-fatal)
    mod_genome_res = None
    if hasattr(config, 'mod_genome') and config.mod_genome:
        try:
            mod_genome_res = validate_genome(config.mod_genome, mod_genome_settings)
            report.write(mod_genome_res, file_type="genome")
        except ValidationError as e:
            logger.warning(f"Optional mod_genome validation failed: {e}")

    # Inter-genome validation — only if both genomes validated successfully
    if mod_genome_res is not None:
        try:
            res = genomexgenome_validation(ref_genome_res, mod_genome_res, genomexgenome_settings)
            report.write(res, file_type="genomexgenome")
        except ValidationError as e:
            logger.warning(f"Inter-genome validation failed: {e}")

    # Validate plasmid genomes (optional — non-fatal)
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

    # Validate reads (required — fatal on failure)
    reads_res = validate_reads(config.reads, reads_settings)
    report.write(reads_res, file_type="read")

    # Add interread validation
    try:
        readxread_res = readxread_validation(reads_res, readxread_settings)
        report.write(readxread_res, file_type="readxread")
    except ValidationError as e:
        logger.warning(f"Inter-read validation failed: {e}")

    # Validate features (optional — non-fatal)
    if hasattr(config, 'ref_feature') and config.ref_feature:
        try:
            res = validate_feature(config.ref_feature, ref_feature_settings)
            report.write(res, file_type="feature")
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
    # Step 4: Generate nextflow_schema.json with validation results
    # ========================================================================
    validation_results = {
        "ref_genome": ref_genome_res,
        "mod_genome": mod_genome_res,
    }
    schema = nf_schema.build_schema(validation_results)
    project_root = Path(__file__).resolve().parent.parent.parent
    nf_schema.write_schema(schema, project_root / "nextflow_schema.json")

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
