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
    validate_reads,
    setup_logging,
    get_logger,
    readxread_validation,
    genomexgenome_validation
)

from pathlib import Path

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
    log_file = output_dir / "validation.log"
    setup_logging(console_level='DEBUG', log_file=log_file)

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
        replace_id_with='chr',
        min_sequence_length=100
    )

    # Settings for modified genome
    mod_genome_settings = GenomeValidator.Settings(
        plasmids_to_one=True,
        main_longest=True,
        coding_type=None,
        output_filename_suffix='mod',
        replace_id_with='chr',
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
    genomexgenome_settings = GenomeXGenomeSettings()

    # Inter read validation settings (using defaults)
    readxread_settings = ReadXReadSettings()

    # ========================================================================
    # Step 3: Run validation using functional API
    # ========================================================================
    report = ValidationReport(output_dir / "report.txt")

    # Validate reference genome
    ref_genome_res = validate_genome(config.ref_genome, ref_genome_settings)
    report.write(ref_genome_res,file_type = "genome")

    # Validate modified genome
    mod_genome_res = validate_genome(config.mod_genome, mod_genome_settings)
    report.write(mod_genome_res,file_type = "genome")

    # Add intergenome validation
    res = genomexgenome_validation(ref_genome_res,mod_genome_res,genomexgenome_settings)
    report.write(res,file_type = "genomexgenome")

    # Validate plasmid genomes (if present in config)
    if hasattr(config, 'ref_plasmid') and config.ref_plasmid:
        res = validate_genome(config.ref_plasmid, plasmid_settings)
        report.write(res,file_type = "genome")

    if hasattr(config, 'mod_plasmid') and config.mod_plasmid:
        res = validate_genome(config.mod_plasmid, plasmid_settings)
        report.write(res,file_type = "genome")

    # Validate reads
    reads_res = validate_reads(config.reads, reads_settings)
    report.write(reads_res,file_type = "read")

    # Add interread validation
    readxread_res = readxread_validation(reads_res,readxread_settings)
    report.write(readxread_res,file_type = "readxread")

    # Validate features
    if config.ref_feature:
        res = validate_feature(config.ref_feature, ref_feature_settings)
        report.write(res,file_type = "feature")

    if config.mod_feature:
        res = validate_feature(config.mod_feature, mod_feature_settings)
        report.write(res,file_type = "feature")

    report.flush(format='text')
    print(f"Log file: {log_file}")
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as e:
        logger = get_logger()
        logger.error(f"✗ Fatal error: {e}")
        logger.debug(traceback.format_exc())
        if len(sys.argv) >= 2:
            log_file = Path(sys.argv[1]).resolve().parent.parent / "valid" / "validation.log"
            print(f"Log file: {log_file}")
        sys.exit(1)
