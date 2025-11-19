#!/usr/bin/env python3

import sys
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

    # Setup logging,
    setup_logging(console_level='DEBUG',log_file=Path("/EFSA_workspace/data/outputs/validation.log"))

    # ========================================================================
    # Step 1: Read and validate config
    # ========================================================================
    config_path = sys.argv[1]
    config = ConfigManager.load(config_path)

    # ========================================================================
    # Step 2: Edit settings for each validator
    # ========================================================================

    # Settings for reference genome
    ref_genome_settings = GenomeValidator.Settings()
    ref_genome_settings = ref_genome_settings.update(
        plasmids_to_one=True,
        main_longest=True,
        coding_type=None,
        output_filename_suffix='ref',
        replace_id_with='chr',
        min_sequence_length=100
    )

    # Settings for modified genome
    mod_genome_settings = GenomeValidator.Settings()
    mod_genome_settings = mod_genome_settings.update(
        plasmids_to_one=True,
        main_longest=True,
        coding_type=None,
        output_filename_suffix='mod',
        replace_id_with='chr',
        min_sequence_length=100,
    )

    # Settings for plasmid genomes (if you have them)
    plasmid_settings = GenomeValidator.Settings()
    plasmid_settings = plasmid_settings.update(
        is_plasmid=True,
        plasmids_to_one=True,
        coding_type=None,
        output_filename_suffix='plasmid'
    )
    
    # Settings for reads
    reads_settings = ReadValidator.Settings()
    reads_settings = reads_settings.update(
        coding_type='gz',
        outdir_by_ngs_type = True,
    )

    # Settings for features
    ref_feature_settings = FeatureValidator.Settings()
    ref_feature_settings = ref_feature_settings.update(
        sort_by_position=False,
        check_coordinates=False,
        replace_id_with='chr',
        coding_type=None,
        output_filename_suffix='ref'
    )

        # Settings for features
    
    mod_feature_settings = FeatureValidator.Settings()
    mod_feature_settings = mod_feature_settings.update(
        sort_by_position=False,
        check_coordinates=False,
        replace_id_with='chr',
        coding_type=None,
        output_filename_suffix='mod'
    )
    
    # Inter genome validation settings
    genomexgenome_settings = GenomeXGenomeSettings()
    genomexgenome_settings.update(
        same_number_of_sequences = True
    )

    # Inter read validation settings
    readxread_settings = ReadXReadSettings()
    readxread_settings.update(
        pair_end_basename = True
    )

    # ========================================================================
    # Step 3: Run validation using functional API
    # ========================================================================
    report = ValidationReport(Path("logs/report.txt"))

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

    report.flush()
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as e:
        print(f"\n✗ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
