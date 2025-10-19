#!/usr/bin/env python3
"""
Simple Usage Script - Exactly matching your desired API
========================================================

This script demonstrates the EXACT usage pattern you proposed:

    # read and valid config
    config_path = sys.argv[1]
    config = ConfigManager.load(config_path)

    # edit settings
    ref_genome_settings = GenomeValidator.Settings()
    ref_genome_settings = ref_genome_settings.update(<my settings>)
    ...same for mod_genome, plasmid_genomes, features and reads

    # run validation
    validate_genome(config.ref_genome, config.output_dir, ref_genome_settings)
    ...same for others

Usage:
    python simple_usage.py config.json
"""

import sys
from validation_pkg import (
    ConfigManager,
    GenomeValidator,
    ReadValidator,
    FeatureValidator,
    validate_genome,
    validate_reads,
    validate_feature,
    setup_logging
)


def main():
    # Check command line arguments
    if len(sys.argv) < 2:
        print("Usage: python simple_usage.py <config_path>")
        print("\nExample:")
        print("  python simple_usage.py config.json")
        return 1

    # ========================================================================
    # Step 1: Read and validate config
    # ========================================================================
    config_path = sys.argv[1]
    config = ConfigManager.load(config_path)

    print(f"Configuration loaded from: {config_path}")
    print(f"Output directory: {config.output_dir}\n")

    # Setup logging
    setup_logging(console_level='INFO')

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
        outdir_by_ngs_type=True
    )

    # Settings for features
    ref_feature_settings = FeatureValidator.Settings()
    ref_feature_settings = ref_feature_settings.update(
        sort_by_position=True,
        check_coordinates=True,
        allow_zero_length=False,
        replace_id_with='chr',
        coding_type=None,
        output_filename_suffix='ref'
    )

        # Settings for features
    mod_feature_settings = FeatureValidator.Settings()
    mod_feature_settings = mod_feature_settings.update(
        sort_by_position=True,
        check_coordinates=True,
        allow_zero_length=False,
        replace_id_with='chr',
        coding_type=None,
        output_filename_suffix='mod'
    )
    # ========================================================================
    # Step 3: Run validation using functional API
    # ========================================================================

    print("="*70)
    print("Starting Validation Workflow")
    print("="*70)

    # Validate reference genome
    if config.ref_genome:
        print(f"\n[1/4] Validating reference genome: {config.ref_genome.filename}")
        stats = validate_genome(config.ref_genome, config.output_dir, ref_genome_settings)

    # Validate modified genome
    if config.mod_genome:
        print(f"\n[2/4] Validating modified genome: {config.mod_genome.filename}")
        stats = validate_genome(config.mod_genome, config.output_dir, mod_genome_settings)

    # Validate plasmid genomes (if present in config)
    if hasattr(config, 'ref_plasmid') and config.ref_plasmid:
        print(f"\n[2.1/4] Validating reference plasmid: {config.ref_plasmid.filename}")
        stats = validate_genome(config.ref_plasmid, config.output_dir, plasmid_settings)

    if hasattr(config, 'mod_plasmid') and config.mod_plasmid:
        print(f"\n[2.2/4] Validating modified plasmid: {config.mod_plasmid.filename}")
        stats = validate_genome(config.mod_plasmid, config.output_dir, plasmid_settings)

    # Validate reads
    if config.reads:
        print(f"\n[3/4] Validating {len(config.reads)} read file(s)...")
        stats_list = validate_reads(config.reads, config.output_dir, reads_settings)

    # Validate features
    if config.ref_feature:
        print(f"\n[4/4] Validating feature file: {config.ref_feature.filename}")
        stats = validate_feature(config.ref_feature, config.output_dir, ref_feature_settings)

    if config.mod_feature:
        print(f"\n[4/4] Validating feature file: {config.mod_feature.filename}")
        stats = validate_feature(config.mod_feature, config.output_dir, mod_feature_settings)

    # ========================================================================
    # Done!
    # ========================================================================

    print("\n" + "="*70)
    print("✓ All validations completed successfully!")
    print("="*70)
    print(f"\nOutput files saved to: {config.output_dir}")

    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as e:
        print(f"\n✗ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
