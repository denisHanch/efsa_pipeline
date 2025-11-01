#!/usr/bin/env python3

import sys
from validation_pkg import (
    ConfigManager,
    GenomeValidator,
    ReadValidator,
    FeatureValidator,
    validate_genome,
    validate_reads,
    validate_feature,
    validate_reads,
    setup_logging
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
    setup_logging(console_level='DEBUG',log_file=Path("logs/validation.log"),report_file=Path("logs/report.txt"))

    # ========================================================================
    # Step 1: Read and validate config
    # ========================================================================
    config_path = sys.argv[1]
    config = ConfigManager.load(config_path)
    print(f"Configuration loaded from: {config_path}")

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
    
    # ========================================================================
    # Step 3: Run validation using functional API
    # ========================================================================
    import time

    print("="*70)
    print("Starting Validation Workflow")
    print("="*70)

    # Validate reference genome
    start_time = time.time()
    print(f"\n[1/4] Validating reference genome: {config.ref_genome.filename}")
    validate_genome(config.ref_genome, ref_genome_settings)
    end_time = time.time()
    print(f"Execution time: {end_time - start_time:.6f} seconds")

    # Validate modified genome
    start_time = time.time()
    print(f"\n[2/4] Validating modified genome: {config.mod_genome.filename}")
    validate_genome(config.mod_genome, mod_genome_settings)
    end_time = time.time()
    print(f"Execution time: {end_time - start_time:.6f} seconds")

    # Validate plasmid genomes (if present in config)
    if hasattr(config, 'ref_plasmid') and config.ref_plasmid:
        start_time = time.time()
        print(f"\n[2.1/4] Validating reference plasmid: {config.ref_plasmid.filename}")
        validate_genome(config.ref_plasmid, plasmid_settings)
        end_time = time.time()
        print(f"Execution time: {end_time - start_time:.6f} seconds")

    if hasattr(config, 'mod_plasmid') and config.mod_plasmid:
        start_time = time.time()
        print(f"\n[2.2/4] Validating modified plasmid: {config.mod_plasmid.filename}")
        validate_genome(config.mod_plasmid, plasmid_settings)
        end_time = time.time()
        print(f"Execution time: {end_time - start_time:.6f} seconds")

    # Validate reads
    start_time = time.time()
    print(f"\n[3/4] Validating {len(config.reads)} read file(s)...")
    validate_reads(config.reads, reads_settings)
    end_time = time.time()
    print(f"Execution time: {end_time - start_time:.6f} seconds")

    # Validate features
    if config.ref_feature:
        start_time = time.time()
        print(f"\n[4/4] Validating feature file: {config.ref_feature.filename}")
        validate_feature(config.ref_feature, ref_feature_settings)
        end_time = time.time()
        print(f"Execution time: {end_time - start_time:.6f} seconds")

    if config.mod_feature:
        start_time = time.time()
        print(f"\n[4/4] Validating feature file: {config.mod_feature.filename}")
        validate_feature(config.mod_feature, mod_feature_settings)
        end_time = time.time()
        print(f"Execution time: {end_time - start_time:.6f} seconds")

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
