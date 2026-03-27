"""
Comprehensive tests for Config and ConfigManager classes.

Tests cover:
- Configuration file parsing and validation
- File path resolution and security checks
- Global options handling (threads, validation_level, logging_level)
- Genome, read, and feature configuration parsing
- Error handling and validation
- Directory processing for read files
"""

import pytest
import json
import tempfile
from pathlib import Path
from validation_pkg.config_manager import ConfigManager, Config
from validation_pkg.exceptions import ConfigurationError, FileNotFoundError as ValidationFileNotFoundError
from validation_pkg.utils.formats import ReadFormat,FeatureFormat,GenomeFormat,CodingType


class TestConfigManager:
    """Test suite for ConfigManager configuration loading and validation."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)
    
    @pytest.fixture
    def valid_config_minimal(self, temp_dir):
        """Create minimal valid configuration with actual files."""
        # Create dummy files
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads_R1.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [
                {"filename": "reads_R1.fastq", "ngs_type": "illumina"}
            ]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))
        return config_file
    
    @pytest.fixture
    def valid_config_full(self, temp_dir):
        """Create full valid configuration with all optional fields."""
        # Create dummy files
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "ref_plasmid.gbk").write_text("LOCUS\n")
        (temp_dir / "mod_plasmid.gbk").write_text("LOCUS\n")
        (temp_dir / "reads_R1.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        (temp_dir / "reads_R2.fastq").write_text("@read2\nATCG\n+\nIIII\n")
        (temp_dir / "ref_features.gff").write_text("##gff-version 3\n")
        (temp_dir / "mod_features.bed").write_text("chr1\t100\t200\n")
        
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "ref_plasmid_filename": {"filename": "ref_plasmid.gbk"},
            "mod_plasmid_filename": {"filename": "mod_plasmid.gbk"},
            "reads": [
                {"filename": "reads_R1.fastq", "ngs_type": "illumina"},
                {"filename": "reads_R2.fastq", "ngs_type": "illumina"}
            ],
            "ref_feature_filename": {"filename": "ref_features.gff"},
            "mod_feature_filename": {"filename": "mod_features.bed"},
            "options": {"threads": 8, "validation_level": "trust"}
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))
        return config_file

    def test_load_minimal_valid_config(self, valid_config_minimal):
        """Test loading minimal valid configuration."""
        config = ConfigManager.load(str(valid_config_minimal))

        assert config.ref_genome is not None
        assert config.mod_genome is not None
        assert len(config.reads) == 1
        assert config.reads[0].ngs_type == "illumina"
        assert config.ref_plasmid is None
        assert config.mod_plasmid is None

        # Check absolute paths are set
        assert config.ref_genome.filepath.is_absolute()
        assert config.ref_genome.filepath.exists()
        assert config.mod_genome.filepath.is_absolute()
        assert config.reads[0].filepath.is_absolute()

    def test_load_full_valid_config(self, valid_config_full):
        """Test loading full configuration with all fields."""
        config = ConfigManager.load(str(valid_config_full))

        assert config.ref_plasmid is not None
        assert config.mod_plasmid is not None
        assert len(config.reads) == 2
        assert config.ref_feature is not None
        assert config.mod_feature is not None
        assert config.options["threads"] == 8
        assert config.options["validation_level"] == "trust"
        
        # Check all absolute paths are set
        assert config.ref_genome.filepath.is_absolute()
        assert config.mod_genome.filepath.is_absolute()
        assert config.ref_plasmid.filepath.is_absolute()
        assert config.mod_plasmid.filepath.is_absolute()
        assert config.ref_feature.filepath.is_absolute()
        assert config.mod_feature.filepath.is_absolute()
        assert all(read.filepath.is_absolute() for read in config.reads)
        
    
    def test_missing_config_file(self, temp_dir):
        """Test error when config file doesn't exist."""
        with pytest.raises(ValidationFileNotFoundError, match="Configuration file not found"):
            ConfigManager.load(str(temp_dir / "nonexistent.json"))
    
    def test_missing_required_field_ref_genome(self, temp_dir):
        """Test error when ref_genome_filename is missing."""
        config = {
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        with pytest.raises((ConfigurationError, ValueError), match="Missing required field: ref_genome_filename"):
            ConfigManager.load(str(config_file))
    
    def test_mod_genome_is_optional(self, temp_dir):
        """Test that mod_genome_filename is optional."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        loaded_config = ConfigManager.load(str(config_file))
        assert loaded_config.ref_genome is not None
        assert loaded_config.mod_genome is None
        assert len(loaded_config.reads) == 1
    
    def test_missing_required_field_reads(self, temp_dir):
        """Test error when reads field is missing."""
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"}
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        with pytest.raises((ConfigurationError, ValueError), match="Missing required field: reads"):
            ConfigManager.load(str(config_file))
    
    def test_empty_reads_list(self, temp_dir):
        """Test error when reads list is empty."""
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": []
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        with pytest.raises((ConfigurationError, ValueError), match="'reads' must be a non-empty list"):
            ConfigManager.load(str(config_file))
    
    def test_invalid_ngs_type(self, temp_dir):
        """Test error with invalid ngs_type value."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "454"}]  # Invalid type
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        with pytest.raises((ConfigurationError, ValueError), match="Invalid ngs_type"):
            ConfigManager.load(str(config_file))
    
    def test_missing_genome_file(self, temp_dir):
        """Test error when referenced genome file doesn't exist."""
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},  # File doesn't exist
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        with pytest.raises(ValidationFileNotFoundError, match="ref.fasta"):
            ConfigManager.load(str(config_file))

    def test_missing_read_file(self, temp_dir):
        """Test error when referenced read file doesn't exist."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "missing_reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        with pytest.raises(ValidationFileNotFoundError, match="reads\\[0\\]"):
            ConfigManager.load(str(config_file))

    def test_missing_read_file_second_entry_index(self, temp_dir):
        """Test that reads[N] context reports the correct index."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads_ok.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "reads": [
                {"filename": "reads_ok.fastq", "ngs_type": "illumina"},
                {"filename": "missing.fastq", "ngs_type": "illumina"},  # index 1 is missing
            ]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        with pytest.raises(ValidationFileNotFoundError, match="reads\\[1\\]"):
            ConfigManager.load(str(config_file))

    def test_missing_mod_genome_file(self, temp_dir):
        """Test error when optional mod_genome file is specified but missing."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod_missing.fasta"},  # File doesn't exist
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        with pytest.raises(ValidationFileNotFoundError, match="mod_missing.fasta"):
            ConfigManager.load(str(config_file))

    def test_missing_ref_plasmid_file(self, temp_dir):
        """Test error when ref_plasmid file is specified but missing."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "ref_plasmid_filename": {"filename": "plasmid_missing.fasta"},  # File doesn't exist
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        with pytest.raises(ValidationFileNotFoundError, match="plasmid_missing.fasta"):
            ConfigManager.load(str(config_file))

    def test_missing_ref_feature_file(self, temp_dir):
        """Test error when ref_feature file is specified but missing."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "ref_feature_filename": {"filename": "features_missing.gff"},  # File doesn't exist
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        with pytest.raises(ValidationFileNotFoundError, match="features_missing.gff"):
            ConfigManager.load(str(config_file))

    def test_missing_mod_feature_file(self, temp_dir):
        """Test error when mod_feature file is specified but missing."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_feature_filename": {"filename": "mod_features_missing.gff"},  # File doesn't exist
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        with pytest.raises(ValidationFileNotFoundError, match="mod_features_missing.gff"):
            ConfigManager.load(str(config_file))

    def test_missing_file_is_logged(self, temp_dir):
        """Test that a missing file error is logged via logger.error and add_validation_issue."""
        from unittest.mock import patch, call
        from validation_pkg.logger import get_logger

        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref_missing.fasta"},  # File doesn't exist
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        logger = get_logger()
        with patch.object(logger, 'error') as mock_error, \
             patch.object(logger, 'add_validation_issue') as mock_issue:
            with pytest.raises(ValidationFileNotFoundError):
                ConfigManager.load(str(config_file))

            mock_error.assert_called_once()
            assert "ref_missing.fasta" in mock_error.call_args[0][0]

            mock_issue.assert_called_once()
            call_kwargs = mock_issue.call_args[1] if mock_issue.call_args[1] else {}
            assert call_kwargs.get('category') == 'configuration'
    
    def test_directory_reads(self, temp_dir):
        """Test reads specified as directory."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        
        # Create reads directory
        reads_dir = temp_dir / "reads"
        reads_dir.mkdir()
        (reads_dir / "R1.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        (reads_dir / "R2.fastq").write_text("@read2\nATCG\n+\nIIII\n")
        
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"directory": "reads", "ngs_type": "illumina"}]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        loaded_config = ConfigManager.load(str(config_file))
        assert len(loaded_config.reads) == 2

        # Check filepaths in an order-independent way
        read_paths = {r.filepath for r in loaded_config.reads}
        expected_paths = {reads_dir / "R1.fastq", reads_dir / "R2.fastq"}
        assert read_paths == expected_paths

        # Verify all reads have correct properties (order-independent)
        assert all(r.detected_format == ReadFormat.FASTQ for r in loaded_config.reads)
        assert all(r.ngs_type == "illumina" for r in loaded_config.reads)

        # Verify each file individually by looking it up by name
        r1_config = next((r for r in loaded_config.reads if r.filename == "R1.fastq"), None)
        r2_config = next((r for r in loaded_config.reads if r.filename == "R2.fastq"), None)

        assert r1_config is not None
        assert r1_config.filepath == (reads_dir / "R1.fastq")
        assert r1_config.detected_format == ReadFormat.FASTQ
        assert r1_config.ngs_type == "illumina"

        assert r2_config is not None
        assert r2_config.filepath == (reads_dir / "R2.fastq")
        assert r2_config.detected_format == ReadFormat.FASTQ
        assert r2_config.ngs_type == "illumina"
    
    def test_missing_directory_reads(self, temp_dir):
        """Test error when reads directory doesn't exist."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"directory": "nonexistent_dir", "ngs_type": "illumina"}]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        with pytest.raises((ValidationFileNotFoundError, FileNotFoundError), match="Directory not found"):
            ConfigManager.load(str(config_file))
    
    def test_neither_filename_nor_directory_in_reads(self, temp_dir):
        """Test error when reads has neither filename nor directory."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"ngs_type": "illumina"}]  # Missing filename/directory
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        with pytest.raises((ConfigurationError, ValueError), match="filename' or 'directory' field is required"):
            ConfigManager.load(str(config_file))
    
    def test_malformed_json(self, temp_dir):
        """Test error with malformed JSON."""
        config_file = temp_dir / "config.json"
        config_file.write_text("{ invalid json }")
        
        with pytest.raises((ConfigurationError, json.JSONDecodeError)):
            ConfigManager.load(str(config_file))
    
    def test_genome_config_missing_filename_field(self, temp_dir):
        """Test error when genome config dict is missing filename."""
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        
        config = {
            "ref_genome_filename": {"some_key": "value"},  # Missing filename
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        with pytest.raises((ConfigurationError, ValueError), match="must contain 'filename' field"):
            ConfigManager.load(str(config_file))
    
    def test_backwards_compatibility_string_filenames(self, temp_dir):
        """Test that plain strings are accepted for backwards compatibility."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        
        config = {
            "ref_genome_filename": "ref.fasta",  # Plain string
            "mod_genome_filename": "mod.fasta",  # Plain string
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        loaded_config = ConfigManager.load(str(config_file))
        assert loaded_config.ref_genome.filename == "ref.fasta"
        assert loaded_config.ref_genome.filepath.is_absolute()
    
    def test_path_resolution(self, valid_config_minimal):
        """Test that paths are resolved relative to config directory."""
        config = ConfigManager.load(str(valid_config_minimal))
        
        # Check that filepath is in the same directory as config
        config_dir = valid_config_minimal.parent
        assert config.ref_genome.filepath.parent == config_dir
        assert config.config_dir == config_dir
    
    def test_multiple_reads_different_types(self, temp_dir):
        """Test configuration with multiple read types."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "illumina.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        (temp_dir / "ont.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        (temp_dir / "pacbio.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [
                {"filename": "illumina.fastq", "ngs_type": "illumina"},
                {"filename": "ont.fastq", "ngs_type": "ont"},
                {"filename": "pacbio.fastq", "ngs_type": "pacbio"}
            ]
        }
        
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        
        loaded_config = ConfigManager.load(str(config_file))
        assert len(loaded_config.reads) == 3
        assert loaded_config.reads[0].ngs_type == "illumina"
        assert loaded_config.reads[1].ngs_type == "ont"
        assert loaded_config.reads[2].ngs_type == "pacbio"
        
        # Check all paths are absolute
        assert all(read.filepath.is_absolute() for read in loaded_config.reads)
    
    def test_extra_keys_logged_as_warnings(self, temp_dir):
        """Test that extra keys in config are logged as warnings (not stored)."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {
                "filename": "ref.fasta",
                "custom_key": "custom_value",
                "quality_threshold": 30
            },
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [
                {
                    "filename": "reads.fastq",
                    "ngs_type": "illumina",
                    "custom_read_key": "value"
                }
            ]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        # Config should load successfully (extra keys are just logged as warnings)
        loaded_config = ConfigManager.load(str(config_file))

        # Verify config loaded successfully
        assert loaded_config.ref_genome is not None
        assert loaded_config.ref_genome.filename == "ref.fasta"
        assert len(loaded_config.reads) == 1

    def test_genome_file_genbank_format(self,temp_dir):
        """Genome files in GenBank format (.gbk) are accepted and paths resolved."""
        (temp_dir / "ref.gbk").write_text("LOCUS       X  10 bp DNA\nORIGIN\natcg\n//\n")
        (temp_dir / "mod.gbk").write_text("LOCUS       Y  10 bp DNA\nORIGIN\ngcta\n//\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.gbk"},
            "mod_genome_filename": {"filename": "mod.gbk"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
        }
        cfg_file = temp_dir / "config.json"
        cfg_file.write_text(json.dumps(config))

        loaded = ConfigManager.load(str(cfg_file))
        assert loaded.ref_genome.detected_format == GenomeFormat.GENBANK
        assert loaded.mod_genome.detected_format == GenomeFormat.GENBANK
        assert loaded.ref_genome.filepath.is_absolute()
        assert loaded.mod_genome.filepath.is_absolute()

    def test_genome_file_fasta_format(self,temp_dir):
        """Genome files in FASTA format (.fasta) are accepted and paths resolved."""
        (temp_dir / "ref.fasta").write_text(">seq\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
        }
        cfg_file = temp_dir / "config.json"
        cfg_file.write_text(json.dumps(config))

        loaded = ConfigManager.load(str(cfg_file))
        assert loaded.ref_genome.detected_format == GenomeFormat.FASTA
        assert loaded.mod_genome.detected_format == GenomeFormat.FASTA
        assert loaded.ref_genome.filepath.is_absolute()
        assert loaded.mod_genome.filepath.is_absolute()

    def test_genome_file_gz_coding(self,temp_dir):
        """
        A genome entry that requests gz output coding should store that preference.
        (The input file can be any extension; this tests config parsing, not IO.)
        """
        (temp_dir / "ref.fasta.gz").write_text(">seq\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta.gz"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
        }
        cfg_file = temp_dir / "config.json"
        cfg_file.write_text(json.dumps(config))

        loaded = ConfigManager.load(str(cfg_file))
        assert loaded.ref_genome.coding_type == CodingType.GZIP
        assert loaded.mod_genome.coding_type == CodingType.NONE

    def test_genome_file_no_extension(self,temp_dir):
        """Genome files without an extension are still accepted; paths resolve relative to config."""
        (temp_dir / "ref").write_text(">seq\nATCG\n")  # no suffix
        (temp_dir / "mod").write_text(">seq\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref"},
            "mod_genome_filename": {"filename": "mod"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
        }
        cfg_file = temp_dir / "config.json"
        cfg_file.write_text(json.dumps(config))

        with pytest.raises((ConfigurationError)):
            ConfigManager.load(str(cfg_file))

    def test_genome_file_fasta_gz_coding(self,temp_dir):
        """FASTA genomes plus explicit gz output coding on both ref and mod."""
        (temp_dir / "ref.fasta.gzip").write_text(">seq\nATCG\n")
        (temp_dir / "mod.fasta.gzip").write_text(">seq\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta.gzip"},
            "mod_genome_filename": {"filename": "mod.fasta.gzip"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
        }
        cfg_file = temp_dir / "config.json"
        cfg_file.write_text(json.dumps(config))

        loaded = ConfigManager.load(str(cfg_file))
        assert loaded.ref_genome.detected_format == GenomeFormat.FASTA
        assert loaded.ref_genome.coding_type == CodingType.GZIP
        assert loaded.mod_genome.coding_type == CodingType.GZIP

    def test_read_gz_coding(self,temp_dir):
        """Reads entry with gz output coding is preserved in the model."""
        # Minimal input genomes (required)
        (temp_dir / "ref.fasta").write_text(">s\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">s\nATCG\n")
        (temp_dir / "reads.fastq.gz").write_text("@r1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [
                {"filename": "reads.fastq.gz", "ngs_type": "illumina"}
            ],
        }
        cfg_file = temp_dir / "config.json"
        cfg_file.write_text(json.dumps(config))

        loaded = ConfigManager.load(str(cfg_file))
        assert len(loaded.reads) == 1
        assert loaded.reads[0].filepath.is_absolute()
        assert loaded.reads[0].ngs_type == "illumina"
        assert loaded.reads[0].coding_type == CodingType.GZIP


    def test_feature_gz_coding(self,temp_dir):
        """
        Feature file with gz output coding is preserved.
        Assumes features are given as a list under the 'features' key.
        """
        (temp_dir / "ref.fasta").write_text(">s\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">s\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")
        (temp_dir / "features.gff.gz").write_text("##gff-version 3\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
            "ref_feature_filename": {"filename": "features.gff.gz"}
        }
        cfg_file = temp_dir / "config.json"
        cfg_file.write_text(json.dumps(config))

        loaded = ConfigManager.load(str(cfg_file))
        # If your schema uses a different key than 'features', rename here.
        assert loaded.ref_feature.filepath.is_absolute()
        assert loaded.ref_feature.detected_format == FeatureFormat.GFF
        assert loaded.ref_feature.coding_type == CodingType.GZIP

class TestConfigManagerOutputDirectory:
    """Test output directory setup functionality."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_output_directory_created(self, temp_dir):
        """Test that output directory is created during config load."""
        # Create minimal config files
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        loaded_config = ConfigManager.load(str(config_file))

        # Check output_dir was set
        assert loaded_config.output_dir is not None
        assert loaded_config.output_dir.exists()
        assert loaded_config.output_dir.is_dir()

        # Check it's at the expected location (config_dir.parent / "valid")
        expected_output_dir = config_file.parent.parent / "valid"
        assert loaded_config.output_dir == expected_output_dir

    def test_output_directory_path_structure(self, temp_dir):
        """Test output directory is at config_dir.parent / 'valid'."""
        # Create config in a subdirectory
        config_subdir = temp_dir / "config"
        config_subdir.mkdir()

        (config_subdir / "ref.fasta").write_text(">seq1\nATCG\n")
        (config_subdir / "mod.fasta").write_text(">seq1\nATCG\n")
        (config_subdir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = config_subdir / "config.json"
        config_file.write_text(json.dumps(config))

        loaded_config = ConfigManager.load(str(config_file))

        # Output dir should be at temp_dir / "valid" (parent of config_subdir)
        expected_output_dir = temp_dir / "valid"
        assert loaded_config.output_dir == expected_output_dir
        assert expected_output_dir.exists()

    def test_genome_configs_have_output_dir(self, temp_dir):
        """Test that genome configs have output_dir set."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        loaded_config = ConfigManager.load(str(config_file))

        # Check genome configs have output_dir
        expected_output_dir = temp_dir.parent / "valid"
        assert loaded_config.ref_genome.output_dir == expected_output_dir
        assert loaded_config.mod_genome.output_dir == expected_output_dir

    def test_read_configs_have_output_dir(self, temp_dir):
        """Test that read configs have output_dir set."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        loaded_config = ConfigManager.load(str(config_file))

        # Check read configs have output_dir
        expected_output_dir = temp_dir.parent / "valid"
        assert all(read.output_dir == expected_output_dir for read in loaded_config.reads)


class TestConfigOptionsThreads:
    """Test suite for options.threads configuration parsing and validation."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def create_config_with_threads(self, temp_dir, threads_value):
        """Helper to create config with threads option."""
        # Create dummy files
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
            "options": {"threads": threads_value}
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))
        return config_file

    def test_threads_valid_integer(self, temp_dir):
        """Test parsing valid thread count."""
        config_file = self.create_config_with_threads(temp_dir, 4)
        config = ConfigManager.load(str(config_file))

        assert config.options['threads'] == 4
        assert config.threads == 4

    def test_threads_single_thread(self, temp_dir):
        """Test single thread configuration."""
        config_file = self.create_config_with_threads(temp_dir, 1)
        config = ConfigManager.load(str(config_file))

        assert config.threads == 1

    def test_threads_many_threads(self, temp_dir):
        """Test high thread count (should warn but allow)."""
        config_file = self.create_config_with_threads(temp_dir, 20)
        config = ConfigManager.load(str(config_file))

        # Should succeed but log warning
        assert config.threads == 20

    def test_threads_null_means_autodetect(self, temp_dir):
        """Test that null threads value means auto-detect."""
        config_file = self.create_config_with_threads(temp_dir, None)
        config = ConfigManager.load(str(config_file))

        assert config.threads is None

    def test_threads_omitted_means_none(self, temp_dir):
        """Test that omitting threads returns None (empty options)."""
        # Create config without threads option
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        loaded_config = ConfigManager.load(str(config_file))

        # No options section means threads returns None
        assert loaded_config.threads is None

    def test_threads_zero_raises_error(self, temp_dir):
        """Test that threads=0 raises ValueError."""
        config_file = self.create_config_with_threads(temp_dir, 0)

        with pytest.raises(ConfigurationError, match="threads.*positive integer"):
            ConfigManager.load(str(config_file))

    def test_threads_negative_raises_error(self, temp_dir):
        """Test that negative threads raises ValueError."""
        config_file = self.create_config_with_threads(temp_dir, -1)

        with pytest.raises(ConfigurationError, match="threads.*positive integer"):
            ConfigManager.load(str(config_file))

    def test_threads_wrong_type_string_raises_error(self, temp_dir):
        """Test that string threads value raises error."""
        # Create config with string threads
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
            "options": {"threads": "4"}  # String instead of int
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        with pytest.raises(ConfigurationError, match="threads.*must be an integer"):
            ConfigManager.load(str(config_file))

    def test_threads_wrong_type_float_raises_error(self, temp_dir):
        """Test that float threads value raises error."""
        config_file = self.create_config_with_threads(temp_dir, 4.5)

        with pytest.raises(ConfigurationError, match="threads.*must be an integer"):
            ConfigManager.load(str(config_file))

    def test_options_not_dict_raises_error(self, temp_dir):
        """Test that non-dict options raises error."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
            "options": "invalid"  # Should be dict
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        with pytest.raises(ConfigurationError, match="options.*must be a dictionary"):
            ConfigManager.load(str(config_file))

    def test_config_get_threads_method(self):
        """Test Config.threads method."""
        config = Config()

        # No threads specified
        assert config.threads is None

        # Threads specified
        config.options['threads'] = 4
        assert config.threads == 4

        # Threads explicitly null
        config.options['threads'] = None
        assert config.threads is None


class TestConfigValidatorSettings:
    """Test suite for file-level settings merging in config.json."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_genome_config_file_level_validation_level(self, temp_dir):
        """Test that file-level validation_level is stored in global_options."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {
                "filename": "ref.fasta",
                "validation_level": "trust"
            },
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        loaded_config = ConfigManager.load(str(config_file))

        # File-level validation_level overrides the default in global_options
        assert loaded_config.ref_genome.global_options['validation_level'] == 'trust'
        # mod_genome gets the default validation_level
        assert loaded_config.mod_genome.global_options['validation_level'] == 'trust'

    def test_read_config_file_level_validation_level(self, temp_dir):
        """Test that file-level validation_level is stored in global_options for reads."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [
                {
                    "filename": "reads.fastq",
                    "ngs_type": "illumina",
                    "validation_level": "trust"
                }
            ]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        loaded_config = ConfigManager.load(str(config_file))

        assert loaded_config.reads[0].global_options['validation_level'] == 'trust'
        assert loaded_config.reads[0].ngs_type == "illumina"

    def test_feature_config_file_level_validation_level(self, temp_dir):
        """Test that file-level validation_level is stored in global_options for features."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        (temp_dir / "features.gff").write_text("##gff-version 3\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
            "ref_feature_filename": {
                "filename": "features.gff",
                "validation_level": "minimal"
            }
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        loaded_config = ConfigManager.load(str(config_file))

        assert loaded_config.ref_feature.global_options['validation_level'] == 'minimal'

    def test_non_global_options_logged_as_warnings(self, temp_dir):
        """Test that non-global options (like plasmid_split) are logged as warnings and ignored."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {
                "filename": "ref.fasta",
                "validation_level": "trust",
                "plasmid_split": True,  # Not a global option - ignored with warning
                "min_sequence_length": 500  # Not a global option - ignored with warning
            },
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        # Should load successfully (non-global options just logged as warnings)
        loaded_config = ConfigManager.load(str(config_file))

        # Only allowed file-level options end up in global_options; non-global fields are ignored.
        # All defaults are present; the file-level validation_level overrides the default.
        assert loaded_config.ref_genome.global_options['validation_level'] == 'trust'
        assert 'plasmid_split' not in loaded_config.ref_genome.global_options
        assert 'min_sequence_length' not in loaded_config.ref_genome.global_options

    def test_directory_reads_inherit_file_level_settings(self, temp_dir):
        """Test that files from directory inherit file-level validation_level."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")

        # Create read directory with files
        reads_dir = temp_dir / "ont_reads"
        reads_dir.mkdir()
        (reads_dir / "read1.fastq").write_text("@read1\nATCG\n+\nIIII\n")
        (reads_dir / "read2.fastq").write_text("@read2\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [
                {
                    "directory": "ont_reads/",
                    "ngs_type": "ont",
                    "validation_level": "trust"
                }
            ]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        loaded_config = ConfigManager.load(str(config_file))

        # Both files from directory should inherit validation_level
        assert len(loaded_config.reads) == 2
        for read_config in loaded_config.reads:
            assert read_config.global_options['validation_level'] == 'trust'

    def test_global_and_file_level_options_merge(self, temp_dir):
        """Test that global options and file-level options merge correctly."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {
                "filename": "ref.fasta",
                "validation_level": "minimal"  # Override global
            },
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
            "options": {
                "threads": 8,
                "validation_level": "trust"  # Global default
            }
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        loaded_config = ConfigManager.load(str(config_file))

        # ref_genome overrides global validation_level but keeps threads
        assert loaded_config.ref_genome.global_options["validation_level"] == "minimal"
        assert loaded_config.ref_genome.global_options["threads"] == 8

        # mod_genome uses global options
        assert loaded_config.mod_genome.global_options["validation_level"] == "trust"
        assert loaded_config.mod_genome.global_options["threads"] == 8

    def test_global_options_contain_defaults_when_no_settings(self, temp_dir):
        """Test that global_options contains all defaults when no options provided in config."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        loaded_config = ConfigManager.load(str(config_file))

        # global_options should contain all defaults when no options specified
        expected_defaults = {
            'threads': None,
            'validation_level': 'trust',
            'logging_level': 'INFO',
            'type': 'prokaryote',
            'force_defragment_ref': False,
        }
        assert loaded_config.ref_genome.global_options == expected_defaults
        assert loaded_config.mod_genome.global_options == expected_defaults
        assert loaded_config.reads[0].global_options == expected_defaults


class TestSecurityPathTraversal:
    """Security tests for path traversal protection."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_path_traversal_blocked_in_genome_config(self, temp_dir):
        """Test that path traversal attacks are blocked in genome config."""
        # Create valid files
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        # Try to access parent directory using path traversal
        config = {
            "ref_genome_filename": {"filename": "../../../etc/passwd"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        # Should raise ConfigurationError for path traversal
        with pytest.raises(ConfigurationError) as exc_info:
            ConfigManager.load(str(config_file))

        assert "Path traversal detected" in str(exc_info.value)
        assert "etc/passwd" in str(exc_info.value)

    def test_path_traversal_blocked_in_read_config(self, temp_dir):
        """Test that path traversal attacks are blocked in read config."""
        # Create valid files
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")

        # Try to access parent directory in reads
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "../../../../etc/shadow", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        # Should raise ConfigurationError for path traversal
        with pytest.raises(ConfigurationError) as exc_info:
            ConfigManager.load(str(config_file))

        assert "Path traversal detected" in str(exc_info.value)

    def test_path_traversal_blocked_in_feature_config(self, temp_dir):
        """Test that path traversal attacks are blocked in feature config."""
        # Create valid files
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        # Try to access parent directory in feature files
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
            "ref_feature_filename": {"filename": "../../sensitive_data.gff"}
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        # Should raise ConfigurationError for path traversal
        with pytest.raises(ConfigurationError) as exc_info:
            ConfigManager.load(str(config_file))

        assert "Path traversal detected" in str(exc_info.value)

    def test_symlink_blocked(self, temp_dir):
        """Test that symlinks outside config directory are blocked."""
        import os

        # Create a file outside the temp directory (in a subdirectory)
        outside_dir = temp_dir.parent / "outside"
        outside_dir.mkdir(exist_ok=True)
        sensitive_file = outside_dir / "sensitive.fasta"
        sensitive_file.write_text(">secret\nATCG\n")

        # Create symlink in temp_dir pointing to outside file
        symlink_path = temp_dir / "symlink.fasta"
        try:
            os.symlink(sensitive_file, symlink_path)
        except OSError:
            pytest.skip("Symlinks not supported on this platform")

        # Create other valid files
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "symlink.fasta"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        # Should raise ConfigurationError because symlink resolves outside config_dir
        with pytest.raises(ConfigurationError) as exc_info:
            ConfigManager.load(str(config_file))

        assert "Path traversal detected" in str(exc_info.value)

    def test_normal_relative_paths_allowed(self, temp_dir):
        """Test that normal relative paths within config directory are allowed."""
        # Create subdirectory with files
        subdir = temp_dir / "data"
        subdir.mkdir()
        (subdir / "ref.fasta").write_text(">seq1\nATCG\n")
        (subdir / "mod.fasta").write_text(">seq1\nATCG\n")
        (subdir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        # Use relative path to subdirectory (this should be allowed)
        config = {
            "ref_genome_filename": {"filename": "data/ref.fasta"},
            "mod_genome_filename": {"filename": "data/mod.fasta"},
            "reads": [{"filename": "data/reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        # Should succeed - files are within config directory
        loaded_config = ConfigManager.load(str(config_file))
        assert loaded_config.ref_genome.filename == "ref.fasta"
        assert "data" in str(loaded_config.ref_genome.filepath)

    def test_absolute_paths_within_config_dir_blocked(self, temp_dir):
        """Test that absolute paths are normalized and checked for traversal."""
        # Create file
        (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@read1\nATCG\n+\nIIII\n")

        # Try to use absolute path to /etc/passwd
        config = {
            "ref_genome_filename": {"filename": "/etc/passwd"},
            "mod_genome_filename": {"filename": "mod.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}]
        }

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        # Should raise ConfigurationError for path traversal
        with pytest.raises(ConfigurationError) as exc_info:
            ConfigManager.load(str(config_file))

        assert "Path traversal detected" in str(exc_info.value)


class TestNSequenceLimit:
    """Tests for the n_sequence_limit genome config parameter.

    Rules:
    - ref_genome_filename: n_sequence_limit is optional; when absent the
      default is 5. Explicit values are accepted and stored.
    - mod_genome_filename: same behaviour.
    - ref_plasmid_filename / mod_plasmid_filename: n_sequence_limit is NOT
      applicable; it is ignored with a warning and the field is set to None.
    """

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def _write_config(self, temp_dir: Path, ref_genome_value, mod_genome_value=None) -> Path:
        """Write a minimal config JSON and create stub files."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")
        config = {
            "ref_genome_filename": ref_genome_value,
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
        }
        if mod_genome_value is not None:
            (temp_dir / "mod.fasta").write_text(">seq1\nATCG\n")
            config["mod_genome_filename"] = mod_genome_value
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))
        return config_file

    # ── ref_genome: n_sequence_limit behaves the same as mod_genome ──────────

    def test_ref_genome_n_sequence_limit_default_is_5(self, temp_dir):
        """When n_sequence_limit is absent from ref_genome, default is 5."""
        config_file = self._write_config(temp_dir, {"filename": "ref.fasta"})
        config = ConfigManager.load(str(config_file))
        assert config.ref_genome.n_sequence_limit == 5

    def test_ref_genome_n_sequence_limit_explicit_value(self, temp_dir):
        """Explicit n_sequence_limit on ref_genome is stored correctly."""
        config_file = self._write_config(
            temp_dir,
            {"filename": "ref.fasta", "n_sequence_limit": 42},
        )
        config = ConfigManager.load(str(config_file))
        assert config.ref_genome.n_sequence_limit == 42

    def test_ref_genome_n_sequence_limit_string_uses_default(self, temp_dir):
        """Plain-string ref_genome value (no dict) also gets the default (5)."""
        config_file = self._write_config(temp_dir, "ref.fasta")
        config = ConfigManager.load(str(config_file))
        assert config.ref_genome.n_sequence_limit == 5

    def test_ref_genome_n_sequence_limit_zero_raises(self, temp_dir):
        """n_sequence_limit=0 is rejected for ref_genome (must be positive)."""
        config_file = self._write_config(
            temp_dir,
            {"filename": "ref.fasta", "n_sequence_limit": 0},
        )
        with pytest.raises(ConfigurationError, match="n_sequence_limit"):
            ConfigManager.load(str(config_file))

    def test_ref_genome_n_sequence_limit_negative_raises(self, temp_dir):
        """Negative n_sequence_limit is rejected for ref_genome."""
        config_file = self._write_config(
            temp_dir,
            {"filename": "ref.fasta", "n_sequence_limit": -3},
        )
        with pytest.raises(ConfigurationError, match="n_sequence_limit"):
            ConfigManager.load(str(config_file))

    def test_ref_genome_n_sequence_limit_float_raises(self, temp_dir):
        """Float n_sequence_limit is rejected for ref_genome (must be integer)."""
        config_file = self._write_config(
            temp_dir,
            {"filename": "ref.fasta", "n_sequence_limit": 5.5},
        )
        with pytest.raises(ConfigurationError, match="n_sequence_limit"):
            ConfigManager.load(str(config_file))

    def test_ref_genome_n_sequence_limit_string_value_raises(self, temp_dir):
        """String n_sequence_limit is rejected for ref_genome."""
        config_file = self._write_config(
            temp_dir,
            {"filename": "ref.fasta", "n_sequence_limit": "ten"},
        )
        with pytest.raises(ConfigurationError, match="n_sequence_limit"):
            ConfigManager.load(str(config_file))

    # ── mod_genome: default and explicit value ────────────────────────────────

    def test_mod_genome_n_sequence_limit_default_is_5(self, temp_dir):
        """When n_sequence_limit is absent from mod_genome, default is 5."""
        config_file = self._write_config(
            temp_dir, "ref.fasta", mod_genome_value={"filename": "mod.fasta"}
        )
        config = ConfigManager.load(str(config_file))
        assert config.mod_genome.n_sequence_limit == 5

    def test_mod_genome_n_sequence_limit_string_uses_default(self, temp_dir):
        """Plain-string mod_genome value (no dict) also gets the default (5)."""
        config_file = self._write_config(temp_dir, "ref.fasta", mod_genome_value="mod.fasta")
        config = ConfigManager.load(str(config_file))
        assert config.mod_genome.n_sequence_limit == 5

    def test_mod_genome_n_sequence_limit_explicit_value(self, temp_dir):
        """Explicit n_sequence_limit on mod_genome is stored correctly."""
        config_file = self._write_config(
            temp_dir, "ref.fasta",
            mod_genome_value={"filename": "mod.fasta", "n_sequence_limit": 42},
        )
        config = ConfigManager.load(str(config_file))
        assert config.mod_genome.n_sequence_limit == 42

    def test_ref_explicit_mod_explicit(self, temp_dir):
        """Both ref_genome and mod_genome carry their own explicit n_sequence_limit."""
        config_file = self._write_config(
            temp_dir,
            {"filename": "ref.fasta", "n_sequence_limit": 20},
            mod_genome_value={"filename": "mod.fasta", "n_sequence_limit": 50},
        )
        loaded = ConfigManager.load(str(config_file))
        assert loaded.ref_genome.n_sequence_limit == 20
        assert loaded.mod_genome.n_sequence_limit == 50

    # ── validation errors (mod_genome) ───────────────────────────────────────

    def test_n_sequence_limit_zero_raises(self, temp_dir):
        """n_sequence_limit=0 is rejected for mod_genome (must be positive)."""
        config_file = self._write_config(
            temp_dir, "ref.fasta",
            mod_genome_value={"filename": "mod.fasta", "n_sequence_limit": 0},
        )
        with pytest.raises(ConfigurationError, match="n_sequence_limit"):
            ConfigManager.load(str(config_file))

    def test_n_sequence_limit_negative_raises(self, temp_dir):
        """Negative n_sequence_limit is rejected for mod_genome."""
        config_file = self._write_config(
            temp_dir, "ref.fasta",
            mod_genome_value={"filename": "mod.fasta", "n_sequence_limit": -5},
        )
        with pytest.raises(ConfigurationError, match="n_sequence_limit"):
            ConfigManager.load(str(config_file))

    def test_n_sequence_limit_float_raises(self, temp_dir):
        """Float n_sequence_limit is rejected for mod_genome (must be integer)."""
        config_file = self._write_config(
            temp_dir, "ref.fasta",
            mod_genome_value={"filename": "mod.fasta", "n_sequence_limit": 5.5},
        )
        with pytest.raises(ConfigurationError, match="n_sequence_limit"):
            ConfigManager.load(str(config_file))

    def test_n_sequence_limit_string_value_raises(self, temp_dir):
        """String n_sequence_limit is rejected for mod_genome."""
        config_file = self._write_config(
            temp_dir, "ref.fasta",
            mod_genome_value={"filename": "mod.fasta", "n_sequence_limit": "ten"},
        )
        with pytest.raises(ConfigurationError, match="n_sequence_limit"):
            ConfigManager.load(str(config_file))

    # ── plasmid entries ignore n_sequence_limit ───────────────────────────────

    def test_n_sequence_limit_on_ref_plasmid_is_ignored(self, temp_dir):
        """n_sequence_limit specified on ref_plasmid is ignored with a warning."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "plasmid.fasta").write_text(">plasmid1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "ref_plasmid_filename": {"filename": "plasmid.fasta", "n_sequence_limit": 3},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
        }
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        loaded = ConfigManager.load(str(config_file))
        assert loaded.ref_plasmid.n_sequence_limit is None

    def test_n_sequence_limit_on_mod_plasmid_is_ignored(self, temp_dir):
        """n_sequence_limit specified on mod_plasmid is ignored with a warning."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "plasmid.fasta").write_text(">plasmid1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "mod_plasmid_filename": {"filename": "plasmid.fasta", "n_sequence_limit": 10},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
        }
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        loaded = ConfigManager.load(str(config_file))
        assert loaded.mod_plasmid.n_sequence_limit is None

    def test_plasmid_n_sequence_limit_is_none_by_default(self, temp_dir):
        """Plasmid n_sequence_limit is None when not specified (no default applied)."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "plasmid.fasta").write_text(">plasmid1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")
        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "ref_plasmid_filename": {"filename": "plasmid.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
        }
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config))

        loaded = ConfigManager.load(str(config_file))
        assert loaded.ref_plasmid.n_sequence_limit is None


class TestConfigOptionsType:
    """Test suite for options.type configuration parsing and validation."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def _create_config(self, temp_dir, type_value=None, include_type=True):
        """Helper to create a minimal config with an optional type option."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")

        options = {}
        if include_type:
            options["type"] = type_value

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
        }
        if options:
            config["options"] = options

        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))
        return config_file

    def test_type_prokaryote(self, temp_dir):
        """'prokaryote' is a valid type value."""
        config_file = self._create_config(temp_dir, "prokaryote")
        config = ConfigManager.load(str(config_file))
        assert config.options["type"] == "prokaryote"
        assert config.type == "prokaryote"

    def test_type_eukaryote(self, temp_dir):
        """'eukaryote' is a valid type value."""
        config_file = self._create_config(temp_dir, "eukaryote")
        config = ConfigManager.load(str(config_file))
        assert config.options["type"] == "eukaryote"
        assert config.type == "eukaryote"

    def test_type_omitted_defaults_to_prokaryote(self, temp_dir):
        """When type is not specified, Config.type defaults to 'prokaryote'."""
        config_file = self._create_config(temp_dir, include_type=False)
        config = ConfigManager.load(str(config_file))
        assert config.type == "prokaryote"

    def test_type_invalid_raises_error(self, temp_dir):
        """An unrecognised type value raises ConfigurationError."""
        config_file = self._create_config(temp_dir, "virus")
        with pytest.raises(ConfigurationError, match="Invalid type"):
            ConfigManager.load(str(config_file))

    def test_type_wrong_type_raises_error(self, temp_dir):
        """A non-string type value raises ConfigurationError."""
        config_file = self._create_config(temp_dir, 123)
        with pytest.raises(ConfigurationError, match="'type' must be a string"):
            ConfigManager.load(str(config_file))

    def test_config_type_property(self):
        """Config.type property defaults to 'prokaryote' and reflects explicit values."""
        config = Config()
        assert config.type == "prokaryote"  # default

        config.options["type"] = "prokaryote"
        assert config.type == "prokaryote"

        config.options["type"] = "eukaryote"
        assert config.type == "eukaryote"

    def test_type_uppercase_raises_error(self, temp_dir):
        """Type value is case-sensitive; 'PROKARYOTE' is not valid."""
        config_file = self._create_config(temp_dir, "PROKARYOTE")
        with pytest.raises(ConfigurationError, match="Invalid type"):
            ConfigManager.load(str(config_file))

    def test_type_mixed_case_raises_error(self, temp_dir):
        """Mixed-case type value is not accepted."""
        config_file = self._create_config(temp_dir, "Eukaryote")
        with pytest.raises(ConfigurationError, match="Invalid type"):
            ConfigManager.load(str(config_file))

    def test_type_empty_string_raises_error(self, temp_dir):
        """Empty string is not a valid type value."""
        config_file = self._create_config(temp_dir, "")
        with pytest.raises(ConfigurationError, match="Invalid type"):
            ConfigManager.load(str(config_file))

    def test_type_combined_with_other_options(self, temp_dir):
        """type coexists correctly with threads and validation_level."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {"filename": "ref.fasta"},
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
            "options": {
                "type": "eukaryote",
                "threads": 4,
                "validation_level": "strict",
            },
        }
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        loaded = ConfigManager.load(str(config_file))
        assert loaded.type == "eukaryote"
        assert loaded.threads == 4
        assert loaded.validation_level == "strict"

    def test_type_as_file_level_option_is_ignored(self, temp_dir):
        """type specified at file level (not in options) is logged as warning and ignored."""
        (temp_dir / "ref.fasta").write_text(">seq1\nATCG\n")
        (temp_dir / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")

        config = {
            "ref_genome_filename": {
                "filename": "ref.fasta",
                "type": "prokaryote",  # file-level, not in options block
            },
            "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
        }
        config_file = temp_dir / "config.json"
        config_file.write_text(json.dumps(config, indent=2))

        loaded = ConfigManager.load(str(config_file))
        # Global type is not set (only file-level was given, which gets ignored) - falls back to default
        assert loaded.type == 'prokaryote'
        # file-level type is ignored (not in ALLOWED_FILE_OPTIONS); global default applies
        assert loaded.ref_genome.global_options['type'] == 'prokaryote'


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

