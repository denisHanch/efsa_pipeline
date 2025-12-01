"""
Tests for validation report generation module.

Tests cover:
- FileValidationRecord and InterFileValidationRecord dataclasses
- ValidationReport construction and data collection
- Text and JSON report generation
- Metadata formatting and handling
- Error handling and edge cases
"""

import pytest
import json
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass, field
from typing import List

from validation_pkg.report import (
    ValidationReport,
    FileValidationRecord,
    InterFileValidationRecord
)


# ============================================================================
# Test Fixtures
# ============================================================================

@pytest.fixture
def temp_report_path(tmp_path):
    """Provide a temporary path for report output."""
    return tmp_path / "test_report.txt"


@pytest.fixture
def sample_genome_metadata():
    """Sample genome validation output metadata (dict format)."""
    return {
        'output_file': '/output/genome.fasta.gz',
        'output_filename': 'genome.fasta.gz',
        'num_sequences': 3,
        'total_length': 5000000,
        'sequence_ids': ['chr1', 'chr2', 'chr3'],
        'sequence_lengths': [2000000, 1800000, 1200000],
        'gc_content': 42.5,
        'validation_level': 'strict'
    }


@pytest.fixture
def sample_read_metadata():
    """Sample read validation output metadata (dict format)."""
    return {
        'output_file': '/output/reads_R1.fastq.gz',
        'output_filename': 'reads_R1.fastq.gz',
        'num_sequences': 1000000,
        'total_length': 150000000,
        'base_name': 'reads',
        'read_number': 1,
        'ngs_type_detected': 'illumina',
        'validation_level': 'trust'
    }


@pytest.fixture
def sample_feature_metadata():
    """Sample feature validation output metadata (dict format)."""
    return {
        'output_file': '/output/features.gff.gz',
        'output_filename': 'features.gff.gz',
        'num_features': 5000,
        'feature_types': ['gene', 'mRNA', 'exon', 'CDS'],
        'sequence_ids': ['chr1', 'chr2'],
        'validation_level': 'strict'
    }


@pytest.fixture
def sample_genome_settings():
    """Sample genome validator settings (as dict with to_dict method)."""
    # Create a mock settings object that has to_dict() method
    class MockSettings:
        def to_dict(self):
            return {
                'validation_level': 'strict',
                'plasmid_split': True,
                'main_longest': True,
                'min_sequence_length': 1000,
                'check_duplicate_ids': True,
                'coding_type': 'gz'
            }
    return MockSettings()


@pytest.fixture
def sample_interfile_result():
    """Sample inter-file validation result."""
    return {
        'passed': True,
        'errors': [],
        'warnings': ['Genome length mismatch: expected 5000000, got 5000100'],
        'metadata': {
            'reference_sequences': 3,
            'modified_sequences': 3,
            'matching_sequences': 3
        }
    }


# ============================================================================
# Test FileValidationRecord
# ============================================================================

class TestFileValidationRecord:
    """Test FileValidationRecord dataclass with simplified API."""

    def test_record_creation_minimal(self):
        """Test creating record with minimal required fields."""
        output_data = {
            'input_file': 'genome.fasta',
            'output_file': '/output/genome.fasta.gz'
        }
        record = FileValidationRecord(
            output_data=output_data,
            validator_type='genome'
        )

        assert record.output_data == output_data
        assert record.validator_type == 'genome'
        assert record.input_settings is None

    def test_record_creation_complete(self, sample_genome_metadata, sample_genome_settings):
        """Test creating record with all fields."""
        record = FileValidationRecord(
            output_data=sample_genome_metadata,
            validator_type='genome',
            input_settings=sample_genome_settings.to_dict()
        )

        assert record.output_data == sample_genome_metadata
        assert record.validator_type == 'genome'
        assert record.input_settings == sample_genome_settings.to_dict()

    def test_record_with_dict_settings(self):
        """Test creating record with dict settings."""
        output_data = {'input_file': 'test.fasta', 'output_file': '/output/test.fasta.gz'}
        settings = {'validation_level': 'strict', 'threads': 4}

        record = FileValidationRecord(
            output_data=output_data,
            validator_type='genome',
            input_settings=settings
        )

        assert record.input_settings == settings


# ============================================================================
# Test InterFileValidationRecord
# ============================================================================

class TestInterFileValidationRecord:
    """Test InterFileValidationRecord dataclass."""

    def test_record_creation_passed(self):
        """Test creating passing inter-file validation record."""
        record = InterFileValidationRecord(
            validation_type='genomexgenome',
            status='PASSED',
            passed=True
        )

        assert record.validation_type == 'genomexgenome'
        assert record.status == 'PASSED'
        assert record.passed is True
        assert record.errors == []
        assert record.warnings == []
        assert record.metadata == {}
        assert isinstance(record.timestamp, str)

    def test_record_creation_failed(self):
        """Test creating failed inter-file validation record."""
        errors = ['Sequence count mismatch', 'ID mismatch: chr1 != contig1']
        warnings = ['Different sequence lengths']
        metadata = {'ref_seqs': 3, 'mod_seqs': 2}

        record = InterFileValidationRecord(
            validation_type='genomexgenome',
            status='FAILED',
            passed=False,
            errors=errors,
            warnings=warnings,
            metadata=metadata
        )

        assert record.validation_type == 'genomexgenome'
        assert record.status == 'FAILED'
        assert record.passed is False
        assert record.errors == errors
        assert record.warnings == warnings
        assert record.metadata == metadata


# ============================================================================
# Test ValidationReport Construction
# ============================================================================

class TestValidationReportConstruction:
    """Test ValidationReport initialization and setup."""

    def test_initialization(self, temp_report_path):
        """Test report initialization."""
        report = ValidationReport(temp_report_path)

        assert report.report_path == temp_report_path
        assert report.file_records == []
        assert report.interfile_records == []
        assert isinstance(report.start_time, datetime)

    def test_parent_directory_creation(self, tmp_path):
        """Test that parent directories are created."""
        nested_path = tmp_path / "reports" / "subdir" / "report.txt"
        report = ValidationReport(nested_path)

        assert nested_path.parent.exists()

    def test_path_conversion(self, tmp_path):
        """Test that string paths are converted to Path objects."""
        path_str = str(tmp_path / "report.txt")
        report = ValidationReport(path_str)

        assert isinstance(report.report_path, Path)
        assert str(report.report_path) == path_str


# ============================================================================
# Test Report Data Collection - File Validation
# ============================================================================

class TestFileValidationRecording:
    """Test adding file validation results to report."""

    def test_write_genome_result_dict(self, temp_report_path, sample_genome_metadata, sample_genome_settings):
        """Test recording genome validation result (dict format)."""
        report = ValidationReport(temp_report_path)

        # Add input_file to metadata
        output_data = sample_genome_metadata.copy()
        output_data['input_file'] = 'genome.fasta'

        report.write(
            output_data,
            file_type='genome',
            settings=sample_genome_settings
        )

        assert len(report.file_records) == 1
        record = report.file_records[0]
        assert record.output_data['input_file'] == 'genome.fasta'
        assert record.validator_type == 'genome'
        assert record.input_settings == sample_genome_settings.to_dict()
        assert record.output_data['output_file'] == sample_genome_metadata['output_file']

    def test_write_read_result_dict(self, temp_report_path, sample_read_metadata):
        """Test recording read validation result (dict format)."""
        report = ValidationReport(temp_report_path)

        output_data = sample_read_metadata.copy()
        output_data['input_file'] = 'reads.fastq'

        report.write(output_data, file_type='read')

        assert len(report.file_records) == 1
        assert report.file_records[0].validator_type == 'read'

    def test_write_feature_result_dict(self, temp_report_path, sample_feature_metadata):
        """Test recording feature validation result (dict format)."""
        report = ValidationReport(temp_report_path)

        output_data = sample_feature_metadata.copy()
        output_data['input_file'] = 'features.gff'

        report.write(output_data, file_type='feature')

        assert len(report.file_records) == 1
        assert report.file_records[0].validator_type == 'feature'

    def test_write_multiple_results(self, temp_report_path, sample_genome_metadata, sample_read_metadata):
        """Test recording multiple file validation results."""
        report = ValidationReport(temp_report_path)

        genome_result = sample_genome_metadata.copy()
        genome_result['input_file'] = 'genome.fasta'

        read_result = sample_read_metadata.copy()
        read_result['input_file'] = 'reads.fastq'

        report.write(genome_result, file_type='genome')
        report.write(read_result, file_type='read')

        assert len(report.file_records) == 2
        assert report.file_records[0].validator_type == 'genome'
        assert report.file_records[1].validator_type == 'read'

    def test_write_result_list(self, temp_report_path, sample_read_metadata):
        """Test recording list of results (multiple read files)."""
        report = ValidationReport(temp_report_path)

        results = [
            {**sample_read_metadata, 'input_file': 'reads1.fastq', 'output_file': '/output/reads1_R1.fastq.gz'},
            {**sample_read_metadata, 'input_file': 'reads2.fastq', 'output_file': '/output/reads2_R1.fastq.gz'}
        ]

        report.write(results, file_type='read')

        assert len(report.file_records) == 2

    def test_write_without_settings(self, temp_report_path, sample_genome_metadata):
        """Test recording result without settings object."""
        report = ValidationReport(temp_report_path)

        output_data = sample_genome_metadata.copy()
        output_data['input_file'] = 'genome.fasta'

        report.write(output_data, file_type='genome')

        assert len(report.file_records) == 1
        assert report.file_records[0].input_settings is None

    def test_write_infers_input_filename(self, temp_report_path, sample_genome_metadata):
        """Test that output_data is stored correctly."""
        report = ValidationReport(temp_report_path)

        output_data = sample_genome_metadata.copy()
        output_data['input_file'] = 'genome.fasta'

        report.write(output_data, file_type='genome')

        assert len(report.file_records) == 1
        assert report.file_records[0].output_data['input_file'] == 'genome.fasta'


# ============================================================================
# Test Report Data Collection - Inter-file Validation
# ============================================================================

class TestInterFileValidationRecording:
    """Test adding inter-file validation results to report."""

    def test_write_genomexgenome_passed(self, temp_report_path):
        """Test recording passing genome inter-file validation."""
        report = ValidationReport(temp_report_path)

        result = {
            'passed': True,
            'errors': [],
            'warnings': [],
            'metadata': {'matching_sequences': 3}
        }

        report.write(result, file_type='genomexgenome')

        assert len(report.interfile_records) == 1
        record = report.interfile_records[0]
        assert record.validation_type == 'genomexgenome'
        assert record.status == 'PASSED'
        assert record.passed is True
        assert record.errors == []

    def test_write_genomexgenome_failed(self, temp_report_path):
        """Test recording failed genome inter-file validation."""
        report = ValidationReport(temp_report_path)

        result = {
            'passed': False,
            'errors': ['Sequence count mismatch: 3 vs 2'],
            'warnings': ['Length difference detected'],
            'metadata': {'ref_seqs': 3, 'mod_seqs': 2}
        }

        report.write(result, file_type='genomexgenome')

        assert len(report.interfile_records) == 1
        record = report.interfile_records[0]
        assert record.validation_type == 'genomexgenome'
        assert record.status == 'FAILED'
        assert record.passed is False
        assert len(record.errors) == 1
        assert len(record.warnings) == 1

    def test_write_readxread(self, temp_report_path):
        """Test recording read inter-file validation."""
        report = ValidationReport(temp_report_path)

        result = {
            'passed': True,
            'errors': [],
            'warnings': [],
            'metadata': {'paired_files': 2}
        }

        report.write(result, file_type='readxread')

        assert len(report.interfile_records) == 1
        assert report.interfile_records[0].validation_type == 'readxread'

    def test_write_featurexgenome(self, temp_report_path):
        """Test recording feature-genome inter-file validation."""
        report = ValidationReport(temp_report_path)

        result = {
            'passed': False,
            'errors': ['Sequence chr4 in features not found in genome'],
            'warnings': [],
            'metadata': {}
        }

        report.write(result, file_type='featurexgenome')

        assert len(report.interfile_records) == 1
        assert report.interfile_records[0].validation_type == 'featurexgenome'
        assert len(report.interfile_records[0].errors) == 1


# ============================================================================
# Test Error Handling
# ============================================================================

class TestErrorHandling:
    """Test error handling in ValidationReport."""

    def test_invalid_file_type(self, temp_report_path):
        """Test that invalid file_type raises ValueError."""
        report = ValidationReport(temp_report_path)

        with pytest.raises(ValueError, match="Unknown file_type"):
            report.write({'output_file': 'test.txt'}, file_type='invalid_type')

    def test_invalid_flush_format(self, temp_report_path):
        """Test that invalid flush format raises ValueError."""
        report = ValidationReport(temp_report_path)

        with pytest.raises(ValueError, match="Unknown format"):
            report.flush(format='xml')


# ============================================================================
# Test Text Report Generation
# ============================================================================

class TestTextReportGeneration:
    """Test text format report generation."""

    def test_flush_text_creates_file(self, temp_report_path):
        """Test that flush creates the report file."""
        report = ValidationReport(temp_report_path)
        report.flush(format='text')

        assert temp_report_path.exists()

    def test_text_report_basic_structure(self, temp_report_path, sample_genome_metadata):
        """Test that text report has expected structure."""
        report = ValidationReport(temp_report_path)

        output_data = sample_genome_metadata.copy()
        output_data['input_file'] = 'genome.fasta'

        report.write(output_data, file_type='genome')

        report.flush(format='text')

        content = temp_report_path.read_text()
        assert "VALIDATION PIPELINE REPORT" in content
        assert "SUMMARY" in content
        assert "Generated:" in content
        assert "Total Duration:" in content

    def test_text_report_file_counts(self, temp_report_path, sample_genome_metadata, sample_read_metadata):
        """Test that summary section shows correct file counts."""
        report = ValidationReport(temp_report_path)

        genome_result = sample_genome_metadata.copy()
        genome_result['input_file'] = 'genome.fasta'

        read_result = sample_read_metadata.copy()
        read_result['input_file'] = 'reads.fastq'

        report.write(genome_result, file_type='genome')
        report.write(read_result, file_type='read')

        report.flush(format='text')

        content = temp_report_path.read_text()
        assert "Files Processed: 2" in content
        assert "Genomes:  1" in content
        assert "Reads:    1" in content

    def test_text_report_interfile_status(self, temp_report_path):
        """Test that inter-file validation status is shown."""
        report = ValidationReport(temp_report_path)

        result = {
            'passed': False,
            'errors': ['Test error'],
            'warnings': [],
            'metadata': {}
        }
        report.write(result, file_type='genomexgenome')

        report.flush(format='text')

        content = temp_report_path.read_text()
        assert "Overall Status: ✗ FAILED" in content
        assert "Inter-file Validations: 1" in content
        assert "Failed:  1" in content

    def test_text_report_empty(self, temp_report_path):
        """Test generating text report with no results."""
        report = ValidationReport(temp_report_path)
        report.flush(format='text')

        content = temp_report_path.read_text()
        assert "VALIDATION PIPELINE REPORT" in content
        assert "Files Processed: 0" in content


# ============================================================================
# Test JSON Report Generation
# ============================================================================

class TestJSONReportGeneration:
    """Test JSON format report generation."""

    def test_flush_json_creates_file(self, temp_report_path):
        """Test that flush creates JSON report file."""
        report = ValidationReport(temp_report_path)
        report.flush(format='json')

        json_path = temp_report_path.with_suffix('.json')
        assert json_path.exists()

    def test_json_report_structure(self, temp_report_path):
        """Test that JSON report has expected structure."""
        report = ValidationReport(temp_report_path)
        report.flush(format='json')

        json_path = temp_report_path.with_suffix('.json')
        data = json.loads(json_path.read_text())

        assert 'report_metadata' in data
        assert 'summary' in data
        assert 'file_validations' in data
        assert 'inter_file_validations' in data

        assert 'generated' in data['report_metadata']
        assert 'duration_seconds' in data['report_metadata']
        assert 'version' in data['report_metadata']

    def test_json_report_file_validations(self, temp_report_path, sample_genome_metadata, sample_genome_settings):
        """Test that file validations are included in JSON."""
        report = ValidationReport(temp_report_path)

        output_data = sample_genome_metadata.copy()
        output_data['input_file'] = 'genome.fasta'

        report.write(
            output_data,
            file_type='genome',
            settings=sample_genome_settings
        )

        report.flush(format='json')

        json_path = temp_report_path.with_suffix('.json')
        data = json.loads(json_path.read_text())

        assert len(data['file_validations']) == 1
        validation = data['file_validations'][0]
        assert validation['validator_type'] == 'genome'
        assert validation['input_settings'] == sample_genome_settings.to_dict()
        assert validation['output_data']['input_file'] == 'genome.fasta'
        assert validation['output_data']['output_file'] == sample_genome_metadata['output_file']

    def test_json_report_interfile_validations(self, temp_report_path):
        """Test that inter-file validations are included in JSON."""
        report = ValidationReport(temp_report_path)

        result = {
            'passed': False,
            'errors': ['Error 1', 'Error 2'],
            'warnings': ['Warning 1'],
            'metadata': {'key': 'value'}
        }
        report.write(result, file_type='genomexgenome')

        report.flush(format='json')

        json_path = temp_report_path.with_suffix('.json')
        data = json.loads(json_path.read_text())

        assert len(data['inter_file_validations']) == 1
        validation = data['inter_file_validations'][0]
        assert validation['validation_type'] == 'genomexgenome'
        assert validation['status'] == 'FAILED'
        assert validation['passed'] is False
        assert len(validation['errors']) == 2
        assert len(validation['warnings']) == 1
        assert validation['metadata'] == {'key': 'value'}

    def test_json_report_valid_json(self, temp_report_path, sample_genome_metadata):
        """Test that generated JSON is valid and parseable."""
        report = ValidationReport(temp_report_path)

        output_data = sample_genome_metadata.copy()
        output_data['input_file'] = 'genome.fasta'

        report.write(output_data, file_type='genome')

        report.flush(format='json')

        json_path = temp_report_path.with_suffix('.json')
        # Should not raise exception
        data = json.loads(json_path.read_text())
        assert isinstance(data, dict)


# ============================================================================
# Test Metadata Handling
# ============================================================================

class TestMetadataHandling:
    """Test metadata conversion and formatting."""

    def test_output_data_stored_correctly(self, temp_report_path):
        """Test that output data is stored correctly in the new format."""
        report = ValidationReport(temp_report_path)

        output_data = {
            'input_file': 'test.fasta',
            'output_file': 'output/test.fasta',
            'num_sequences': 10,
            'optional_field': None
        }

        report.write(output_data, file_type='genome')

        assert len(report.file_records) == 1
        assert report.file_records[0].output_data == output_data

    def test_dict_format_with_settings(self, temp_report_path):
        """Test dict format with settings object."""
        report = ValidationReport(temp_report_path)

        output_data = {
            'input_file': 'genome.fasta',
            'output_file': '/output/genome.fasta.gz',
            'num_sequences': 3
        }

        settings_dict = {'validation_level': 'strict', 'threads': 4}

        report.write(output_data, file_type='genome', settings={'to_dict': lambda: settings_dict})

        assert len(report.file_records) == 1
        # Settings should be converted via to_dict if callable
        # In the actual implementation, write() handles both dict and object with to_dict()


# ============================================================================
# Test Edge Cases
# ============================================================================

class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_empty_metadata(self, temp_report_path):
        """Test handling of empty metadata."""
        report = ValidationReport(temp_report_path)

        output_data = {
            'input_file': 'test.fasta',
            'output_file': '/output/test.fasta'
        }

        report.write(output_data, file_type='genome')

        assert len(report.file_records) == 1
        assert report.file_records[0].output_data == output_data

    def test_missing_output_file(self, temp_report_path):
        """Test handling result without output_file."""
        report = ValidationReport(temp_report_path)

        output_data = {
            'input_file': 'test.fasta',
            'num_sequences': 10
        }

        report.write(output_data, file_type='genome')

        assert len(report.file_records) == 1
        # Should handle gracefully

    def test_unicode_in_filenames(self, temp_report_path):
        """Test handling of Unicode characters in filenames."""
        report = ValidationReport(temp_report_path)

        output_data = {
            'input_file': 'génome_données.fasta',
            'output_file': '/output/génome_données.fasta'
        }

        report.write(output_data, file_type='genome')

        assert len(report.file_records) == 1
        assert report.file_records[0].output_data['input_file'] == 'génome_données.fasta'

    def test_very_long_error_messages(self, temp_report_path):
        """Test handling of very long error messages."""
        report = ValidationReport(temp_report_path)

        long_error = 'Error: ' + 'x' * 10000
        result = {
            'passed': False,
            'errors': [long_error],
            'warnings': [],
            'metadata': {}
        }

        report.write(result, file_type='genomexgenome')

        assert len(report.interfile_records) == 1
        assert len(report.interfile_records[0].errors[0]) > 10000

    def test_large_number_of_records(self, temp_report_path, sample_genome_metadata):
        """Test handling large number of validation records."""
        report = ValidationReport(temp_report_path)

        # Add 100 file records
        for i in range(100):
            output_data = sample_genome_metadata.copy()
            output_data['input_file'] = f'file{i}.fasta'
            output_data['output_file'] = f'/output/file{i}.fasta'

            report.write(output_data, file_type='genome')

        assert len(report.file_records) == 100

        # Should be able to flush without issues
        report.flush(format='text')
        assert temp_report_path.exists()
