"""
Basic tests for the logging functionality.
Tests core features without over-testing Python's logging module.
"""

import pytest
import tempfile
from pathlib import Path
from validation_pkg.logger import ValidationLogger, setup_logging, get_logger


class TestLogger:
    """Test suite for ValidationLogger."""
    
    @pytest.fixture(autouse=True)
    def reset_logger(self):
        """Reset logger state before each test."""
        logger = get_logger()
        logger.clear_issues()
        # Reset logger instance (structlog doesn't use handlers in the same way)
        logger.logger = None
        yield
        logger.clear_issues()
    
    def test_logger_singleton(self):
        """Test that logger follows singleton pattern."""
        logger1 = get_logger()
        logger2 = get_logger()
        assert logger1 is logger2
    
    def test_setup_logging_creates_log_file(self):
        """Test that setup_logging creates log file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = Path(tmpdir) / "test.log"
            
            logger = setup_logging(log_file=log_file)
            logger.info("Test message")
            
            assert log_file.exists()
            content = log_file.read_text()
            assert "Test message" in content
    
    def test_add_validation_issue_stores_issue(self):
        """Test that validation issues are stored correctly."""
        logger = get_logger()
        logger.setup()
        
        logger.add_validation_issue(
            level='ERROR',
            category='genome',
            message='Test error message',
            details={'file': 'test.fasta', 'line': 42}
        )
        
        assert len(logger.validation_issues) == 1
        issue = logger.validation_issues[0]
        assert issue['level'] == 'ERROR'
        assert issue['category'] == 'genome'
        assert issue['message'] == 'Test error message'
        assert issue['details']['file'] == 'test.fasta'
    
    def test_add_multiple_validation_issues(self):
        """Test adding multiple issues."""
        logger = get_logger()
        logger.setup()
        
        logger.add_validation_issue('ERROR', 'genome', 'Error 1')
        logger.add_validation_issue('WARNING', 'feature', 'Warning 1')
        logger.add_validation_issue('ERROR', 'read', 'Error 2')
        
        assert len(logger.validation_issues) == 3
    
    def test_get_summary_counts_issues(self):
        """Test that get_summary counts issues correctly."""
        logger = get_logger()
        logger.setup()
        
        logger.add_validation_issue('ERROR', 'genome', 'Error 1')
        logger.add_validation_issue('ERROR', 'genome', 'Error 2')
        logger.add_validation_issue('WARNING', 'feature', 'Warning 1')
        
        summary = logger.get_summary()
        
        assert summary['total_issues'] == 3
        assert summary['errors'] == 2
        assert summary['warnings'] == 1
        assert summary['passed'] == False
    
    def test_get_summary_passed_when_no_errors(self):
        """Test that passed=True when there are no errors."""
        logger = get_logger()
        logger.setup()
        
        logger.add_validation_issue('WARNING', 'feature', 'Warning only')
        
        summary = logger.get_summary()
        
        assert summary['errors'] == 0
        assert summary['passed'] == True
    
    def test_get_summary_empty(self):
        """Test get_summary with no issues."""
        logger = get_logger()
        logger.setup()
        
        summary = logger.get_summary()
        
        assert summary['total_issues'] == 0
        assert summary['errors'] == 0
        assert summary['warnings'] == 0
        assert summary['passed'] == True
    
    def test_generate_report_creates_file(self):
        """Test that generate_report creates report file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_file = Path(tmpdir) / "report.txt"
            
            logger = setup_logging(report_file=report_file)
            logger.add_validation_issue('ERROR', 'genome', 'Test error')
            
            report_content = logger.generate_report()
            
            assert report_file.exists()
            assert "BIOINFORMATICS VALIDATION REPORT" in report_content
            assert "Test error" in report_content
    
    def test_generate_report_contains_all_issues(self):
        """Test that report contains all logged issues."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_file = Path(tmpdir) / "report.txt"
            
            logger = setup_logging(report_file=report_file)
            logger.add_validation_issue('ERROR', 'genome', 'Error message')
            logger.add_validation_issue('WARNING', 'feature', 'Warning message')
            
            report_content = logger.generate_report()
            
            assert "Error message" in report_content
            assert "Warning message" in report_content
            assert "[ERROR]" in report_content
            assert "[WARNING]" in report_content
    
    def test_generate_report_shows_details(self):
        """Test that report includes issue details."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_file = Path(tmpdir) / "report.txt"
            
            logger = setup_logging(report_file=report_file)
            logger.add_validation_issue(
                level='ERROR',
                category='genome',
                message='Test error',
                details={'file': 'test.fasta', 'line': 42}
            )
            
            report_content = logger.generate_report()
            
            assert "test.fasta" in report_content
            assert "42" in report_content
    
    def test_generate_report_no_issues(self):
        """Test report when no issues found."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_file = Path(tmpdir) / "report.txt"
            
            logger = setup_logging(report_file=report_file)
            report_content = logger.generate_report()
            
            assert "No issues found" in report_content
            assert "all validations passed" in report_content.lower()
    
    def test_clear_issues(self):
        """Test that clear_issues removes all issues."""
        logger = get_logger()
        logger.setup()
        
        logger.add_validation_issue('ERROR', 'genome', 'Error 1')
        logger.add_validation_issue('WARNING', 'feature', 'Warning 1')
        
        assert len(logger.validation_issues) == 2
        
        logger.clear_issues()
        
        assert len(logger.validation_issues) == 0
    
    def test_log_file_directory_creation(self):
        """Test that log file parent directories are created."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = Path(tmpdir) / "nested" / "dir" / "test.log"
            
            logger = setup_logging(log_file=log_file)
            logger.info("Test")
            
            assert log_file.exists()
            assert log_file.parent.exists()
    
    def test_report_file_directory_creation(self):
        """Test that report file parent directories are created."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_file = Path(tmpdir) / "nested" / "dir" / "report.txt"
            
            logger = setup_logging(report_file=report_file)
            logger.generate_report()
            
            assert report_file.exists()
            assert report_file.parent.exists()


class TestParallelValidationLogging:
    """Test suite for parallel validation logging output."""

    @pytest.fixture(autouse=True)
    def reset_logger(self):
        """Reset logger state before each test."""
        logger = get_logger()
        logger.clear_issues()
        logger.logger = None
        yield
        logger.clear_issues()

    def _create_test_fastq(self, path, num_reads=1000):
        """Create a test FASTQ file"""
        with open(path, 'w') as f:
            for i in range(num_reads):
                f.write(f"@read{i}\n")
                f.write("ATCGATCGATCGATCG\n")
                f.write("+\n")
                f.write("IIIIIIIIIIIIIIII\n")

    def test_sequential_validation_logging(self, capsys):
        """Test that sequential validation logs correctly"""
        from validation_pkg.validators.read_validator import ReadValidator
        from validation_pkg.config_manager import ReadConfig
        from validation_pkg.utils.formats import CodingType, ReadFormat

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            fastq_file = tmpdir / "test.fastq"
            self._create_test_fastq(fastq_file, num_reads=100)

            # Create logger
            logger = setup_logging()

            config = ReadConfig(
                filename=fastq_file.name,
                filepath=fastq_file,
                ngs_type="illumina",
                coding_type=CodingType.NONE,
                detected_format=ReadFormat.FASTQ,
                output_dir=tmpdir,
                global_options={"threads": 1, "validation_level": "strict"}
            )

            settings = ReadValidator.Settings(check_invalid_chars=True)
            validator = ReadValidator(config, settings)

            # Validate (should use sequential)
            validator._parse_file()
            validator._validate_sequences()

            # Check output
            captured = capsys.readouterr()
            assert "Sequential validation" in captured.out or "Validating" in captured.err

    def test_parallel_validation_logging(self, capsys):
        """Test that parallel validation logs parallelism state"""
        from validation_pkg.validators.read_validator import ReadValidator
        from validation_pkg.config_manager import ReadConfig
        from validation_pkg.utils.formats import CodingType, ReadFormat

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            fastq_file = tmpdir / "test.fastq"
            self._create_test_fastq(fastq_file, num_reads=5000)

            # Create logger
            logger = setup_logging()

            config = ReadConfig(
                filename=fastq_file.name,
                filepath=fastq_file,
                ngs_type="illumina",
                coding_type=CodingType.NONE,
                detected_format=ReadFormat.FASTQ,
                output_dir=tmpdir,
                global_options={"threads": 4, "validation_level": "strict"}
            )

            settings = ReadValidator.Settings(check_invalid_chars=True)
            validator = ReadValidator(config, settings)

            # Validate (should use parallel)
            validator._parse_file()
            validator._validate_sequences()

            # Check output mentions parallel validation
            captured = capsys.readouterr()
            output = captured.out + captured.err
            assert "Parallel validation enabled" in output or "parallel" in output.lower()

    def test_parallel_logging_shows_worker_count(self, capsys):
        """Test that parallel logging shows worker count"""
        from validation_pkg.validators.read_validator import ReadValidator
        from validation_pkg.config_manager import ReadConfig
        from validation_pkg.utils.formats import CodingType, ReadFormat

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            fastq_file = tmpdir / "test.fastq"
            self._create_test_fastq(fastq_file, num_reads=2000)

            # Create logger
            logger = setup_logging()

            threads = 8
            config = ReadConfig(
                filename=fastq_file.name,
                filepath=fastq_file,
                ngs_type="illumina",
                coding_type=CodingType.NONE,
                detected_format=ReadFormat.FASTQ,
                output_dir=tmpdir,
                global_options={"threads": threads, "validation_level": "strict"}
            )

            settings = ReadValidator.Settings(check_invalid_chars=True)
            validator = ReadValidator(config, settings)

            # Validate
            validator._parse_file()
            validator._validate_sequences()

            # Check output shows worker count
            captured = capsys.readouterr()
            output = captured.out + captured.err
            assert f"{threads} worker" in output or f"threads={threads}" in output

    def test_trust_mode_shows_sequential(self, capsys):
        """Test that trust mode always shows sequential validation"""
        from validation_pkg.validators.read_validator import ReadValidator
        from validation_pkg.config_manager import ReadConfig
        from validation_pkg.utils.formats import CodingType, ReadFormat

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            fastq_file = tmpdir / "test.fastq"
            self._create_test_fastq(fastq_file, num_reads=100)

            # Create logger
            logger = setup_logging()

            config = ReadConfig(
                filename=fastq_file.name,
                filepath=fastq_file,
                ngs_type="illumina",
                coding_type=CodingType.NONE,
                detected_format=ReadFormat.FASTQ,
                output_dir=tmpdir,
                global_options={"threads": 8, "validation_level": "trust"}  # threads=8 but trust mode
            )

            settings = ReadValidator.Settings(check_invalid_chars=True)
            validator = ReadValidator(config, settings)

            # Validate (should use sequential despite threads=8)
            validator._parse_file()
            validator._validate_sequences()

            # Check output shows sequential
            captured = capsys.readouterr()
            output = captured.out + captured.err
            assert "sequential" in output.lower() or "Trust mode" in output


class TestConditionalProcessID:
    """Test that process_id only appears during parallel validation."""

    @pytest.fixture(autouse=True)
    def reset_logger(self):
        """Reset logger state before each test."""
        logger = get_logger()
        logger.clear_issues()
        logger.logger = None
        # Ensure parallel logging is disabled
        logger.disable_parallel_logging()
        yield
        logger.clear_issues()
        logger.disable_parallel_logging()

    def test_enable_disable_parallel_logging(self):
        """Test enable/disable parallel logging methods"""
        logger = setup_logging()

        # Initially disabled
        assert not logger.is_parallel_logging_enabled()

        # Enable
        logger.enable_parallel_logging()
        assert logger.is_parallel_logging_enabled()

        # Disable
        logger.disable_parallel_logging()
        assert not logger.is_parallel_logging_enabled()

    def test_sequential_validation_no_process_id_in_logs(self, tmp_path):
        """Test that sequential validation does NOT include process_id in logs"""
        from validation_pkg.validators.read_validator import ReadValidator
        from validation_pkg.config_manager import ReadConfig
        from validation_pkg.utils.formats import CodingType, ReadFormat

        fastq_file = tmp_path / "test.fastq"
        with open(fastq_file, 'w') as f:
            for i in range(100):
                f.write(f"@read{i}\n")
                f.write("ATCGATCGATCGATCG\n")
                f.write("+\n")
                f.write("IIIIIIIIIIIIIIII\n")

        # Create logger with file output to check JSON logs
        log_file = tmp_path / "test.log"
        logger = setup_logging(log_file=log_file)

        config = ReadConfig(
            filename=fastq_file.name,
            filepath=fastq_file,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=tmp_path,
            global_options={"threads": 1, "validation_level": "strict"}
        )

        settings = ReadValidator.Settings(check_invalid_chars=True)
        validator = ReadValidator(config, settings)

        # Validate (sequential mode)
        validator._parse_file()
        validator._validate_sequences()

        # Check log file - should NOT contain process_id
        log_content = log_file.read_text()

        # Log file uses JSON format, so check for the keys
        assert '"process_id"' not in log_content or log_content.count('"process_id"') == 0
    def test_parallel_validation_has_process_id_in_logs(self, tmp_path):
        """Test that parallel validation DOES include process_id in logs"""
        from validation_pkg.validators.read_validator import ReadValidator
        from validation_pkg.config_manager import ReadConfig
        from validation_pkg.utils.formats import CodingType, ReadFormat

        fastq_file = tmp_path / "test.fastq"
        with open(fastq_file, 'w') as f:
            for i in range(5000):
                f.write(f"@read{i}\n")
                f.write("ATCGATCGATCGATCG\n")
                f.write("+\n")
                f.write("IIIIIIIIIIIIIIII\n")

        # Create logger with file output to check JSON logs
        log_file = tmp_path / "test.log"
        logger = setup_logging(log_file=log_file)

        config = ReadConfig(
            filename=fastq_file.name,
            filepath=fastq_file,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=tmp_path,
            global_options={"threads": 4, "validation_level": "strict"}
        )

        settings = ReadValidator.Settings(check_invalid_chars=True)
        validator = ReadValidator(config, settings)

        # Validate (parallel mode)
        validator._parse_file()
        validator._validate_sequences()

        # Check log file - SHOULD contain process_id during parallel section
        log_content = log_file.read_text()

        # Log file uses JSON format, so check for the keys
        assert '"process_id"' in log_content
    def test_parallel_logging_cleanup_after_validation(self, tmp_path):
        """Test that parallel logging is properly disabled after validation completes"""
        from validation_pkg.validators.read_validator import ReadValidator
        from validation_pkg.config_manager import ReadConfig
        from validation_pkg.utils.formats import CodingType, ReadFormat

        fastq_file = tmp_path / "test.fastq"
        with open(fastq_file, 'w') as f:
            for i in range(2000):
                f.write(f"@read{i}\n")
                f.write("ATCGATCGATCGATCG\n")
                f.write("+\n")
                f.write("IIIIIIIIIIIIIIII\n")

        logger = setup_logging()

        # Verify it's disabled at start
        assert not logger.is_parallel_logging_enabled()

        config = ReadConfig(
            filename=fastq_file.name,
            filepath=fastq_file,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=tmp_path,
            global_options={"threads": 4, "validation_level": "strict"}
        )

        settings = ReadValidator.Settings(check_invalid_chars=True)
        validator = ReadValidator(config, settings)

        # Validate (parallel mode)
        validator._parse_file()
        validator._validate_sequences()

        # After validation, parallel logging should be disabled
        assert not logger.is_parallel_logging_enabled()

        # Log something after validation
        logger.info("Post-validation log")
        # This message should NOT have process_id

    def test_trust_mode_no_process_id_in_logs(self, tmp_path):
        """Test that trust mode (always sequential) does NOT include process_id"""
        from validation_pkg.validators.read_validator import ReadValidator
        from validation_pkg.config_manager import ReadConfig
        from validation_pkg.utils.formats import CodingType, ReadFormat

        fastq_file = tmp_path / "test.fastq"
        with open(fastq_file, 'w') as f:
            for i in range(1000):
                f.write(f"@read{i}\n")
                f.write("ATCGATCGATCGATCG\n")
                f.write("+\n")
                f.write("IIIIIIIIIIIIIIII\n")

        # Create logger with file output
        log_file = tmp_path / "test.log"
        logger = setup_logging(log_file=log_file)

        config = ReadConfig(
            filename=fastq_file.name,
            filepath=fastq_file,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=tmp_path,
            global_options={"threads": 8, "validation_level": "trust"}  # threads=8 but trust mode
        )

        settings = ReadValidator.Settings(check_invalid_chars=True)
        validator = ReadValidator(config, settings)

        # Validate (trust mode - always sequential)
        validator._parse_file()
        validator._validate_sequences()

        # Check log file - should NOT contain process_id
        log_content = log_file.read_text()

        # Trust mode is sequential, so no process_id
        assert '"process_id"' not in log_content or log_content.count('"process_id"') == 0
    def test_full_integration_sequential_vs_parallel(self, tmp_path):
        """Integration test: Compare sequential vs parallel logging in the same test"""
        from validation_pkg.validators.read_validator import ReadValidator
        from validation_pkg.config_manager import ReadConfig
        from validation_pkg.utils.formats import CodingType, ReadFormat

        fastq_file = tmp_path / "test.fastq"
        with open(fastq_file, 'w') as f:
            for i in range(5000):
                f.write(f"@read{i}\n")
                f.write("ATCGATCGATCGATCG\n")
                f.write("+\n")
                f.write("IIIIIIIIIIIIIIII\n")

        # Test 1: Sequential validation
        log_file_seq = tmp_path / "sequential.log"
        logger = setup_logging(log_file=log_file_seq)

        config_seq = ReadConfig(
            filename=fastq_file.name,
            filepath=fastq_file,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=tmp_path,
            global_options={"threads": 1, "validation_level": "strict"}
        )

        validator_seq = ReadValidator(config_seq, ReadValidator.Settings(check_invalid_chars=True))
        validator_seq._parse_file()
        validator_seq._validate_sequences()

        # Test 2: Parallel validation (clear and setup new logger)
        logger.clear_issues()
        logger.logger = None
        log_file_par = tmp_path / "parallel.log"
        logger = setup_logging(log_file=log_file_par)

        config_par = ReadConfig(
            filename=fastq_file.name,
            filepath=fastq_file,
            ngs_type="illumina",
            coding_type=CodingType.NONE,
            detected_format=ReadFormat.FASTQ,
            output_dir=tmp_path,
            global_options={"threads": 4, "validation_level": "strict"}
        )

        validator_par = ReadValidator(config_par, ReadValidator.Settings(check_invalid_chars=True))
        validator_par._parse_file()
        validator_par._validate_sequences()

        # Compare log files
        seq_log = log_file_seq.read_text()
        par_log = log_file_par.read_text()

        # Sequential should NOT have process_id
        assert '"process_id"' not in seq_log or seq_log.count('"process_id"') == 0
        # Parallel SHOULD have process_id
        assert '"process_id"' in par_log
        # Verify parallel logging was cleaned up
        assert not logger.is_parallel_logging_enabled()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])