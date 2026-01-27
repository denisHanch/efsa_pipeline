"""
Tests for the logging functionality.

Tests cover:
- Logger initialization and configuration
- Log level management
- Timer functionality for performance tracking
- Validation issue tracking
- File timing summary
- Console and file output formatting
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
                basename="test",
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
                basename="test",
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
                basename="test",
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
                basename="test",
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


class TestTimingMethods:
    """Test suite for timing functionality."""

    @pytest.fixture(autouse=True)
    def reset_logger(self):
        """Reset logger state before each test."""
        logger = get_logger()
        logger.clear_issues()
        logger.logger = None
        # Clear timers
        with logger._timers_lock:
            logger._timers.clear()
        yield
        logger.clear_issues()
        with logger._timers_lock:
            logger._timers.clear()

    def test_start_stop_timer_basic(self):
        """Test basic timer functionality."""
        import time
        logger = get_logger()
        logger.setup()

        logger.start_timer("test_operation")
        time.sleep(0.1)  # Sleep for 100ms
        elapsed = logger.stop_timer("test_operation")

        # Should be approximately 0.1 seconds (allow some tolerance)
        assert 0.09 < elapsed < 0.15, f"Expected ~0.1s, got {elapsed}s"

    def test_stop_timer_without_start_raises_error(self):
        """Test that stopping a non-existent timer raises KeyError."""
        logger = get_logger()
        logger.setup()

        with pytest.raises(KeyError, match="Timer 'nonexistent' was never started"):
            logger.stop_timer("nonexistent")

    def test_multiple_timers_independent(self):
        """Test that multiple named timers work independently."""
        import time
        logger = get_logger()
        logger.setup()

        # Start two timers
        logger.start_timer("timer1")
        time.sleep(0.05)
        logger.start_timer("timer2")
        time.sleep(0.05)

        # Stop timer1 (should be ~0.1s)
        elapsed1 = logger.stop_timer("timer1")
        # Stop timer2 (should be ~0.05s)
        elapsed2 = logger.stop_timer("timer2")

        assert 0.09 < elapsed1 < 0.15, f"Timer1: expected ~0.1s, got {elapsed1}s"
        assert 0.04 < elapsed2 < 0.08, f"Timer2: expected ~0.05s, got {elapsed2}s"

    def test_get_timers_returns_copy(self):
        """Test that get_timers returns a copy of timers dict."""
        import time
        logger = get_logger()
        logger.setup()

        logger.start_timer("test")
        time.sleep(0.05)
        logger.stop_timer("test")

        timers = logger.get_timers()
        assert "test" in timers
        assert isinstance(timers, dict)

        # Modifying the copy should not affect internal state
        timers["test"] = 999.0
        timers_again = logger.get_timers()
        assert timers_again["test"] != 999.0


class TestFileAutoIncrement:
    """Test suite for log and report file auto-increment functionality."""

    @pytest.fixture(autouse=True)
    def reset_logger(self):
        """Reset logger state before each test."""
        logger = get_logger()
        logger.clear_issues()
        logger.logger = None
        with logger._timers_lock:
            logger._timers.clear()
        yield
        logger.clear_issues()
        with logger._timers_lock:
            logger._timers.clear()

    def test_log_file_auto_increment(self):
        """Test that log files auto-increment instead of overwriting."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = Path(tmpdir) / "validation.log"

            # First setup - should create validation.log
            logger1 = get_logger()
            logger1.setup(log_file=log_file)
            logger1.info("First run")

            assert log_file.exists()
            first_content = log_file.read_text()
            assert "First run" in first_content

            # Reset logger
            logger1.logger = None
            logger1.clear_issues()

            # Second setup - should create validation_001.log
            logger2 = get_logger()
            logger2.setup(log_file=log_file)
            logger2.info("Second run")

            log_file_001 = Path(tmpdir) / "validation_001.log"
            assert log_file_001.exists()
            second_content = log_file_001.read_text()
            assert "Second run" in second_content

            # Original file should still have first content
            first_content_after = log_file.read_text()
            assert "First run" in first_content_after
            assert "Second run" not in first_content_after


if __name__ == "__main__":
    pytest.main([__file__, "-v"])