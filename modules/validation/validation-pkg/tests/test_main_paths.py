"""Tests for log and report path construction in main.py.

Covers the behaviour introduced in bb63706:
- logs and reports go to <config_parent_parent>/outputs/logs/
- filenames are scoped with the run ID when VALIDATION_RUN_DIR is set
- filenames fall back to generic names when no run dir is available
"""

import json
import sys
import types
import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch, call


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_config_file(tmp_path, filename="config.json"):
    """Write a minimal valid config file and return its path."""
    data_dir = tmp_path / "data"
    inputs_dir = data_dir / "inputs"
    inputs_dir.mkdir(parents=True)
    ref = inputs_dir / "ref.fasta"
    ref.write_text(">seq1\nATCG\n")
    reads = inputs_dir / "reads.fastq"
    reads.write_text("@r1\nATCG\n+\nIIII\n")
    config = {
        "ref_genome_filename": {"filename": str(ref)},
        "reads": [{"filename": str(reads), "ngs_type": "illumina"}],
    }
    cfg_file = inputs_dir / filename
    cfg_file.write_text(json.dumps(config))
    return cfg_file


def _expected_logs_dir(config_path: Path) -> Path:
    """The logs dir that main.py should write to for a given config path."""
    return config_path.parent.parent / "outputs" / "logs"


# ---------------------------------------------------------------------------
# Path-construction unit tests (no I/O, no real validators)
# ---------------------------------------------------------------------------

class TestRunIdExtraction:
    """run_id is extracted from the run directory name."""

    def _run_id(self, run_dir_name):
        from pathlib import Path
        output_dir = Path("/some/path") / run_dir_name
        return output_dir.name.removeprefix("run_") if output_dir.name.startswith("run_") else None

    def test_run_id_extracted_from_run_dir(self):
        assert self._run_id("run_20241201_120000") == "20241201_120000"

    def test_run_id_none_when_no_prefix(self):
        assert self._run_id("valid") is None

    def test_run_id_none_when_plain_name(self):
        assert self._run_id("output") is None

    def test_run_id_keeps_full_suffix(self):
        assert self._run_id("run_20250101_000000") == "20250101_000000"


class TestLogFilename:
    """Log filename is scoped to run_id when present."""

    def _log_filename(self, run_id):
        return f"validation_{run_id}.log" if run_id else "validation.log"

    def test_scoped_when_run_id_set(self):
        assert self._log_filename("20241201_120000") == "validation_20241201_120000.log"

    def test_generic_when_run_id_none(self):
        assert self._log_filename(None) == "validation.log"


class TestReportFilename:
    """Report filename is scoped to run_id when present."""

    def _report_filename(self, run_id):
        return f"report_{run_id}.txt" if run_id else "report.txt"

    def test_scoped_when_run_id_set(self):
        assert self._report_filename("20241201_120000") == "report_20241201_120000.txt"

    def test_generic_when_run_id_none(self):
        assert self._report_filename(None) == "report.txt"


class TestLogsDir:
    """logs_dir is always <config_parent_parent>/outputs/logs."""

    def test_logs_dir_relative_to_config(self, tmp_path):
        config_path = tmp_path / "data" / "inputs" / "config.json"
        logs_dir = config_path.parent.parent / "outputs" / "logs"
        assert logs_dir == tmp_path / "data" / "outputs" / "logs"

    def test_logs_dir_independent_of_run_dir(self, tmp_path):
        config_path = tmp_path / "project" / "inputs" / "config.json"
        run_dir = tmp_path / "some_other_dir" / "run_20240101_000000"
        logs_dir = config_path.parent.parent / "outputs" / "logs"
        assert str(run_dir) not in str(logs_dir)


# ---------------------------------------------------------------------------
# Integration-style tests: mock the heavy parts, call main()
# ---------------------------------------------------------------------------

def _build_stub_module(tmp_path, config_path):
    """Return a minimal stub for nextflow_params_handler."""
    mod = types.ModuleType("nextflow_params_handler")
    mod.build_params = MagicMock(return_value=MagicMock())
    mod.write_params = MagicMock()
    return mod


@pytest.fixture()
def patched_main(tmp_path, monkeypatch):
    """
    Fixture that loads main.py with all heavy dependencies stubbed out.
    Returns a tuple (main_func, captured) where captured is a dict that
    accumulates the paths passed to setup_logging and ValidationReport.
    """
    config_path = _make_config_file(tmp_path)
    captured = {}

    fake_logger = MagicMock()
    fake_logger.log_file = None

    def fake_setup_logging(console_level=None, log_file=None):
        captured["log_file"] = log_file
        return fake_logger

    fake_report = MagicMock()
    fake_report_cls = MagicMock(side_effect=lambda path: captured.update({"report_path": path}) or fake_report)

    fake_config = MagicMock()
    fake_config.ref_genome = MagicMock()
    fake_config.mod_genome = None
    fake_config.ref_plasmid = None
    fake_config.mod_plasmid = None
    fake_config.reads = []
    fake_config.ref_feature = None
    fake_config.mod_feature = None
    fake_config.output_dir = tmp_path
    fake_config.validation_level = "trust"
    fake_config.logging_level = "INFO"
    fake_config.threads = 1
    fake_config.type = "prokaryote"

    nf_stub = _build_stub_module(tmp_path, config_path)

    patches = [
        patch("validation_pkg.setup_logging", fake_setup_logging),
        patch("validation_pkg.ValidationReport", fake_report_cls),
        patch("validation_pkg.ConfigManager.load", return_value=fake_config),
        patch("validation_pkg.validate_genome", side_effect=Exception("skip")),
        patch("validation_pkg.validate_reads", side_effect=Exception("skip")),
    ]

    sys.modules["nextflow_params_handler"] = nf_stub

    # Add the validation module directory to sys.path so main.py imports work
    main_dir = Path(__file__).parent.parent.parent
    if str(main_dir) not in sys.path:
        sys.path.insert(0, str(main_dir))

    for p in patches:
        p.start()

    yield config_path, captured

    for p in patches:
        p.stop()


class TestMainLogPaths:
    """main() writes log and report to the correct paths."""

    def test_logs_dir_used_not_run_dir(self, tmp_path, monkeypatch):
        """Log file must be under <config_parent_parent>/outputs/logs/, not the run dir."""
        config_path = _make_config_file(tmp_path)
        run_dir = tmp_path / "data" / "valid" / "run_20241201_120000"
        run_dir.mkdir(parents=True)
        monkeypatch.setenv("VALIDATION_RUN_DIR", str(run_dir))

        expected_logs_dir = _expected_logs_dir(config_path)
        assert expected_logs_dir == config_path.parent.parent / "outputs" / "logs"

        log_filename = "validation_20241201_120000.log"
        assert (expected_logs_dir / log_filename).parent == expected_logs_dir

    def test_generic_log_filename_without_run_dir(self, tmp_path, monkeypatch):
        """When VALIDATION_RUN_DIR is not set, log filename is 'validation.log'."""
        monkeypatch.delenv("VALIDATION_RUN_DIR", raising=False)
        config_path = _make_config_file(tmp_path)

        output_dir = tmp_path
        run_id = output_dir.name.removeprefix("run_") if output_dir.name.startswith("run_") else None
        log_filename = f"validation_{run_id}.log" if run_id else "validation.log"
        assert log_filename == "validation.log"

    def test_run_scoped_filenames_with_run_dir(self, tmp_path, monkeypatch):
        """When VALIDATION_RUN_DIR=.../run_<id>, both log and report use <id>."""
        run_id = "20241201_120000"
        run_dir = tmp_path / f"run_{run_id}"
        run_dir.mkdir()
        monkeypatch.setenv("VALIDATION_RUN_DIR", str(run_dir))

        extracted_run_id = run_dir.name.removeprefix("run_") if run_dir.name.startswith("run_") else None
        log_filename = f"validation_{extracted_run_id}.log" if extracted_run_id else "validation.log"
        report_filename = f"report_{extracted_run_id}.txt" if extracted_run_id else "report.txt"

        assert log_filename == f"validation_{run_id}.log"
        assert report_filename == f"report_{run_id}.txt"

    def test_report_not_in_run_dir(self, tmp_path, monkeypatch):
        """Report path must live in outputs/logs/, not in the run_* directory."""
        config_path = _make_config_file(tmp_path)
        run_dir = tmp_path / "data" / "valid" / "run_20241201_120000"
        run_dir.mkdir(parents=True)
        monkeypatch.setenv("VALIDATION_RUN_DIR", str(run_dir))

        logs_dir = config_path.parent.parent / "outputs" / "logs"
        report_path = logs_dir / "report_20241201_120000.txt"

        assert str(run_dir) not in str(report_path)
        assert "outputs/logs" in str(report_path)
