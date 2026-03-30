"""Tests for global option defaults and CLI option fallback in ConfigManager."""

import pytest
import json
import tempfile
from pathlib import Path

from validation_pkg.config_manager import ConfigManager, Config
from validation_pkg.exceptions import ConfigurationError


# ---------------------------------------------------------------------------
# Shared helper
# ---------------------------------------------------------------------------

def _write_config(directory: Path, options: dict = None) -> Path:
    """Write a minimal config.json and return its path."""
    (directory / "ref.fasta").write_text(">seq1\nATCG\n")
    (directory / "reads.fastq").write_text("@r1\nATCG\n+\nIIII\n")
    data = {
        "ref_genome_filename": {"filename": "ref.fasta"},
        "reads": [{"filename": "reads.fastq", "ngs_type": "illumina"}],
    }
    if options is not None:
        data["options"] = options
    config_file = directory / "config.json"
    config_file.write_text(json.dumps(data))
    return config_file


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

class TestParameterDefaults:
    """All 5 global options default to the correct values when omitted from config.json."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def config(self, temp_dir):
        return ConfigManager.load(str(_write_config(temp_dir)))

    def test_threads_default_is_none(self, config):
        assert config.threads is None

    def test_validation_level_default_is_trust(self, config):
        assert config.validation_level == 'trust'

    def test_logging_level_default_is_info(self, config):
        assert config.logging_level == 'INFO'

    def test_type_default_is_prokaryote(self, config):
        assert config.type == 'prokaryote'

    def test_force_defragment_ref_default_is_false(self, config):
        assert config.force_defragment_ref is False

    def test_all_defaults_stored_in_options_dict(self, config):
        """Defaults are in config.options, not only accessible via properties."""
        assert config.options['threads'] is None
        assert config.options['validation_level'] == 'trust'
        assert config.options['logging_level'] == 'INFO'
        assert config.options['type'] == 'prokaryote'
        assert config.options['force_defragment_ref'] is False

    def test_defaults_propagated_to_ref_genome_global_options(self, config):
        gopt = config.ref_genome.global_options
        assert gopt['threads'] is None
        assert gopt['validation_level'] == 'trust'
        assert gopt['logging_level'] == 'INFO'
        assert gopt['type'] == 'prokaryote'
        assert gopt['force_defragment_ref'] is False

    def test_defaults_propagated_to_read_global_options(self, config):
        gopt = config.reads[0].global_options
        assert gopt['threads'] is None
        assert gopt['validation_level'] == 'trust'
        assert gopt['logging_level'] == 'INFO'
        assert gopt['type'] == 'prokaryote'
        assert gopt['force_defragment_ref'] is False


# ---------------------------------------------------------------------------
# force_defragment_ref
# ---------------------------------------------------------------------------

class TestForceDefragmentRef:
    """Tests for force_defragment_ref option parsing from config.json."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_omitted_defaults_to_false(self, temp_dir):
        config = ConfigManager.load(str(_write_config(temp_dir)))
        assert config.force_defragment_ref is False

    def test_set_true(self, temp_dir):
        config = ConfigManager.load(str(_write_config(temp_dir, {"force_defragment_ref": True})))
        assert config.force_defragment_ref is True

    def test_set_false_explicitly(self, temp_dir):
        config = ConfigManager.load(str(_write_config(temp_dir, {"force_defragment_ref": False})))
        assert config.force_defragment_ref is False

    def test_string_true_raises_error(self, temp_dir):
        with pytest.raises(ConfigurationError, match="must be a boolean"):
            ConfigManager.load(str(_write_config(temp_dir, {"force_defragment_ref": "true"})))

    def test_integer_raises_error(self, temp_dir):
        with pytest.raises(ConfigurationError, match="must be a boolean"):
            ConfigManager.load(str(_write_config(temp_dir, {"force_defragment_ref": 1})))

    def test_propagated_to_file_global_options(self, temp_dir):
        config = ConfigManager.load(str(_write_config(temp_dir, {"force_defragment_ref": True})))
        assert config.ref_genome.global_options['force_defragment_ref'] is True


# ---------------------------------------------------------------------------
# CLI options fallback
# ---------------------------------------------------------------------------

class TestCliOptionsFallback:
    """CLI options (cli_options arg) apply when a key is absent from config.json."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_cli_threads_applied(self, temp_dir):
        config = ConfigManager.load(str(_write_config(temp_dir)), cli_options={'threads': 4})
        assert config.threads == 4

    def test_cli_threads_none_means_autodetect(self, temp_dir):
        config = ConfigManager.load(str(_write_config(temp_dir)), cli_options={'threads': None})
        assert config.threads is None

    def test_cli_validation_level_applied(self, temp_dir):
        config = ConfigManager.load(str(_write_config(temp_dir)), cli_options={'validation_level': 'strict'})
        assert config.validation_level == 'strict'

    def test_cli_logging_level_applied(self, temp_dir):
        config = ConfigManager.load(str(_write_config(temp_dir)), cli_options={'logging_level': 'DEBUG'})
        assert config.logging_level == 'DEBUG'

    def test_cli_type_applied(self, temp_dir):
        config = ConfigManager.load(str(_write_config(temp_dir)), cli_options={'type': 'eukaryote'})
        assert config.type == 'eukaryote'

    def test_cli_force_defragment_ref_applied(self, temp_dir):
        config = ConfigManager.load(str(_write_config(temp_dir)), cli_options={'force_defragment_ref': True})
        assert config.force_defragment_ref is True

    def test_cli_options_propagated_to_file_global_options(self, temp_dir):
        config = ConfigManager.load(
            str(_write_config(temp_dir)),
            cli_options={'validation_level': 'minimal', 'threads': 2}
        )
        assert config.ref_genome.global_options['validation_level'] == 'minimal'
        assert config.ref_genome.global_options['threads'] == 2
        assert config.reads[0].global_options['validation_level'] == 'minimal'

    def test_none_cli_options_uses_defaults(self, temp_dir):
        config = ConfigManager.load(str(_write_config(temp_dir)), cli_options=None)
        assert config.validation_level == 'trust'
        assert config.threads is None

    def test_empty_cli_options_uses_defaults(self, temp_dir):
        config = ConfigManager.load(str(_write_config(temp_dir)), cli_options={})
        assert config.validation_level == 'trust'
        assert config.threads is None


# ---------------------------------------------------------------------------
# Priority: config.json > cli_options > defaults
# ---------------------------------------------------------------------------

class TestOptionsPriority:
    """Priority order: config.json beats CLI; CLI beats built-in defaults."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    # --- config.json beats CLI ---

    def test_config_beats_cli_threads(self, temp_dir):
        config = ConfigManager.load(
            str(_write_config(temp_dir, {"threads": 8})),
            cli_options={'threads': 2}
        )
        assert config.threads == 8

    def test_config_beats_cli_validation_level(self, temp_dir):
        config = ConfigManager.load(
            str(_write_config(temp_dir, {"validation_level": "strict"})),
            cli_options={'validation_level': 'minimal'}
        )
        assert config.validation_level == 'strict'

    def test_config_beats_cli_logging_level(self, temp_dir):
        config = ConfigManager.load(
            str(_write_config(temp_dir, {"logging_level": "WARNING"})),
            cli_options={'logging_level': 'DEBUG'}
        )
        assert config.logging_level == 'WARNING'

    def test_config_beats_cli_type(self, temp_dir):
        config = ConfigManager.load(
            str(_write_config(temp_dir, {"type": "eukaryote"})),
            cli_options={'type': 'prokaryote'}
        )
        assert config.type == 'eukaryote'

    def test_config_beats_cli_force_defragment_ref(self, temp_dir):
        config = ConfigManager.load(
            str(_write_config(temp_dir, {"force_defragment_ref": False})),
            cli_options={'force_defragment_ref': True}
        )
        assert config.force_defragment_ref is False

    # --- CLI beats default ---

    def test_cli_beats_default_threads(self, temp_dir):
        config = ConfigManager.load(str(_write_config(temp_dir)), cli_options={'threads': 6})
        assert config.threads == 6

    def test_cli_beats_default_validation_level(self, temp_dir):
        config = ConfigManager.load(str(_write_config(temp_dir)), cli_options={'validation_level': 'minimal'})
        assert config.validation_level == 'minimal'

    def test_cli_beats_default_type(self, temp_dir):
        config = ConfigManager.load(str(_write_config(temp_dir)), cli_options={'type': 'eukaryote'})
        assert config.type == 'eukaryote'

    def test_cli_beats_default_force_defragment_ref(self, temp_dir):
        config = ConfigManager.load(str(_write_config(temp_dir)), cli_options={'force_defragment_ref': True})
        assert config.force_defragment_ref is True

    # --- Mixed: each option resolved from the correct source ---

    def test_mixed_sources(self, temp_dir):
        """Some keys from config.json, some from CLI, rest from defaults."""
        config = ConfigManager.load(
            str(_write_config(temp_dir, {"threads": 8, "type": "eukaryote"})),
            cli_options={
                'threads': 2,                # ignored — config.json wins
                'validation_level': 'minimal',  # CLI fallback
            }
        )
        assert config.threads == 8              # from config.json
        assert config.type == 'eukaryote'       # from config.json
        assert config.validation_level == 'minimal'   # from CLI
        assert config.logging_level == 'INFO'   # from default
        assert config.force_defragment_ref is False   # from default

    def test_mixed_sources_reflected_in_file_global_options(self, temp_dir):
        """Priority resolution is visible in file-level global_options."""
        config = ConfigManager.load(
            str(_write_config(temp_dir, {"threads": 8, "type": "eukaryote"})),
            cli_options={'validation_level': 'minimal'}
        )
        gopt = config.ref_genome.global_options
        assert gopt['threads'] == 8
        assert gopt['type'] == 'eukaryote'
        assert gopt['validation_level'] == 'minimal'
        assert gopt['logging_level'] == 'INFO'
        assert gopt['force_defragment_ref'] is False


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
