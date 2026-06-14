import importlib.util
import math
import sys
from pathlib import Path
from types import SimpleNamespace


def load_create_sv_output(monkeypatch):
    def isna(value):
        return value is None or (isinstance(value, float) and math.isnan(value))

    monkeypatch.setitem(sys.modules, "numpy", SimpleNamespace(generic=(), nan=float("nan")))
    monkeypatch.setitem(sys.modules, "pandas", SimpleNamespace(isna=isna))
    monkeypatch.setitem(sys.modules, "structlog", SimpleNamespace())

    script = Path(__file__).resolve().parents[1] / "create_sv_output.py"
    spec = importlib.util.spec_from_file_location("create_sv_output", script)
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, "create_sv_output", module)
    spec.loader.exec_module(module)
    return module


def test_delly_supporting_reads_uses_split_read_support_for_single_end_like_rows(monkeypatch):
    create_sv_output = load_create_sv_output(monkeypatch)

    assert create_sv_output._resolve_delly_supporting_reads(
        {"supporting_reads": "0", "PE": "0", "SR": "1"}
    ) == 1


def test_delly_supporting_reads_counts_paired_end_support_as_one_fragment(monkeypatch):
    create_sv_output = load_create_sv_output(monkeypatch)

    assert create_sv_output._resolve_delly_supporting_reads(
        {"supporting_reads": "1", "PE": "1", "SR": "0"}
    ) == 1


def test_delly_supporting_reads_uses_variant_junction_read_support(monkeypatch):
    create_sv_output = load_create_sv_output(monkeypatch)

    assert create_sv_output._resolve_delly_supporting_reads(
        {"supporting_reads": "0", "PE": "0", "SR": "0", "DV": "0", "RV": "1"}
    ) == 1


def test_delly_supporting_reads_prefers_format_evidence_over_info_evidence(monkeypatch):
    create_sv_output = load_create_sv_output(monkeypatch)

    assert create_sv_output._resolve_delly_supporting_reads(
        {"supporting_reads": "12", "PE": "10", "SR": "5", "DV": "4", "RV": "1"}
    ) == 5


def test_delly_supporting_reads_falls_back_to_legacy_supporting_reads(monkeypatch):
    create_sv_output = load_create_sv_output(monkeypatch)

    assert create_sv_output._resolve_delly_supporting_reads(
        {"supporting_reads": "7"}
    ) == 7
