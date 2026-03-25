"""Unsupported workaround: merge a highly fragmented reference into one sequence."""

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from validation_pkg.config_manager import GenomeConfig
from validation_pkg.logger import get_logger
from validation_pkg.utils import file_handler


def defragment_reference(
    genome_config: GenomeConfig,
) -> tuple[Path, Path]:
    """Merge all sequences in genome_config into a single FASTA sequence.

    This is an UNSUPPORTED workaround for highly fragmented references.
    Results produced from the merged file are NOT guaranteed to be meaningful.

    Returns:
        (merged_fasta_path, tsv_path) — both written next to the original file.
    """
    logger = get_logger()
    _warn_before_start(genome_config, logger)

    tsv_rows, merged_seq = _consume_sequences(genome_config, logger)

    merged_fasta_path = _write_fasta(genome_config, merged_seq, logger)
    tsv_path = _write_tsv(genome_config, tsv_rows, logger)

    logger.warning(
        "Defragmentation complete. "
        "This merged file is a pipeline workaround only. "
        "Do NOT draw biological conclusions from results produced with it."
    )

    return merged_fasta_path, tsv_path


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _warn_before_start(genome_config: GenomeConfig, logger) -> None:
    logger.warning("=" * 70)
    logger.warning("DEFRAGMENTATION: UNSUPPORTED OPERATION — READ CAREFULLY")
    logger.warning("=" * 70)
    logger.warning(f"Input file  : {genome_config.filepath}")
    logger.warning(
        f"Format      : {genome_config.detected_format}  |  Compression: {genome_config.coding_type}"
    )
    logger.warning(
        "All contig boundaries and original sequence IDs will be LOST "
        "in the merged output. The join order TSV is the only record "
        "of what was concatenated and in which order."
    )
    logger.warning(
        "Do NOT use results produced from this merged reference for "
        "regulatory submissions or biological conclusions without "
        "independent expert review."
    )
    logger.warning("=" * 70)


def _consume_sequences(
    genome_config: GenomeConfig,
    logger,
) -> tuple[list[tuple[str, int, int]], Seq]:
    """Read all sequences with BioPython, log each one, return TSV rows + merged Seq."""
    tsv_rows: list[tuple[str, int, int]] = []
    parts: list[str] = []
    offset = 0

    fmt = genome_config.detected_format.to_biopython()
    with file_handler.open_file_with_coding_type(
        genome_config.filepath, genome_config.coding_type, "rt"
    ) as handle:
        for i, record in enumerate(SeqIO.parse(handle, fmt), start=1):
            seq_str = str(record.seq)
            seq_len = len(seq_str)
            logger.warning(f"  [{i}] Consuming contig '{record.id}' ({seq_len} bp)")
            tsv_rows.append((record.id, seq_len, offset + 1))
            parts.append(seq_str)
            offset += seq_len

    if not tsv_rows:
        raise ValueError(
            f"No sequences found in {genome_config.filepath}. "
            "Cannot defragment an empty file."
        )

    logger.warning(f"Consumed {len(tsv_rows)} contig(s), total merged length: {offset} bp")
    return tsv_rows, Seq("".join(parts))


def _write_fasta(genome_config: GenomeConfig, merged_seq: Seq, logger) -> Path:
    record = SeqRecord(
        merged_seq,
        id=f"merged_{genome_config.basename}",
        description="force-defragmented — unsupported workaround",
    )
    out_path = (
        genome_config.filepath.parent
        / f"{genome_config.basename}_defragmented.fasta"
    )
    with open(out_path, "w") as fh:
        SeqIO.write(record, fh, "fasta")
    logger.warning(f"Merged FASTA written: {out_path}")
    return out_path


def _write_tsv(
    genome_config: GenomeConfig,
    tsv_rows: list[tuple[str, int, int]],
    logger,
) -> Path:
    # data/inputs/../outputs/tables  →  data/outputs/tables
    tsv_dir = genome_config.filepath.parent.parent / "outputs" / "tables" / "tsv"
    tsv_dir.mkdir(parents=True, exist_ok=True)
    out_path = tsv_dir / f"{genome_config.basename}_defragmented_join_order.tsv"
    with open(out_path, "w") as fh:
        fh.write("seq_id\tlength\tstart\n")
        for seq_id, length, start in tsv_rows:
            fh.write(f"{seq_id}\t{length}\t{start}\n")
    logger.warning(f"Join order TSV written: {out_path}")
    return out_path
