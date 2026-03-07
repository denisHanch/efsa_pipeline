"""Inter-file validation for genome files."""

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Any
from Bio import SeqIO
from ..utils.base_settings import BaseSettings
from ..exceptions import InterFileValidationError
from ..logger import get_logger
from ..utils.file_handler import check_tool_available, open_compressed_writer
from ..utils.path_utils import strip_all_extensions
from ..utils.formats import CodingType


def _get_metadata_field(metadata_obj, field_name, default=None):
    """Get field from either OutputMetadata object or dict."""
    if hasattr(metadata_obj, field_name):
        # OutputMetadata object - use attribute access
        return getattr(metadata_obj, field_name, default)
    else:
        # Dictionary - use dict access
        return metadata_obj.get(field_name, default)


@dataclass
class GenomeXGenomeSettings(BaseSettings):
    """Settings for genome-to-genome inter-file validation."""
    same_number_of_sequences: bool = True
    same_sequence_ids: bool = False
    same_sequence_lengths: bool = False
    characterize: bool = True

    def __post_init__(self):
        if self.same_sequence_lengths and not self.same_sequence_ids:
            raise ValueError(
                "same_sequence_lengths requires same_sequence_ids=True "
                "(cannot compare lengths without matching IDs)"
            )


def genomexgenome_validation(
    ref_genome_result,  # OutputMetadata or Dict[str, Any]
    mod_genome_result,  # OutputMetadata or Dict[str, Any]
    settings: Optional[GenomeXGenomeSettings] = None
) -> Dict[str, Any]:
    """Validate consistency between two genome files (reference vs modified).

    When settings.characterize is True (default), also runs minimap2 to align
    the modified genome against the reference and separates contigs (mapped)
    from plasmids (unmapped), writing each to output FASTA files.
    Characterization failures are non-fatal and recorded as warnings.
    """
    settings = settings or GenomeXGenomeSettings()
    logger = get_logger()

    warnings = []
    errors = []

    logger.info("Running inter-file validation: genome-to-genome consistency")

    # Extract metadata from OutputMetadata objects
    ref_meta = ref_genome_result
    mod_meta = mod_genome_result

    # Initialize result metadata
    metadata = {
        'ref_num_sequences': _get_metadata_field(ref_meta, 'num_sequences', 0),
        'mod_num_sequences': _get_metadata_field(mod_meta, 'num_sequences', 0),
        'common_sequence_ids': [],
        'ref_only_ids': [],
        'mod_only_ids': [],
        'length_mismatches': {}
    }

    # Check 1: Same number of sequences
    if settings.same_number_of_sequences:
        ref_count = _get_metadata_field(ref_meta, 'num_sequences', 0)
        mod_count = _get_metadata_field(mod_meta, 'num_sequences', 0)

        if ref_count != mod_count:
            error_msg = (
                f"Genome sequence count mismatch: "
                f"reference has {ref_count} sequence(s), "
                f"modified has {mod_count} sequence(s)"
            )
            errors.append(error_msg)
            logger.error(error_msg)
        else:
            logger.debug(f"Sequence count match: {ref_count} sequence(s)")

    # Check 2: Same sequence IDs
    if settings.same_sequence_ids:
        ref_ids = set(_get_metadata_field(ref_meta, 'sequence_ids', []))
        mod_ids = set(_get_metadata_field(mod_meta, 'sequence_ids', []))

        common_ids = ref_ids & mod_ids
        ref_only = ref_ids - mod_ids
        mod_only = mod_ids - ref_ids

        metadata['common_sequence_ids'] = sorted(common_ids)
        metadata['ref_only_ids'] = sorted(ref_only)
        metadata['mod_only_ids'] = sorted(mod_only)

        if ref_only or mod_only:
            error_parts = []
            if ref_only:
                error_parts.append(f"reference-only: {sorted(ref_only)}")
            if mod_only:
                error_parts.append(f"modified-only: {sorted(mod_only)}")

            error_msg = f"Genome sequence ID mismatch: {', '.join(error_parts)}"
            errors.append(error_msg)
            logger.error(error_msg)
        else:
            logger.debug(f"Sequence IDs match: {len(common_ids)} common ID(s)")

    # Check 3: Same sequence lengths for common IDs
    if settings.same_sequence_lengths and common_ids:
        ref_lengths = _get_metadata_field(ref_meta, 'sequence_lengths', {})
        mod_lengths = _get_metadata_field(mod_meta, 'sequence_lengths', {})
        for seq_id in common_ids:

            ref_len = ref_lengths.get(seq_id)
            mod_len = mod_lengths.get(seq_id)
            if ref_len is None or mod_len is None:
                warning_msg = f"Missing length information for sequence '{seq_id}'"
                warnings.append(warning_msg)
                logger.warning(warning_msg)
                continue

            if ref_len != mod_len:
                metadata['length_mismatches'][seq_id] = {
                    'ref_length': ref_len,
                    'mod_length': mod_len,
                    'difference': mod_len - ref_len
                }

                error_msg = (
                    f"Sequence length mismatch for '{seq_id}': "
                    f"reference={ref_len} bp, modified={mod_len} bp "
                    f"(diff={mod_len - ref_len:+d} bp)"
                )
                errors.append(error_msg)
                logger.error(error_msg)

        if not metadata['length_mismatches']:
            logger.debug(f"Sequence lengths match for {len(common_ids)} sequence(s)")

    # Characterize contigs vs plasmids via minimap2 (non-fatal)
    if settings.characterize:
        logger.info("Running genome characterization: contigs vs plasmids via minimap2")
        try:
            _characterize_into_metadata(ref_genome_result, mod_genome_result, metadata, logger)
        except Exception as e:
            error_msg = f"Genome characterization skipped: {e}"
            errors.append(error_msg)
            logger.error(error_msg)
            metadata.update({
                'contigs_found': None,
                'plasmids_found': None,
                'contig_files': [],
                'plasmid_file': None,
            })

    # Determine if validation passed (no ERROR-level issues)
    passed = len(errors) == 0
    
    result = {
        'passed': passed,
        'warnings': warnings,
        'errors': errors,
        'metadata': metadata
    }

    if passed:
        logger.info(f"✓ Genome inter-file validation passed")
    else:
        logger.error(f"✗ Genome inter-file validation failed: {len(errors)} error(s), {len(warnings)} warning(s)")
        error_summary = '\n  - '.join([''] + errors)
        raise InterFileValidationError(f"Genome inter-file validation failed:{error_summary}")
    return result


def _parse_paf_best_hits(paf_output: str) -> Dict[str, Dict]:
    """Parse PAF output and return the best-scoring hit per query sequence.

    PAF columns (0-indexed):
      0  query name      1  query length
      2  q_start         3  q_end
      4  strand (+/-)    5  ref name
      6  ref length      7  r_start
      8  r_end           9  residue matches
      10 alignment block length   11 mapping quality

    For each query name the hit with the largest alignment block length
    (col 10) is selected.  This correctly handles multi-copy elements such
    as ribosomal RNA that produce many hits: only the highest-scoring
    reference assignment is kept.

    Returns
    -------
    Dict mapping query_name -> {
        'ref_name'     : str,   # name of best-matching reference sequence (expected to be only one reference = should be same all the time)
        'strand'       : str,   # '+' or '-'
        'alignment_len': int,   # alignment block length of the winning hit
    }
    """
    best_hits: Dict[str, Dict] = {}
    # Track (q_start, r_start) per (query_name, ref_name) pair to detect ambiguous alignments
    pair_starts: Dict[tuple, tuple] = {}
    for line in paf_output.splitlines():
        if not line.strip():
            continue
        cols = line.split('\t')
        if len(cols) < 12:
            continue
        query_name = cols[0]
        strand = cols[4]
        ref_name = cols[5]
        try:
            q_start = int(cols[2])
            r_start = int(cols[7])
            alignment_block_len = int(cols[10])
        except ValueError:
            continue

        pair_key = (query_name, ref_name)
        if pair_key in pair_starts:
            seen_q_start, seen_r_start = pair_starts[pair_key]
            if q_start != seen_q_start:
                raise InterFileValidationError(
                    f"Ambiguous minimap2 alignment: query '{query_name}' maps to "
                    f"reference '{ref_name}' with conflicting q_start values "
                    f"({seen_q_start} and {q_start})"
                )
            if r_start != seen_r_start:
                raise InterFileValidationError(
                    f"Ambiguous minimap2 alignment: query '{query_name}' maps to "
                    f"reference '{ref_name}' with conflicting r_start values "
                    f"({seen_r_start} and {r_start})"
                )
        else:
            pair_starts[pair_key] = (q_start, r_start)

        current = best_hits.get(query_name)
        if current is None or alignment_block_len > current['alignment_len']:
            best_hits[query_name] = {
                'ref_name': ref_name,
                'strand': strand,
                'alignment_len': alignment_block_len,
            }
    return best_hits


def _characterize_into_metadata(ref_genome_result, mod_genome_result, metadata, logger):
    """Run minimap2 alignment and populate characterization keys into metadata.

    Raises an exception on any failure so the caller can demote it to a warning.
    """
    if not check_tool_available('minimap2'):
        raise InterFileValidationError(
            "minimap2 is required for genome characterization but was not found. "
            "Please install minimap2 and ensure it is on your PATH."
        )

    ref_path = _get_metadata_field(ref_genome_result, 'output_file')
    mod_path = _get_metadata_field(mod_genome_result, 'output_file')
    mod_filename = _get_metadata_field(mod_genome_result, 'output_filename')

    logger.debug(f"Running minimap2: ref={ref_path} query={mod_path}")
    try:
        proc = subprocess.run(
            ['minimap2', '-x', 'asm5', str(ref_path), str(mod_path)],
            capture_output=True,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        raise InterFileValidationError(
            f"minimap2 failed (exit code {e.returncode}): {e.stderr.strip()}"
        )

    # Parse PAF output: for each query keep only the best-scoring hit so that
    # multi-copy elements (e.g. ribosomal RNA) are resolved to a single
    # reference assignment, and orientation is captured.
    best_hits = _parse_paf_best_hits(proc.stdout)
    mapped_ids = set(best_hits.keys())
    
    logger.debug(f"minimap2 mapped {len(mapped_ids)} sequence(s) to reference")

    # Split sequences into contigs (mapped) and plasmids (unmapped)
    contig_seqs: List = []
    plasmid_seqs: List = []
    for record in SeqIO.parse(str(mod_path), 'fasta'):
        if record.id in mapped_ids:
            contig_seqs.append(record)
        else:
            plasmid_seqs.append(record)

    output_dir = Path(mod_path).parent
    base_name = strip_all_extensions(mod_filename)

    # Write each contig to its own file
    contig_files: List[str] = []
    for i, seq in enumerate(contig_seqs):
        hit = best_hits[seq.id]
        new_id = seq.id.rstrip('0123456789')
        if hit['strand'] == '-':
            out_seq = seq.reverse_complement(id=new_id, description='')
        else:
            out_seq = seq[:]
            out_seq.id = new_id
            out_seq.description = ''
        contig_path = output_dir / f"{base_name}_contig_{i}.fasta"
        with open_compressed_writer(contig_path, CodingType.NONE) as handle:
            SeqIO.write([out_seq], handle, 'fasta')
        contig_files.append(str(contig_path))
        logger.debug(
            f"Contig {i}: {new_id} ({len(out_seq.seq)} bp) "
            f"→ {hit['ref_name']} [{hit['strand']}] → {contig_path.name}"
        )

    # Write all plasmids to a single file
    plasmid_file: Optional[str] = None
    if plasmid_seqs:
        plasmid_path = output_dir / f"{base_name}_plasmids.fasta"
        with open_compressed_writer(plasmid_path, CodingType.NONE) as handle:
            SeqIO.write(plasmid_seqs, handle, 'fasta')
        plasmid_file = str(plasmid_path)
        logger.debug(f"Plasmids ({len(plasmid_seqs)}) → {plasmid_path.name}")

    logger.info(
        f"Genome characterization complete: "
        f"{len(contig_seqs)} contig(s), {len(plasmid_seqs)} plasmid(s)"
    )

    metadata.update({
        'contigs_found': len(contig_seqs),
        'plasmids_found': len(plasmid_seqs),
        'contig_files': contig_files,
        'plasmid_file': plasmid_file,
        'contig_orientations': {seq.id: best_hits[seq.id]['strand'] for seq in contig_seqs},
        'contig_ref_names': {seq.id: best_hits[seq.id]['ref_name'] for seq in contig_seqs},
    })
