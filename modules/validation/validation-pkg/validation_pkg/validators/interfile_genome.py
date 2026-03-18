"""Inter-file validation for genome files."""

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Any
from Bio import SeqIO
from ..utils.base_settings import BaseSettings
from ..exceptions import InterFileValidationError, ValidationError
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
    # Concatenation settings (used when characterize=True)
    concatenate: bool = False
    min_mapping_quality: int = 0      # PAF col 11; 0 = accept all
    min_query_coverage: float = 0.0   # alignment_len / q_len; 0.0 = accept all
    min_identity: float = 0.0         # residue_matches / alignment_len; 0.0 = accept all
    max_ref_overlap: int = 0          # max allowed overlap (bp) between adjacent queries on ref

    def __post_init__(self):
        if self.same_sequence_lengths and not self.same_sequence_ids:
            raise ValidationError(
                "same_sequence_lengths requires same_sequence_ids=True "
                "(cannot compare lengths without matching IDs)"
            )
        if not (0.0 <= self.min_query_coverage <= 1.0):
            raise ValidationError("min_query_coverage must be between 0.0 and 1.0")
        if not (0.0 <= self.min_identity <= 1.0):
            raise ValidationError("min_identity must be between 0.0 and 1.0")
        if self.min_mapping_quality < 0:
            raise ValidationError("min_mapping_quality must be >= 0")
        if self.max_ref_overlap < 0:
            raise ValidationError("max_ref_overlap must be >= 0")


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
            _characterize_into_metadata(ref_genome_result, mod_genome_result, metadata, logger, settings)
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


def _parse_paf_best_hits(paf_source: str) -> Dict[str, Dict[str, Dict]]:
    """Parse a PAF file (or raw PAF string) and return the best-scoring hit per (query, ref) pair.

    PAF columns (0-indexed):
      0  query name      1  query length
      2  q_start         3  q_end
      4  strand (+/-)    5  ref name
      6  ref length      7  r_start
      8  r_end           9  residue matches
      10 alignment block length   11 mapping quality

    For each (query_name, ref_name) pair the hit with the largest alignment
    block length (col 10) is selected.  A query that maps to multiple distinct
    reference sequences produces one entry per ref_name, enabling a unique
    contig file to be written for every such mapping.

    Parameters
    ----------
    paf_source : str
        Either a file path to a PAF file, or a raw PAF string (containing newlines).

    Returns
    -------
    Dict mapping query_name -> ref_name -> {
        'strand'        : str,   # '+' or '-'
        'alignment_len' : int,   # alignment block length (col 10)
        'q_len'         : int,   # query length (col 1)
        'q_end'         : int,   # query alignment end (col 3)
        'r_start'       : int,   # reference alignment start (col 7)
        'r_end'         : int,   # reference alignment end (col 8)
        'residue_matches': int,  # matching residues (col 9)
        'mapq'          : int,   # mapping quality (col 11)
    }
    """
    if '\t' in paf_source or paf_source == '':
        lines = paf_source.splitlines()
    else:
        with open(paf_source, 'r') as fh:
            lines = fh.readlines()

    best_hits: Dict[str, Dict[str, Dict]] = {}
    for line in lines:
        if not line.strip():
            continue
        cols = line.split('\t')
        if len(cols) < 12:
            continue
        query_name = cols[0]
        strand = cols[4]
        ref_name = cols[5]
        try:
            alignment_block_len = int(cols[10])
            q_len   = int(cols[1])
            q_end   = int(cols[3])
            r_start = int(cols[7])
            r_end   = int(cols[8])
            residue_matches = int(cols[9])
            mapq    = int(cols[11])
        except ValueError:
            continue

        ref_hits = best_hits.setdefault(query_name, {})
        current = ref_hits.get(ref_name)
        if current is None or alignment_block_len > current['alignment_len']:
            ref_hits[ref_name] = {
                'strand':          strand,
                'alignment_len':   alignment_block_len,
                'q_len':           q_len,
                'q_end':           q_end,
                'r_start':         r_start,
                'r_end':           r_end,
                'residue_matches': residue_matches,
                'mapq':            mapq,
            }
    return best_hits


def _filter_and_group_hits(
    best_hits: Dict[str, Dict[str, Dict]],
    settings: GenomeXGenomeSettings,
    logger,
) -> tuple:
    """Filter PAF hits by quality thresholds and group passing queries by reference.

    (query_name, ref_name) pairs failing per-hit thresholds, or belonging to a
    group whose adjacent reference mappings overlap more than max_ref_overlap,
    are returned in individual_ids and written as standalone contig files by
    the caller.

    Returns
    -------
    groups : Dict[ref_name, List[query_name]]
        Queries that passed all checks, sorted by r_start within each group.
    individual_ids : set[Tuple[str, str]]
        (query_name, ref_name) pairs that failed a threshold or overlap check.
    """
    passed: Dict[tuple, Dict] = {}
    individual_ids: set = set()

    for query_name, ref_hits in best_hits.items():
        for ref_name, hit in ref_hits.items():
            if hit['mapq'] < settings.min_mapping_quality:
                logger.debug(
                    f"  {query_name} → {ref_name}: MAPQ {hit['mapq']} < "
                    f"{settings.min_mapping_quality}, writing individually"
                )
                individual_ids.add((query_name, ref_name))
                continue

            query_coverage = hit['alignment_len'] / hit['q_len'] if hit['q_len'] > 0 else 0.0
            if query_coverage < settings.min_query_coverage:
                logger.debug(
                    f"  {query_name} → {ref_name}: query_coverage {query_coverage:.3f} < "
                    f"{settings.min_query_coverage}, writing individually"
                )
                individual_ids.add((query_name, ref_name))
                continue

            identity = (
                hit['residue_matches'] / hit['alignment_len']
                if hit['alignment_len'] > 0 else 0.0
            )
            if identity < settings.min_identity:
                logger.debug(
                    f"  {query_name} → {ref_name}: identity {identity:.3f} < "
                    f"{settings.min_identity}, writing individually"
                )
                individual_ids.add((query_name, ref_name))
                continue

            passed[(query_name, ref_name)] = hit

    # Group by ref_name, sort each group by r_start
    groups: Dict[str, List] = {}
    for (query_name, ref_name), hit in passed.items():
        groups.setdefault(ref_name, []).append(query_name)
    for ref_name in groups:
        groups[ref_name].sort(key=lambda qn: passed[(qn, ref_name)]['r_start'])

    # Per-group overlap check: demote entire group if any adjacent pair overlaps too much
    clean_groups: Dict[str, List] = {}
    for ref_name, query_names in groups.items():
        has_overlap = False
        for i in range(len(query_names) - 1):
            q_a, q_b = query_names[i], query_names[i + 1]
            overlap = max(0, passed[(q_a, ref_name)]['r_end'] - passed[(q_b, ref_name)]['r_start'])
            if overlap > settings.max_ref_overlap:
                logger.warning(
                    f"Concat group for ref '{ref_name}' has reference overlap "
                    f"{overlap} bp between '{q_a}' and '{q_b}' "
                    f"(max_ref_overlap={settings.max_ref_overlap}). "
                    "Writing sequences individually."
                )
                has_overlap = True
                break
        if has_overlap:
            for qn in query_names:
                individual_ids.add((qn, ref_name))
        else:
            clean_groups[ref_name] = query_names

    return clean_groups, individual_ids


def _characterize_into_metadata(ref_genome_result, mod_genome_result, metadata, logger, settings: GenomeXGenomeSettings = None):
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

    # Parse PAF file: for each (query, ref) pair keep the best-scoring hit so
    # that a query mapping to multiple reference sequences produces one contig
    # file per reference mapping.
    best_hits = _parse_paf_best_hits(str(paf_path))
    mapped_ids = set(best_hits.keys())

    logger.debug(f"minimap2 mapped {len(mapped_ids)} sequence(s) to reference")

    # Determine grouping for concatenation mode
    if settings is not None and settings.concatenate:
        groups, individual_ids = _filter_and_group_hits(best_hits, settings, logger)
        concat_pairs = {(qn, ref) for ref, qns in groups.items() for qn in qns}
        logger.debug(
            f"Concatenation: {len(groups)} group(s), "
            f"{len(individual_ids)} (query,ref) pair(s) written individually"
        )
    else:
        groups = {}
        concat_pairs: set = set()

    # All (query_name, ref_name) pairs from best_hits
    all_pairs = {(qn, ref) for qn, ref_hits in best_hits.items() for ref in ref_hits}
    individual_pairs = all_pairs - concat_pairs

    # Split sequences: mapped → contigs, unmapped → plasmids
    contig_seqs: List = []
    plasmid_seqs: List = []
    for record in SeqIO.parse(str(mod_path), 'fasta'):
        if record.id in mapped_ids:
            contig_seqs.append(record)
        else:
            plasmid_seqs.append(record)

    seq_by_id = {seq.id: seq for seq in contig_seqs}
    contig_files: List[str] = []
    contig_index = 0
    concat_groups_metadata: Dict[str, Any] = {}

    # --- Concatenated groups ---
    for ref_name, query_names in groups.items():
        oriented_parts = []
        for qn in query_names:
            hit = best_hits[qn][ref_name]
            seq = seq_by_id[qn]
            new_id = seq.id.rstrip('0123456789')
            if hit['strand'] == '-':
                part = seq.reverse_complement(id=new_id, description='')
            else:
                part = seq[:]
                part.id = new_id
                part.description = ''
            oriented_parts.append(part)

        concat_seq = oriented_parts[0]
        for part in oriented_parts[1:]:
            concat_seq = concat_seq + part
        concat_seq.id = ref_name.rstrip('0123456789')
        concat_seq.description = f"concatenated_{len(query_names)}_queries"

        contig_path = output_dir / f"{base_name}_contig_{contig_index}.fasta"
        with open_compressed_writer(contig_path, CodingType.NONE) as handle:
            SeqIO.write([concat_seq], handle, 'fasta')
        contig_files.append(str(contig_path))
        concat_groups_metadata[ref_name] = {
            'query_names': query_names,
            'contig_file': str(contig_path),
        }
        logger.debug(
            f"Concat contig {contig_index}: {len(query_names)} queries → "
            f"{ref_name} ({len(concat_seq.seq)} bp) → {contig_path.name}"
        )
        contig_index += 1

    # --- Individual contigs: one file per (query_name, ref_name) pair ---
    for query_name, ref_name in sorted(individual_pairs):
        hit = best_hits[query_name][ref_name]
        seq = seq_by_id[query_name]
        new_id = seq.id.rstrip('0123456789')
        if hit['strand'] == '-':
            out_seq = seq.reverse_complement(id=new_id, description='')
        else:
            out_seq = seq[:]
            out_seq.id = new_id
            out_seq.description = ''
        contig_path = output_dir / f"{base_name}_contig_{contig_index}.fasta"
        with open_compressed_writer(contig_path, CodingType.NONE) as handle:
            SeqIO.write([out_seq], handle, 'fasta')
        contig_files.append(str(contig_path))
        logger.debug(
            f"Contig {contig_index}: {new_id} ({len(out_seq.seq)} bp) "
            f"→ {ref_name} [{hit['strand']}] → {contig_path.name}"
        )
        contig_index += 1

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
        f"{len(contig_files)} contig file(s) "
        f"({len(groups)} concatenated group(s), {len(individual_pairs)} individual (query,ref) pair(s)), "
        f"{len(plasmid_seqs)} plasmid(s)"
    )

    # Build orientation/ref-name metadata: nested for groups, tuple-keyed for individual pairs
    contig_orientations: Dict[str, Any] = {
        ref_name: {qn: best_hits[qn][ref_name]['strand'] for qn in qns}
        for ref_name, qns in groups.items()
    }
    contig_orientations.update(
        {(qn, ref): best_hits[qn][ref]['strand'] for qn, ref in individual_pairs}
    )
    contig_ref_names: Dict[str, Any] = {ref_name: ref_name for ref_name in groups}
    contig_ref_names.update(
        {(qn, ref): ref for qn, ref in individual_pairs}
    )

    metadata.update({
        'contigs_found': len(contig_files),
        'plasmids_found': len(plasmid_seqs),
        'contig_files': contig_files,
        'plasmid_file': plasmid_file,
        'contig_orientations': contig_orientations,
        'contig_ref_names': contig_ref_names,
        'concat_groups': concat_groups_metadata or None,
    })
