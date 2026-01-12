"""Validation report generation module."""

from dataclasses import dataclass, field
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional, Any, Union
import json
from validation_pkg.utils.path_utils import get_incremented_path


@dataclass
class FileValidationRecord:
    """Complete record of a single file's validation."""
    # Input information
    output_data: dict
    validator_type: str  # "genome", "read", "feature"
    input_settings: Optional[dict] = None  # Serialized Settings object


@dataclass
class InterFileValidationRecord:
    """Record of inter-file validation check."""
    validation_type: str  # "genomexgenome", "readxread", etc.
    status: str  # "PASSED", "FAILED"
    passed: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())


class ValidationReport:
    """Comprehensive validation report builder."""

    def __init__(self, report_path: Path):
        """Initialize validation report."""

        self.report_path = Path(report_path)
        self.report_path.parent.mkdir(parents=True, exist_ok=True)

        # Auto-increment if file exists
        self.report_path = get_incremented_path(self.report_path)

        # Storage for file validation records
        self.file_records: List[FileValidationRecord] = []

        # Storage for inter-file validation records
        self.interfile_records: List[InterFileValidationRecord] = []

        # Track start time
        self.start_time = datetime.now()

    @staticmethod
    def _has_value(data: Dict[str, Any], key: str) -> bool:
        """Check if a key exists in dictionary and its value is not None."""
        return key in data and data[key] is not None

    def write(
        self,
        result: Union[Any, List[Any]],
        file_type: str,
        settings: Optional[Any] = None
    ) -> None:
        """Add a validation result to the report."""
        # Handle file validation results
        if file_type in ("genome", "read", "feature"):
            # Handle both single result and list of results
            results_list = result if isinstance(result, list) else [result]

            if settings and not isinstance(settings, dict):
                settings = settings.to_dict()

            for single_result in results_list:
                if not isinstance(single_result, dict):
                    single_result = single_result.to_dict()
    
                # Create file record
                record = FileValidationRecord(
                    output_data=single_result,
                    validator_type=file_type,
                    input_settings= settings if settings else None
                )

                self.file_records.append(record)

        # Handle inter-file validation results
        elif file_type in ("genomexgenome", "readxread", "featurexgenome"):
            passed = result.get('passed', False)
            status = "PASSED" if passed else "FAILED"

            record = InterFileValidationRecord(
                validation_type=file_type,
                status=status,
                passed=passed,
                errors=result.get('errors', []),
                warnings=result.get('warnings', []),
                metadata=result.get('metadata', {})
            )

            self.interfile_records.append(record)
        else:
            raise ValueError(f"Unknown file_type: {file_type}")

    def flush(self, format: str = "text") -> None:
        """Generate and write the final report to file."""
        if format == "text":
            self._flush_text()
        elif format == "json":
            self._flush_json()
        else:
            raise ValueError(f"Unknown format: {format}. Use 'text' or 'json'")

    def _flush_text(self) -> None:
        """Generate text format report."""
        end_time = datetime.now()
        duration = (end_time - self.start_time).total_seconds()

        lines = []

        # Header with nice formatting
        lines.append("=" * 100)
        lines.append("")
        lines.append("  VALIDATION PIPELINE REPORT")
        lines.append("")
        lines.append("=" * 100)
        lines.append("")
        lines.append(f"  Generated:       {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append(f"  Total Duration:  {duration:.2f}s")
        lines.append("")

        # Summary
        lines.extend(self._generate_summary_section())

        # File-specific results
        if self.file_records:
            lines.extend(self._generate_file_results_section())

        # Inter-file validation results
        if self.interfile_records:
            lines.extend(self._generate_interfile_section())

        # Write to file
        report_text = '\n'.join(lines)
        self.report_path.write_text(report_text, encoding='utf-8')

        print(f"\n✓ Report written to: {self.report_path}")

    def _flush_json(self) -> None:
        """Generate JSON format report."""
        end_time = datetime.now()
        duration = (end_time - self.start_time).total_seconds()

        report_data = {
            "report_metadata": {
                "generated": end_time.isoformat(),
                "duration_seconds": duration,
                "version": "1.0.0"
            },
            "summary": self._get_summary_data(),
            "file_validations": [
                {
                    "validator_type": rec.validator_type,
                    "input_settings": rec.input_settings,
                    "output_data": rec.output_data
                }
                for rec in self.file_records
            ],
            "inter_file_validations": [
                {
                    "validation_type": rec.validation_type,
                    "status": rec.status,
                    "passed": rec.passed,
                    "errors": rec.errors,
                    "warnings": rec.warnings,
                    "metadata": rec.metadata,
                    "timestamp": rec.timestamp
                }
                for rec in self.interfile_records
            ]
        }

        json_path = self.report_path.with_suffix('.json')
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(report_data, f, indent=2, ensure_ascii=False)

        print(f"\n✓ JSON report written to: {json_path}")

    def _generate_summary_section(self) -> List[str]:
        """Generate summary section."""
        lines = []

        # Count files by type
        genome_count = sum(1 for r in self.file_records if r.validator_type == "genome")
        read_count = sum(1 for r in self.file_records if r.validator_type == "read")
        feature_count = sum(1 for r in self.file_records if r.validator_type == "feature")
        total_files = len(self.file_records)

        # Count inter-file validations
        interfile_passed = sum(1 for r in self.interfile_records if r.passed)
        interfile_failed = len(self.interfile_records) - interfile_passed

        # Overall status
        overall_status = "✓ PASSED" if interfile_failed == 0 else "✗ FAILED"

        lines.append("SUMMARY")
        lines.append("=" * 100)
        lines.append("")
        lines.append(f"  Overall Status: {overall_status}")
        lines.append("")
        lines.append(f"  Files Processed: {total_files}")
        if genome_count > 0:
            lines.append(f"    ├─ Genomes:  {genome_count}")
        if read_count > 0:
            lines.append(f"    ├─ Reads:    {read_count}")
        if feature_count > 0:
            lines.append(f"    └─ Features: {feature_count}")
        lines.append("")

        if self.interfile_records:
            lines.append(f"  Inter-file Validations: {len(self.interfile_records)}")
            if interfile_passed > 0:
                lines.append(f"    ├─ Passed:  {interfile_passed} ✓")
            if interfile_failed > 0:
                lines.append(f"    └─ Failed:  {interfile_failed} ✗")
            lines.append("")

        return lines

    def _get_summary_data(self) -> Dict[str, Any]:
        """Get summary data."""
        genome_count = sum(1 for r in self.file_records if r.validator_type == "genome")
        read_count = sum(1 for r in self.file_records if r.validator_type == "read")
        feature_count = sum(1 for r in self.file_records if r.validator_type == "feature")

        interfile_passed = sum(1 for r in self.interfile_records if r.passed)
        interfile_failed = len(self.interfile_records) - interfile_passed

        return {
            "total_files": len(self.file_records),
            "genome_files": genome_count,
            "read_files": read_count,
            "feature_files": feature_count,
            "interfile_validations": len(self.interfile_records),
            "interfile_passed": interfile_passed,
            "interfile_failed": interfile_failed,
            "overall_status": "PASSED" if interfile_failed == 0 else "FAILED"
        }

    def _generate_file_results_section(self) -> List[str]:
        """Generate file results section."""
        lines = []
        lines.append("=" * 100)
        lines.append("FILE VALIDATION RESULTS".center(100))
        lines.append("=" * 100)
        lines.append("")

        for idx, record in enumerate(self.file_records, 1):
            lines.extend(self._format_file_record(idx, record))

        return lines

    def _format_file_record(self, idx: int, record: FileValidationRecord) -> List[str]:
        """Format file validation record."""
        lines = []

        # Header with validator type
        validator_label = record.validator_type.upper()
        lines.append(f"[{idx}] {validator_label} FILE")
        lines.append("")

        # Extract key info from output_data for header
        input_file = record.output_data.get('input_file', 'Unknown')
        output_file = record.output_data.get('output_file', 'Unknown')
        elapsed_time = record.output_data.get('elapsed_time')

        lines.append(f"  Input:  {input_file}")
        lines.append(f"  Output: {output_file}")
        if elapsed_time is not None:
            lines.append(f"  Time:   {elapsed_time:.2f}s")
        lines.append("")

        # Statistics section (validator-specific formatting)
        lines.append("  Statistics:")
        lines.extend(self._format_validator_data(record.output_data, record.validator_type, indent="    ", verbose_labels=False))
        lines.append("")

        # Settings used (collapsed for readability)
        if record.input_settings:
            # Only show non-default/interesting settings
            interesting_settings = {
                k: v for k, v in record.input_settings.items()
                if k in ['validation_level', 'threads', 'plasmid_split', 'sort_by_position',
                         'check_invalid_chars', 'min_sequence_length', 'coding_type']
                and v not in [None, False, '', []]
            }

            if interesting_settings:
                lines.append("  Settings:")
                for key, value in sorted(interesting_settings.items()):
                    lines.append(f"    {key}: {value}")
                lines.append("")

        lines.append("-" * 100)
        lines.append("")

        return lines

    def _format_validator_data(
        self,
        data: Dict[str, Any],
        validator_type: str,
        indent: str = "    ",
        verbose_labels: bool = False
    ) -> List[str]:
        """
        Format validator data (statistics or metadata) based on validator type.

        Args:
            data: Dictionary containing validation data
            validator_type: Type of validator ("genome", "read", "feature")
            indent: Indentation string for formatting (default: 4 spaces)
            verbose_labels: Use verbose label style (default: False for concise)

        Returns:
            List of formatted strings
        """
        lines = []

        if validator_type == "genome":
            # Genome-specific data
            if self._has_value(data, 'num_sequences'):
                lines.append(f"{indent}Sequences: {data['num_sequences']:,}")

            if self._has_value(data, 'sequence_ids'):
                seq_ids = data['sequence_ids']
                if verbose_labels:
                    # Metadata style - show IDs sub-indented
                    if len(seq_ids) <= 5:
                        lines.append(f"{indent}  IDs: {', '.join(seq_ids)}")
                    else:
                        lines.append(f"{indent}  IDs: {', '.join(seq_ids[:5])}, ... ({len(seq_ids)} total)")
                else:
                    # Statistics style - inline
                    if len(seq_ids) <= 3:
                        lines.append(f"{indent}Sequence IDs: {', '.join(seq_ids)}")
                    else:
                        lines.append(f"{indent}Sequence IDs: {seq_ids[0]}, {seq_ids[1]}, ... (+{len(seq_ids)-2} more)")

            if self._has_value(data, 'sequence_lengths'):
                seq_lengths = data['sequence_lengths']
                if isinstance(seq_lengths, dict):
                    total_len = sum(seq_lengths.values())
                else:
                    total_len = sum(seq_lengths)
                lines.append(f"{indent}Total Length: {total_len:,} bp")

        elif validator_type == "read":
            # Read-specific data
            if self._has_value(data, 'num_reads'):
                lines.append(f"{indent}Reads: {data['num_reads']:,}")

            if self._has_value(data, 'ngs_type_detected'):
                lines.append(f"{indent}NGS Type: {data['ngs_type_detected']}")

            if self._has_value(data, 'base_name') and self._has_value(data, 'read_number'):
                lines.append(f"{indent}Paired-End: R{data['read_number']} (base: {data['base_name']})")

            if self._has_value(data, 'total_bases'):
                label = "Total Bases" if not verbose_labels else "Total Bases"
                lines.append(f"{indent}{label}: {data['total_bases']:,} bp")

            if self._has_value(data, 'mean_read_length'):
                label = "Mean Length" if not verbose_labels else "Mean Read Length"
                lines.append(f"{indent}{label}: {data['mean_read_length']:.1f} bp")

            if self._has_value(data, 'n50'):
                lines.append(f"{indent}N50: {data['n50']:,} bp")

            if verbose_labels:
                # Metadata style - separate longest/shortest
                if self._has_value(data, 'longest_read_length'):
                    lines.append(f"{indent}Longest Read: {data['longest_read_length']:,} bp")
                if self._has_value(data, 'shortest_read_length'):
                    lines.append(f"{indent}Shortest Read: {data['shortest_read_length']:,} bp")
            else:
                # Statistics style - combined range
                if self._has_value(data, 'longest_read_length') and self._has_value(data, 'shortest_read_length'):
                    lines.append(f"{indent}Length Range: {data['shortest_read_length']:,} - {data['longest_read_length']:,} bp")

        elif validator_type == "feature":
            # Feature-specific data
            if self._has_value(data, 'num_features'):
                lines.append(f"{indent}Features: {data['num_features']:,}")

            if self._has_value(data, 'feature_types'):
                types = data['feature_types']
                if verbose_labels:
                    # Metadata style - inline
                    lines.append(f"{indent}  Types: {', '.join(types)}")
                else:
                    # Statistics style - with truncation
                    if len(types) <= 5:
                        lines.append(f"{indent}Types: {', '.join(types)}")
                    else:
                        lines.append(f"{indent}Types: {', '.join(list(types)[:5])}, ... (+{len(types)-5} more)")

            if self._has_value(data, 'sequence_ids'):
                seq_ids = data['sequence_ids']
                if verbose_labels:
                    # Metadata style
                    if len(seq_ids) <= 5:
                        lines.append(f"{indent}  Sequences: {', '.join(seq_ids)}")
                    else:
                        lines.append(f"{indent}  Sequences: {', '.join(seq_ids[:5])}, ... ({len(seq_ids)} total)")
                else:
                    # Statistics style
                    if len(seq_ids) <= 3:
                        lines.append(f"{indent}Sequences: {', '.join(seq_ids)}")
                    else:
                        lines.append(f"{indent}Sequences: {seq_ids[0]}, {seq_ids[1]}, ... (+{len(seq_ids)-2} more)")

        return lines

    def _generate_interfile_section(self) -> List[str]:
        """Generate inter-file validation section."""
        lines = []
        lines.append("=" * 100)
        lines.append("INTER-FILE VALIDATION RESULTS".center(100))
        lines.append("=" * 100)
        lines.append("")

        for idx, record in enumerate(self.interfile_records, 1):
            lines.extend(self._format_interfile_record(idx, record))

        return lines

    def _format_interfile_record(self, idx: int, record: InterFileValidationRecord) -> List[str]:
        """Format inter-file validation record."""
        lines = []

        # Title with nice formatting
        if record.validation_type == "genomexgenome":
            title = "Genome ↔ Genome Validation"
        elif record.validation_type == "readxread":
            title = "Read ↔ Read Validation (Paired-End)"
        elif record.validation_type == "featurexgenome":
            title = "Feature ↔ Genome Validation"
        else:
            title = f"Inter-file Validation ({record.validation_type})"

        lines.append(f"[{idx}] {title}")
        lines.append("")

        # Status with visual indicator
        status_symbol = "✓" if record.passed else "✗"
        status_text = "PASSED" if record.passed else "FAILED"
        lines.append(f"  Status: {status_symbol} {status_text}")
        lines.append("")

        # Errors (with indentation and clear formatting)
        if record.errors:
            lines.append(f"  Errors ({len(record.errors)}):")
            for i, error in enumerate(record.errors, 1):
                # Multi-line errors with proper indentation
                error_lines = error.split('\n')
                lines.append(f"    {i}. {error_lines[0]}")
                for extra_line in error_lines[1:]:
                    lines.append(f"       {extra_line}")
            lines.append("")

        # Warnings (with indentation and clear formatting)
        if record.warnings:
            lines.append(f"  Warnings ({len(record.warnings)}):")
            for i, warning in enumerate(record.warnings, 1):
                warning_lines = warning.split('\n')
                lines.append(f"    {i}. {warning_lines[0]}")
                for extra_line in warning_lines[1:]:
                    lines.append(f"       {extra_line}")
            lines.append("")

        # Details/Metadata (formatted nicely)
        if record.metadata:
            lines.append("  Details:")
            for key, value in sorted(record.metadata.items()):
                # Format lists nicely
                if isinstance(value, list):
                    if len(value) == 0:
                        lines.append(f"    {key}: (none)")
                    elif len(value) <= 3:
                        lines.append(f"    {key}: {', '.join(str(v) for v in value)}")
                    else:
                        lines.append(f"    {key}: {len(value)} items ({', '.join(str(v) for v in value[:2])}, ...)")
                # Format dicts
                elif isinstance(value, dict):
                    if len(str(value)) > 80:
                        lines.append(f"    {key}: {len(value)} entries")
                    else:
                        lines.append(f"    {key}: {value}")
                # Format large numbers
                elif isinstance(value, int) and value > 9999:
                    lines.append(f"    {key}: {value:,}")
                else:
                    lines.append(f"    {key}: {value}")
            lines.append("")

        lines.append("-" * 100)
        lines.append("")
        return lines
