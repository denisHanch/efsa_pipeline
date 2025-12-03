"""
Validation report generation module.

Provides comprehensive reporting for validation results including
input settings, output metadata, and inter-file validation results.
"""

from dataclasses import dataclass, field
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional, Any, Union
import json


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
    """
    Comprehensive validation report builder.

    Collects results from validators and inter-file validation,
    tracks input settings and output metadata, then generates a
    detailed report when flush() is called.
    """

    def __init__(self, report_path: Path):
        """
        Initialize validation report.

        Args:
            report_path: Path where the report will be written
        """
        from validation_pkg.utils.file_handler import get_incremented_path

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
        """
        Check if a key exists in dictionary and its value is not None.

        Args:
            data: Dictionary to check
            key: Key name to look for

        Returns:
            True if key exists and value is not None, False otherwise
        """
        return key in data and data[key] is not None

    def write(
        self,
        result: Union[Any, List[Any]],
        file_type: str,
        settings: Optional[Any] = None
    ) -> None:
        """
        Add a validation result to the report.

        Args:
            result: Result dict from validator or inter-file validation
            settings: Settings object used for validation (for file validators)
            file_type: Type of result:
                - "genome": GenomeValidator result
                - "read": ReadValidator result (can be single or list)
                - "feature": FeatureValidator result
                - "genomexgenome": Genome inter-file validation
                - "readxread": Read inter-file validation
                - "featurexgenome": Feature-genome inter-file validation
        """
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
        """
        Generate and write the final report to file.

        Args:
            format: Report format - "text" (default) or "json"

        Creates a comprehensive report with:
        - Summary statistics
        - Per-file validation results with settings
        - Inter-file validation results
        - Issue listings
        """
        if format == "text":
            self._flush_text()
        elif format == "json":
            self._flush_json()
        else:
            raise ValueError(f"Unknown format: {format}. Use 'text' or 'json'")

    def _flush_text(self) -> None:
        """Generate and write text format report."""
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
        """Generate and write JSON format report."""
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
        """Generate summary section of report."""
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
        """Get summary data for JSON export."""
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
        """Generate file-specific results section."""
        lines = []
        lines.append("=" * 100)
        lines.append("FILE VALIDATION RESULTS".center(100))
        lines.append("=" * 100)
        lines.append("")

        for idx, record in enumerate(self.file_records, 1):
            lines.extend(self._format_file_record(idx, record))

        return lines

    def _format_file_record(self, idx: int, record: FileValidationRecord) -> List[str]:
        """Format a single file validation record."""
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
        lines.extend(self._format_statistics(record.output_data, record.validator_type))
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

    def _format_statistics(self, output_data: Dict[str, Any], validator_type: str) -> List[str]:
        """Format statistics from output_data based on validator type."""
        lines = []

        if validator_type == "genome":
            # Genome-specific statistics
            if self._has_value(output_data, 'num_sequences'):
                lines.append(f"    Sequences: {output_data['num_sequences']:,}")

            if self._has_value(output_data, 'sequence_lengths'):
                seq_lengths = output_data['sequence_lengths']
                if isinstance(seq_lengths, dict):
                    total_len = sum(seq_lengths.values())
                else:
                    total_len = sum(seq_lengths)
                lines.append(f"    Total Length: {total_len:,} bp")

            if self._has_value(output_data, 'sequence_ids'):
                seq_ids = output_data['sequence_ids']
                if len(seq_ids) <= 3:
                    lines.append(f"    Sequence IDs: {', '.join(seq_ids)}")
                else:
                    lines.append(f"    Sequence IDs: {seq_ids[0]}, {seq_ids[1]}, ... (+{len(seq_ids)-2} more)")

        elif validator_type == "read":
            # Read-specific statistics
            if self._has_value(output_data, 'num_reads'):
                lines.append(f"    Reads: {output_data['num_reads']:,}")

            if self._has_value(output_data, 'total_bases'):
                lines.append(f"    Total Bases: {output_data['total_bases']:,} bp")

            if self._has_value(output_data, 'mean_read_length'):
                lines.append(f"    Mean Length: {output_data['mean_read_length']:.1f} bp")

            if self._has_value(output_data, 'n50'):
                lines.append(f"    N50: {output_data['n50']:,} bp")

            if self._has_value(output_data, 'longest_read_length') and self._has_value(output_data, 'shortest_read_length'):
                lines.append(f"    Length Range: {output_data['shortest_read_length']:,} - {output_data['longest_read_length']:,} bp")

            if self._has_value(output_data, 'ngs_type_detected'):
                ngs_type = output_data['ngs_type_detected']
                lines.append(f"    NGS Type: {ngs_type}")

            if self._has_value(output_data, 'base_name') and self._has_value(output_data, 'read_number'):
                base_name = output_data['base_name']
                read_num = output_data['read_number']
                lines.append(f"    Paired-End: R{read_num} (base: {base_name})")

        elif validator_type == "feature":
            # Feature-specific statistics
            if self._has_value(output_data, 'num_features'):
                lines.append(f"    Features: {output_data['num_features']:,}")

            if self._has_value(output_data, 'feature_types'):
                types = output_data['feature_types']
                if len(types) <= 5:
                    lines.append(f"    Types: {', '.join(types)}")
                else:
                    lines.append(f"    Types: {', '.join(list(types)[:5])}, ... (+{len(types)-5} more)")

            if self._has_value(output_data, 'sequence_ids'):
                seq_ids = output_data['sequence_ids']
                if len(seq_ids) <= 3:
                    lines.append(f"    Sequences: {', '.join(seq_ids)}")
                else:
                    lines.append(f"    Sequences: {seq_ids[0]}, {seq_ids[1]}, ... (+{len(seq_ids)-2} more)")

        return lines

    def _format_metadata(self, metadata: Dict[str, Any], validator_type: str) -> List[str]:
        """Format metadata based on validator type."""
        lines = []

        if validator_type == "genome":
            # Genome-specific metadata
            if self._has_value(metadata, 'num_sequences'):
                lines.append(f"  Sequences: {metadata['num_sequences']}")

            if self._has_value(metadata, 'sequence_ids'):
                seq_ids = metadata['sequence_ids']
                if len(seq_ids) <= 5:
                    lines.append(f"    IDs: {', '.join(seq_ids)}")
                else:
                    lines.append(f"    IDs: {', '.join(seq_ids[:5])}, ... ({len(seq_ids)} total)")

            if self._has_value(metadata, 'sequence_lengths'):
                # Handle both list and dict formats
                seq_lengths = metadata['sequence_lengths']
                if isinstance(seq_lengths, dict):
                    total_len = sum(seq_lengths.values())
                else:  # Assume list
                    total_len = sum(seq_lengths)
                lines.append(f"  Total Length: {total_len:,} bp")

        elif validator_type == "read":
            # Read-specific metadata
            if self._has_value(metadata, 'num_reads'):
                lines.append(f"  Reads: {metadata['num_reads']:,}")

            if self._has_value(metadata, 'ngs_type_detected'):
                lines.append(f"  NGS Type: {metadata['ngs_type_detected']}")

            if self._has_value(metadata, 'base_name') and self._has_value(metadata, 'read_number'):
                lines.append(f"  Paired-End: R{metadata['read_number']} (base: {metadata['base_name']})")

            # Read statistics
            if self._has_value(metadata, 'total_bases'):
                lines.append(f"  Total Bases: {metadata['total_bases']:,} bp")

            if self._has_value(metadata, 'mean_read_length'):
                lines.append(f"  Mean Read Length: {metadata['mean_read_length']:.1f} bp")

            if self._has_value(metadata, 'n50'):
                lines.append(f"  N50: {metadata['n50']:,} bp")

            if self._has_value(metadata, 'longest_read_length'):
                lines.append(f"  Longest Read: {metadata['longest_read_length']:,} bp")

            if self._has_value(metadata, 'shortest_read_length'):
                lines.append(f"  Shortest Read: {metadata['shortest_read_length']:,} bp")

        elif validator_type == "feature":
            # Feature-specific metadata
            if self._has_value(metadata, 'num_features'):
                lines.append(f"  Features: {metadata['num_features']}")

            if self._has_value(metadata, 'feature_types'):
                types = metadata['feature_types']
                lines.append(f"    Types: {', '.join(types)}")

            if self._has_value(metadata, 'sequence_ids'):
                seq_ids = metadata['sequence_ids']
                if len(seq_ids) <= 5:
                    lines.append(f"    Sequences: {', '.join(seq_ids)}")
                else:
                    lines.append(f"    Sequences: {', '.join(seq_ids[:5])}, ... ({len(seq_ids)} total)")

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
        """Format a single inter-file validation record."""
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
