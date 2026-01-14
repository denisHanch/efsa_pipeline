"""Validation report generation module."""

from dataclasses import dataclass, field
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional, Any, Union
import json
from validation_pkg.utils.path_utils import get_incremented_path
from validation_pkg.validators.genome_validator import GenomeOutputMetadata
from validation_pkg.validators.read_validator import ReadOutputMetadata
from validation_pkg.validators.feature_validator import FeatureOutputMetadata


@dataclass
class FileValidationRecord:
    """Complete record of a single file's validation."""
    # Input information
    output_data: dict
    validator_type: str  # "genome", "read", "feature"
    input_settings: Optional[dict] = None  # Serialized Settings object

    def _get_metadata(self):
        """Get OutputMetadata object from output_data."""
        from dataclasses import fields

        if self.validator_type == "genome":
            metadata_class = GenomeOutputMetadata
        elif self.validator_type == "read":
            metadata_class = ReadOutputMetadata
        elif self.validator_type == "feature":
            metadata_class = FeatureOutputMetadata
        else:
            raise ValueError(f"Unknown validator_type: {self.validator_type}")

        # Filter output_data to only include fields that exist in the metadata class
        valid_fields = {f.name for f in fields(metadata_class)}
        filtered_data = {k: v for k, v in self.output_data.items() if k in valid_fields}

        return metadata_class(**filtered_data)

    def format_title(self, idx: int) -> list[str]:
        """Format the title section for this file validation."""
        lines = []

        validator_label = self.validator_type.upper()

        # Check if this is a plasmid file
        is_plasmid = (self.input_settings
                     and self.input_settings.get('is_plasmid') == True)

        if is_plasmid:
            lines.append(f"[{idx}] {validator_label} FILE (Plasmid)")
        else:
            lines.append(f"[{idx}] {validator_label} FILE")
        lines.append("")

        return lines

    def format_common_fields(self, indent: str = "  ") -> list[str]:
        """Format common fields (input, output, time)."""
        metadata = self._get_metadata()
        lines = metadata.format_common_fields(indent=indent)
        lines.append("")
        return lines

    def format_statistics(self, indent: str = "  ") -> list[str]:
        """Format validator-specific statistics."""
        lines = []
        lines.append(f"{indent}Statistics:")

        metadata = self._get_metadata()
        # Statistics use double indent (4 spaces by default)
        stat_indent = indent + "  "
        lines.extend(metadata.format_statistics(indent=stat_indent, input_settings=self.input_settings))
        lines.append("")

        return lines

    def format_settings(self, indent: str = "  ") -> list[str]:
        """Format settings section (only interesting/non-default settings)."""
        lines = []

        if self.input_settings:
            # Only show non-default/interesting settings
            interesting_settings = {
                k: v for k, v in self.input_settings.items()
                if k in ['validation_level', 'threads', 'plasmid_split', 'sort_by_position',
                         'check_invalid_chars', 'min_sequence_length', 'coding_type']
                and v not in [None, False, '', []]
            }

            if interesting_settings:
                lines.append(f"{indent}Settings:")
                settings_indent = indent + "  "
                for key, value in sorted(interesting_settings.items()):
                    lines.append(f"{settings_indent}{key}: {value}")
                lines.append("")

        return lines


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

    def format_title(self, idx: int) -> list[str]:
        """Format the title section for this inter-file validation."""
        lines = []

        # Title with nice formatting
        if self.validation_type == "genomexgenome":
            title = "Genome ↔ Genome Validation"
        elif self.validation_type == "readxread":
            title = "Read ↔ Read Validation (Paired-End)"
        elif self.validation_type == "featurexgenome":
            title = "Feature ↔ Genome Validation"
        else:
            title = f"Inter-file Validation ({self.validation_type})"

        lines.append(f"[{idx}] {title}")
        lines.append("")

        return lines

    def format_status(self, indent: str = "  ") -> list[str]:
        """Format the status section (passed/failed indicator)."""
        lines = []

        status_symbol = "✓" if self.passed else "✗"
        status_text = "PASSED" if self.passed else "FAILED"
        lines.append(f"{indent}Status: {status_symbol} {status_text}")
        lines.append("")

        return lines

    def format_errors(self, indent: str = "  ") -> list[str]:
        """Format the errors section."""
        lines = []

        if self.errors:
            lines.append(f"{indent}Errors ({len(self.errors)}):")
            for i, error in enumerate(self.errors, 1):
                # Multi-line errors with proper indentation
                error_lines = error.split('\n')
                lines.append(f"{indent}  {i}. {error_lines[0]}")
                for extra_line in error_lines[1:]:
                    lines.append(f"{indent}     {extra_line}")
            lines.append("")

        return lines

    def format_warnings(self, indent: str = "  ") -> list[str]:
        """Format the warnings section."""
        lines = []

        if self.warnings:
            lines.append(f"{indent}Warnings ({len(self.warnings)}):")
            for i, warning in enumerate(self.warnings, 1):
                warning_lines = warning.split('\n')
                lines.append(f"{indent}  {i}. {warning_lines[0]}")
                for extra_line in warning_lines[1:]:
                    lines.append(f"{indent}     {extra_line}")
            lines.append("")

        return lines

    def format_metadata(self, indent: str = "  ") -> list[str]:
        """Format the metadata/details section."""
        lines = []

        if self.metadata:
            lines.append(f"{indent}Details:")
            for key, value in sorted(self.metadata.items()):
                # Format lists nicely
                if isinstance(value, list):
                    if len(value) == 0:
                        lines.append(f"{indent}  {key}: (none)")
                    elif len(value) <= 3:
                        lines.append(f"{indent}  {key}: {', '.join(str(v) for v in value)}")
                    else:
                        lines.append(f"{indent}  {key}: {len(value)} items ({', '.join(str(v) for v in value[:2])}, ...)")
                # Format dicts
                elif isinstance(value, dict):
                    if len(str(value)) > 80:
                        lines.append(f"{indent}  {key}: {len(value)} entries")
                    else:
                        lines.append(f"{indent}  {key}: {value}")
                # Format large numbers
                elif isinstance(value, int) and value > 9999:
                    lines.append(f"{indent}  {key}: {value:,}")
                else:
                    lines.append(f"{indent}  {key}: {value}")
            lines.append("")

        return lines


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

        # Count files by type, separating plasmids from genomes
        plasmid_count = sum(1 for r in self.file_records
                           if r.validator_type == "genome"
                           and r.input_settings
                           and r.input_settings.get('is_plasmid') == True)

        genome_count = sum(1 for r in self.file_records
                          if r.validator_type == "genome"
                          and (not r.input_settings or not r.input_settings.get('is_plasmid')))

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

        # Build file type list dynamically to determine tree characters
        file_types = []
        if genome_count > 0:
            file_types.append(('Genomes', genome_count))
        if plasmid_count > 0:
            file_types.append(('Plasmids', plasmid_count))
        if read_count > 0:
            file_types.append(('Reads', read_count))
        if feature_count > 0:
            file_types.append(('Features', feature_count))

        # Display file types with proper tree characters
        for i, (label, count) in enumerate(file_types):
            is_last = (i == len(file_types) - 1)
            tree_char = "└─" if is_last else "├─"
            lines.append(f"    {tree_char} {label}: {count:>2}")

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

        # Count plasmid files separately from genomes
        plasmid_count = sum(1 for r in self.file_records
                           if r.validator_type == "genome"
                           and r.input_settings
                           and r.input_settings.get('is_plasmid') == True)

        genome_count = sum(1 for r in self.file_records
                          if r.validator_type == "genome"
                          and (not r.input_settings or not r.input_settings.get('is_plasmid')))

        read_count = sum(1 for r in self.file_records if r.validator_type == "read")
        feature_count = sum(1 for r in self.file_records if r.validator_type == "feature")

        interfile_passed = sum(1 for r in self.interfile_records if r.passed)
        interfile_failed = len(self.interfile_records) - interfile_passed

        return {
            "total_files": len(self.file_records),
            "genome_files": genome_count,
            "plasmid_files": plasmid_count,
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

        # Title
        lines.extend(record.format_title(idx))

        # Common fields (input, output, time)
        lines.extend(record.format_common_fields(indent="  "))

        # Statistics
        lines.extend(record.format_statistics(indent="  "))

        # Settings
        lines.extend(record.format_settings(indent="  "))

        lines.append("-" * 100)
        lines.append("")

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

        # Title
        lines.extend(record.format_title(idx))

        # Status
        lines.extend(record.format_status(indent="  "))

        # Errors
        lines.extend(record.format_errors(indent="  "))

        # Warnings
        lines.extend(record.format_warnings(indent="  "))

        # Metadata/Details
        lines.extend(record.format_metadata(indent="  "))

        lines.append("-" * 100)
        lines.append("")
        return lines
