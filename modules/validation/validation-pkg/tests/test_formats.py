"""
Unit tests for the formats module.

Tests format enums including:
- CodingType enum (compression types)
- GenomeFormat enum
- ReadFormat enum
- FeatureFormat enum
- _missing_ method for flexible input
- to_biopython() conversion
"""

import pytest
from pathlib import Path

from validation_pkg.utils.formats import (
    CodingType,
    GenomeFormat,
    ReadFormat,
    FeatureFormat
)


class TestCodingType:
    """Tests for CodingType enum."""

    def test_enum_values(self):
        """Test that enum values are correct."""
        assert CodingType.GZIP.value == "gzip"
        assert CodingType.BZIP2.value == "bzip2"
        assert CodingType.NONE.value == "none"

    def test_direct_enum_access(self):
        """Test direct enum member access."""
        assert CodingType.GZIP == CodingType.GZIP
        assert CodingType.BZIP2 == CodingType.BZIP2
        assert CodingType.NONE == CodingType.NONE

    def test_missing_with_short_extensions(self):
        """Test _missing_ method with short extensions."""
        assert CodingType('gz') == CodingType.GZIP
        assert CodingType('bz2') == CodingType.BZIP2
        assert CodingType('gzip') == CodingType.GZIP
        assert CodingType('bzip2') == CodingType.BZIP2

    def test_missing_with_dot_extensions(self):
        """Test _missing_ method with dot-prefixed extensions."""
        assert CodingType('.gz') == CodingType.GZIP
        assert CodingType('.bz2') == CodingType.BZIP2
        assert CodingType('.gzip') == CodingType.GZIP
        assert CodingType('.bzip2') == CodingType.BZIP2

    def test_missing_with_filenames(self):
        """Test _missing_ method with full filenames."""
        assert CodingType('genome.fasta.gz') == CodingType.GZIP
        assert CodingType('reads.fastq.bz2') == CodingType.BZIP2
        assert CodingType('features.gff.gz') == CodingType.GZIP

    def test_missing_with_plain_files(self):
        """Test _missing_ method with uncompressed filenames."""
        assert CodingType('genome.fasta') == CodingType.NONE
        assert CodingType('reads.fastq') == CodingType.NONE
        assert CodingType('file.txt') == CodingType.NONE

    def test_missing_case_insensitive(self):
        """Test that _missing_ is case-insensitive."""
        assert CodingType('GZ') == CodingType.GZIP
        assert CodingType('BZ2') == CodingType.BZIP2
        assert CodingType('GZIP') == CodingType.GZIP
        assert CodingType('.GZ') == CodingType.GZIP
        assert CodingType('genome.FASTA.GZ') == CodingType.GZIP

    def test_missing_with_empty_string(self):
        """Test _missing_ with empty string defaults to NONE."""
        assert CodingType('') == CodingType.NONE

    def test_missing_with_invalid_value(self):
        """Test _missing_ with invalid value defaults to NONE."""
        assert CodingType('invalid') == CodingType.NONE
        assert CodingType('xyz') == CodingType.NONE

    def test_missing_with_path_object(self):
        """Test _missing_ with Path objects."""
        assert CodingType(Path('genome.fasta.gz')) == CodingType.GZIP
        assert CodingType(Path('reads.fastq.bz2')) == CodingType.BZIP2

    def test_missing_strips_whitespace(self):
        """Test that _missing_ strips whitespace."""
        assert CodingType('  gz  ') == CodingType.GZIP
        assert CodingType('  .bz2  ') == CodingType.BZIP2


class TestGenomeFormat:
    """Tests for GenomeFormat enum."""

    def test_enum_values(self):
        """Test that enum values are correct."""
        assert GenomeFormat.FASTA.value == "fasta"
        assert GenomeFormat.GENBANK.value == "genbank"

    def test_to_biopython(self):
        """Test to_biopython() method."""
        assert GenomeFormat.FASTA.to_biopython() == "fasta"
        assert GenomeFormat.GENBANK.to_biopython() == "genbank"

    def test_missing_with_format_names(self):
        """Test _missing_ with format names."""
        assert GenomeFormat('fasta') == GenomeFormat.FASTA
        assert GenomeFormat('genbank') == GenomeFormat.GENBANK

    def test_missing_with_fasta_extensions(self):
        """Test _missing_ with FASTA extensions."""
        assert GenomeFormat('fa') == GenomeFormat.FASTA
        assert GenomeFormat('fasta') == GenomeFormat.FASTA
        assert GenomeFormat('fna') == GenomeFormat.FASTA
        assert GenomeFormat('faa') == GenomeFormat.FASTA
        assert GenomeFormat('.fa') == GenomeFormat.FASTA
        assert GenomeFormat('.fasta') == GenomeFormat.FASTA

    def test_missing_with_genbank_extensions(self):
        """Test _missing_ with GenBank extensions."""
        assert GenomeFormat('gb') == GenomeFormat.GENBANK
        assert GenomeFormat('gbk') == GenomeFormat.GENBANK
        assert GenomeFormat('genbank') == GenomeFormat.GENBANK
        assert GenomeFormat('.gb') == GenomeFormat.GENBANK
        assert GenomeFormat('.gbk') == GenomeFormat.GENBANK

    def test_missing_with_filenames(self):
        """Test _missing_ with full filenames."""
        assert GenomeFormat('genome.fasta') == GenomeFormat.FASTA
        assert GenomeFormat('genome.fa') == GenomeFormat.FASTA
        assert GenomeFormat('genome.gb') == GenomeFormat.GENBANK
        assert GenomeFormat('genome.gbk') == GenomeFormat.GENBANK

    def test_missing_with_compressed_filenames(self):
        """Test _missing_ with compressed filenames.

        Note: Format enums extract the last extension only, so compressed files
        (genome.fasta.gz) will fail since .gz is not a valid format. Users should
        use file_handler utilities (get_base_filename) to strip compression first.
        """
        with pytest.raises(ValueError):
            GenomeFormat('genome.fasta.gz')

        with pytest.raises(ValueError):
            GenomeFormat('genome.gb.bz2')

    def test_missing_case_insensitive(self):
        """Test that _missing_ is case-insensitive."""
        assert GenomeFormat('FASTA') == GenomeFormat.FASTA
        assert GenomeFormat('FA') == GenomeFormat.FASTA
        assert GenomeFormat('GENBANK') == GenomeFormat.GENBANK
        assert GenomeFormat('GB') == GenomeFormat.GENBANK

    def test_missing_with_invalid_value_raises_error(self):
        """Test that invalid values raise ValueError."""
        with pytest.raises(ValueError) as exc_info:
            GenomeFormat('invalid')
        assert "is not a valid GenomeFormat" in str(exc_info.value)

        with pytest.raises(ValueError):
            GenomeFormat('txt')

        with pytest.raises(ValueError):
            GenomeFormat('fastq')

    def test_missing_strips_whitespace(self):
        """Test that _missing_ strips whitespace."""
        assert GenomeFormat('  fasta  ') == GenomeFormat.FASTA
        assert GenomeFormat('  .gb  ') == GenomeFormat.GENBANK


class TestReadFormat:
    """Tests for ReadFormat enum."""

    def test_enum_values(self):
        """Test that enum values are correct."""
        assert ReadFormat.FASTQ.value == "fastq"
        assert ReadFormat.BAM.value == "bam"

    def test_to_biopython(self):
        """Test to_biopython() method."""
        assert ReadFormat.FASTQ.to_biopython() == "fastq"
        assert ReadFormat.BAM.to_biopython() == "bam"

    def test_missing_with_format_names(self):
        """Test _missing_ with format names."""
        assert ReadFormat('fastq') == ReadFormat.FASTQ
        assert ReadFormat('bam') == ReadFormat.BAM

    def test_missing_with_fastq_extensions(self):
        """Test _missing_ with FASTQ extensions."""
        assert ReadFormat('fq') == ReadFormat.FASTQ
        assert ReadFormat('fastq') == ReadFormat.FASTQ
        assert ReadFormat('.fq') == ReadFormat.FASTQ
        assert ReadFormat('.fastq') == ReadFormat.FASTQ

    def test_missing_with_bam_extensions(self):
        """Test _missing_ with BAM extensions."""
        assert ReadFormat('bam') == ReadFormat.BAM
        assert ReadFormat('.bam') == ReadFormat.BAM

    def test_missing_with_filenames(self):
        """Test _missing_ with full filenames."""
        assert ReadFormat('reads.fastq') == ReadFormat.FASTQ
        assert ReadFormat('reads.fq') == ReadFormat.FASTQ
        assert ReadFormat('aligned.bam') == ReadFormat.BAM

    def test_missing_with_compressed_filenames(self):
        """Test _missing_ with compressed filenames.

        Note: Format enums extract the last extension only, so compressed files
        (reads.fastq.gz) will fail since .gz is not a valid format. Users should
        use file_handler utilities (get_base_filename) to strip compression first.
        """
        with pytest.raises(ValueError):
            ReadFormat('reads.fastq.gz')

        with pytest.raises(ValueError):
            ReadFormat('reads.fq.bz2')

    def test_missing_case_insensitive(self):
        """Test that _missing_ is case-insensitive."""
        assert ReadFormat('FASTQ') == ReadFormat.FASTQ
        assert ReadFormat('FQ') == ReadFormat.FASTQ
        assert ReadFormat('BAM') == ReadFormat.BAM

    def test_missing_with_invalid_value_raises_error(self):
        """Test that invalid values raise ValueError."""
        with pytest.raises(ValueError) as exc_info:
            ReadFormat('invalid')
        assert "is not a valid ReadFormat" in str(exc_info.value)

        with pytest.raises(ValueError):
            ReadFormat('fasta')

        with pytest.raises(ValueError):
            ReadFormat('txt')

    def test_missing_strips_whitespace(self):
        """Test that _missing_ strips whitespace."""
        assert ReadFormat('  fastq  ') == ReadFormat.FASTQ
        assert ReadFormat('  .bam  ') == ReadFormat.BAM


class TestFeatureFormat:
    """Tests for FeatureFormat enum."""

    def test_enum_values(self):
        """Test that enum values are correct."""
        assert FeatureFormat.GFF.value == "gff"
        assert FeatureFormat.BED.value == "bed"

    def test_to_biopython(self):
        """Test to_biopython() method."""
        assert FeatureFormat.GFF.to_biopython() == "gff"
        assert FeatureFormat.BED.to_biopython() == "bed"

    def test_missing_with_format_names(self):
        """Test _missing_ with format names."""
        assert FeatureFormat('gff') == FeatureFormat.GFF
        assert FeatureFormat('bed') == FeatureFormat.BED

    def test_missing_with_gff_extensions(self):
        """Test _missing_ with GFF extensions."""
        assert FeatureFormat('gff') == FeatureFormat.GFF
        assert FeatureFormat('gff3') == FeatureFormat.GFF
        assert FeatureFormat('gff2') == FeatureFormat.GTF
        assert FeatureFormat('gtf') == FeatureFormat.GTF
        assert FeatureFormat('.gff') == FeatureFormat.GFF
        assert FeatureFormat('.gff3') == FeatureFormat.GFF
        assert FeatureFormat('.gtf') == FeatureFormat.GTF

    def test_missing_with_bed_extensions(self):
        """Test _missing_ with BED extensions."""
        assert FeatureFormat('bed') == FeatureFormat.BED
        assert FeatureFormat('.bed') == FeatureFormat.BED

    def test_missing_with_filenames(self):
        """Test _missing_ with full filenames."""
        assert FeatureFormat('features.gff') == FeatureFormat.GFF
        assert FeatureFormat('features.gff3') == FeatureFormat.GFF
        assert FeatureFormat('genes.gtf') == FeatureFormat.GTF
        assert FeatureFormat('regions.bed') == FeatureFormat.BED

    def test_missing_with_compressed_filenames(self):
        """Test _missing_ with compressed filenames.

        Note: Format enums extract the last extension only, so compressed files
        (features.gff.gz) will fail since .gz is not a valid format. Users should
        use file_handler utilities (get_base_filename) to strip compression first.
        """
        with pytest.raises(ValueError):
            FeatureFormat('features.gff.gz')

        with pytest.raises(ValueError):
            FeatureFormat('genes.gtf.bz2')

        with pytest.raises(ValueError):
            FeatureFormat('regions.bed.gz')

    def test_missing_case_insensitive(self):
        """Test that _missing_ is case-insensitive."""
        assert FeatureFormat('GFF') == FeatureFormat.GFF
        assert FeatureFormat('GTF') == FeatureFormat.GTF
        assert FeatureFormat('BED') == FeatureFormat.BED

    def test_missing_with_invalid_value_raises_error(self):
        """Test that invalid values raise ValueError."""
        with pytest.raises(ValueError) as exc_info:
            FeatureFormat('invalid')
        assert "is not a valid FeatureFormat" in str(exc_info.value)

        with pytest.raises(ValueError):
            FeatureFormat('fasta')

        with pytest.raises(ValueError):
            FeatureFormat('txt')

    def test_missing_strips_whitespace(self):
        """Test that _missing_ strips whitespace."""
        assert FeatureFormat('  gff  ') == FeatureFormat.GFF
        assert FeatureFormat('  .bed  ') == FeatureFormat.BED


class TestEnumComparison:
    """Test comparison operations between enum members."""

    def test_equality(self):
        """Test that enum members compare equal correctly."""
        assert CodingType.GZIP == CodingType.GZIP
        assert GenomeFormat.FASTA == GenomeFormat.FASTA
        assert ReadFormat.FASTQ == ReadFormat.FASTQ
        assert FeatureFormat.GFF == FeatureFormat.GFF

    def test_inequality(self):
        """Test that different enum members are not equal."""
        assert CodingType.GZIP != CodingType.BZIP2
        assert GenomeFormat.FASTA != GenomeFormat.GENBANK
        assert ReadFormat.FASTQ != ReadFormat.BAM
        assert FeatureFormat.GFF != FeatureFormat.BED

    def test_identity(self):
        """Test that enum members are singletons."""
        assert CodingType('gz') is CodingType.GZIP
        assert GenomeFormat('fasta') is GenomeFormat.FASTA
        assert ReadFormat('fastq') is ReadFormat.FASTQ
        assert FeatureFormat('gff') is FeatureFormat.GFF


class TestEnumIteration:
    """Test iteration over enum members."""

    def test_coding_type_iteration(self):
        """Test iterating over CodingType members."""
        members = list(CodingType)
        assert len(members) == 3
        assert CodingType.GZIP in members
        assert CodingType.BZIP2 in members
        assert CodingType.NONE in members

    def test_genome_format_iteration(self):
        """Test iterating over GenomeFormat members."""
        members = list(GenomeFormat)
        assert len(members) == 2
        assert GenomeFormat.FASTA in members
        assert GenomeFormat.GENBANK in members

    def test_read_format_iteration(self):
        """Test iterating over ReadFormat members."""
        members = list(ReadFormat)
        assert len(members) == 2
        assert ReadFormat.FASTQ in members
        assert ReadFormat.BAM in members

    def test_feature_format_iteration(self):
        """Test iterating over FeatureFormat members."""
        members = list(FeatureFormat)
        assert len(members) == 3
        assert FeatureFormat.GFF in members
        assert FeatureFormat.GTF in members
        assert FeatureFormat.BED in members


class TestEnumStringRepresentation:
    """Test string representation of enum members."""

    def test_str(self):
        """Test __str__ method."""
        assert str(CodingType.GZIP) == "CodingType.GZIP"
        assert str(GenomeFormat.FASTA) == "GenomeFormat.FASTA"
        assert str(ReadFormat.FASTQ) == "ReadFormat.FASTQ"
        assert str(FeatureFormat.GFF) == "FeatureFormat.GFF"

    def test_repr(self):
        """Test __repr__ method."""
        assert repr(CodingType.GZIP) == "<CodingType.GZIP: 'gzip'>"
        assert repr(GenomeFormat.FASTA) == "<GenomeFormat.FASTA: 'fasta'>"
        assert repr(ReadFormat.FASTQ) == "<ReadFormat.FASTQ: 'fastq'>"
        assert repr(FeatureFormat.GFF) == "<FeatureFormat.GFF: 'gff'>"

    def test_value_access(self):
        """Test accessing the value attribute."""
        assert CodingType.GZIP.value == "gzip"
        assert GenomeFormat.FASTA.value == "fasta"
        assert ReadFormat.FASTQ.value == "fastq"
        assert FeatureFormat.GFF.value == "gff"
