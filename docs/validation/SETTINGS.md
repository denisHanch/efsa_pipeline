# Validation Settings

The validation module supports three validation levels to balance thoroughness and performance.

## Validation Levels

### Comparison Table

| Level | Parsing | Validation | Edits | Output | Speed | Use Case |
|-------|---------|------------|-------|--------|-------|----------|
| **strict** (default) | All data | All data | All applied | BioPython write | Slowest | Structure validation sequence by sequence, statistics gathering |
| **trust** | All data (genome)<br>First record only (reads) | First sequence only | All applied (genome)<br>None (reads) | BioPython write (genome)<br>File copy (reads) | Fast | Trust data, adapt file coding, name and location |
| **minimal** | None | None | None | File copy | Fastest | Rename and move files to meet the requirements |

## Level Details

### Strict Mode (Default)
- Validates every sequence in the file
- Performs comprehensive quality checks
- Generates detailed statistics
- Recommended for first-time data processing
- Ensures highest data quality

**Use when:**
- Processing data for the first time
- Data quality is uncertain
- Detailed statistics are needed

### Trust Mode
- Validates only the first sequence
- Assumes remaining data is consistent
- Faster processing for large files
- Still performs format conversions

**Use when:**
- Data has been pre-validated
- Files are from trusted sources
- You need faster processing times

### Minimal Mode
- No validation performed
- Only renames and moves files
- Fastest processing
- Use with caution

**Use when:**
- Files are already in correct format
- You're absolutely certain of data quality
- Only file organization is needed

## Configuration

Set the validation level in your `config.json`:

```json
{
  "ref_genome_filename": {
    "filename": "reference.fasta",
    "validation_level": "strict",
    "threads": 8
  }
}
```

Or set globally in options:

```json
{
  "options": {
    "validation_level": "trust",
    "threads": 8
  }
}
```

## Performance Impact

| Mode | Relative Speed | Resource Usage |
|------|---------------|----------------|
| minimal | 1x (baseline) | Minimal |
| trust | 10-15x slower | Moderate |
| strict | 100x+ slower | High |

## See Also

- [Performance Tips](PERFORMANCE.md) - Optimization strategies
- [Configuration Guide](CONFIG_GUIDE.md) - Full configuration options
