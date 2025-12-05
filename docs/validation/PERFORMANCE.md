# Performance Tips

Optimize your validation performance with these recommendations.

## Quick Optimization Checklist

1. **Install parallel compression tools**
2. **Use appropriate validation levels**
3. **Configure thread counts**

## Parallel Compression Tools

Install these tools for faster compression/decompression:

```bash
sudo apt-get install pigz pbzip2
```

**Performance gains:**
- `pigz` - parallel gzip (2-4x faster)
- `pbzip2` - parallel bzip2 (3-6x faster)

## Validation Level Selection

Choose the right validation level for your use case:

| Scenario | Recommended Level | Speed Gain |
|----------|------------------|------------|
| First-time data processing | `strict` | Baseline |
| Pre-validated data | `trust` | 10-15x faster |
| Files already in correct format | `minimal` | 100x+ faster |

Example configuration:

```json
{
  "ref_genome_filename": {
    "filename": "reference.fasta",
    "validation_level": "trust",
    "threads": 8
  }
}
```

See [Validation Settings](SETTINGS.md) for detailed information on validation levels.

## Thread Configuration

### Per-File Thread Settings

Configure threads for each file:

```json
{
  "ref_genome_filename": {
    "filename": "reference.fasta",
    "validation_level": "strict",
    "threads": 16
  }
}
```

### Global Thread Settings

Or set globally for all files:

```json
{
  "options": {
    "threads": 16,
    "validation_level": "trust"
  }
}
```

**Thread performance:**
- Strict mode: 3-7x faster with 8+ threads
- Trust mode: Minimal benefit from threading
- Recommended: Use 8-16 threads for strict mode

## Performance Summary

| Optimization | Performance Gain | Effort |
|--------------|-----------------|--------|
| validation_level='trust' | 10-15x faster | Low (config change) |
| validation_level='minimal' | 100x+ faster | Low (config change) |
| threads=16 (strict mode) | 3-7x faster | Low (config change) |
| parallel compression tools | 2-6x faster | Medium (installation) |

## Best Practices

1. **For production:** Use `strict` mode first, then `trust` for subsequent runs
2. **For development:** Use `trust` or `minimal` modes
3. **For large datasets:** Always use parallel compression tools and maximum threads
4. **Monitor resources:** Check CPU and memory usage with `htop`

## Logging Performance

Reduce logging verbosity for slight performance improvements:

```json
{
  "options": {
    "logging_level": "WARNING"
  }
}
```

Logging levels (from most to least verbose):
- `DEBUG` - Most detailed
- `INFO` - Standard
- `WARNING` - Warnings only
- `ERROR` - Errors only

## See Also

- [Validation Settings](SETTINGS.md) - Detailed validation level information
- [Configuration Guide](CONFIG_GUIDE.md) - Complete configuration options
