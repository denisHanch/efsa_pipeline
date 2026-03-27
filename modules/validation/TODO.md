## In progress / partial
- when loading required params in config fails, the validation must fails — `ConfigManager.load()` returns 1 on failure ✓, but required validation failures (ref_genome, reads) still exit 0; needs `sys.exit(1)` when `fatal_errors` is non-empty
- unify upper case for validation_level, debug etc. — `logging_level` is `.upper()`-ed ✓, but `validation_level` is not `.lower()`-ed before the `in VALID_LEVELS` check

## Remaining
- add nextflow parameter `overall_status` that returns whether a fatal failure (for required files) appeared — `NextflowParams` dataclass has no such field yet
- optimalizuj kopirovani souboru — currently plain `shutil.copy2`; consider `os.sendfile` or chunked streaming for large genomic files
- update minimap2 version
