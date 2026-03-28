#   How to run

Load your data and configuration file into `/EFSA_workspace/data/inputs/` (description of config.json bellow)

Run validation as 
```
validate
```

or use specific config path with
```
validate --config=path/to/config.json
```

output of the validation is on `/EFSA_workspace/data/valid/` 

All outputs go to `data/valid/`. Previous runs are preserved via auto-incremented filenames (`validation_001.log`, `report_1.txt`, etc.).

#   What to export
-   log: `data/valid/validation.log`
-   report: `data/valid/report.txt`

# Configuration File — Specification & Guide
- **[docs/validation/CONFIG_GUIDE.md](../../docs/validation/CONFIG_GUIDE.md)**

## Key option notes
- `validation_level`, `type`, `ngs_type`, and `logging_level` are **case-insensitive** — `"Trust"`, `"TRUST"`, and `"trust"` are all accepted.
