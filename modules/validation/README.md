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

**!! WARNING !! - Every new run of validation rewrite the previous one**

#   What to export
-   log: `/EFSA_workspace/data/valid/validation.log`
-   report: `/EFSA_workspace/data/valid/report.txt`

# Documentation
- **[read here](validation-pkg/docs/INDEX.md)**