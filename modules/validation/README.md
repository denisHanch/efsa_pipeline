#   How to run

Load your data and configuration file into `/EFSA_workspace/data/inputs/` (description of config.json bellow)

Run validation as 
```
python3 /EFSA_workspace/modules/validation/main.py /EFSA_workspace/data/inputs/config.json > stdout.log
```

output of the validation is on `/EFSA_workspace/data/valid/` 

#   What to export
logs are on `.logs/` and `./stdout.log`
Please share them with me.

# Configuration File — Specification & Guide
- **[validation-pkg/docs/CONFIG_GUIDE.md](validation-pkg/docs/CONFIG_GUIDE.md)**
