# Input Validation Overview

> **Important!**
>
> When container is built please follow the steps to preprocess the data with a validation package.

The input validation module preprocesses and verifies all input data to ensure it meets the required format and structure before the Nextflow pipeline is executed.

## Purpose

The validation module ensures that:

- All input files are in the correct format
- Files are properly structured and can be parsed
- Data meets quality standards
- Files are converted to standardized formats for pipeline processing

## How It Works

The validation process:

1. Reads the configuration file from `data/inputs/config.json`
2. Validates each input file according to its type (genome, reads, features)
3. Converts files to standardized formats
4. Outputs validated files to `data/valid/`
5. Generates logs and reports in `logs/`

## Running Validation

```bash
python3 ./modules/validation/main.py ./data/inputs/config.json
```

## Related Documentation

- [Supported File Formats](FILE_FORMATS.md)
- [Configuration Guide](CONFIG_GUIDE.md)
- [Validation Settings](SETTINGS.md)
- [Performance Tips](PERFORMANCE.md)

## Output

After successful validation:

- Validated files are placed in `data/valid/`
- Log file created in `./logs/validation_ID.log`
- Report file created in `./logs/report_ID.txt`
