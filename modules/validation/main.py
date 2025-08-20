#!/usr/bin/env python3
import validation.validate as v
from validation.utils import *
import json
import sys


#   TODO  
#       jak rozlisit jestli je long nebo short read ????? 
#       print funkce trochu automatizovat, aby se neopakovalo self.__class__.__name__, __name__ 
#       cleanup doesnt work

def check_cfg(cfg: dict):
    """
    Validate that all mandatory arguments exist in the configuration.

    Args:
        cfg (dict): Configuration dictionary parsed from JSON.

    Raises:
        SystemExit: If any mandatory argument is missing.
    """
    try:
        print(f"{__name__}:check_cfg: +INFO: Checking mandatory arguments.")
        mandatory_arguments = [
            "ref_genome_filename",
            "mod_genome_filename",
            "reads" #   TODO: check also read filename and extensions?
        ]
        
        missing = [arg for arg in mandatory_arguments if arg not in cfg]
        
        if missing:
            raise Exception(f"Missing mandatory arguments: {', '.join(missing)}")
        
        if len(cfg["reads"]) == 0:
            raise Exception(f"Missing Reads")
        
        if cfg["ref_genome_filename"] == cfg["mod_genome_filename"]:
            raise Exception("Same ref and mod genome file")
        print(f"{__name__}:check_cfg: |-DONE\n")
    except Exception as e:
        print(f"    {__name__}:check_cfg: |--ERROR: {e}\n")
        raise


def validate_genome(filepath, ext, suff):
    try:
        print(f"{__name__}:validate_genome: +INFO: Validating {suff} genome on {filepath}")

        if check_file_format(ext,FASTA_FORMATS):
            mod_validator = v.FASTA_Validator(filepath, "FASTA" ,outsuffix=suff)        
        elif check_file_format(ext,GBK_FORMATS):
            mod_validator = v.GBK_Validator(filepath, "FASTA" ,outsuffix=suff)        
        else:
            raise Exception(f"Unknown file format {ext}, expected fasta in {FASTA_FORMATS} or gene bank in {GBK_FORMATS}")
        status = mod_validator.run()
        print(f"{__name__}:validate_genome: |-DONE: {str(status)}\n")
    except Exception as e:
        print(f"    {__name__}:validate_genome: |--ERROR: {e}\n")
        raise


def validate_reads(filepath, ext, suff):
    try:
        print(f"{__name__}:validate_reads +INFO: Validating Reads on {filepath}.")
        if check_file_format(ext,FASTQ_FORMATS):
            mod_validator = v.FASTQ_Validator(filepath, "FASTQ" ,outsuffix=suff)        
        elif check_file_format(ext,BAM_FORMATS):
            mod_validator = v.BAM_Validator(filepath, "FASTQ" ,outsuffix=suff)        
        else:
            raise Exception(f"Unknown file format {ext}, expected fastq in {FASTQ_FORMATS} or bam in {BAM_FORMATS}")
        status = mod_validator.run()
        print(f"{__name__} |-DONE: {str(status)}\n")
    except Exception as e:
        print(f"    {__name__}:validate_reads: |--ERROR: {e}\n")
        raise


def validate_feature(filepath, ext, suff):
    try:
        print(f"{__name__}:validate_feature: +INFO: Validating Features on {filepath}.")
        if check_file_format(ext,GTF_FORMATS):
            mod_validator = v.GTF_Validator(filepath, "GTF" ,outsuffix=suff)        
        elif check_file_format(ext,GFF_FORMATS):
            mod_validator = v.GFF_Validator(filepath, "GTF" ,outsuffix=suff)        
        else:
            raise Exception(f"Unknown file format {ext}, expected gtf in {GTF_FORMATS} or gff in {GFF_FORMATS}")
        status = mod_validator.run()
        print(f"{__name__}:validate_feature: |-DONE: {str(status)}\n")
    except Exception as e:
        print(f"    {__name__}:validate_feature: |--ERROR: {e}\n")
        raise

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   MAIN
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

if __name__ == '__main__':

    try:
        # Arguments: <config.json>
        if len(sys.argv) != 2:
            raise BaseException(f"Invalid number of arguments. Usage: python {sys.argv[0]} <config_file.json>")

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   Config file mng

        cfg_path = sys.argv[1]

        #   Read config file
        print(f"{__name__} +INFO: Reading config_file on {cfg_path}")
        cfg = json.load(open(cfg_path))
        print(f"{__name__} |-DONE\n")

        base_path = os.path.dirname(cfg_path)

        #   Check expected number of arguments in config file 
        check_cfg(cfg)


        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   Run validation on every expected file
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   REFERENCE GENOME
        #       File with reference genome (in FASTA or geneBank format), keyed as "ref_genome_..."
        filepath = os.path.join(base_path,cfg["ref_genome_filename"])
        filext = '.'.join(cfg["ref_genome_filename"].split('.')[1:]).upper()
        validate_genome(filepath,filext,"ref")

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   MODIFIED GENOME
        #       File with modified genome (in FASTA or geneBank format), keyed as "mod_genome_..."
        filepath = os.path.join(base_path,cfg["mod_genome_filename"])
        filext = '.'.join(cfg["mod_genome_filename"].split('.')[1:]).upper()
        validate_genome(filepath,filext,"mod")

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   READS
        #       Files with reads (in FASTQ or BAM format), keyed as "reads_..."
        for filename in cfg["reads"]:
            filepath = os.path.join(base_path,filename)
            filext = '.'.join(filename.split('.')[1:]).upper()
            validate_reads(filepath,filext,"mod")

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   FEATURES
        #       GTF or GFF files with additional information about reference/modified genomes
        try:
            validate_feature(os.path.join(base_path,cfg["ref_feature_filename"], cfg["ref_feature_extension"].upper(),"ref"))
        except KeyError: #  TODO?? proc
            print(f"{__name__} |-WARNING: no feature file for reference genome \n")
        
        try:
            validate_feature(os.path.join(base_path,cfg["mod_feature_filename"], cfg["mod_feature_extension"].upper(),"mod"))
        except KeyError:
            print(f"{__name__} |-WARNING: no feature file for modified genome \n")
        
    except Exception as e:
        print(f"{__name__} |-ERROR: {e}\n")
        



