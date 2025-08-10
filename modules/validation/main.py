#!/usr/bin/env python3
import validation.validate as v
import validation.format as r
from validation.utils import *

#   TODO  
#       Validation of fasta raise Warning about deprecated comments - ??
#       temporary file extension in config file?       
#       test for duplicate filepath? (ref & mod same file?)
#       

#   //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//

if __name__ == '__main__':
    try:
        # Arguments: <config.json>
        if len(sys.argv) != 2:
            error_message = (
                f"Invalid number of arguments. "
                f"Usage: python {sys.argv[0]} <config_file.json>"
            )
            error("00",error_message)

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        cfg_path = sys.argv[1]
        print(f"INFO: Reading config_file on {cfg_path}")

        cfg = read_cfg(cfg_path)
        print(f"    -   DONE: Reading config_file was sucessful")

        base_path = os.path.dirname(cfg_path)
        print(base_path)


        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        #   Check expected number of arguments in config file 
        print("INFO: Checking mandatory arguments.")
        check_cfg(cfg)
        print(f"    -   DONE: All mandatory arguments are present")

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        #   Create tmp directory if not present
        print("INFO: Checking tmp directory.")
        tmp_dir = create_tmp()
        print(f"    -   DONE: created on {tmp_dir}")

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   Run validation on every expected file
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        #   REFERENCE GENOME
        #       File with reference genome (in FASTA or geneBank format), keyed as ref_genome_filepath
        print("INFO: Validating Reference genome.")
        ref_genome_filepath = os.path.join(base_path,cfg["ref_genome_filename"])
        ref_genome_ext = cfg["ref_genome_extension"].lower()    #   expected one of fasta or gbk

        expected_fasta = (ref_genome_ext.lower() == "fasta") 
        expected_gbk = (ref_genome_ext.lower() == "gbk") 

        if not (expected_fasta ^ expected_gbk):
            error_message = (
                    f"Invalid reference genome format "
                    f"Expected one of: fasta, gbk "
                )
            error("03",error_message)

        #   validate fasta format
        if expected_fasta:
            if not v.is_fasta(ref_genome_filepath):
                error("06","Reference genome is not in FASTA format as expected")

        #   validate gbk format
        if expected_gbk:
            if not v.is_gbk(ref_genome_filepath):
                error("06","Reference genome is not in GBK format as expected")
            else:
                #   reformat to FASTA
                #   replace filepath
                # ref_genome_filepath = r.gbk_to_fasta(ref_genome_filepath)
                TODO()

        #   rename and save normalized reference fasta file
        save_file(ref_genome_filepath, tmp_dir+"/ref_genome.fasta")
        print(f"    -   DONE: Saved on {tmp_dir}ref_genome.fasta")
        
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        
        #   MODIFIED GENOME
        #       File with modificated genome (in FASTA or geneBank format), keyed as mod_genome_filepath
        print("INFO: Validating Modified genome.")
        mod_genome_filepath = os.path.join(base_path,cfg["mod_genome_filename"])
        mod_genome_ext = cfg["mod_genome_extension"].lower()    #   expected one of fasta or gbk

        expected_fasta = (mod_genome_ext == "fasta") 
        expected_gbk = (mod_genome_ext == "gbk") 

        if not (expected_fasta ^ expected_gbk):
            error_message = (
                    f"Invalid modified genome format "
                    f"Expected one of: fasta, gbk "
                )
            error("03",error_message)

        #   validate fasta format
        if expected_fasta:
            if not v.is_fasta(mod_genome_filepath):
                error("06","Modified genome is not in FASTA format as expected")

        #   validate gbk format
        if expected_gbk:
            if not v.is_gbk(mod_genome_filepath):
                error("06","Modified genome is not in GBK format as expected")
            else:
                #   reformat to FASTA
                #   replace filepath
                # mod_genome_filepath = r.gbk_to_fasta(mod_genome_filepath)
                TODO()

        #   rename and save normalized reference fasta file
        save_file(mod_genome_filepath, tmp_dir+"/mod_genome.fasta")
        print(f"    -   DONE: Saved on {tmp_dir}mod_genome.fasta")
        
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        #   READS
        #       File with reference genome (in FASTA or geneBank format), keyed as mod_reads_filepath
        print("INFO: Validating Reads of modified genome.")

        mod_reads_filepath = os.path.join(base_path,cfg["mod_reads_filename"])
        mod_reads_ext = cfg["mod_reads_extension"].lower()  #   expected one of fastq or bam

        expected_fastq = (mod_reads_ext == "fastq") 
        expected_bam = (mod_reads_ext == "bam") 

        if not (expected_fastq ^ expected_bam):
            error_message = (
                    f"Invalid reads format "
                    f"Expected one of: fastq, bam "
                )
            error("03",error_message)

        #   validate fasta format
        if expected_fastq:
            if not v.is_fastq(mod_reads_filepath):
                error("06","Modified reads is not in FASTQ format as expected")

        #   validate gbk format
        if expected_bam:
            if not v.is_bam(mod_reads_filepath):
                error("06","Modified reads is not in BAM format as expected")
            else:
                #   reformat to FASTA
                #   replace filepath
                # mod_reads_filepath = r.gbk_to_fasta(mod_reads_filepath)
                TODO()

        #   rename and save normalized reference fasta file
        save_file(mod_reads_filepath, tmp_dir+"/mod_reads.fastq")
        print(f"    -   DONE: Saved on {tmp_dir}mod_reads.fastq")
        

    except Exception as e:
        print(e)
        error("10","Unhandled main exception")




