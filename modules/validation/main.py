#!/usr/bin/env python3
import validation.validate as v
import validation.format as r
from validation.utils import *

#   TODO  
#       zatim pridat gtf bez validace, ale aby prosel
#       temporary file extension in config file? testovat pro fasta, fa, fna   
#       test for duplicate filepath? (ref & mod same file?)
#       jak rozlisit jestli je long nebo short read ?????

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
        tmp_dir = create_tmp(base_path)
        print(f"    -   DONE: created on {tmp_dir}")

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        #   Run validation on every expected file
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        #   REFERENCE GENOME
        #       File with reference genome (in FASTA or geneBank format), keyed as ref_genome_filepath
        print("INFO: Validating Reference genome.")
        ref_genome_filepath = os.path.join(base_path,cfg["ref_genome_filename"])
        ref_genome_ext = cfg["ref_genome_extension"].lower()    #   expected one of fasta or gbk

        #   rename, save reference fasta file
        save_file(ref_genome_filepath, f"{tmp_dir}/ref.{ref_genome_ext}")
        ref_genome_filepath = f"{tmp_dir}/ref.{ref_genome_ext}"

        #   decode gzip
        if ref_genome_ext.endswith(".gz"):
            r.gz_decode(ref_genome_filepath,ref_genome_filepath[:-3])
            ref_genome_filepath = ref_genome_filepath[:-3]
            ref_genome_ext = ref_genome_ext[:-3] 

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

        #   rename, save and edit normalized reference fasta file
        v.edit_fasta(ref_genome_filepath)

        print(f"    -   DONE: Saved on {tmp_dir}/ref.{ref_genome_ext}")
        
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        
        #   MODIFIED GENOME
        #       File with modificated genome (in FASTA or geneBank format), keyed as mod_genome_filepath
        print("INFO: Validating Modified genome.")
        mod_genome_filepath = os.path.join(base_path,cfg["mod_genome_filename"])
        mod_genome_ext = cfg["mod_genome_extension"].lower()    #   expected one of fasta or gbk

        #   rename, save modified fasta file
        save_file(mod_genome_filepath, f"{tmp_dir}/mod.{mod_genome_ext}")
        mod_genome_filepath = f"{tmp_dir}/mod.{mod_genome_ext}"

        #   decode gzip
        if mod_genome_ext.endswith(".gz"):
            r.gz_decode(mod_genome_filepath,mod_genome_filepath[:-3])
            mod_genome_filepath = mod_genome_filepath[:-3]
            mod_genome_ext = mod_genome_ext[:-3] 

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

        #   edit normalized reference fasta file
        v.edit_fasta(mod_genome_filepath)

        print(f"    -   DONE: Saved on {tmp_dir}/mod.{mod_genome_ext}")
        
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        #   READS
        #       File with reference genome (in FASTA or geneBank format), keyed as mod_reads_filepath
        print("INFO: Validating Reads of modified genome.")

        for reads_file in cfg["reads"]:
            
            reads_filepath = os.path.join(base_path,reads_file["reads_filename"])
            reads_ext = reads_file["reads_extension"].lower()  #   expected one of fastq or bam

            #   rename, save reference fasta file
            save_file(reads_filepath, f'{tmp_dir}/{reads_file["reads_filename"]}')
            reads_filepath = f'{tmp_dir}/{reads_file["reads_filename"]}'

            #   decode gzip
            if reads_ext.endswith(".gz"):
                r.gz_decode(reads_filepath,reads_filepath[:-3])
                reads_filepath = reads_filepath[:-3]
                reads_ext = reads_ext[:-3] 

            expected_fastq = (reads_ext == "fastq") 
            expected_bam = (reads_ext == "bam") 

            if not (expected_fastq ^ expected_bam):
                error_message = (
                        f"Invalid reads format "
                        f"Expected one of: fastq, bam "
                    )
                error("03",error_message)

            #   validate fasta format
            if expected_fastq:
                if not v.is_fastq(reads_filepath):
                    error("06","Modified reads is not in FASTQ format as expected")

            #   validate gbk format
            if expected_bam:
                if not v.is_bam(reads_filepath):
                    error("06","Modified reads is not in BAM format as expected")
                else:
                    #   reformat to FASTA
                    #   replace filepath
                    # reads_filepath = r.gbk_to_fasta(reads_filepath)
                    TODO()


            #   encode gzip
            r.gz_encode(reads_filepath,reads_filepath+".gz")
            reads_filepath = reads_filepath+".gz"
            reads_ext = reads_ext + ".gz"

            print(f'    -   DONE: Saved on {reads_filepath}')

    except Exception as e:
        print("INFO: Cleaning tmp files")
        os.removedirs(tmp_dir)
        error("10",f"Unhandled main exception {e}")



