#!/usr/bin/env python3
import sys 

from validation.utils import file_exists, TODO,error
from Bio import SeqIO

#   TODO - make it as Validator class?


#   TODO - return some statistics - like number of sequences?
@file_exists
def is_fasta(filepath: str, **kwargs) -> bool:
    try:
        with open(filepath,'r') as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            return len(records) > 0
    except Exception as e:
        return False

@file_exists
def is_gbk(filepath: str) -> bool:
    TODO()

#   TODO - return some statistics - like number of reads, quality?
@file_exists
def is_fastq(filepath: str) -> bool:
    try:
        with open(filepath) as handle:
            records = list(SeqIO.parse(handle, "fastq"))
            return len(records) > 0
    except Exception:
        return False

@file_exists
def is_gff(filepath: str) -> bool:
    TODO()

@file_exists
def is_gtf(filepath: str) -> bool:
    TODO()

@file_exists
def is_bam(filepath: str) -> bool:
    TODO()

@file_exists
def edit_fasta(filepath: str, **kwargs: dict):
    """_summary_
    Check number of sequences and keep only the first one. All the others are saved into their own files with suffix given by the order of the sequences in the original file.
    The "chr1" prefix is added to the first sequence name
    Args:
        filepath (str): filepath to the fasta file
    """
    records = None
    try:
        with open(filepath,'r') as handle:
            records = list(SeqIO.parse(handle, "fasta"))

        if len(records) > 1:
            print("WARNING: more than one sequence found, first will be taken as reference, the rest will be suppressed")
            for i,record in enumerate(records):
                if i == 0:
                    #   first seq = ref
                    with open(filepath,'w') as f:
                        f.write(">chr1 " + record.name+'\n')
                        f.write(str(record.seq)+'\n')
                else:
                    p = filepath[:-6]+f"_{i}.fasta"
                    with open(p,'w') as f:
                        f.write(">"+record.name+'\n')
                        f.write(str(record.seq)+'\n')
    except Exception as e:
        error("07",f"Edit Fasta error: {e}")