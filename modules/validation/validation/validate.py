#!/usr/bin/env python3
import sys 

from validation.utils import file_exists, TODO
from Bio import SeqIO


#   TODO - return some statistics - like number of sequences?
@file_exists
def is_fasta(filepath: str) -> bool:
    try:
        with open(filepath) as handle:
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