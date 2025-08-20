#!/usr/bin/env python3
import os
import shutil
import subprocess


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   GLOBALS
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


#   the firts in tuple will be standard for naming all of such files
GZIP_EXTS = ("GZ")
BZIP_EXTS = ("BZ2")
ENCODING_EXTS = GZIP_EXTS+BZIP_EXTS
FASTA_FORMATS = ("FASTA","FA","FNA")
GBK_FORMATS = ("GBK","GENE_BANK","GENE_BANK")
FASTQ_FORMATS = ("FASTQ")
BAM_FORMATS = ("BAM")
GTF_FORMATS = ("GTF")
GFF_FORMATS = ("GFF")

def left_join_with_coding(what):
    return what + tuple(f"{f}.{e}" for f in what for e in ENCODING_EXTS)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   SYSTEM FILE MANAGEMENT
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def split_filepath(filepath:str) -> list:
    filepath = os.path.abspath(filepath)
    basename = os.path.basename(filepath)
    basedir = os.path.dirname(filepath)
    tmp = basename.split('.')
    if len(tmp) > 3 or len(tmp) < 2:
        raise Exception(f"Unexpected file extension, allowed one for file format and one for coding type, got {len(tmp)}")
    if len(tmp) == 3:
        return basedir,tmp[0],tmp[1].upper(),tmp[2].upper()    
    if len(tmp) == 2:
        return basedir,tmp[0],tmp[1].upper(), None    


def create_tmpdir(basepath = os.getcwd()):
    tmpdir = os.path.join(basepath,"tmp")
    os.makedirs(tmpdir,exist_ok=True)
    return tmpdir


def remove_dir(dirpath):
    """
    Remove a directory and all its contents.
    If the directory does not exist, do nothing.
    If removal fails, exit with error.
    """
    shutil.rmtree(dirpath)


def file_exists(filepath:str):
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"file {filepath} not found")


def move_file(src, dst):
    """
    Move a file from src to dst using subprocess and handle errors.
    """
    subprocess.run(['mv', src, dst], check=True)


def copy_file(src, dst):
    """
    Copy a file from src to dst using subprocess and handle errors.
    """
    subprocess.run(['cp', src, dst], check=True)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#   OTHER
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def check_coding_format(coding_format):
    if coding_format not in ENCODING_EXTS:
        raise Exception(f"Unknown coding format")


def check_file_format(format, allowed):
    if format not in allowed:
        return False
    return True


def TODO(pref = None):
    print(f"{pref}WARNING: This part is not implemented yet")