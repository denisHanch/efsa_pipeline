#!/usr/bin/env python3
from validation.utils import file_exists, TODO, error
import gzip
import shutil
import os

@file_exists
def gbk_to_fasta(filepath):
    TODO()


@file_exists
def bam_to_fastq(filepath):
    TODO()

@file_exists
def gz_decode(filepath: str, outpath: str, replace:bool = True) -> str:
    """Decompress a .gz file and save it without the .gz extension.

    Args:
        filepath (str): Path to the .gz file.

    Returns:
        str: Path to the decompressed file.
    """
    try:
        with gzip.open(filepath, 'rb') as f_in:
            with open(outpath, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        if replace:
            os.remove(filepath)

    except Exception as e:
        error("05",f"Unable to decode gzip file {filepath}")

@file_exists
def gz_encode(filepath: str, outpath: str, compresslevel: int = 6, replace: bool = True) -> str:
    """
    Compress a file with gzip.

    Args:
        filepath: Path to the source file to compress.
        compresslevel: gzip compression level (0-9). 6 is default/good balance.

    Returns:
        The path to the created .gz file.
    """
    try:
        # Prepare atomic write
        tmppath = outpath + ".tmp"

        # Compress streaming
        with open(filepath, "rb") as f_in:
            with gzip.open(tmppath, "wb", compresslevel=compresslevel) as f_out:
                shutil.copyfileobj(f_in, f_out)

        # Atomic move into place
        os.replace(tmppath, outpath)

        if replace:
            os.remove(filepath)
    except Exception as e:
        error("05",f"Unable to decode gzip file {filepath}")

