#!/usr/bin/env python3
import os
import subprocess
import validation.utils as u


def remove_coding(coding_type,**kwargs):
    if coding_type in u.GZIP_EXTS:
        gz_decode(**kwargs)
    elif coding_type in u.BZIP_EXTS:
        bz_decode(**kwargs)
    

def add_coding(coding_type,**kwargs):
    if coding_type in u.GZIP_EXTS:
        gz_encode(**kwargs)
    elif coding_type in u.BZIP_EXTS:
        bz_encode(**kwargs)


def gz_decode(filepath: str, outpath: str, replace: bool = False):
    """
    Decompress a .gz file using subprocess and save it to outpath.
    """
    try:
        with open(outpath, 'wb') as f_out:
            subprocess.run(['gzip', '-dc', filepath], stdout=f_out, check=True)
        if replace:
            os.remove(filepath)
    except Exception as e:
        raise RuntimeError(f"Failed to decompress {filepath}: {e}")


def gz_encode(filepath: str, outpath: str, compresslevel: int = 6, replace: bool = False):
    """
    Compress a file with gzip using subprocess.
    """
    try:
        with open(filepath, 'rb') as f_in, open(outpath, 'wb') as f_out:
            subprocess.run(['gzip', f'-{compresslevel}'], stdin=f_in, stdout=f_out, check=True)
        if replace:
            os.remove(filepath)
    except Exception as e:
        raise RuntimeError(f"Failed to compress {filepath}: {e}")


def bz_decode(filepath: str, outpath: str, replace: bool = False):
    """
    Decompress a .bz2 file using subprocess and save it to outpath.
    """
    try:
        with open(outpath, 'wb') as f_out:
            subprocess.run(['bzip2', '-dc', filepath], stdout=f_out, check=True)
        if replace:
            os.remove(filepath)
    except Exception as e:
        raise RuntimeError(f"Failed to decompress {filepath}: {e}")


def bz_encode(filepath: str, outpath: str, compresslevel: int = 9, replace: bool = False):
    """
    Compress a file with bzip2 using subprocess.
    """
    try:
        with open(filepath, 'rb') as f_in, open(outpath, 'wb') as f_out:
            subprocess.run(['bzip2', f'-{compresslevel}'], stdin=f_in, stdout=f_out, check=True)
        if replace:
            os.remove(filepath)
    except Exception as e:
        raise RuntimeError(f"Failed to compress {filepath}: {e}")