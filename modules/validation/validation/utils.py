#!/usr/bin/env python3
import os
import sys
import shutil
import json

#   //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
#   ERROR HANDLING

ERROR_CODES = {
    "00": "Invalid number of arguments",
    "01": "Missing arguments",
    "02": "",
    "03": "Unknown format",
    "04": "",
    "05": "File management failed",
    "06": "File validation failed",
    "07": "",
    "08": "",
    "09": "",
    "10": "Unknown Exception",

}

def error(code, msg=None):
    print(f"ERROR_{code}: {msg or ERROR_CODES.get(code, 'Unknown error')}", file=sys.stderr)
    sys.exit(1)

def success():
    print("SUCCESS")
    sys.exit(0)

#   //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
#   SYSTEM FILES MANAGEMENT

def check_cfg(cfg: dict):
    """
    Validate that all mandatory arguments exist in the configuration.

    Args:
        cfg (dict): Configuration dictionary parsed from JSON.

    Raises:
        SystemExit: If any mandatory argument is missing.
    """
    mandatory_arguments = [
        "ref_genome_filename",
        "mod_genome_filename",
        "mod_reads_filename",
        "ref_genome_extension", # data format specification - temporary?
        "mod_genome_extension", # data format specification - temporary?
        "mod_reads_extension"   # data format specification - temporary?
    ]
    
    missing = [arg for arg in mandatory_arguments if arg not in cfg]
    
    if missing:
        error("01",f"Missing mandatory arguments: {', '.join(missing)}")

def read_cfg(cfg_path: str) -> dict:
    """
    Load a JSON configuration file.

    Args:
        cfg_path (str): Path to the JSON config file.

    Returns:
        dict: Dictionary of input arguments
    
    Raises:
        FileNotFoundError: If the file does not exist.
        json.JSONDecodeError: If the file is not valid JSON.
    """
    try:
        cfg = json.load(open(cfg_path))
        return cfg
    except (FileNotFoundError, json.JSONDecodeError) as e:
        error("05",f"Reading configuration file failed: {e}")    

def create_tmp():
    try:
        cwd = os.getcwd()
        tmp_dir = os.path.join(cwd,"tmp")
        os.makedirs(tmp_dir,exist_ok=True)
        return tmp_dir
    except (PermissionError, OSError) as e:
        error("05",f"Creating tmp directory failed: {e}")

def save_file(what, on):
    try:
        # Copy file to destination
        shutil.copy2(what, on)
    except Exception as e:
        error("05", f"Failed to copy file from {what} to {on}: {e}")

#   //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
#   OTHER

def file_exists(func):
    def inner(filepath, *args, **kwargs):
        if not os.path.isfile(filepath):
            error("05",f"Unable to open file on path {filepath}")
        abs_path = os.path.abspath(filepath)
        return func(abs_path, *args, **kwargs)
    return inner

def TODO():
    print("WARNING: This part is not implemented yet")
