#!/usr/bin/env python3
import sys

# Dummy validation script, prints TEST MODE
if __name__ == '__main__':
    # Arguments: <input_file> <output_file>
    print("TEST MODE")
    with open(sys.argv[2], 'w') as out:
        out.write("TEST MODE\n")