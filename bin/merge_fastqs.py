#!/usr/bin/env python3

import sys

def interleave_sample_pair(r1_path, r2_path, out_handle):
    """Write a single sample's reads to the output interleaved file"""
    with open(r1_path, 'r') as r1, open(r2_path, 'r') as r2:
        while True:
            r1_lines = [r1.readline() for _ in range(4)]
            r2_lines = [r2.readline() for _ in range(4)]

            # Stop at end of file
            if not r1_lines[0] or not r2_lines[0]:
                break

            # Write interleaved reads
            out_handle.writelines(r1_lines)
            out_handle.writelines(r2_lines)

def interleave_multiple_samples(pairs, out_path):
    """Interleave all sample pairs into a single output file"""
    with open(out_path, 'w') as out:
        for r1_path, r2_path in pairs:
            interleave_sample_pair(r1_path, r2_path, out)

if __name__ == "__main__":
    if len(sys.argv) < 4 or (len(sys.argv[1:]) - 1) % 2 != 0:
        print(f"Usage: {sys.argv[0]} <output.interleaved.fastq> <R1_1> <R2_1> [<R1_2> <R2_2> ...]")
        sys.exit(1)

    out_path = sys.argv[1]
    input_files = sys.argv[2:]
    pairs = list(zip(input_files[::2], input_files[1::2]))

    interleave_multiple_samples(pairs, out_path)
    print(f"Interleaved FASTQ written to {out_path}")

