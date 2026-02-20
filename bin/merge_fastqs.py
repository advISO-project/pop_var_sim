#!/usr/bin/env python3

#!/usr/bin/env python3
import argparse
from pathlib import Path
from collections import defaultdict

def interleave_pair(f1, f2, out_handle):
    """Interleave reads from two FASTQ files."""
    with open(f1, 'r') as fh1, open(f2, 'r') as fh2:
        while True:
            # Read one FASTQ record (4 lines) from each file
            r1 = [fh1.readline() for _ in range(4)]
            r2 = [fh2.readline() for _ in range(4)]

            if not r1[0] or not r2[0]:
                # End of either file
                break

            # Write interleaved
            out_handle.writelines(r1)
            out_handle.writelines(r2)

def write_single(fq, out_handle):
    """Write single FASTQ file as-is."""
    with open(fq, 'r') as fh:
        for line in fh:
            out_handle.write(line)

def main():
    parser = argparse.ArgumentParser(
        description="Interleave paired FASTQs and add unpaired files"
    )
    parser.add_argument("-o", "--output", required=True, help="Output file name")
    parser.add_argument("inputs", nargs='+', help="Input FASTQ files")

    args = parser.parse_args()

    # Group files by prefix
    pairs = defaultdict(dict)  # prefix -> {'1': file, '2': file}
    singles = []

    for path in args.inputs:
        fname = Path(path).name
        if fname.endswith('_1.fq'):
            prefix = fname.rsplit('_1.fq', 1)[0]
            pairs[prefix]['1'] = path
        elif fname.endswith('_2.fq'):
            prefix = fname.rsplit('_2.fq', 1)[0]
            pairs[prefix]['2'] = path
        else:
            singles.append(path)

    # Files that couldn't form pairs are considered singles
    for prefix, files in pairs.items():
        if '1' not in files or '2' not in files:
            for f in files.values():
                singles.append(f)
            pairs[prefix] = None  # mark as invalid pair

    # Open output file
    with open(args.output, 'w') as out_handle:
        # Interleave valid pairs
        for prefix, files in pairs.items():
            if files is not None:
                interleave_pair(files['1'], files['2'], out_handle)

        # Add unpaired files
        for fq in singles:
            write_single(fq, out_handle)

if __name__ == "__main__":
    main()

