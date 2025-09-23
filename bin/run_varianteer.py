#!/usr/bin/env python3
import os, argparse
from varianteer import Varianteer

def get_args():
    parser = argparse.ArgumentParser(prog="Varianteer: produce reference sequences with custom mutations.")
    parser.add_argument("-s", "--sample_id", type = str, required = True, help = "Sample ID to label outputs with.")
    parser.add_argument("-f", "--fasta_file", type = str, required = True, help = "The base reference sequence in fasta format.")
    parser.add_argument("-vf", "--vcf_file", type = str, required = False, help = "The set of variants to apply.")
    parser.add_argument("-b", "--bed_file", type = str, required = False, help = "The bed file containing amplicon start and end coordinates.")
    parser.add_argument("-cs", "--chunk_size", type = int, required = False, default=5000, help = "Size of chunks into which to break down large input sequences. [Default 5000]")
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    sample_id = args.sample_id
    fasta = os.path.abspath(args.fasta_file)
    if args.vcf_file:
        vcf = os.path.abspath(args.vcf_file)
    else:
        vcf = None
    if args.bed_file:
        bed = os.path.abspath(args.bed_file)
    else:
        bed = None
    chunk_size = args.chunk_size

    if not fasta or (file_not_found := not(os.path.exists(fasta))):
        if file_not_found:
            raise FileNotFoundError(f"Specified fasta file not found! Expected fasta file at path: {fasta}")
        raise ValueError("Fasta file is required!")

    if not vcf and not bed:
        print("Nothing to do!")
        return

    run_obj = Varianteer(sample_id=sample_id, fasta_file=fasta, vcf_file=vcf, bed_file=bed)
    if bed:
        run_obj.process_amplicons()
        run_obj.write_outputs(outfile_name=sample_id, mode="amplicon")
    else:
        run_obj.process_whole_genome_chunks(chunk_size)
        run_obj.write_outputs(outfile_name=sample_id, mode="wgs")

if __name__ == "__main__":
    main()


