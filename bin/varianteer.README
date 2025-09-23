# Varianteer: Create Custom Refernce Sequences

## Quick Start

```bash
run_varianteer.py \
    -s label \
    -f path/to/fasta_file \
    -vf path/to/vcf_file \
    -b $bed_file
```
---
### Parameters
- `-s`, `--sample_id` `[REQUIRED]`: Sample ID to label the outputs with.
- `-f`, `--fasta_file` `[REQUIRED]`: Path to the fasta file that will form the basis of the modified reference sequence.
- `-vf`, `--vcf_file` `[OPTIONAL]`: Path to the VCF-like TSV file detailing the variants to be applied to the base sequence.
- `-b`, `--bed_file` `[OPTIONAL]`: Path to the BED TSV file defining the start and end coordinates of amplicons to be extracted from the base sequence.

---
### Outputs
If a BED file is provided, then the modified reference sequence file will contain amplicon sequences. Each amplicon sequence in the output file is named as: `>{chromosome_name}_{start_coord}_{end_coord}`. Two output files are produced:
- `{sample_id}_amp.fa`: Amplicon sequences, with no mutations applied
- `{sample_id}_amp_mutated.fa`: Amplicon sequences, with mutations applied (when mutations are provided via the variants file; otherwise this file has the same contents as the first.)

If a BED file is **not** provided, then the modified sequence file will contian the same sequences as the base reference, but with mutations applied. In this case, one output file is created:
- `{sample_id}_mutated.fa`: The input fasta file, with mutations applied

> NB: In the latter case, if a variants file is **not** provided, then nothing is done as the output would be identical to the input fasta file.

---
### Example Variants File
This is a VCF-like tab-separated file. An example is shown below (all columns shown below are **required**).
> COORDINATES IN THIS FILE ARE EXPECTED TO BE 1-BASED
```tsv
#CHROM	  POS	REF	ALT
chrom_1	  22	T	A
chrom_1	  138	C	A
chrom_1	  232	G	T
chrom_1	  302	A	G
chrom_1	  366	T	C
```

### Example BED File
This is a tab-separated BED file. An example is shown below (columns `#CHROM`, `START`, and `END` are **required**. All other columns are optional and will not be used).
> COORDINATES IN THIS FILE ARE EXPECTED TO BE 0-BASED
```tsv
#CHROM	  START	    END	     DESC
chrom_1	  145388	145662	 chrom_1-145388-145662-amplicon
chrom_1	  162867	163115	 chrom_1-162867-163115-amplicon
chrom_1	  181512	181761	 chrom_1-181512-181761-amplicon
chrom_1	  455794	456054	 chrom_1-455794-456054-amplicon
chrom_1	  528859	529104	 chrom_1-528859-529104-amplicon
chrom_1	  535965	536239	 chrom_1-535965-536239-amplicon
```
---

# TODO

## Features
- [ ] Add support for non-SNP variants

## Development
- [ ] Add docstrings to functions
- [ ] Add unit tests
- [ ] Pull code out into own repo