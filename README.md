# Simulate Amplicon and WGS Reads From Mutated or Unmodified References

## Quick Start

```bash
nextflow run path/to/main.nf \
    --haplotype_manifest path/to/haplotype_manifest.csv \
    --sample_design_file path/to/sample_manifest.json \
    --run_name "RUN_NAME" \
    -profile docker \
    -resume
```
---

## Overview

The overall workflow involves the following steps:
1. Parse a haplotype definition manifest to create modified references that will be used for simulation (referred to as simulation references in this README)
2. Parse the sample design file, which defines sample IDs and the combinations of haplotypes to be used to create that sample
3. Using the modified references created in (1) and sample information from (2), simulate per-sample-per-haplotype fastq files; and merge files on sample IDs to create a single pair of fastq files per sample.
4. Run fastqc on the simulated data, and collate these outputs using multiqc
5. Additionally create a `{RUN_NAME}.manifest.csv` file, each row of which contains the following information: sample_id, path to simulated R1 fastq file, path to simulated R2 fastq file


## Inputs

### Haplotype Manifest
This is a CSV file containing the definitions of the haplotypes to be simulated. An example is shown below.

```csv
haplotype   , base_fasta , variants_file  , bed_file

haplotype_1 , ref_1.fa   , variants_1.vcf , regions_1.bed
haplotype_2 , ref_2.fa   , variants_2.vcf , regions_2.bed
haplotype_3 , ref_3.fa   , variants_3.vcf , regions_3.bed
haplotype_4 , ref_4.fa   , variants_4.vcf ,
haplotype_5 , ref_5.fa   ,                , regions_5.bed
```

> NB: The values in the columns `base_fasta`, `variants_file`, and `bed_file` are expected to be paths to the respective files - the most robust option here is to use absolute paths.

#### Columns
---
- `haplotype` `[REQUIRED, UNIQUE VALUES EXPECTED]`: An alphanumeric string used as the label for this particular simulation reference
- `base_fasta` `[REQUIRED, NON-UNIQUE VALUES ALLOWED]`: Path to the reference fasta which will be the basis of the simulation reference
- `variants_file` `[OPTIONAL, NON-UNIQUE VALUES ALLOWED]`: Path to a (VCF-like) TSV file containing variant information (1-based coordinates) [OPTIONAL]
- `bed_file` `[OPTIONAL, NON-UNIQUE VALUES ALLOWED]`: Path to Path to a (BED-like) TSV file containing amplicon start/end coordinates (0-based coordinates)

All column names are expected and required in the header. However, leaving values in `variants_file` and `bed_file` is allowed. If neither is specified, the simulation reference is the same as the base fasta file. If only a variants file is specified, the simulation reference will be a "mutated" whole-genome. If only a BED file is specified, the simulation reference will be a collection of amplicons with no "mutations".

### Sample Design File
This is a JSON file containing details of the samples that need to be simulated. An example is shown below.

```json
[
    {
        "sample_id": "sample_1", // sample_id
        "genotypes": ["haplotype_1", "haplotype_2"], // haplotypes needed in this sample
        "proportions": [0.7, 0.3], // proportions of each haplotype specified above
        "num_reads": 100000 // total number of reads to simulate in this sample
    },
    {
        "sample_id": "sample_2",
        "genotypes": ["haplotype_2", "haplotype_3"],
        "proportions": [0.6, 0.4],
        "num_reads": 1000000
    }
]
```

#### Fields
---
- `sample_id` `[REQUIRED, UNIQUE VALUES EXPECTED]`: What to name the sample output files
- `genotypes` `[REQUIRED, NON-UNIQUE VALUES ALLOWED]`: The list of haplotypes to include in this sample
- `proportions` `[REQUIRED, NON-UNIQUE VALUES ALLOWED]`: A corresponding list of fractions defining the proportions of each haplotype specified in `genotypes`
- `num_reads` `[REQUIRED, NON-UNIQUE VALUES ALLOWED]`: The total number of reads to simulate for this sample. The number of reads in the output fastq files will be approximately this value.

> NB: The haplotype names listed in the `genotypes` field **MUST** match those used in the `haplotype` column in the haplotype manifest.  
>
> The number of items in the list `genotypes` **MUST** equal the number of items in the corresponding list `proportions`.

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
- [ ] Add support for defining amplicons using primer-sets instead of coordinates
- [ ] Add support to enable the simulation of optional host-read contamination

## Development
- [ ] Review and clean up code
- [ ] Add param-check function code
- [ ] Add tests