import os, csv

def load_file(file_path:str, file_format:str, parse_mode:str = None):
    full_path = _check_file_exists(file_path)
    match file_format:
        case "fa" | "fasta":
            if not parse_mode:
                file_small_enough = False
                parse_mode = "iterator" if not (file_small_enough:= _file_size_small(full_path)) else "list"
                if file_small_enough:
                    print("Input fasta file size is <= 100MB, using parse_mode 'list'...")
            match parse_mode:
                case "iter" | "iterator":
                    return _parse_fasta(full_path)
                case "list":
                    return list(_parse_fasta(full_path))
                case _:
                    raise ValueError(f"For fasta file, invalid parse_mode provided: {parse_mode}. Must be 'iter' OR 'iterator'/'list'. ")

        case "bed":
            return _load_bed_file(full_path)

        case "vcf":
            return _load_vcf_file(full_path)
        case _:
            raise ValueError(f"For fasta file, invalid file_format provided: {file_format}. Must be 'fa' OR 'fasta'/'bed'/'vcf'. ")

def _parse_fasta(fasta_file):
    header = None
    seq = []
    with open(fasta_file, "r") as fastafile_handle:
        for line in fastafile_handle:
            line = line.strip()
            ## empty line
            if not line:
                continue
            ## header line
            if line.startswith(">"):
                ## if header pre-exists here, this line is a new header, yield record
                if header:
                    rec_id = header.split(" ")[0]
                    desc = " ".join(header.split(" ")[1:])
                    yield {"id": rec_id, "desc": desc, "seq": "".join(seq)}
                header = line[1:]
                seq = []
            ## seq line
            else:
                seq.append(line.upper())
        ## last record
        if header:
            rec_id = header.split(" ")[0]
            desc = " ".join(header.split(" ")[1:])
            yield {"id": rec_id, "desc": desc, "seq": "".join(seq)}

def _load_bed_file(bed_file):
    amplicon_coords = {}
    count = 0
    with open(bed_file, "r") as bedfile_handle:
        reader = csv.DictReader(bedfile_handle, delimiter="\t")
        for row in reader:
            chrom = row["#CHROM"]
            # BED is half-open, so adjust start 
            start = row["START"]
            end = row["END"]
            row.update({"amplicon_id": f"{chrom}_{start}_{end}"})
            if chrom in amplicon_coords.keys():
                amplicon_coords[chrom].append(row)
                count += 1
            else:
                amplicon_coords[chrom] = [row]
                count += 1

    print(f"Found {count} amplicons in the BED file provided.")
    return amplicon_coords

def _load_vcf_file(vcf_file):
    variant_coords = {}
    variant_data = {}
    with open(vcf_file, "r") as vcffile_handle:
        reader = csv.DictReader(vcffile_handle, delimiter="\t")
        found = 0
        coords_loaded = 0
        vars_loaded = 0
        for row in reader:
            found += 1
            chrom = row["#CHROM"]
            ## Adjust VCF coordinates: these are 1-indexed
            pos = int(row["POS"])
            if chrom in variant_coords.keys():
                variant_coords[chrom].append(pos)
                coords_loaded += 1
            else:
                variant_coords[chrom] = [pos]
                coords_loaded += 1

            if chrom in variant_data.keys():
                variant_data[chrom].update({pos: (pos, row["REF"], row["ALT"])})
                vars_loaded += 1
            else:
                variant_data[chrom] = {pos: (pos, row["REF"], row["ALT"])}
                vars_loaded += 1

    print(f"Found {found} variants in file, loaded {coords_loaded} co-ordinates and {vars_loaded} variants.")
    return variant_coords, variant_data

def _check_file_exists(file_path:str):
    file_abspath = os.path.abspath(file_path)
    if not os.path.exists(file_abspath):
        raise FileNotFoundError(f"{file_abspath} does not exist!")
    return file_abspath

def _file_size_small(file_path:str): #, size_limit:int = 100):
    size_limit = 100
    file_size_mb = os.path.getsize(file_path) / (1024 * 1024)
    return file_size_mb <= size_limit
