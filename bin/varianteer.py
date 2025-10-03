from data_loader import load_file
from multiprocessing import Pool

class Varianteer():
    def __init__(self, sample_id:str, fasta_file:str, vcf_file:str, bed_file:str):
        self.sample_id = sample_id
        self.genome = load_file(file_path=fasta_file, file_format="fasta")

        self.variant_coord, self.variants = {}, {}
        if vcf_file:
            self.variant_coord, self.variants = load_file(file_path=vcf_file, file_format="vcf")

        if bed_file:
            self.regions = load_file(file_path=bed_file, file_format="bed")
            self.mode = "amplicon"
        else:
            self.regions = None
            self.mode = "wgs"

        self.mutated_data = {}

    def process_amplicons(self):
        for record in self.genome:
            rec_id = record["id"]
            if not rec_id in self.regions.keys():
                continue
            seq = record["seq"]
            relevant_coords = [(int(c["START"]), int(c["END"]), c["amplicon_id"]) for c in self.regions[rec_id]]
            relevant_variants = self.variants.get(rec_id, {})

            self.amplicons = self._collect_amplicons_using_bed(seq, relevant_coords, relevant_variants)
            for amp_id, amp_data in self.amplicons.items():
                mutated = self.apply_mutations(amp_id, amp_data["seq"], amp_data["variants"])
                self.mutated_data.update(mutated)

    def _collect_amplicons_using_bed(self, seq: str, coords: list, variants: dict):
        amplicons = {}
        for (amp_start, amp_end, amp_id) in coords:
            variant_subset = [i for i in [variants.get(i, ()) for i in range(amp_start, amp_end+1)] if i]
            if not variant_subset:
                amplicons[amp_id] = {"start": amp_start, "end": amp_end, "seq": seq[amp_start:amp_end+1], "variants": []}
            else:
                variant_subset_corrected = self._offset_variant_coords(variant_subset, amp_start)
                amplicons[amp_id] = {"start": amp_start, "end": amp_end, "seq": seq[amp_start:amp_end+1], "variants": variant_subset_corrected}

        return amplicons

    def _offset_variant_coords(self, variants:list, amp_start: int):
        # print(f"Original {variants = }")
        corrected_variants_list = []
        for (pos, ref, alt) in variants:
            corrected_variants_list.append((pos-amp_start, ref, alt))
        # print(f"{corrected_variants_list = }")
        return corrected_variants_list

    def process_whole_genome_naive(self):
        for record in self.genome:
            rec_id = record["id"]
            seq = record["seq"]
            relevant_variants = [data for pos, data in self.variants[rec_id].items()]
            self.mutated_data.update(self.apply_mutations(rec_id, seq, relevant_variants))

    def process_whole_genome_pool(self):
        items = []
        for record in self.genome:
            rec_id = record["id"]
            seq = record["seq"]
            relevant_variants = [data for pos, data in self.variants[rec_id].items()]
            items.append((rec_id, seq, relevant_variants))

        with Pool(4) as pool:
            results = pool.starmap(self.apply_mutations, items)
        for r in results:
            self.mutated_data.update(r)

    def process_whole_genome_chunks(self, chunk_size):
        for record in self.genome:
            rec_id = record["id"]
            seq = record["seq"]
            rec_vars = self.variants.get(rec_id, {})
            ## if no variants on this record, do nothing, update mutated_data, skip to next
            if not rec_vars:
                self.mutated_data[rec_id] = {"original": seq, "mutated": seq}
                continue

            if len(seq) > chunk_size:
                chunk_mutated_data = {}
                seq_chunk_idxs = [i for i in range(0, len(seq), chunk_size)]

                seq_chunk_idxs.append(len(seq))
                seq_chunk_coords = [(seq_chunk_idxs[i], seq_chunk_idxs[i+1]) for i in range(len(seq_chunk_idxs)-1)]

                counter = 0
                for (seq_chunk_start, seq_chunk_end) in seq_chunk_coords:
                    chunk_id = f"{rec_id}.{counter}"
                    chunk_variant_all = [rec_vars.get(i, ()) for i in range(seq_chunk_start, seq_chunk_end)]
                    chunk_variant_subset = [i for i in chunk_variant_all if i]
                    corrected_chunk_variant_subset = self._offset_variant_coords(chunk_variant_subset, seq_chunk_start)
                    seq_chunk = seq[seq_chunk_start:seq_chunk_end]
                    chunk_mutated_data.update(self.apply_mutations(chunk_id, seq_chunk, corrected_chunk_variant_subset))
                    counter += 1

                consolidated = {"original": "", "mutated": ""}
                for k, v in chunk_mutated_data.items():
                    consolidated["original"] += v["original"]
                    consolidated["mutated"] += v["mutated"]

                self.mutated_data[rec_id] = consolidated

            else:
                relevant_variants = [data for pos, data in self.variants[rec_id].items()]
                self.mutated_data.update(self.apply_mutations(rec_id, seq, relevant_variants))

    def apply_mutations(self, seq_id: str, seq: str, variants: list):
        mutated_data = {}
        if not variants:
            print(f"No Variants for {seq_id}!!!!")
            mutated_data[seq_id] = {"original": seq, "mutated": seq}
            return mutated_data
        b_seq = bytearray(seq.encode("utf-8"))
        for pos, ref, alt in variants:
            try:
                if seq[pos] != ref:
                    print(f"WARNING: in given sequence at position {pos}, expected {ref} but found {seq[pos]}")

                b_seq[pos] = alt.encode("utf-8")[0]

                mutated_data[seq_id] = {"original": seq, "mutated": b_seq.decode("utf-8")}
            except IndexError as ie:
                raise ie

        return mutated_data

    def write_outputs(self, outfile_name: str, mode: str):
        match mode.lower():
            case "wgs":
                with open(f"{outfile_name}_mutated.fa", "w") as outfile_handle:
                    for k, v in self.mutated_data.items():
                        outfile_handle.write(f">{k}\n")
                        outfile_handle.write(f"{v['mutated']}\n")

                outfile_handle.close()
                print(f"Wrote fasta file with mutations to: {outfile_name}_mutated.fa")

            case "amplicon":
                orig_amp_file = open(f"{outfile_name}_amp.fa", "w")
                mut_amp_file = open(f"{outfile_name}_amp_mutated.fa", "w")
                for k, v in self.mutated_data.items():
                    orig_amp_file.write(f">{k}\n")
                    orig_amp_file.write(f"{v['original']}\n")

                    mut_amp_file.write(f">{k}\n")
                    mut_amp_file.write(f"{v['mutated']}\n")

                print(f"Wrote fasta file with non-mutated amplicons to: {outfile_name}_amp.fa")
                print(f"Wrote fasta file with mutated amplicons to: {outfile_name}_amp_mutated.fa")
