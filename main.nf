#!/usr/bin/env nextflow
include {run_varianteer} from "./modules/run_varianteer.nf"
include {run_art} from "./modules/run_simulation.nf"
include {merge_files; write_output_manifest; run_fastqc; run_multiqc} from "./modules/prepare_outputs.nf"
nextflow.enable.dsl = 2

workflow {

    check_params()

    // [haplotype_label, base_fasta, [variants_file, bed_file], is_amplicon, amplicon_count]
    input_ch = parse_manifest(params.haplotype_manifest)

    // GENERATE SIMREFS
    run_varianteer(input_ch)

    // [halotype_mutated.fa, haplotype_label, is_amplicon, num_amplicons]
    haplotypes_ch = run_varianteer.out.haplotype_ch

    // [sample_id, haplotype_label, fraction_reads]
    sample_design_ch = parse_sample_design(params.sample_design_file)

    // [sample_id, haplotype_label, haplo_fasta, amplicon_fraction_reads, is_amplicon]
    // OR
    // [sample_id, haplotype_label, haplo_fasta, fraction_reads, is_amplicon]
    simulation_in_ch = sample_design_ch
                            .combine(haplotypes_ch, by:1)
                            .map{ haplotype_label, sample_id, num_reads, haplo_fasta, is_amplicon, num_amplicons ->
                                    if (is_amplicon) {
                                        def amplicon_fraction_reads = Math.round(num_reads/num_amplicons)
                                        tuple(sample_id, haplotype_label, haplo_fasta, amplicon_fraction_reads, is_amplicon)
                                    } else {
                                        tuple(sample_id, haplotype_label, haplo_fasta, Math.round(num_reads), is_amplicon)
                                    }
                                }

    run_art(simulation_in_ch) // run_art.out is like [sample_1, sample_1.haplo1_1.fq, sample_1.haplo1_2.fq]

    // collapse run_art.out so it is like [ [sample_1, [sample_1.haplo1_1.fq, sample_1.haplo2_1.fq], [sample_1.haplo1_2.fq, sample_1.haplo2_2.fq]] ]
    grouped_ch = run_art.out.map { sample_id, fq1, fq2 ->
                                tuple([sample_id, '1', fq1], [sample_id, '2', fq2])
                            }
                        .flatMap()
                        .groupTuple(by: [0, 1])
                        .groupTuple(by: [0])
                        .map {sample_id, __, file_list ->
                                tuple(sample_id, file_list[0], file_list[1])
                            }

    merge_files(grouped_ch)

    merge_files.out
        .collect(flat:false, sort:true)
        .set{ manifest_lines_ch }

    merge_files.out
            .map { it ->
                    tuple(it[1], it[2])
                }
            .collect()
            .set{fastqc_ch}


    write_output_manifest(params.run_name, manifest_lines_ch)
    run_fastqc(fastqc_ch)

    run_multiqc(params.run_name, run_fastqc.out.fastqc_zips)

}

def check_params() {
    return 0
}

def parse_manifest(manifest_path) {
    def rows = Channel
                        .fromPath(manifest_path)
                        .splitCsv(header: true, sep: ',')

    def haplotype_manifest_ch = rows.map { row ->
                                            def haplotype_label = row.haplotype.trim().toLowerCase().toString()
                                            def base_fasta = file(row.base_fasta)
                                            def variants_file_val = row.variants_file
                                            def bed_file_val = row.bed_file

                                            def variants_file = null
                                            def bed_file = null
                                            def amplicon_count = null
                                            def is_amplicon = false

                                            if (variants_file_val == ""){
                                                    variants_file = null
                                                } else {
                                                    variants_file = file(row.variants_file)
                                                }

                                            if (bed_file_val == ""){
                                                    bed_file = null
                                                } else {
                                                    is_amplicon = true
                                                    bed_file = file(row.bed_file)
                                                    amplicon_count = bed_file.readLines().size() - 1
                                                }

                                            tuple(haplotype_label, base_fasta, [variants_file, bed_file], is_amplicon, amplicon_count)
                                    }


    return haplotype_manifest_ch
}

def parse_sample_design(sample_design_file) {
    def json_slurper = new groovy.json.JsonSlurper()
    def json_data = json_slurper.parse(new File(sample_design_file))
    def flattened_design_data = []

    json_data.each { entry ->
                        def sample_id = entry.sample_id
                        def haplotypes = entry.genotypes
                        def proportions = entry.proportions
                        def num_reads = entry.num_reads

                    haplotypes.eachWithIndex {
                                                haplotype, i ->
                                                    def hap_proportion = proportions[i]
                                                    def fraction_reads = num_reads * hap_proportion

                                                flattened_design_data << [sample_id, haplotype.trim().toLowerCase().toString(), fraction_reads]
                                            }
                    }

    return Channel.from(flattened_design_data) // [[sample_id, haplotype1, haplotype_fraction_total_reads]]]
}
