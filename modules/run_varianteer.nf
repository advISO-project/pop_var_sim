nextflow.enable.dsl = 2
process run_varianteer {
    /*
    * Runs Varianteer to simulate variants in a given haplotype.
    */

    tag { haplotype_label }
    input:
        tuple val(haplotype_label), path(fasta_file), val(amp_and_var_files), val(is_amplicon), val(num_amplicons)

    output:
        tuple path("${haplotype_label}*_mutated.fa"), val(haplotype_label), val(is_amplicon), val(num_amplicons), emit: haplotype_ch


    script:
        def vcf_file = null
        def bed_file = null
        def args = []
        if (amp_and_var_files[0] != null) {
            vcf_file = amp_and_var_files[0]
            args << "-vf $vcf_file"
        }
        if (amp_and_var_files[1] != null) {
            bed_file = amp_and_var_files[1]
            args << "-b $bed_file"
        }
        def arg_str = args.join(" ")

        // If no VCF or BED file is provided, just copy the original FASTA to the output with the new name
        if (vcf_file == null && bed_file == null)
            """
            cp $fasta_file "${haplotype_label}_mutated.fa"
            """
        // Otherwise, run Varianteer with the provided arguments
        else
            """
            run_varianteer.py \
                -s $haplotype_label \
                -f $fasta_file \
                $arg_str
            """

    stub:
        """
        touch ${haplotype_label}_mutated.fa
        echo ${haplotype_label} > ${haplotype_label}_mutated.fa
        """
}