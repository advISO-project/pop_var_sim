
process run_art {
    /*
    * Simulates Illumina reads from a given genome using ART.
    * Supports both amplicon and whole-genome sequencing simulations.
    * Outputs paired-end FASTQ files for each sample and haplotype.
    */
    tag { sample_id }

    input:
        tuple val(sample_id), val(haplotype_label), path(genome), val(read_count), val(is_amplicon)
        val read_length
        val seq_system
        val paired_end
        val std_dev
        val mean_frag_len
        val random_seed

    output:
        tuple val(sample_id), path("*.fq")

    script:
        def args = task.ext.args ?: ''
        def paired_flag = paired_end ? "-p" : ""
        if (is_amplicon)
            """
            art_illumina \
                -c $read_count \
                -i $genome \
                -o "${sample_id}.${haplotype_label}_" \
                -d "${haplotype_label}" \
                -na \
                -amp \
                -l $read_length \
                -ss $seq_system \
                $paired_flag \
                -s $std_dev \
                -m $mean_frag_len \
                -rs $random_seed \
                $args
            """
        else
            """
            art_illumina \
                -c $read_count \
                -i $genome \
                -o "${sample_id}.${haplotype_label}_" \
                -d "${haplotype_label}" \
                -na \
                -l $read_length \
                -ss $seq_system \
                $paired_flag \
                -s $std_dev \
                -m $mean_frag_len \
                -rs $random_seed \
                $args
            """

    stub:
        """
        touch $sample_id.${haplotype_label}_1.fq
        touch $sample_id.${haplotype_label}_2.fq

        echo '$sample_id || $haplotype_label || $read_count reads' > $sample_id.${haplotype_label}_1.fq
        echo '$sample_id || $haplotype_label || $read_count reads' > $sample_id.${haplotype_label}_2.fq
        """
}