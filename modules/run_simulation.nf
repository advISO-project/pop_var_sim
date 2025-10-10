nextflow.enable.dsl = 2
process run_art {
    /*
    * Simulates Illumina reads from a given genome using ART.
    * Supports both amplicon and whole-genome sequencing simulations.
    * Outputs paired-end FASTQ files for each sample and haplotype.
    */
    tag { sample_id }

    input:
        tuple val(sample_id), val(haplotype_label), path(genome), val(read_count), val(is_amplicon)

    output:
        tuple val(sample_id), path("$sample_id.${haplotype_label}_1.fq"), path("$sample_id.${haplotype_label}_2.fq")

    script:
        def args = task.ext.args ?: ''

        if (is_amplicon)
            """
            art_illumina \
                -c $read_count \
                -i $genome \
                -o "${sample_id}.${haplotype_label}_" \
                -amp \
                $args
            """
        else
            """
            art_illumina \
                -c $read_count \
                -i $genome \
                -o ${sample_id}.${haplotype_label}_ \
                $args
            """

    stub:
        """
        touch $sample_id.${haplotype_label}_1.fq
        touch $sample_id.${haplotype_label}_2.fq

        echo '$sample_id || $haplotype_label || $read_count reads' > $sample_id.${haplotype_label}_1.fq
        echo '$sample_id || $haplotype_label || $read_count reads' > $sample_id.${haplotype_label}_2.fq

        #gzip $sample_id.${haplotype_label}_1.fq
        #gzip $sample_id.${haplotype_label}_2.fq

        """
}