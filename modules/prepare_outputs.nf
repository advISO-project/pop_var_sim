nextflow.enable.dsl = 2
process merge_files {
    input:
        tuple val(sample_id), path(fastq1), path(fastq2)

    output:
        tuple val(sample_id), path("${sample_id}_1.fq.gz"), path("${sample_id}_2.fq.gz")

    script:
        """
        cat $fastq1 > ${sample_id}_1.fq
        cat $fastq2 > ${sample_id}_2.fq

        gzip ${sample_id}_1.fq
        gzip ${sample_id}_2.fq
        """

    stub:
        """
        cat $fastq1 > ${sample_id}_1.fq.gz
        cat $fastq2 > ${sample_id}_2.fq.gz
        """
}

process write_output_manifest {
    input:
        val(run_name)
        val(list_of_lines)

    output:
        path("${run_name}.manifest.csv")

    script:
        def header = "sample_id,fastq_1,fastq2"
        def csvFile = file("out.csv")
            csvFile.withWriter { writer ->
                writer.writeLine(header)
                list_of_lines
                                .sort{ a, b -> a[0] <=> b[0] }
                                .each { row ->
                                        def row_sample_id = row[0]
                                        def f1 = (row[1]).name
                                        def f2 = (row[2]).name
                                        def path1 = file("${params.outdir}/fastqs/$f1")
                                        def path2 = file("${params.outdir}/fastqs/$f2")
                                        def line = "$row_sample_id,$path1,$path2"
                                        writer.writeLine(line)
                                        }
                                    }
        """
        mv $csvFile ${run_name}.manifest.csv
        """

    stub:
        """
        touch ${run_name}.manifest.csv
        """
}

process run_fastqc {
    input:
        val(fq_list)

    output:
        path "*_fastqc.zip", emit: fastqc_zips
        path "*_fastqc.html", emit: fastqc_htmls

    script:
        """
        fastqc -o . ${fq_list.join(' ')}
        """
}

process run_multiqc {
    input:
        val(run_name)
        path(fastqc_zips), stageAs: "?/*"

    output:
        path "${run_name}.multiqc_report.html", emit: multiqc_html

    script:
        """
        multiqc .
        mv multiqc_report.html ${run_name}.multiqc_report.html
        """
}