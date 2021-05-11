process BWAMEM {
    /**
    * Maps trimmed paired fastq using BWA (http://bio-bwa.sourceforge.net/)
    * Uses samtools to convert to BAM, sort and index sorted BAM (http://www.htslib.org/doc/samtools.html)
    * @input
    * @output
    */

    tag {"Mapping_${seq}"}
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*.sorted.bam", mode: 'copy'

    cpus 4

    input:
        tuple(path(seq),path(ref))
        path(indexed_reference)

    output:
        path("*.sorted.bam"), emit: bam

    script:
      """
      bwa mem -t ${task.cpus} ${ref} ${seq} | \
      samtools sort -o ${seq.baseName}.sorted.bam
      """
}
