/*
 * Filter BAM file
 */
process BAM_FRAGMENT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::deeptools=3.5.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0' :
        'quay.io/biocontainers/deeptools:3.5.1--py_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
 
    output:
    tuple val(meta),  path("*.fragment_median.txt"), emit: fragment_median

    path "*.fragmentSize.png" , emit: fragment_png
    path "*.fragmentSize.txt" , emit: fragment_txt
    path "versions.yml"           , emit: versions

    script:
    def prefix           = task.ext.prefix ?: "${meta.id}"
    """   
    # get insert size stat and automatically select those with at most three times the median
    bamPEFragmentSize \\
            --bamfiles $bam \\
            --histogram ${prefix}.fragmentSize.png \\
            --numberOfProcessors ${task.cpus} \\
            --plotTitle "Fragment size of PE data" \\
            --maxFragmentLength 0 \\
            --table ${prefix}.fragmentSize.txt    

    # get the correct colum
    col=\$(awk -F'\\t' -v var='Frag. Len. Median' '{ for(i=1;i<=NF;i++) { if(\$i == var) { printf(i) } } exit 0 }' ${prefix}.fragmentSize.txt)

    # We consider three times the median
    filename=\$(basename $bam)
    identifier=\${filename%%.*}
    frag_median=\$(awk -v var="^\$identifier" \\
                        -v col="\$col" \\
                        '\$0~var{print \$col}' \\
                        ${prefix}.fragmentSize.txt)

    echo \$frag_median > ${prefix}.fragment_median.txt 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(echo \$(deeptools --version 2>&1) | sed 's/^.*deeptools //; s/Part .*\$//')
    END_VERSIONS
    """
    
}
