process COUNT_NORM {
    tag "$meta.id"
    label 'process_medium'

    // (Bio)conda packages have intentionally not been pinned to a specific version
    // This was to avoid the pipeline failing due to package conflicts whilst creating the environment when using -profile conda
    conda (params.enable_conda ? "conda-forge::r-base bioconda::bioconductor-deseq2 bioconda::bioconductor-biocparallel bioconda::bioconductor-tximport bioconda::bioconductor-complexheatmap conda-forge::r-optparse conda-forge::r-ggplot2 conda-forge::r-rcolorbrewer conda-forge::r-pheatmap" : null)
    container 'docker://fgualdr/envnorm_fix'

    input:
    tuple val(meta), path(counts)
    path deseq2_pca_header
    path deseq2_clustering_header

    output:

    path "normalization"        , optional:true, emit: noamlization
    path "scaling_dat.txt" , optional:true, emit: noamlization_txt
    path "*.pdf"                , optional:true, emit: pdf
    path "*.RData"              , optional:true, emit: rdata
    path "*.rds"                , optional:true, emit: rds
    path "*pca.vals.txt"        , optional:true, emit: pca_txt
    path "*pca.vals_mqc.tsv"    , optional:true, emit: pca_multiqc
    path "*sample.dists.txt"    , optional:true, emit: dists_txt
    path "*sample.dists_mqc.tsv", optional:true, emit: dists_multiqc
    path "*.log"                , optional:true, emit: log
    path "size_factors"         , optional:true, emit: size_factors
    path "versions.yml"         , emit: versions

    script:
    def args      = task.ext.args ?: ''
    def peak_type = params.narrow_peak ? 'narrowPeak' : 'broadPeak'
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def norm = params.normalize ? 'TRUE' : 'FALSE'
    def sigma_times = params.sigma_times ? params.sigma_times : '1'
    def n_pop = params.n_pop ? params.n_pop : 'NULL'

    """ 
    general_normalizer.r \\
        --count_file $counts \\
        --outdir ./ \\
        --norm $norm \\
        --sigma_times $sigma_times \\
        --n_pop $n_pop \\
        --outprefix $prefix \\
        --cores $task.cpus \\
        $args

    sed 's/deseq2_pca/deseq2_pca_${task.index}/g' <$deseq2_pca_header >tmp.txt
    sed -i -e 's/DESeq2 /${meta.id} DESeq2 /g' tmp.txt
    cat tmp.txt ${prefix}.pca.vals.txt > ${prefix}.pca.vals_mqc.tsv

    sed 's/deseq2_clustering/deseq2_clustering_${task.index}/g' <$deseq2_clustering_header >tmp.txt
    sed -i -e 's/DESeq2 /${meta.id} DESeq2 /g' tmp.txt
    cat tmp.txt ${prefix}.sample.dists.txt > ${prefix}.sample.dists_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}
