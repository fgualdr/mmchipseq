/*
 * Filter BAM file
 */
process BAM_EM_PREP {
    tag "$meta.id"
    label 'process_samtools_sort'

    conda (params.enable_conda ? "bedtools=2.31.0 samtools=1.17 coreutils=9.3 gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-9d3a458f6420e5712103ae2af82c94d26d63f059:60b54b43045e8cf39ba307fd683c69d4c57240ce-0':
        'quay.io/biocontainers//mulled-v2-9d3a458f6420e5712103ae2af82c94d26d63f059:60b54b43045e8cf39ba307fd683c69d4c57240ce-0' }"

    input:
    tuple val(meta), path(bam)
    path chrom_sizes

    output:
    tuple val(meta), path("*.read_target_match.txt"), emit: read_target_match
    tuple val(meta), path("*.labeled.nsort.bam"), emit: bam_label
    tuple val(meta), path("*.multi.hotspots.bed"), emit: mm_hotspots
    path "versions.yml"           , emit: versions

    script:
    def prefix           = task.ext.prefix ?: "${meta.id}"

    // added to control memory usage on HPC

    def overlap = params.label_overlap ? params.label_overlap : 0.25
    // convert to megabytes the memory available
    def avail_mem = task.memory ? task.memory.mega : false
    // sort memory is the available memory minus 1GB divided by the number of cpus
    def sort_memory = avail_mem ? ((avail_mem - (1000*task.cpus)) /task.cpus).intValue() : 2000
    def dist_r = params.inser_size ? params.inser_size.toInteger() : 1000

    """
    # Add the HI:i bit to each
    samtools view -h $bam | awk 'BEGIN{OFS="\\t"} \$1~/^@/{print;next} {for(i=12;i<=NF;i++){if(\$i~/^HI:i:/){ \$1=\$1"_"\$i}};print}' | samtools view -bS - > ${prefix}.labeled.bam
    samtools sort -@ $task.cpus -m ${sort_memory}M -n -o ${prefix}.labeled.nsort.bam -T ${prefix}.labeled.nsort ${prefix}.labeled.bam
    bedtools bamtobed -i ${prefix}.labeled.nsort.bam -bedpe > ${prefix}.labeled.nsort.bedpe 

    #################################################################################
    #   GENERATE THE APPROPRIATE FILES FOR EM                                       #
    #################################################################################

    # Generate the files required for the EM step:
    # Separate multi from single - this is done on the BAM first then on the BEDPE
    samtools view -h ${prefix}.labeled.nsort.bam | grep -v 'NH:i:1'  | samtools view -bS - > ${prefix}.labeled.multi.bam
    bedtools bamtobed -i ${prefix}.labeled.multi.bam -bedpe > ${prefix}.labeled.nsort.multi.bedpe  

    # Hotspots
    cut -f 1,2,6 ${prefix}.labeled.nsort.multi.bedpe  > ${prefix}.multi.bed
    sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.multi.bed > ${prefix}.multi.sort.bed
    bedtools merge -d $dist_r -i ${prefix}.multi.sort.bed > ${prefix}.multi.sort.merged.bed
    awk 'BEGIN {FS="\\t"; OFS="\t"} {print \$0"\\t""id_"NR}' ${prefix}.multi.sort.merged.bed > ${prefix}.multi.hotspots.temp.bed
    sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.multi.hotspots.temp.bed > ${prefix}.multi.hotspots.bed
    
    # find intersections - at least 50% of the read has to be within default!
    cut -f 1,2,6,7 ${prefix}.labeled.nsort.bedpe  > ${prefix}.labeled.nsort.bed
    sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.labeled.nsort.bed > ${prefix}.labeled.sort.bed
    bedtools intersect \\
            -a ${prefix}.labeled.sort.bed \\
            -b ${prefix}.multi.hotspots.bed \\
            -sorted \\
            -wo -loj > ${prefix}.nsort.match.bed
    
    # Replace the no match "." with "id_nomatch"
    awk 'BEGIN{FS=OFS="\\t"} \$8=="."{\$8="id_nomatch"}1' ${prefix}.nsort.match.bed > ${prefix}.nsort.match.fix.bed

    #################################
    #   OUTPUT DATA                 #
    #################################

    # Generate the read_target table for the EM which is composed of read_full_id, read_id_ori, target_id
    awk -F'\\t' 'BEGIN{OFS="\\t"} {id=\$4;sub("_HI:.*", "", \$4); print id, \$4, \$8}' ${prefix}.nsort.match.fix.bed > ${prefix}.read_target.txt
    # Separate the ${prefix}.read_target.txt in two with the target_id  == 'id_nomatch' or != 'id_nomatch' which is found of the 3rd column 
    awk -F'\\t' 'BEGIN{OFS="\\t"} \$3=="id_nomatch"{print \$1, \$2, \$3}' ${prefix}.read_target.txt > ${prefix}.read_target_nomatch.txt
    awk -F'\\t' 'BEGIN{OFS="\\t"} \$3!="id_nomatch"{print \$1, \$2, \$3}' ${prefix}.read_target.txt > ${prefix}.read_target_match.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
