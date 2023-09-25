/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    aligners       : [  'star' ]
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowChipseq.initialise(params, log, valid_params)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.multiqc_config,
    params.fasta,
    params.gtf, params.gff, params.gene_bed,
    params.star_index,
    params.blacklist,
    params.bamtools_filter_pe_config, params.bamtools_filter_se_config //,params.rerpmsk
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
//if (params.rerpmsk) { ch_rerpmsk = file(params.rerpmsk) } else { exit 1, 'rerpmsk must be provided!' }

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

// JSON files required by BAMTools for alignment filtering
ch_bamtools_filter_se_config = file(params.bamtools_filter_se_config, checkIfExists: true)
ch_bamtools_filter_pe_config = file(params.bamtools_filter_pe_config, checkIfExists: true)

// Header files for MultiQC
ch_spp_nsc_header           = file("$projectDir/assets/multiqc/spp_nsc_header.txt", checkIfExists: true)
ch_spp_rsc_header           = file("$projectDir/assets/multiqc/spp_rsc_header.txt", checkIfExists: true)
ch_spp_correlation_header   = file("$projectDir/assets/multiqc/spp_correlation_header.txt", checkIfExists: true)
ch_peak_count_header        = file("$projectDir/assets/multiqc/peak_count_header.txt", checkIfExists: true)
ch_frip_score_header        = file("$projectDir/assets/multiqc/frip_score_header.txt", checkIfExists: true)
ch_peak_annotation_header   = file("$projectDir/assets/multiqc/peak_annotation_header.txt", checkIfExists: true)
ch_deseq2_pca_header        = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_deseq2_clustering_header = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)

ch_with_inputs = params.with_inputs.toBoolean()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { FRIP_SCORE                          } from '../modules/local/frip_score'
include { PLOT_MACS2_QC                       } from '../modules/local/plot_macs2_qc'
include { PLOT_HOMER_ANNOTATEPEAKS            } from '../modules/local/plot_homer_annotatepeaks'
include { MACS2_CONSENSUS                     } from '../modules/local/macs2_consensus'
include { ANNOTATE_BOOLEAN_PEAKS              } from '../modules/local/annotate_boolean_peaks'
include { COUNT_NORM                          } from '../modules/local/count_normalization'
include { IGV                                 } from '../modules/local/igv'
include { MULTIQC                             } from '../modules/local/multiqc'
include { MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS } from '../modules/local/multiqc_custom_phantompeakqualtools'
include { MULTIQC_CUSTOM_PEAKS                } from '../modules/local/multiqc_custom_peaks'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK         } from '../subworkflows/local/input_check'
include { PREPARE_GENOME      } from '../subworkflows/local/prepare_genome'
include { BAM_FILTER_EM } from '../subworkflows/local/bam_filter_em'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { PICARD_MERGESAMFILES          } from '../modules/nf-core/modules/picard/mergesamfiles/main'
include { PICARD_COLLECTMULTIPLEMETRICS } from '../modules/nf-core/modules/picard/collectmultiplemetrics/main'
include { PHANTOMPEAKQUALTOOLS          } from '../modules/nf-core/modules/phantompeakqualtools/main'
include { DEEPTOOLS_BIGWIG              } from '../modules/local/deeptools_bw' 
include { DEEPTOOLS_BIGWIG_NORM         } from '../modules/local/deeptools_bw_norm' 
include { DEEPTOOLS_COMPUTEMATRIX       } from '../modules/nf-core/modules/deeptools/computematrix/main'
include { DEEPTOOLS_PLOTPROFILE         } from '../modules/nf-core/modules/deeptools/plotprofile/main'
include { DEEPTOOLS_PLOTHEATMAP         } from '../modules/nf-core/modules/deeptools/plotheatmap/main'
include { DEEPTOOLS_PLOTFINGERPRINT     } from '../modules/nf-core/modules/deeptools/plotfingerprint/main'
include { KHMER_UNIQUEKMERS             } from '../modules/nf-core/modules/khmer/uniquekmers/main'
include { MACS2_CALLPEAK as MACS2_CALLPEAK_SINGLE          } from '../modules/nf-core/modules/macs2/callpeak/main'
include { MACS2_CALLPEAK as MACS2_CALLPEAK_MERGED          } from '../modules/nf-core/modules/macs2/callpeak/main'
include { SUBREAD_FEATURECOUNTS         } from '../modules/nf-core/modules/subread/featurecounts/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_MACS2     } from '../modules/nf-core/modules/homer/annotatepeaks/main'
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_CONSENSUS } from '../modules/nf-core/modules/homer/annotatepeaks/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//

include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main'
include { FASTQ_FASTQC_UMITOOLS_FASTP      } from '../subworkflows/nf-core/fastq_fastqc_umitools_fastp/main'

include { FASTQC_TRIMGALORE      } from '../subworkflows/nf-core/fastqc_trimgalore'
include { ALIGN_STAR             } from '../subworkflows/nf-core/align_star'
include { MARK_DUPLICATES_PICARD } from '../subworkflows/nf-core/mark_duplicates_picard'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CHIPSEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    PREPARE_GENOME (
        params.aligner
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input),
        params.seq_center
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with TrimGalore!
    //
    ch_filtered_reads      = Channel.empty()
    ch_fastqc_raw_multiqc  = Channel.empty()
    ch_fastqc_trim_multiqc = Channel.empty()
    ch_trim_log_multiqc    = Channel.empty()
    ch_trim_read_count     = Channel.empty()
    if (params.trimmer == 'trimgalore') {
        FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
            INPUT_CHECK.out.reads,
            params.skip_fastqc || params.skip_qc,
            false,
            false,
            params.skip_trimming,
            0
        )
        ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
        ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
        ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip
        ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log
        ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)
    }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with fastp
    //
    if (params.trimmer == 'fastp') {
        FASTQ_FASTQC_UMITOOLS_FASTP (
            INPUT_CHECK.out.reads,
            params.skip_fastqc || params.skip_qc,
            false,
            false,
            0,
            params.skip_trimming,
            [],
            params.save_trimmed,
            params.save_trimmed,
            params.min_trimmed_reads
        )
        ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
        ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip
        ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip
        ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json
        ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_read_count
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)
    }

    //
    // SUBWORKFLOW: Alignment with STAR & BAM QC
    //
    if (params.aligner == 'star') {
        ALIGN_STAR (
            ch_filtered_reads,
            PREPARE_GENOME.out.star_index
        )
        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_genome_bam_index  = ALIGN_STAR.out.bai
        
        ch_samtools_stats    = ALIGN_STAR.out.stats
        ch_samtools_flagstat = ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = ALIGN_STAR.out.idxstats
        ch_star_multiqc      = ALIGN_STAR.out.log_final

        ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)
    }

    //
    // MODULE: Merge resequenced BAM files - 
    // Sligtly off.. would be better to identify before hand which ones are to be merged.. this can be done from the sample_sheet
    // This is a point of collection to is stops before proceeding:
    // It removes the "T" and merges all of them later will remove the "R[0-9]" bit
    //
    ch_genome_bam
        .map {
            meta, bam ->
                def meta_clone = meta.clone()
                meta_clone.remove('read_group')
                meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
                [ meta_clone, bam ]  
        }
        .groupTuple(by: [0])
        .map { 
            it ->
                [ it[0], it[1].flatten() ] 
        }
        .set { ch_sort_bam }

    PICARD_MERGESAMFILES (
        ch_sort_bam
    )
    ch_versions = ch_versions.mix(PICARD_MERGESAMFILES.out.versions.first().ifEmpty(null))

    //
    // SUBWORKFLOW: Mark duplicates & filter BAM files after merging
    //
    MARK_DUPLICATES_PICARD (
        PICARD_MERGESAMFILES.out.bam
    )
    ch_versions = ch_versions.mix(MARK_DUPLICATES_PICARD.out.versions)

    //
    // SUBWORKFLOW: Filter BAM file with BamTools 
    //
    
    BAM_FILTER_EM (
        MARK_DUPLICATES_PICARD.out.bam.join(MARK_DUPLICATES_PICARD.out.bai, by: [0]),
        PREPARE_GENOME.out.filtered_bed.first(),
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.chrom_sizes,
        ch_bamtools_filter_se_config,
        ch_bamtools_filter_pe_config
    )
    ch_versions = ch_versions.mix(BAM_FILTER_EM.out.versions.first().ifEmpty(null))

    //
    // MODULE: Picard post alignment QC
    //
    ch_picardcollectmultiplemetrics_multiqc = Channel.empty()
    if (!params.skip_picard_metrics) {
        PICARD_COLLECTMULTIPLEMETRICS (
            BAM_FILTER_EM.out.bam,
            PREPARE_GENOME.out.fasta,
            []
        )
        ch_picardcollectmultiplemetrics_multiqc = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
    }

    //
    // MODULE: Phantompeaktools strand cross-correlation and QC metrics
    //
    PHANTOMPEAKQUALTOOLS (
        BAM_FILTER_EM.out.bam
    )
    ch_versions = ch_versions.mix(PHANTOMPEAKQUALTOOLS.out.versions.first())

    //
    // MODULE: MultiQC custom content for Phantompeaktools
    //
    MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS (
        PHANTOMPEAKQUALTOOLS.out.spp.join(PHANTOMPEAKQUALTOOLS.out.rdata, by: [0]),
        ch_spp_nsc_header,
        ch_spp_rsc_header,
        ch_spp_correlation_header
    )


    //
    // Create channels: [ meta, [ ip_bam, control_bam ] [ ip_bai, control_bai ] ]
    // Differently from standard nf-core chipseq we can evaluate the possibility to run the chip-seq w/o inputs
    // This needs to be assessed on the fly i.e. check if there ar
    
    BAM_FILTER_EM
        .out
        .bam
        .join(BAM_FILTER_EM.out.bai, by: [0])
        .set { ch_genome_bam_bai }

    //
    // MODULE: deepTools plotFingerprint thi will assess in sample only
    //
    ch_deeptoolsplotfingerprint_multiqc = Channel.empty()
    if (!params.skip_plot_fingerprint ) {
        DEEPTOOLS_PLOTFINGERPRINT (
            ch_genome_bam_bai
        ) 
        ch_deeptoolsplotfingerprint_multiqc = DEEPTOOLS_PLOTFINGERPRINT.out.matrix
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTFINGERPRINT.out.versions.first())
    }
    
    if(!ch_with_inputs){

       println "The value of ch_with_inputs set to w-o input: ${ch_with_inputs}"
        // Create channels: [ meta, ip_bam, ([] for control_bam) ]
        ch_genome_bam_bai
            .map {
                meta, bam, bai -> 
                    !meta.is_input ? [ meta , bam, [] ] : null
            }
            .set { ch_ip_control_bam }
        
        // w/o inputs we simply merge all bams by antibody: from meta,bam,bai 
        ch_ip_control_bam
            .map {
                meta, bam1, bam2 ->
                def new_meta = meta.clone()
                new_meta.id =  meta.antibody
                [new_meta, bam1, bam2]
            }
            .groupTuple(by: 0)
            .map {
                meta, bam1, bam2 ->
                    [ meta , bam1, [] ]
            }
            .set { ch_antibody_bam }
            
            ch_antibody_bam.view()
        
    }else{ 
        println "The value of ch_with_inputs set to with the input: ${ch_with_inputs}"

        ch_genome_bam_bai
        .combine(ch_genome_bam_bai)
        .map { 
            meta1, bam1, bai1, meta2, bam2, bai2 ->
                !meta1.is_input && meta1.which_input == meta2.id ? [ meta1, [ bam1 ], [ bam2 ] ] : null
        }
        .set { ch_ip_control_bam } 

        // w inputs we simply merge all bams by antibody: from meta,bam,bai - we combine the samples and the inputs:
        // we start from the paired and we group it:
        ch_ip_control_bam
            .map {
                meta, bam1, bam2 ->
                def new_meta = meta.clone()
                new_meta.id =  meta.antibody
                [new_meta, bam1, bam2]
            }
            .groupTuple(by: 0)
            .map {
                meta, bam1, bam2 ->
                    [ meta , bam1, bam2 ]
            }
            .set { ch_antibody_bam }

    }
    // 
    // MODULE: Calculute genome size with khmer
    //
    ch_macs_gsize                     = Channel.empty()
    ch_custompeaks_frip_multiqc       = Channel.empty()
    ch_custompeaks_count_multiqc      = Channel.empty()
    ch_plothomerannotatepeaks_multiqc = Channel.empty()
    ch_subreadfeaturecounts_multiqc   = Channel.empty()
    ch_macs_gsize = params.macs_gsize
     
    if (!params.macs_gsize) {
        KHMER_UNIQUEKMERS (
            PREPARE_GENOME.out.fasta,
            params.read_length
        )
        ch_macs_gsize = KHMER_UNIQUEKMERS.out.kmers.map { it.text.trim() }
    }

    //
    // MODULE: Call peaks with MACS2
    //
    MACS2_CALLPEAK_SINGLE (
         ch_ip_control_bam,
         ch_macs_gsize
    )
    ch_versions = ch_versions.mix(MACS2_CALLPEAK_SINGLE.out.versions.first())

    //
    // Filter out samples with 0 MACS2 peaks called
    //
    MACS2_CALLPEAK_SINGLE
        .out
        .peak
        .filter { meta, peaks -> peaks.size() > 0 }
        .set { ch_macs2_peaks }

    // If is narrow we call high conf summits by merging all BAMS:

    MACS2_CALLPEAK_MERGED(
        ch_antibody_bam,
        ch_macs_gsize
    )

    // Create channels: [ meta, ip_bam, peaks ]
    ch_ip_control_bam
        .join(ch_macs2_peaks, by: [0])
        .map { 
            it -> 
                [ it[0], it[1], it[3] ] 
        }
        .set { ch_ip_bam_peaks }


    //
    // MODULE: Calculate FRiP score
    //
    FRIP_SCORE (
        ch_ip_bam_peaks
    )
    ch_versions = ch_versions.mix(FRIP_SCORE.out.versions.first())

    // Create channels: [ meta, peaks, frip ]
    ch_ip_bam_peaks
        .join(FRIP_SCORE.out.txt, by: [0])
        .map { 
            it -> 
                [ it[0], it[2], it[3] ] 
        }
        .set { ch_ip_peaks_frip }

    //
    // MODULE: FRiP score custom content for MultiQC
    //
    MULTIQC_CUSTOM_PEAKS (
        ch_ip_peaks_frip,
        ch_peak_count_header,
        ch_frip_score_header
    )
    ch_custompeaks_frip_multiqc  = MULTIQC_CUSTOM_PEAKS.out.frip
    ch_custompeaks_count_multiqc = MULTIQC_CUSTOM_PEAKS.out.count

    if (!params.skip_peak_annotation) {
        //
        // MODULE: Annotate peaks with MACS2
        //
        HOMER_ANNOTATEPEAKS_MACS2 (
            ch_macs2_peaks,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gtf
        )
        ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS_MACS2.out.versions.first())

        if (!params.skip_peak_qc) {
            //
            // MODULE: MACS2 QC plots with R
            //
            PLOT_MACS2_QC (
                ch_macs2_peaks.collect{it[1]}
            )
            ch_versions = ch_versions.mix(PLOT_MACS2_QC.out.versions)

            //
            // MODULE: Peak annotation QC plots with R
            //
            PLOT_HOMER_ANNOTATEPEAKS (
                HOMER_ANNOTATEPEAKS_MACS2.out.txt.collect{it[1]},
                ch_peak_annotation_header,
                "_peaks.annotatePeaks.txt"
            )
            ch_plothomerannotatepeaks_multiqc = PLOT_HOMER_ANNOTATEPEAKS.out.tsv
            ch_versions = ch_versions.mix(PLOT_HOMER_ANNOTATEPEAKS.out.versions)
        }
    }

    //
    //  Consensus peaks analysis
    //  Here the aim is to generate a global Consensus and a "By_Condition" consensus
    //  Consider selecting by IDR score as best ENCODE practice:
    //  

    ch_macs2_consensus_bed_lib   = Channel.empty()
    ch_macs2_consensus_txt_lib   = Channel.empty()
    ch_deseq2_pca_multiqc        = Channel.empty()
    ch_deseq2_clustering_multiqc = Channel.empty()

    // It makes by default a consensus - this is used to quantify and compute scaling FACTORS:

    // Create channels: [ meta , [ peaks ] ]
    // Where meta = [ id:antibody, multiple_groups:true/false, replicates_exist:true/false ]

    ch_macs2_peaks
        .map { 
            meta, peak -> 
                [ meta.antibody, meta.id.split('_')[0..-2].join('_'), peak ] 
        }
        .groupTuple()
        .map {
            antibody, groups, peaks ->
                [
                    antibody,
                    groups.groupBy().collectEntries { [(it.key) : it.value.size()] },
                    peaks
                ] 
        }
        .map {
            antibody, groups, peaks ->
                def meta_new = [:]
                meta_new.id = antibody
                meta_new.multiple_groups = groups.size() > 1
                meta_new.replicates_exist = groups.max { groups.value }.value > 1
                [ meta_new, peaks ] 
        }
        .set { ch_antibody_peaks }
    
    ch_antibody_peaks.view()
    //
    //  MODULE: Generate consensus peaks across samples
    //  Consider modifying this: i.e. Merge by condition - using IDR score any peak coming out of this will be a true potential peak
    //  Final Get all By_condition peak and perform a final merge - with min_overlap to consider equality across condition 
    //  A final summit has to be computed running MACS2 on all BAMs by antibody - create a channel with sample vs inputs or samples alone and run MACS2
    //

    MACS2_CONSENSUS ( 
        ch_antibody_peaks
    )
    ch_macs2_consensus_bed_lib = MACS2_CONSENSUS.out.bed
    ch_macs2_consensus_txt_lib = MACS2_CONSENSUS.out.txt
    ch_versions = ch_versions.mix(MACS2_CONSENSUS.out.versions)

    if (!params.skip_peak_annotation) {
        //
        // MODULE: Annotate consensus peaks
        //
        HOMER_ANNOTATEPEAKS_CONSENSUS (
            MACS2_CONSENSUS.out.bed,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gtf
        )
        ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS_CONSENSUS.out.versions)
        //
        // MODULE: Add boolean fields to annotated consensus peaks to aid filtering
        //
        ANNOTATE_BOOLEAN_PEAKS (
            MACS2_CONSENSUS.out.boolean_txt.join(HOMER_ANNOTATEPEAKS_CONSENSUS.out.txt, by: [0]),
        )
        ch_versions = ch_versions.mix(ANNOTATE_BOOLEAN_PEAKS.out.versions)
    }

    // Create channels: [ antibody, [ ip_bams ] ]
    ch_ip_control_bam
        .map { 
            meta, ip_bam, control_bam ->
                [ meta.antibody, ip_bam ]
        }
        .groupTuple()
        .set { ch_antibody_bams }
    

    // Create channels: [ meta, [ ip_bams ], saf ]
    MACS2_CONSENSUS
        .out
        .saf
        .map { 
            meta, saf -> 
                [ meta.id, meta, saf ] 
        }
        .join(ch_antibody_bams)
        .map {
            antibody, meta, saf, bams ->
                [ meta, bams.flatten().sort(), saf ]
        }
        .set { ch_saf_bams }

    //
    // MODULE: Quantify peaks across samples with featureCounts
    //
    SUBREAD_FEATURECOUNTS (
        ch_saf_bams
    )
    ch_subreadfeaturecounts_multiqc = SUBREAD_FEATURECOUNTS.out.summary
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

    //
    // Normalize samples with Gualdrini et al. 2016 Method compute scaling factor and generate a channel - [meta , bam, scaling]
    // Here we compute the scaling factor used to normalize the counts and the bigWigs
    // Generate counts as SUBREAD so norm counts rounded for DESEQ2 to work
    // out has: counts, scalings, 
    //
    // MODULE: Compute normalization and quality plots:
    //
    COUNT_NORM (
        SUBREAD_FEATURECOUNTS.out.counts,
        ch_deseq2_pca_header,
        ch_deseq2_clustering_header
    )
    ch_deseq2_pca_multiqc        = COUNT_NORM.out.pca_multiqc
    ch_deseq2_clustering_multiqc = COUNT_NORM.out.dists_multiqc
    
    COUNT_NORM
        .out
        .noamlization_txt
        .splitCsv ( header:true, sep:'\t' )
        .map { row -> 
            def id = row.Sample_ID
            def value = row.scaling
            [ id, value ]
        }
        .set { ch_size_factors }

    // Assemble the channel 
    // Given a tab separated matrix with the first column : Sample_id, Scaling_factor convert the matrix to a channel with [Sample_id, Scaling_factor] pairs
    // Consider that the first line is the header - in principle Sample_id must match the meta.id from BAM_FILTER_EM.out.bam    

    ch_genome_bam_bai
        .combine(ch_size_factors)
        .map { 
            meta1, bam1, bai1, id2, scaling2 ->
                meta1.id == id2 ? [ meta1, bam1, bai1 ,scaling2] : null
        }
        .set { ch_bam_bai_scale } 
        
    ch_bam_bai_scale.view()

    ch_deeptoolsplotprofile_multiqc = Channel.empty()
    // DEEPTOOLS_BIGWIG_NORM.out.bigwig && !DEEPTOOLS_BIGWIG_NORM.out.bigwig.ifEmpty([])
    // execute the below only if params.normalize is false i.e. !params.normalize OR we have one sample OR ch_size_factors is empty
    //
    // Scale to depth of sequencing using Deeptools:
    // 
    DEEPTOOLS_BIGWIG (
        ch_genome_bam_bai
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BIGWIG.out.versions.first())
    ch_big_wig = DEEPTOOLS_BIGWIG.out.bigwig

    if ( params.normalize ) {
        //
        // MODULE: BedGraph coverage tracks.
        //
        DEEPTOOLS_BIGWIG_NORM (
            ch_bam_bai_scale // join with the created channel of scalings - Use deeptool insrtead make directly wiggles
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_BIGWIG_NORM.out.versions.first())
        ch_big_wig = DEEPTOOLS_BIGWIG_NORM.out.bigwig
    } 
    
    if (!params.skip_plot_profile ) {

        // Add an if so that if DEEPTOOLS_BIGWIG_NORM.out.bigwig is empty it will use DEEPTOOLS_BIGWIG.out.bigwig
        // MODULE: deepTools matrix generation for plotting
        //
            
        DEEPTOOLS_COMPUTEMATRIX (
            ch_big_wig,
            PREPARE_GENOME.out.gene_bed
        )

        ch_versions = ch_versions.mix(DEEPTOOLS_COMPUTEMATRIX.out.versions.first())

        //
        // MODULE: deepTools profile plots
        //
        DEEPTOOLS_PLOTPROFILE (
            DEEPTOOLS_COMPUTEMATRIX.out.matrix
        )
        ch_deeptoolsplotprofile_multiqc = DEEPTOOLS_PLOTPROFILE.out.table
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTPROFILE.out.versions.first())

        //
        // MODULE: deepTools heatmaps
        //
        DEEPTOOLS_PLOTHEATMAP (
            DEEPTOOLS_COMPUTEMATRIX.out.matrix
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTHEATMAP.out.versions.first())
    } 

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowChipseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),

            ch_fastqc_raw_multiqc.collect{it[1]}.ifEmpty([]),
            ch_fastqc_trim_multiqc.collect{it[1]}.ifEmpty([]),
            ch_trim_log_multiqc.collect{it[1]}.ifEmpty([]),

            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),

            MARK_DUPLICATES_PICARD.out.stats.collect{it[1]}.ifEmpty([]),
            MARK_DUPLICATES_PICARD.out.flagstat.collect{it[1]}.ifEmpty([]),
            MARK_DUPLICATES_PICARD.out.idxstats.collect{it[1]}.ifEmpty([]),
            MARK_DUPLICATES_PICARD.out.metrics.collect{it[1]}.ifEmpty([]),

            BAM_FILTER_EM.out.stats.collect{it[1]}.ifEmpty([]),
            BAM_FILTER_EM.out.flagstat.collect{it[1]}.ifEmpty([]),
            BAM_FILTER_EM.out.idxstats.collect{it[1]}.ifEmpty([]),
            ch_picardcollectmultiplemetrics_multiqc.collect{it[1]}.ifEmpty([]),
    
            ch_deeptoolsplotprofile_multiqc.collect{it[1]}.ifEmpty([]),
            ch_deeptoolsplotfingerprint_multiqc.collect{it[1]}.ifEmpty([]),
    
            PHANTOMPEAKQUALTOOLS.out.spp.collect{it[1]}.ifEmpty([]),
            MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.nsc.collect{it[1]}.ifEmpty([]),
            MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.rsc.collect{it[1]}.ifEmpty([]),
            MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.correlation.collect{it[1]}.ifEmpty([]),

            ch_custompeaks_frip_multiqc.collect{it[1]}.ifEmpty([]),
            ch_custompeaks_count_multiqc.collect{it[1]}.ifEmpty([]),
            ch_plothomerannotatepeaks_multiqc.collect().ifEmpty([]),
            ch_subreadfeaturecounts_multiqc.collect{it[1]}.ifEmpty([]),

            ch_deseq2_pca_multiqc.collect().ifEmpty([]),
            ch_deseq2_clustering_multiqc.collect().ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
