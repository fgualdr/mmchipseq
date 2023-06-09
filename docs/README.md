# fgualdr/nf_chipmm: Documentation

The fgualdr/mmchip is a nextflow pipeline modification of the nf-core/chipseq pipeline which incorporate MultiMapping Reads placement via Bayesian-EM algorithm

The Pipe at the moment works with PE libraries and at first maps reads using star allowing a set of maximum multi mappers per reads.

The code can be lounched like so:


mmchip=<PATH to THE mmchip folder>
sample_sheet=<PATH to THE sample sheet this has to be a csv file with the following format: rid,sid,sample,replicate,path,lanes,is_input,which_input,antibody>
work_dir=<PATH to THE working directory>

fasta=<PATH to THE fasta file>
gtf=<PATH to THE gtf file>
star_index=<PATH to THE star index if available>

read_length=<READ LENGTH>

nextflow -q run $mmchip/main.nf \
            --input $sample_sheet \
            --outdir $work_dir \
            --fasta $fasta \
            --gtf $gtf \
            --star_index $star_index \
            --publish_dir_mode symlink \
            --read_length $read_length \
            --em_eps 0.0000000001 \
            --em_iter 500 \
            --macs_fdr 0.00001 \
            --outfiltermultimapnmax 50 \
            --outsammultnmax 50 \
            --winanchormultimapnmax 50 \
            --with_inputs false \
            --keep_dups true \
            --keep_multi_map true \
            --min_reps_consensus 2 \
            --macs_model false \
            --narrow_peak true \
            --normalize false \
            -profile <the profile to be used Docker or Singularity> 
