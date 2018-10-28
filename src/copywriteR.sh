#!/bin/env sh

# params
OUTPUT_DIR=~/projects/baf_complex/copywriter/processed
BIG_CONTROL=/home/arendeiro/cll-time_course/data/merged/nonVDJ.merged.sorted.bam
SMALL_CONTROL=/home/arendeiro/cll-time_course/data/merged/nonVDJ.merged.sorted.subsample.bam
SMALL_CONTROL_SIZE=30000000
GENOME=hg19
CAPTURE=~/projects/baf_complex/results/baf_complex_peak_set.bed
EXTENDED_CAPTURE=~/projects/baf_complex/results/baf_complex_peak_set.extended_5kb.bed
RESOLUTIONS=(1000kb 100kb 20kb 10kb)

# prepare
mkdir -p $OUTPUT_DIR
bedtools slop -b 5000 -i $CAPTURE -g ~/resources/genomes/${GENOME}/${GENOME}.chrom.sizes > $CAPTURE
samtools view -@ 8 -b -s 0.0685 $BIG_CONTROL > $SMALL_CONTROL
samtools index $SMALL_CONTROL

# run
find ~/projects/baf_complex/data/ -type f -name "ATAC-seq*filtered.bam" | while read F ; do

    for RESOLUTION in ${RESOLUTIONS[@]}; do

        SAMPLE_NAME=`basename ${F/.trimmed.bowtie2.filtered.bam/}`;
        ARGS=($F $SMALL_CONTROL $GENOME $CAPTURE $RESOLUTION $OUTPUT_DIR $SAMPLE_NAME)
        RUN_NAME=${SAMPLE_NAME}_${RESOLUTION}

        if [ -f ${OUTPUT_DIR}/${RUN_NAME}/CNAprofiles/log2_read_counts.igv ]; then
            continue
        fi

        echo $RUN_NAME

        sbatch -p longq -c 8 --mem 64000 --time 8-00:00:00 -J CopywriteR_${RUN_NAME} -o ${OUTPUT_DIR}/${RUN_NAME}.log --x11 ~/jobs/copywriteR.job.R ${ARGS[@]}
        # sleep 5
    done

    # sleep 3

done

# proceed with copywriteR.py
