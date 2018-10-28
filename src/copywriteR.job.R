#!/bin/env Rscript

library(CopywriteR)
library(BiocParallel)

# Input files
args = commandArgs(trailingOnly=TRUE)

# args = c(
# "/home/arendeiro/projects/cll-timetotreat/data/ATAC-seq_MUW_CLL098_20170103_sorted/mapped/ATAC-seq_MUW_CLL098_20170103_sorted.trimmed.bowtie2.filtered.bam",
# "/home/arendeiro/cll-time_course/data/merged/nonVDJ.merged.sorted.bam",
# "hg19",
# "/home/arendeiro/projects/cll-timetotreat/results/cll-timetotreat_peak_set.extended_5kb.bed",
# "1000kb",
# "/home/arendeiro/projects/cll-timetotreat/copywriter",
# "ATAC-seq_MUW_CLL098_20170103_sorted")
bam_file = args[1]
control_file = args[2]
genome = args[3]
capture.regions.file = args[4]
resolution = args[5]
output_dir = args[6]
sample_name = args[7]

dir.create(output_dir)

output_dir = paste0(output_dir, "/", sample_name, "_", resolution)
dir.create(output_dir)

print(c(args, output_dir))

# # Control files
# control_name = list("nonVDJ")
# control_file = file.path("/home/arendeiro/cll-time_course/data/merged", paste0(control_name, ".merged.sorted.bam"))

# Params
sample.control <- data.frame(c(bam_file, control_file), control_file)
# capture.regions.file <- "/home/arendeiro/cll-time_course/results/cll-time_course_peak_set.extended_5kb.bed"
bp.param <- SnowParam(workers=8, type="SOCK")

# Run
CopywriteR(
    sample.control=sample.control,
    destination.folder=output_dir,
    reference.folder=file.path("~/copywriter_genome_files", paste(genome, resolution, "chr", sep="_")),
    bp.param=bp.param,
    capture.regions.file=capture.regions.file,
    keep.intermediary.files=FALSE)
plotCNA(destination.folder=output_dir, y.min=-4, , y.max=4)


# Load segmentation as plain text CSV
load(file.path(output_dir, "CNAprofiles", "segment.Rdata"))

write.csv(segment.CNA.object$data, file.path(output_dir, "data.csv"), row.names=FALSE, sep=",")
write.csv(segment.CNA.object$output, file.path(output_dir, "output.csv"), row.names=FALSE, sep=",")
write.csv(segment.CNA.object$segRows, file.path(output_dir, "segRows.csv"), row.names=FALSE, sep=",")


# # RUN LIKE THIS:
# # params
# OUTPUT_DIR=~/projects/cll-timetotreat/copywriter
# CONTROL=/home/arendeiro/cll-time_course/data/merged/nonVDJ.merged.sorted.bam
# GENOME=hg19
# CAPTURE=~/projects/cll-timetotreat/results/cll-timetotreat_peak_set.extended_5kb.bed
# RESOLUTIONS=(1000kb 100kb 20kb 10kb)

# # prepare
# mkdir -p $OUTPUT_DIR
# bedtools slop -b 5000 -i ~/projects/cll-timetotreat/results/cll-timetotreat_peak_set.bed -g ~/resources/genomes/hg19/hg19.chrom.sizes > $CAPTURE

# # run
# find ~/projects/cll-timetotreat/data/ -type f -name "*CLL*filtered.bam" | while read F ; do

# for RESOLUTION in ${RESOLUTIONS[@]}; do

# SAMPLE_NAME=`basename ${F/.trimmed.bowtie2.filtered.bam/}`;
# ARGS=($F $CONTROL $GENOME $CAPTURE $RESOLUTION $OUTPUT_DIR $SAMPLE_NAME)
# RUN_NAME=${SAMPLE_NAME}_${RESOLUTION}
# echo $RUN_NAME

# sbatch -p longq -c 8 --mem 24000 --time 8-00:00:00 -J CopywriteR_${RUN_NAME} -o ${OUTPUT_DIR}/${RUN_NAME}.log --x11 ~/jobs/copywriteR.job.R ${ARGS[@]}

# done
# sleep 12
# done

# # /home/arendeiro/projects/cll-timetotreat/data/ATAC-seq_MUW_CLL098_20170103_sorted/mapped/ATAC-seq_MUW_CLL098_20170103_sorted.trimmed.bowtie2.filtered.bam
