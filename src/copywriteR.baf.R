library(CopywriteR)
library(BiocParallel)

# Input files
args = commandArgs(trailingOnly=TRUE)

sample_name = args[1]
resolution = args[2]

bam_file = file.path("/home/arendeiro/baf-kubicek/data/merged", paste0(sample_name, ".merged.sorted.bam"))
output_dir = paste0("/home/arendeiro/baf-kubicek/results/copywriter/", resolution, "_", sample_name)
dir.create(output_dir)

# Control files
control_name = list("ATAC-seq_HAP1_WT")
control_file = file.path("/home/arendeiro/baf-kubicek/data/merged", paste0(control_name, ".merged.sorted.bam"))

# Params
sample.control <- data.frame(c(bam_file, control_file), control_file)
capture.regions.file <- "/home/arendeiro/baf-kubicek/results/baf-kubicek_peak_set.extended_5kb.bed"
bp.param <- SnowParam(workers=2, type="SOCK")

# Run
CopywriteR(
    sample.control=sample.control,
    destination.folder=output_dir,
    reference.folder=file.path("~/copywriter_genome_files", paste0("hg19_", resolution, "_chr")),
    bp.param=bp.param,
    capture.regions.file=capture.regions.file,
    keep.intermediary.files=TRUE)
plotCNA(destination.folder=output_dir, y.min=-4, , y.max=4)


# Load segmentation as plain text CSV
load(file.path(output_dir, "CNAprofiles", "segment.Rdata"))

write.csv(segment.CNA.object$data, file.path(output_dir, "data.csv"), row.names=FALSE, sep=",")
write.csv(segment.CNA.object$output, file.path(output_dir, "output.csv"), row.names=FALSE, sep=",")
write.csv(segment.CNA.object$segRows, file.path(output_dir, "segRows.csv"), row.names=FALSE, sep=",")


# RUN LIKE THIS:
find projects/baf-kubicek/data/merged/ -type f -name "*.bam"  | while read f ; do
NAME=`basename ${f/.merged.sorted.bam/}`;
echo $NAME;

# sbatch -J CopywriteR_1000kb_${NAME} -o projects/baf-kubicek/results/copywriter/1000kb_${NAME}.log ~/jobs/copywriteR.baf.sh $NAME 1000kb
# sbatch -J CopywriteR_100kb_${NAME} -o projects/baf-kubicek/results/copywriter/100kb_${NAME}.log ~/jobs/copywriteR.baf.sh $NAME 100kb
# sbatch -J CopywriteR_20kb_${NAME} -o projects/baf-kubicek/results/copywriter/20kb_${NAME}.log ~/jobs/copywriteR.baf.sh $NAME 20kb
sbatch -J CopywriteR_10kb_${NAME} -o projects/baf-kubicek/results/copywriter/10kb_${NAME}.log ~/jobs/copywriteR.baf.sh $NAME 10kb
sleep 60

done
