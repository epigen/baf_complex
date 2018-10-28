.DEFAULT_GOAL := analysis

requirements:
	pip install -r requirements.txt

process:
	looper run metadata/project_config.yaml

analysis:
	# python -u src/analysis.cnv.py
	python -u src/analysis.atacseq.py
	python -u src/analysis.rnaseq.py
	python -u src/analysis.chipseq.py
	python -u src/analysis.cross_data_type.py

analysis_job:
	mkdir -p log
	TIMESTAMP=`date +"%Y%m%d-%H%M%S"`
	sbatch -p shortq -c 12 --mem 80000 \
	-J baf_complex.cnv-analysis \
	-o log/${TIMESTAMP}.baf_complex.cnv-analysis.log \
	--wrap "python -u src/analysis.cnv.py"
	sbatch -p shortq -c 12 --mem 80000 \
	-J baf_complex.atacseq-analysis \
	-o log/${TIMESTAMP}.baf_complex.atacseq-analysis.log \
	--wrap "python -u src/analysis.atacseq.py"
	sbatch -p shortq -c 12 --mem 80000 \
	-J baf_complex.rnaseq-analysis \
	-o log/${TIMESTAMP}.baf_complex.rnaseq-analysis.log \
	--wrap "python -u src/analysis.rnaseq.py"
	sbatch -p shortq -c 12 --mem 80000 \
	-J baf_complex.chipseq-analysis \
	-o log/${TIMESTAMP}.baf_complex.chipseq-analysis.log \
	--wrap "python -u src/analysis.chipseq.py"
	sbatch -p shortq -c 12 --mem 80000 \
	-J baf_complex.cross_data_type-analysis \
	-o log/${TIMESTAMP}.baf_complex.cross_data_type-analysis.log \
	--wrap "python -u src/analysis.cross_data_type.py"

all: requirements process analysis

.PHONY: requirements process analysis all
