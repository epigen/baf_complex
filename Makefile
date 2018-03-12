.DEFAULT_GOAL := analysis

requirements:
	pip install -r requirements.txt

# process project's data
# with looper/pypiper/pipelines:
# see https://github.com/epigen/looper
# see https://github.com/epigen/pypiper
# see https://github.com/epigen/open_pipelines
process:
	looper run metadata/project_config.yaml

analysis:
	python -u src/analysis.py

analysis_job:
	mkdir -p log
	TIMESTAMP=`date +"%Y%m%d-%H%M%S"`
	sbatch -p shortq -c 12 --mem 80000 -J baf_complex.analysis -o log/${TIMESTAMP}.baf_complex.analysis.log --wrap "python -u src/analysis.py"

all: requirements process analysis

.PHONY: requirements process analysis all
