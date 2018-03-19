<!--[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.231352.svg)](https://doi.org/10.5281/zenodo.231352)-->

### Systematic functional characterization of BAF mutations yields novel intra-complex synthetic lethalities

Sandra Schick, André F. Rendeiro, Bernd Boidol, Kathrin Runggatscher, Melanie Hinkel, Peter Májek, Thomas Penz, Katja Parapatics, Christian Schmidl, Anna Ringler, Guido Boehmelt, Mark  Petronczki, André Mueller, Christoph Bock, Stefan Kubicek

**Paper**: [http://dx.doi.org/](http://dx.doi.org/)

**Website**: [http://www.medical-epigenomics.org/papers/schick2018/](http://http://www.medical-epigenomics.org/papers/schick2018/)

This repository contains scripts used in the analysis of the data in the paper.

<br>

#### Analysis

In the [paper website](http://http://www.medical-epigenomics.org/papers/schick2018/) you can find most of the output of the whole analysis.

Here are a few steps needed to reproduce it:

1. Clone the repository: `git clone git@github.com:epigen/baf_complex.git`
2. Install required software for the analysis:`make requirements` or `pip install -r requirements.txt`

If you wish to reproduce the processing of the raw data (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108390), run these steps:

2. Download the data localy from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108390).
3. Prepare [Looper](https://github.com/epigen/looper) configuration files similar to [these](metadata/project_config.yaml) that fit your local system.
4. Run samples through the pipeline: `looper run metadata/project_config_file.yaml`
5. Get external files (genome annotations mostly): `make external_files` or use the files in the [paper website](http://http://www.medical-epigenomics.org/papers/schick2018/) (`external` folder).
6. Run the analysis: `make analysis`

Additionaly, processed (bigWig and narrowPeak files together with a gene expression matrix) are available from [GEO with accession number GSE108390](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108390).
