# BAF and PBAF structure

## Description

Try to infer the structure of the BAF and PBAF complexes from combined 
knockout and pulldown assays, using a genetic algorithm to find a structure 
fitting well to the observed data.


## Content

Graphs were generated on a computing cluster using the *.jl* files.    
Jupyter notebooks are also provided to depict the analyses performed locally on a laptop.


## Input files

Several external files may be needed to run the scripts:

* *ARID1A-data.csv*, *ARID1A-pval.csv*, *BRG1-data.csv* and *BRG1-pval.csv* are the averaged IP-MS abundances and their corresponding FDR, and can be downloaded from the [paper website](http://www.medical-epigenomics.org/papers/schick2018/).
* All fasta files were downloaded from the [BioMart](https://www.ensembl.org/biomart/martview/) as described in the paper.
* The file *BAF_genefamily.tsv* was downloaded from [HUGO](https://www.genenames.org) for the Gene Family *BAF complex*.


## Running information

Julia Version 0.6.2

### Laptop configuration

Platform Info:

 - OS: macOS (x86_64-apple-darwin14.5.0)
 - CPU: Intel(R) Core(TM) i7-7660U CPU @ 2.50GHz
 - WORD_SIZE: 64
 - BLAS: libopenblas (USE64BITINT DYNAMIC_ARCH NO_AFFINITY Prescott)
 - LAPACK: libopenblas64_
 - LIBM: libopenlibm
 - LLVM: libLLVM-3.9.1 (ORCJIT, broadwell)

23 required packages:

 - Bio                           0.4.7
 - BioSequences                  0.8.3
 - CSV                           0.2.2
 - Clustering                    0.9.1
 - DataFrames                    0.11.5
 - Distributions                 0.15.0
 - DocOpt                        0.3.0
 - EzXML                         0.6.2
 - GeneticAlgorithms             0.0.3
 - GraphIO                       0.2.0
 - GraphPlot                     0.2.0
 - Graphs                        0.9.0
 - IJulia                        1.7.0
 - IterTools                     0.2.1
 - Iterators                     0.3.1
 - JLD                           0.8.3
 - LightGraphs                   0.12.0
 - NetworkViz                    0.0.2
 - Plotly                        0.1.1
 - PlotlyJS                      0.10.1
 - Plots                         0.15.0
 - ProfileView                   0.3.0
 - Stats                         0.1.0

93 additional packages:

 - AutoHashEquals                0.2.0
 - Automa                        0.6.0
 - BGZFStreams                   0.1.4
 - BinDeps                       0.8.7
 - BinaryProvider                0.2.5
 - BioCore                       1.3.0
 - BioSymbols                    1.2.0
 - Blink                         0.6.2
 - Blosc                         0.3.0
 - BufferedStreams               0.4.0
 - Cairo                         0.5.1
 - Calculus                      0.2.2
 - CategoricalArrays             0.3.5
 - CodecZlib                     0.4.2
 - Codecs                        0.4.0
 - ColorTypes                    0.6.7
 - Colors                        0.8.2
 - Combinatorics                 0.5.0
 - CommonSubexpressions          0.0.1
 - Compat                        0.61.0
 - Compose                       0.5.5
 - Conda                         0.7.1
 - Contour                       0.4.0
 - DataStreams                   0.3.4
 - DataStructures                0.7.4
 - DiffResults                   0.0.3
 - DiffRules                     0.0.3
 - Distances                     0.5.0
 - DocStringExtensions           0.4.3
 - FileIO                        0.7.0
 - FixedPointNumbers             0.4.6
 - ForwardDiff                   0.7.3
 - FunctionalCollections         0.3.2
 - GeometryTypes                 0.4.4
 - Graphics                      0.2.0
 - Gtk                           0.13.1
 - GtkReactive                   0.4.0
 - HDF5                          0.8.8
 - Hiccup                        0.1.1
 - Homebrew                      0.6.2
 - HttpCommon                    0.4.0
 - HttpParser                    0.3.1
 - HttpServer                    0.3.1
 - IndexableBitVectors           0.1.2
 - IntervalSets                  0.2.0
 - IntervalTrees                 0.1.0
 - JSON                          0.17.1
 - LaTeXStrings                  0.3.0
 - Lazy                          0.12.0
 - LegacyStrings                 0.3.0
 - Libz                          0.2.4
 - LightXML                      0.6.0
 - MacroTools                    0.4.0
 - MbedTLS                       0.5.8
 - Measures                      0.1.0
 - Missings                      0.2.7
 - Mustache                      0.3.1
 - Mux                           0.2.3
 - NaNMath                       0.3.1
 - NamedTuples                   4.0.0
 - NearestNeighbors              0.3.0
 - Nullables                     0.0.4
 - PDMats                        0.8.0
 - ParserCombinator              1.7.11
 - Patchwork                     0.4.0
 - PlotThemes                    0.2.0
 - PlotUtils                     0.4.4
 - PlotlyBase                    0.1.1
 - Polynomials                   0.2.1
 - QuadGK                        0.2.0
 - Reactive                      0.6.0
 - RecipesBase                   0.2.3
 - Reexport                      0.1.0
 - Requests                      0.5.1
 - Requires                      0.4.3
 - Rmath                         0.3.2
 - Roots                         0.5.0
 - RoundingIntegers              0.0.3
 - SHA                           0.5.6
 - Showoff                       0.1.1
 - SimpleTraits                  0.6.0
 - SortingAlgorithms             0.2.0
 - SpecialFunctions              0.3.8
 - StaticArrays                  0.7.0
 - StatsBase                     0.20.1
 - StatsFuns                     0.5.0
 - ThreeJS                       0.3.0
 - TranscodingStreams            0.5.1
 - Twiddle                       0.4.0
 - URIParser                     0.3.1
 - WeakRefStrings                0.4.3
 - WebSockets                    0.4.0
 - ZMQ                           0.5.1

### Cluster configuration

Platform Info:

 -   OS: Linux (x86_64-pc-linux-gnu)
 -   CPU: Intel(R) Xeon(R) CPU E5-2650 0 @ 2.00GHz
 -   WORD_SIZE: 64
 -   BLAS: libopenblas (USE64BITINT DYNAMIC_ARCH NO_AFFINITY Sandybridge)
 -   LAPACK: libopenblas64_
 -   LIBM: libopenlibm
 -   LLVM: libLLVM-3.9.1 (ORCJIT, sandybridge)

12 required packages:

 - BioAlignments                 0.2.0
 - BioSequences                  0.8.3
 - CSV                           0.2.2
 - DataFrames                    0.11.5
 - Distributions                 0.15.0
 - DocOpt                        0.3.0
 - GraphIO                       0.2.0
 - HDF5                          0.8.8
 - HttpParser                    0.3.1
 - LightGraphs                   0.12.0
 - Plotly                        0.1.1
 - StatsBase                     0.20.1

63 additional packages:

 - AutoHashEquals                0.2.0
 - Automa                        0.6.0
 - BGZFStreams                   0.1.4
 - BinDeps                       0.8.7
 - BinaryProvider                0.2.5
 - BioCore                       1.3.0
 - BioSymbols                    1.2.0
 - Blink                         0.6.2
 - Blosc                         0.3.0
 - BufferedStreams               0.4.0
 - Calculus                      0.2.2
 - CategoricalArrays             0.3.5
 - CodecZlib                     0.4.2
 - Codecs                        0.4.0
 - ColorTypes                    0.6.7
 - Combinatorics                 0.5.0
 - Compat                        0.60.0
 - DataStreams                   0.3.4
 - DataStructures                0.7.4
 - DocStringExtensions           0.4.3
 - EzXML                         0.6.2
 - FileIO                        0.7.0
 - FixedPointNumbers             0.4.6
 - GenomicFeatures               0.2.1
 - Hiccup                        0.1.1
 - HttpCommon                    0.4.0
 - HttpServer                    0.3.1
 - IndexableBitVectors           0.1.2
 - IntervalTrees                 0.4.1
 - IterTools                     0.2.1
 - JLD                           0.8.3
 - JSON                          0.17.1
 - LaTeXStrings                  0.3.0
 - Lazy                          0.12.0
 - LegacyStrings                 0.3.0
 - Libz                          0.2.4
 - MacroTools                    0.4.0
 - MbedTLS                       0.5.8
 - Missings                      0.2.7
 - Mustache                      0.3.1
 - Mux                           0.2.3
 - NamedTuples                   4.0.0
 - Nullables                     0.0.4
 - PDMats                        0.8.0
 - ParserCombinator              1.7.11
 - PlotlyBase                    0.1.1
 - PlotlyJS                      0.10.1
 - Polynomials                   0.2.1
 - QuadGK                        0.2.0
 - Reexport                      0.1.0
 - Requests                      0.5.1
 - Requires                      0.4.3
 - Rmath                         0.3.2
 - SHA                           0.5.6
 - SimpleTraits                  0.6.0
 - SortingAlgorithms             0.2.0
 - SpecialFunctions              0.3.8
 - StatsFuns                     0.5.0
 - TranscodingStreams            0.5.1
 - Twiddle                       0.4.0
 - URIParser                     0.3.1
 - WeakRefStrings                0.4.3
 - WebSockets                    0.4.0
