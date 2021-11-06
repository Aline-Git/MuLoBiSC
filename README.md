# MuLoBiSC
Multi-distance long-read binning (MuLoBi) with silhouette coefficient optimization (SC). It is the source code of the binning algorithm presented in the paper " Disentangle genus microdiversity within a complex microbial community by using a multi-distance long-read binning method: Example of Candidatus Accumulibacter" by Aline Adler, Simon Poirier, Marco Pagni, Julien Maillard and Christof Holliger (in revision). The code is given as it is for the sake of transparency, but we do not intend to maintain it on the long way. It has been tested on linux (ubuntu 16.04 and 20.04) and mac (Big Sur 11.6). 


## MuLoBi
Multi-distance long-read binning. MuLoBi was designed for binning contigs with coverage information from long-read data and 2 different short-read data. It can be adapted for other design by modifying the python script 'coverage_distance_computer'. 

MuLoBi takes as inuput the following arguments :

-i : a file in fasta format containing the sequences of the contigs<br/>
-ci : a short-read median coverage tab separated table : contig_name{tab}coverageA{tab}coverageB<br/>
-cp : a long-read median coverage tab separated table : contig _name{tab}coverage<br/>
-b : a BUSCO tab separated table with a line for each contig and the BUSCO status of each BUSCO gene on each column ('Missing', 'Fragmented', 'Complete', Duplicated')<br/>
-s : distance threshold (e.g. 0.035)<br/>
-p : a list of ponderations for distaces tetranucleotide,coverage short-read,coverage long-read, busco, (e.g. 0.07,0.014,0.004,0.4)<br/>
-n : list of normalization factors for distances tetranucleotide, coverage short-read, coverage long-read (e.g median distances)<br/>
-o : output folder name<br/>
-m : 'True' for the output of distance matrices.

An exemple of command line is provided in the test demo 'launch.sh'.



## Silhouette coefficient optimization (SC)
Silhouette coefficient optimisation will refine the binning of the contigs by increasing the consistency of the global clustering. It also provides a metric, ranging for -1.0 to 1.0, on how well a contig is classified within a bin.

The script takes as input the 'average_distance_table.tsv' and the 'contig_bin_correspondance_with_length.tsv' from MuLoBi output. 

Simon Poirier wrote the silhouette coefficient optimization with the participation of Marco Pagni and Aline Adler.


## Prerequisites
MuLoBi is a python3 script, it requires packages Biophython and argparse. They can be installed for example, via pip.
 

```
pip3 install Biopython
pip3 install argparse
```

SC_binning_optimization.R is a R script, it requires the packages "dplyr", "factoextra", "cluster" and "dendextend". If needed, these packages can be installed in R and may require the installation of other dependencies, depending on the packages already installed. 

```
install.packages("dplyr")
install.packages("factoextra")
install.packages("cluster")
install.packages("dendextend")
```


## Test demo
To test MuLoBiSC on a small dataset. Clone the MuLoBiSC directory and launch the script 'launch.sh'.

```
git clone https://github.com/Aline-Git/MuLoBiSC.git
cd MuLoBiSC
./launch.sh
```

The resulting binning after silhouette coefficient optimisation will be in MuLoBiSC/tmp/test_SC_binning_optimization/binning_round5_SC_sort.csv.

