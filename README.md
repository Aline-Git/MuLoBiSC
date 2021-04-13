# MuLoBiSC
Multi-distance long-read binning with silhouette score optimization

## MuLoBI
Multi-distance long-read binning. MuLoBi was designed for binning contigs with coverage information from long-read data and 2 different short-read data. It can be adapted for other design by modifying the python script 'coverage_distance_computer'. 
MuLoBi takes as inuput (see 'example_test_files' for example of input files):

-i : a file in fasta format containing the sequences of the contigs

-ci : a short-read median coverage tab separated table : contig_name{tab}coverageA{tab}coverageB

-cp : a long-read median coverage tab separated table : contig _name{tab}coverage

-b : a BUSCO tab separated table with a line for each contig and the BUSCO status of each BUSCO gene on each column ('Missing', 'Fragmented', 'Complete', Duplicated').

-s : distance threshold (e.g. 0.035)

-p : a list of ponderations for distaces tetranucleotide,coverage short-read,coverage long-read, busco, (e.g. 0.07,0.014,0.004,0.4) 

-n : list of normalization factors for distances tetranucleotide, coverage short-read, coverage long-read (e.g median distances)

-o : output folder name

-m : 'True' for the output of distance matrices.


To run the binning :
```
MuLoBi/main.py -f example_test_files/configuration_test_file.txt
```

## Silhouette score optimization
Silhouette score optimisation will to refine the binning of the contigs by increasing the consistency of the global clustering. It also provides a metric, ranging for -1.0 to 1.0, on how well a contig is classified within a bin.

The script takes as input the 'average_distance_table.tsv' and the 'contig_bin_correspondance_with_length.tsv' from MuLoBi output. The working directory should be changed line 8. the input files should be placed in a a 'input_files' folder in the working directory. An empty 'output_files' should be placed in the working directory.

Simon Poirier wrote the silhouette score optimization with the participation of Marco Pagni and Aline Adler
