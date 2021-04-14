#!/usr/bin/env python
# -*- coding: utf-8 -*-

"This program is a hierarchical binning based on distance between contigs. It takes in input the fasta file of the contigs, a busco table, an illumina coverage table and a Pacbio coverage table "

##########################################################################
import argparse
import os
import gc
from Bio import SeqIO
from config_file_reader import ConfigFileReader
from table_reader_printer import tableReaderPrinter
from contig import Contig
from contig_collection import ContigCollection
from bin_collection import BinCollection
from tetranu_distance_computer_v4 import tetranuDistanceComputer
from average_distance import averageDistance
##########################################################################

""" The arguments are parsed in the main and are dispatched to the appropriate function/object."""

# command line #
##############################
# /home/aline/teletravail/Scripts_Genomic/new_home_binning/home_binning_v9_for_hybrid_binning/main.py -f <configuration_file>

############################## inputs ####################################

parser = argparse.ArgumentParser(description = 'Iterative home binning based on multiple distances')

# change the parser if a configuration file is provided
parser.add_argument('-f', '--config', help = 'configuration file', required = False)
args = parser.parse_args()

if args.config != None:
    parser = ConfigFileReader(args.config)


# arguments

parser.add_argument('-i', '--input_fasta', help = 'the input fasta file containing the sequences of the contigs', required = True)
parser.add_argument('-ci', '--coverage_illumina', help = 'illumina coverage table for each contig', required = True)
parser.add_argument('-cp', '--coverage_pacbio', help = 'pacbio coverage table for each contig', required = True)
parser.add_argument('-b', '--busco_table', help = 'busco table for each contig', required = True)
parser.add_argument('-e', '--test', help = 'option set to True for the test mode', required = False)

parser.add_argument('-s', '--distance_threshold', help = 'set the maximum distance between a contig and a bin', required = True)
parser.add_argument('-p', '--ponderation', help = 'list of ponderations for distance (4mer,cov_illumina,cov_pacbio,busco)', required = True)
parser.add_argument('-n', '--normalization', help = 'list of normalization factors for distances (4mer,cov_illumina,cov_pacbio)', required = True)
parser.add_argument('-o', '--output_folder', help = 'output folder will contain a summary table of the binning and a fasta file for each bin', required = False)
parser.add_argument('-m', '--output_matrices', help = 'this option output the different distance matrices', required = False)
##############################  arguments  ##############################

if args.config != None:
    parser.read()

args = parser.parse_args()

test_mode = args.test

#input filename of the fasta file containing all the contigs of the assembly
input_filename = args.input_fasta

# output_folder
output_folder_name = args.output_folder

# coverage tables
coverage_illumina_table = args.coverage_illumina
coverage_pacbio_table = args.coverage_pacbio

# busco table
busco_table = args.busco_table

# option

# parameters
distance_threshold = float(args.distance_threshold)
ponderation_list = args.ponderation.split(',')
normalization_factor_list = args.normalization.split(',')


###########################################################################
### 0. Test if the files paths are correct                              ###
###########################################################################
# test output_folder
output_summary_table = open(output_folder_name + '/bin_summary_table.tsv','w')
output_summary_table.close()

# test coverage input
input_test = open(coverage_illumina_table,'r')
input_test.close()

# test coverage pacbio input
input_test = open(coverage_pacbio_table ,'r')
input_test.close()

# test busco input
input_test = open(busco_table,'r')
input_test.close()


###########################################################################
### 1. Initialization of the contig and bin collections                 ###
###########################################################################

contig_collection = ContigCollection()

###########################################################################
### 2. construct the initial contig collection (keep only the > 2000bp) ###
###########################################################################

print('\n***Collecting contigs from assembly***')

# open the input file as a handle
with open(input_filename, "r") as input_contig_handle:
    for record in SeqIO.parse(input_contig_handle, "fasta"):
        
        # if the sequence is > 2kb create a contig and add it in the contig collection
        if len(record.seq) > 2000 :
            contig = Contig(record)
            contig_collection.addContig(contig)
            # print(contig.record.id)

# sort the contig_collection by size from the biggest to the smallest
contig_collection.sortContigsBySize()


print('\n***Initialize average computer***')

# initialise the distance computer and collect information from coverage file and busco file
# 5 empty tables (as an empty dictionnary) are initiated, the illumina, pacbio and busco table 
# are filled from input

average_computer = averageDistance(coverage_illumina_table,coverage_pacbio_table,busco_table,ponderation_list,contig_collection)
#average_computer.busco_distance_computer.complete_dictionary(contig_collection)

###########################################################################
### 3. compute the distance matrices between all the individual contigs ###
###                                                                     ###
###########################################################################


print('\n*** Compute distance tables between contigs ***')
   
# a. compute the non normalized distances between all the contig pairs for every distance, create an    average matrix from all these distances
# the list of all the contigs is the list of all the values of contig_collection.dict sizes : 5xnxn

### STEP1 is in average_computer after the computation of tetranucleotide distance ###
if args.output_matrices == 'T' or args.output_matrices == 'True' :
    average_computer.compute_distance(list(contig_collection.dict.values()),list(contig_collection.dict.values()))

    # print the not normalized tables in tsv output_file
    print('\nprint the non_normalized intermediate tables in output files')

    tetranu_table_printer = tableReaderPrinter()
    tetranu_table_printer.output_in_file(output_folder_name + '/nn_tetranu_distance_table.tsv',average_computer.tetranu_distance_computer.distance_table) 


    coverage_table_printer = tableReaderPrinter()
    coverage_table_printer.output_in_file(output_folder_name + '/nn_coverage_distance_table.tsv',average_computer.coverage_distance_computer.distance_table) 


    pacbio_table_printer = tableReaderPrinter()
    pacbio_table_printer.output_in_file(output_folder_name + '/nn_pacbio_distance_table.tsv',average_computer.pacbio_distance_computer.distance_table) 


    busco_table_printer = tableReaderPrinter() 
    busco_table_printer.output_in_file(output_folder_name + '/nn_busco_distance_table.tsv',average_computer.busco_distance_computer.distance_table) 


    print('collecting garbage')
    gc.collect() 


    # 4. normalize the distance that have to be normalized, output the normalized tables
    average_computer.normalize_distance(normalization_factor_list)
        
    
    print('\nprint the normalized intermediate tables in output files')

    # print the tables in output_files
     
    # print the tetranucleotide distance table
    tetranu_table_printer = tableReaderPrinter()
    tetranu_table_printer.output_in_file(output_folder_name + '/tetranu_distance_table.tsv',average_computer.tetranu_distance_computer.distance_table) 
    tetranu_table_printer = None

    # print the short-read coverage distance table
    coverage_table_printer = tableReaderPrinter()
    coverage_table_printer.output_in_file(output_folder_name + '/coverage_distance_table.tsv',average_computer.coverage_distance_computer.distance_table) 
    coverage_table_printer = None

    # print the long-read coverage distance table
    pacbio_table_printer = tableReaderPrinter()
    pacbio_table_printer.output_in_file(output_folder_name + '/pacbio_distance_table.tsv',average_computer.pacbio_distance_computer.distance_table) 
    pacbio_table_printer = None

    # print the BUSCO distance table
    busco_table_printer = tableReaderPrinter() 
    busco_table_printer.output_in_file(output_folder_name + '/busco_distance_table.tsv',average_computer.busco_distance_computer.distance_table) 
    busco_table_printer = None

    #print('collecting garbage')
    gc.collect() 

    print('\nprint the average distance table in output files')

    #  output the average distance table
    average_table_printer = tableReaderPrinter()
    average_table_printer.output_in_file(output_folder_name + '/average_distance_table.tsv',average_computer.average_distance_table) 
  

###########################################################################
### 5. iteratively construct the bins                                   ###
###########################################################################


print('\n*** starting binning ***')

#use the first contig to create a bin, remove this contig from sorted list
bin_collection = BinCollection()
bin_collection.create_new_bin(contig_collection.dict[contig_collection.sorted_name_list.pop(0)]) 
compte = 2

# go through the contigs in decreasing order of length
for contig_name in contig_collection.sorted_name_list :

    print('\nbinning conting ', contig_name)

    compte += 1    # At the first pass the count is 3

    #get the contig corresponding to the contig name
    contig = contig_collection.dict[contig_name]

    # input_list1 contains the current contig only
    input_list1 = []
    input_list1 = [contig]

    # input_list2 contains the bins existing so far
    input_list2 = []
    input_list2 = list(bin_collection.dict.values())   # the names only : list(bin_collection.dict.keys())

    # compute the distance between the current contig and the existing bins normalize the distances
    # will create 5 vectors of distance of maximal size equal to the number of bis
    average_computer.compute_distance(input_list1,input_list2)

    binning_distance_table = average_computer.normalize_distance(normalization_factor_list)

    # get the bin_name corresponding to the lowest value (distance between a bin and the current contig) 
    selected_bin_name = min(binning_distance_table[contig_name], key=binning_distance_table[contig_name].get) 

    # print information about the contig and the chosen bin
    #print('distance_threshold : ',distance_threshold)
    print('distance between ',contig_name, ' and ', selected_bin_name, ' : ',binning_distance_table[contig_name][selected_bin_name])
    #print('distances all ', binning_distance_table[contig_name]) 

    # if the distance is lower than the threshold, the contig is added to the selected bin
    if binning_distance_table[contig_name][selected_bin_name] <= distance_threshold :
        
        #get the number of the selected bin and add the contig to this bin
        selected_bin_number = selected_bin_name.replace('bin_','')
        bin_collection.add_contig_to_bin(contig,selected_bin_number) 

        print('contig ',contig.record.id,' is added to bin bin_',selected_bin_number)

    # else it is used to create a new bin :
    else :
        print('contig ', contig.record.id, ' is used to form a new bin ; bin number ',bin_collection.current_number)
        bin_collection.create_new_bin(contig)

    #print('collecting garbage')
    gc.collect()     

###########################################################################
### 6. output the bin_summary_table                                     ###
###########################################################################


output_summary_table = open(output_folder_name + '/bin_summary_table.tsv','w')
output_contig_bin = open(output_folder_name + '/contig_bin_correspondance_with_length.tsv','w')
output_tag_list = open(output_folder_name + '/tagged_contig_lists.tsv','w')

#print header in bin_summary_table.tsv
header = 'bin' + '\t' + 'number_of_contigs' + '\t' + 'total_length' + '\t' + 'contig_list' + '\n'
output_summary_table.write(header)

header2 = 'contig' + '\t' + 'bin' + '\t' + 'contig_len' +  '\n'
output_contig_bin.write(header2)

# go through the bins of the bin collection
for bin_number in bin_collection.dict :
    
    # there is no list with only the names of the contigs so we create one
    contig_name_list = []

    # create a fasta file per bin
    bin_filename = output_folder_name + '/bin_'+ str(bin_number)+ '.fasta'
    bin_fasta_file = open(bin_filename,'w')
     
    # the current bin is stored in our_bin
    our_bin = bin_collection.dict[bin_number]

    # go through the contigs of our_bin, write the contigs in the fasta of the bin add the list to the name 
    for contig in our_bin.contig_list :
        SeqIO.write(contig.record,bin_fasta_file,"fasta")
        contig_name_list.append(contig.record.id)

        # write the correspondance between the bin and the contig fro each contig in the correspondance table
        output_line2 = contig.record.id + '\t' + str(our_bin.number) + '\t' + str(len(contig.record.seq))+'\n'
        output_contig_bin.write(output_line2)

    # collect and print summary information in the summary table
    output_line = 'bin_'+ str(our_bin.number) + '\t' + str(len(our_bin.contig_list)) + '\t' + str(our_bin.length) + '\t' + str(contig_name_list) + '\n' 

    output_summary_table.write(output_line)

    bin_fasta_file.close()


# output the tag file
for our_bin in list(bin_collection.dict.values())  :
    our_bin.evaluate_N50_contigs()

header = 'contig_name' + '\t' + 'no_busco' + '\t' +  'bigger_N50' + '\n'

output_tag_list.write(header)

for contig_name in contig_collection.dict :
    output_line = contig_name + '\t'
    if 'no_busco' in contig_collection.dict[contig_name].tags:
        output_line += '1' + '\t'
    else : 
        output_line += '0' + '\t'

    if 'bigger_N50' in contig_collection.dict[contig_name].tags:
        output_line += '1' + '\n'
    else : 
        output_line += '0' + '\n' 
   
    output_tag_list.write(output_line)

output_tag_list.close()
output_summary_table.close()
output_contig_bin.close()


# to implement  
if test_mode :
    pass 









