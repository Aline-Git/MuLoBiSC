#!/usr/bin/env python
# -*- coding: utf-8 -*-

" This module contains tetranuDistanceComputer, daugther of distanceComputer"


##########################################################################
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import math
import time
from table_reader_printer import tableReaderPrinter
from distance_computer import distanceComputer
from our_bin import Bin
##########################################################################


class tetranuDistanceComputer(distanceComputer):

    """
class tetranuDistanceComputer: this object compute a tetranucleotide distance between two sets (lists) of DNA sequences
    """

    # ------------------------------------------------------------------ #
    # Constructors/Destructors                                           #
    # ------------------------------------------------------------------ #

    def __init__(self):
        """__init__: the initialization of  tetranuDistanceComputer is the same as the mother distanceComputer with an empty dictionnary self.distance_table as member """

        distanceComputer.__init__(self)
    


    def __del__(self):
        """__del__: not implemented """
        pass

    # ------------------------------------------------------------------ #
    # Methods                                                            #
    # ------------------------------------------------------------------ #

    # public:

    def compute_distance(self,input_list1,input_list2,ponderation_list=None):
        """compute_distance : compute the euclidean distance table between the median 4mer frequency vectors of two lists of objects (contigs or bins)"""

        # reinitialize the distance table at the begining of each computation
        self.distance_table = {} 

        # compute the mean frequency vector of each item of the two lists
        # here we only take the union of the two input lists
        for item in set(input_list1) | set(input_list2) :
        
            # find the name of the item
	    # if the object is a bin
            if isinstance(item, Bin):
                item_name = 'bin_' + str(item.number)

	    # if the object is a contig
            else : 
                item_name = item.record.id

            # check if the median frequency vector of the object (contig or bin) is already computed
            # and if not compute it
            if item.d_4mers_median == {}:
                print('\ncompute median vector of item ', item_name)
                item.d_4mers_median = self.compute_4mer_median_vector(item)

        # compute the tetranucleotide euclidean distance between each items of the two lists
        for item1 in input_list1 :
	    # get the name of the item if is a bin 
            if isinstance(item1, Bin):
                item_name1 = 'bin_' + str(item1.number)

	    # get the name of the item if is a contig
            else : 
                item_name1 = item1.record.id

            # create a dictionnary for each entry of the first list
            self.distance_table[item_name1] = {}

            # go through each item of the second list and get their name
            for item2 in input_list2 :
                if isinstance(item2, Bin):
                    item_name2 = 'bin_' + str(item2.number) 
                else : 
                    item_name2 = item2.record.id


                # initiate distance between the 2 vectors as zero
                if item_name1 == item_name2 :
                    self.distance_table[item_name1][item_name2]= 0.0

                # compute the distance only if items are different

                # initiate distance between the 2 vectors as zero
                if item_name1 != item_name2 :

                    # check if it is already computed or not
                    distance_computed = False
       
                    # if the reverse entry already exist then distance_computed is set to true
                    try : 
                        self.distance_table[item_name2][item_name1]
                        distance_computed = True

                    except : 
                        distance_computed = False

                    # if the distance was already computed, we dont add another entry 
                    if not distance_computed :
                        self.distance_table[item_name1][item_name2]= 0.0
                    # if distance is not computed yet, compute it 
                        for oligonucl in item1.d_4mers_median:

                            self.distance_table[item_name1][item_name2] += (item2.d_4mers_median[oligonucl]-item1.d_4mers_median[oligonucl])**2
                        self.distance_table[item_name1][item_name2]=math.sqrt(self.distance_table[item_name1][item_name2])

        return self.distance_table





    def compute_4mer_median_vector(self,item):
        """compute_4mer_median_vector : calcule le vecteur median de frequences de tetranucleotides d'un object (contigs ou bin) et le stocke dans un dictionnaire item.d_4mers_median """

        # initialize a dictionnary with the 136 4mer as keys and values 0
        d_4mers_median = self.initialise_freq_vector()

        # dna sequences and name the item will be stored as item_sequence and item_name
        item_sequence = Seq("")
        item_name = ''

	# get their values if the object is a bin
        if isinstance(item, Bin):
            item_name = 'bin_' + str(item.number)
            for contig in item.contig_list :
                item_sequence += contig.record.seq

	# get their values if the object is a contig
        else : 
            item_name = item.record.id
            item_sequence = item.record.seq

        # initialization of the parameters used to compute the frequency vector, a 2000 bp
	# sliding window and a resolution of 200 bp
        sw_len = 2000
        resolution = 200

        # check if there is already a d_4mers_matrix associated to the object, if not compute it
        if item.d_4mers_matrix == {}:
       
            # compute the number of sliding windows
            n= int(math.floor((len(item_sequence) - sw_len)/resolution))
            print("number of sliding windows :", n)

	    # initiate a double dictionnary to store the frequency of each sliding window
            item.d_4mers_matrix = self.initialise_matrix()
            print('calculate 4 mer frequency of item ', item_name)

	    # compute and store the frequency on the sliding window in the matrix
            for j in range(0, n+1):
                if j%1000 == 0 :
                    print('iteration ',j , ' sur ' , n)

                # initiate the frequencies of all oligonucleotides and positions to zero
                for oligonucl in item.d_4mers_matrix:
                    item.d_4mers_matrix[oligonucl][j*resolution] = 0

                # go through the sliding window and compute the frequency of each oligonucleotide
                for i in range(j*resolution,(j*resolution)+sw_len-3):
                    # if the oligonucleotide equals its reverse complement, then add one to it
                    if str(item_sequence[i:i+4].reverse_complement())== str(item_sequence[i:i+4]):
                        item.d_4mers_matrix[str(item_sequence[i:i+4]).upper()][j*resolution]+=1/sw_len

                    # if the oligonucleotide differs from its reverse complement, then add one to one of it
                    else :
                        try :
                            item.d_4mers_matrix[str(item_sequence[i:i+4]).upper()][j*resolution]+=1/sw_len
                        except : 
                            item.d_4mers_matrix[str(item_sequence[i:i+4].reverse_complement()).upper()][j*resolution]+=1/sw_len

        # compute the median over all the sliding windows of a contig
        for oligonucl in item.d_4mers_matrix :
            d_4mers_median[oligonucl] = statistics.median(list(item.d_4mers_matrix[oligonucl].values()))

        return d_4mers_median
     

    def initialise_freq_vector(self):
        """ initialise_freq_vector : returns an empty dictionary with 4mers as keys with reverse-complement 4mers considered as only one 4mer and zeros as values """

        d_4_mers = {}
        alphabet = ['A','C','T','G']
        for nucleotide1 in alphabet:
            for nucleotide2 in alphabet:
                for nucleotide3 in alphabet:
                    for nucleotide4 in alphabet:
                        oligonucl = nucleotide1+nucleotide2+nucleotide3+nucleotide4		
                        if not self.complement(oligonucl) in d_4_mers:
                            d_4_mers[oligonucl]= 0
        return d_4_mers


    def initialise_matrix(self):
        """ initialise_matrix : returns an empty double dictionary with 4mers as keys with reverse-complement 4mers considered as only one 4mer and zeros as values """
        d_4mers_matrix = {}
        alphabet = ['A','C','T','G']
        for nucleotide1 in alphabet:
            for nucleotide2 in alphabet:
                for nucleotide3 in alphabet:
                    for nucleotide4 in alphabet:
                        oligonucl = nucleotide1+nucleotide2+nucleotide3+nucleotide4		
                        if not self.complement(oligonucl) in d_4mers_matrix:
                            d_4mers_matrix[oligonucl]= {}
        return d_4mers_matrix



########## complement return the reverse complement sequence of a given DNA sequence ##############

    def complement(self,dna_sequence):

        rc_dna_sequence = ''
        i=len(dna_sequence)-1
        while i >= 0:
            rc_dna_sequence = rc_dna_sequence + self.complement_nu(dna_sequence[i])
            i=i-1
    
        return rc_dna_sequence
    


    def complement_nu(self,nu):

        if nu == 'A': 
            return 'T'
        if nu == 'T': 
            return 'A'    
        if nu == 'G': 
            return 'C'    
        if nu == 'C': 
            return 'G'
##########################################################################


if __name__ == '__main__':
    pass
