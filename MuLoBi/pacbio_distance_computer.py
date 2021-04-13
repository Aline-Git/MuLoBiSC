#!/usr/bin/env python
# -*- coding: utf-8 -*-

" This module contains the class pacbioDistanceComputer, daugther of distanceComputer"


##########################################################################
from Bio import SeqIO
from Bio.Seq import Seq
import math
import time
from contig import Contig
from table_reader_printer import tableReaderPrinter
from distance_computer import distanceComputer
from our_bin import Bin
##########################################################################


class pacbioDistanceComputer(distanceComputer):

    """
class coverageDistanceComputer: this object compute the distance between two sets (lists) of DNA sequences
    """

    # ------------------------------------------------------------------ #
    # Constructors/Destructors                                           #
    # ------------------------------------------------------------------ #

    def __init__(self):
        """__init__: the initialization of coverageDistanceComputer inherits the self.distance_table = {} from the mother distanceComputer """

        distanceComputer.__init__(self)

        # dictionnary to store coverage information
        self.input_pacbio = {}



    def __del__(self):
        """__del__: not implemented """
        pass

    # ------------------------------------------------------------------ #
    # Methods                                                            #
    # ------------------------------------------------------------------ #

    # public:

    def compute_distance(self,input_list1,input_list2,ponderation_list=None):
        """compute_distance : to describe """

        # reinitialize the distance table at the begining of each computation
        self.distance_table = {}

        # check that the contigs in the lists have a non empty coverage_dict.
        for item in list(set(input_list1 + input_list2)) :

            if isinstance(item, Contig):
                item_name = item.record.id # to remove

                # if the item dont have a pacbio coverage feature associated yet, fill it.
                if item.pacbio_dict == {}:
                    item.pacbio_dict = self.input_pacbio[item_name] # check if the name correspond exaclty, else correct 

        
        # get the name of the items from input_list1 and initialize the distance table dictionnary, whether it is a bin or a contig
        for item1 in input_list1 :
            if isinstance(item1, Bin):
                item_name1 = 'bin_' + str(item1.number) 
            else : 
                item_name1 = item1.record.id

            self.distance_table[item_name1] = {}

            # go through the items of input_list2 and get their names whether it is a bin or a contig
            for item2 in input_list2 :
                if isinstance(item2, Bin):
                    item_name2 = 'bin_' + str(item2.number)
                else : 
                    item_name2 = item2.record.id
                
                # compute euclidian distance 
                if item_name1 == item_name2 :                
                    self.distance_table[item_name1][item_name2]= 0.0

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
                        for feature in item1.pacbio_dict:
                            self.distance_table[item_name1][item_name2] += (float(item1.pacbio_dict[feature])-float(item2.pacbio_dict[feature]))**2
                        self.distance_table[item_name1][item_name2]=math.sqrt(self.distance_table[item_name1][item_name2])

        return self.distance_table



    def get_feature(self,filename, contig_collection):
        """ get_feature : store the pacbio coverage from the file given in input """

        print('get pacbio coverage feature form ', filename)

        # this function reads a coverage table provided in input and stores it in a double dictionnary self.input_pacbio
        table_reader = tableReaderPrinter()
        self.input_pacbio = table_reader.read_file(filename,'\t')

        # attribute a coverage to every contig 
        for contig_name in contig_collection.dict : 
            
            contig_collection.dict[contig_name].pacbio_dict = self.input_pacbio[contig_name] 

        return 0
     




##########################################################################


if __name__ == '__main__':
    pass

