#!/usr/bin/env python
# -*- coding: utf-8 -*-

" module containing coverageDistanceComputer, daugther of distanceComputer"


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


class coverageDistanceComputer(distanceComputer):

    """
class coverageDistanceComputer: this object computes the distance between two sets (lists) of DNA sequences given an input coverage table
    """

    # ------------------------------------------------------------------ #
    # Constructors/Destructors                                           #
    # ------------------------------------------------------------------ #

    def __init__(self):
        """__init__: the initialization of coverageDistanceComputer is the same as the mother distanceComputer """

        distanceComputer.__init__(self) 

        self.coverage_info_stored = False

        # dictionnary to store coverage information
        self.input_coverage = {}




    def __del__(self):
        """__del__: not implemented """
        pass

    # ------------------------------------------------------------------ #
    # Methods                                                            #
    # ------------------------------------------------------------------ #

    # public:

    def compute_distance(self,input_list1,input_list2,ponderation_list=None):
        """compute_distance : compute the euclidean distances between each pairs of item the first from input_list1, the second in input_list2."""
         
        # reinitialize the distance table at the begining of each computation
        self.distance_table = {}

        # check that all the items in the lists have a non empty coverage_dict.
        for item in list(set(input_list1 + input_list2)) :
        
            # find the name of the item
	    # if the object is a contig :
            if isinstance(item, Contig):
                item_name = item.record.id 

                # if the item dont have a coverage feature associated yet, fill it. If removed, the coverage
                if item.coverage_dict == {}:
                    item.coverage_dict = self.input_coverage[item_name] # check if the name correspond exaclty, else correct 

        # get the name of the items from input_list1 and initialize the distance table dictionnary
        for item1 in input_list1 :

	    # if the object is a bin :
            if isinstance(item1, Bin):
                item_name1 = 'bin_' + str(item1.number)

 
	    # if it is a contig
            else : 
                item_name1 = item1.record.id

            self.distance_table[item_name1] = {}

            # go through the items of input_list2 and get their names whether it is a bin or a contig
            for item2 in input_list2 :
                if isinstance(item2, Bin):
                    item_name2 = 'bin_' + str(item2.number) 
                else : 
                    item_name2 = item2.record.id

                # if it is the same, distance is zero
                if item_name1 == item_name2 :                
                    self.distance_table[item_name1][item_name2]= 0.0

                # compute euclidian distance only if not alredy computed for this pair
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
                        # for coverage A and coverage B, add the square of the difference
                        for feature in item1.coverage_dict:
                            
                            self.distance_table[item_name1][item_name2] += (float(item1.coverage_dict[feature])-float(item2.coverage_dict[feature]))**2
                        self.distance_table[item_name1][item_name2]=math.sqrt(self.distance_table[item_name1][item_name2])

        return self.distance_table




    def get_feature(self,filename,contig_collection):
        """ get_feature : store the coverage from the file given in input """


        # this function reads a coverage table provided in input in has the form of a dictionnary
        table_reader = tableReaderPrinter()
        self.input_coverage = table_reader.read_file(filename,'\t')

        # attribute a coverage to every contig 
        for contig_name in contig_collection.dict : 
            
            contig_collection.dict[contig_name].coverage_dict = self.input_coverage[contig_name] 

        # uncomment to check short-read coverage dictionary
        #print('input coverage : ',self.input_coverage)

        return 0
     


##########################################################################


if __name__ == '__main__':
    pass

