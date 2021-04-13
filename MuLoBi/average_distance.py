#!/usr/bin/env python
# -*- coding: utf-8 -*-

" This module contains the class averageDistance"


##########################################################################
from our_bin import Bin
from distance_computer import distanceComputer
from tetranu_distance_computer_v4 import tetranuDistanceComputer
from coverage_distance_computer import coverageDistanceComputer
from pacbio_distance_computer import pacbioDistanceComputer
from busco_distance_computer import buscoDistanceComputer
import json
##########################################################################

class averageDistance(object):

    """
class averageDistance: launches the four different distance computers, have the functions to compute and normalize the different distance matrix and compute the average matrix.
    """

    # ------------------------------------------------------------------ #
    # Constructors/Destructors                                           #
    # ------------------------------------------------------------------ #

    def __init__(self,coverage_illumina_table,coverage_pacbio_table,busco_table,ponderation_list,contig_collection):
        """__init__: the initialization of  tetranuDistanceComputer is the same as the mother distanceComputer"""

        # initialize a tetranu_distance_computer
        self.tetranu_distance_computer = tetranuDistanceComputer()

        # initialize a coverage_distance_computer and store features from the illumina coverage input file size :nx2
        self.coverage_distance_computer = coverageDistanceComputer()
        print('\nget feature from illumina coverage table')
        self.coverage_distance_computer.get_feature(coverage_illumina_table,contig_collection)

        # initialize a coverage_distance_computer and store features from the illumina coverage input file size :n
        self.pacbio_distance_computer = pacbioDistanceComputer()
        print('\nget feature from pacbio coverage table')
        self.pacbio_distance_computer.get_feature(coverage_pacbio_table,contig_collection)

        # initiate a busco_distance_computer and store features from the busco input file, size nx number of busco
        self.busco_distance_computer = buscoDistanceComputer()
        print('\nget feature from busco table')
        self.busco_distance_computer.get_feature(busco_table,contig_collection)
        self.busco_distance_computer.complete_dictionary(contig_collection)
      
        # create an empty average distance tables and stores the ponderation list given in init
        self.average_distance_table = {}
        self.ponderation_list = ponderation_list


       
    def __del__(self):
        """__del__: not implemented """
        pass

    # ------------------------------------------------------------------ #
    # Methods                                                            #
    # ------------------------------------------------------------------ #

    # public:



    def compute_distance(self,item_list1,item_list2):
        """compute_distance : call the four distance computers to compute the four distance table between each item of the two lists, then compute the average matrix form these distance tables and the poderation list """

	# compute for each distance a distance_table which is a double dictionnary containing the distances between the items of list1 and list2 as values and the item names as keys

        print('\ncompute tetranucleotide distance')
        tetranu_distance = self.tetranu_distance_computer.compute_distance(item_list1,item_list2)


        print('\ncompute coverage distance')
        coverage_distance = self.coverage_distance_computer.compute_distance(item_list1,item_list2)

        print('\ncompute pacbio distance')
        pacbio_distance = self.pacbio_distance_computer.compute_distance(item_list1,item_list2)

        print('\ncompute busco distance')
        busco_distance = self.busco_distance_computer.compute_distance(item_list1,item_list2)
       

        # compute the average distance for each pair
        for item1 in item_list1 :

	# get the name of the object if the object is a bin or a contig :
            if isinstance(item1, Bin):
                item_name1 = 'bin_' + str(item1.number)
            else : 
                item_name1= item1.record.id

            # create an entry for item1 and go through item from item_list2
            self.average_distance_table[item_name1] = {}
            for item2 in item_list2 :

	        # get the name of the object if the object is a bin or a contig :
                if isinstance(item2, Bin):
                     item_name2 = 'bin_' + str(item2.number)
                else : 
                    item_name2 = item2.record.id
                
                # if item1 = item2 the distance is 0
                if item_name1 == item_name2 :
                    self.average_distance_table[item_name1][item_name2]= 0.0
          
                # else create an entry only if it does not exist for the pair
 
                else : 

                    # check if it is already computed or not
                    distance_computed = False
       
                    # if the reverse entry already exist then distance_computed is set to true
                    try : 
                        self.average_distance_table[item_name2][item_name1]
                        distance_computed = True

                    except : 
                        distance_computed = False

                    # if the distance was already computed, we dont add another entry 
                    if not distance_computed :

                        # compute the average entry in the average_distance_table for this pair of items
                        self.average_distance_table[item_name1][item_name2] = float(self.ponderation_list[0])*self.tetranu_distance_computer.distance_table[item_name1][item_name2] + float(self.ponderation_list[1])*self.coverage_distance_computer.distance_table[item_name1][item_name2] + float(self.ponderation_list[2])*self.pacbio_distance_computer.distance_table[item_name1][item_name2] + float(self.ponderation_list[3])*self.busco_distance_computer.distance_table[item_name1][item_name2]

        return self.average_distance_table



    def normalize_distance(self,normalization_factor_list):
        """normalize_distance : compute the weighted average distance table from the four normalized distances and their ponderations."""
        print('\nnormalize distances')
        self.tetranu_distance_computer.normalize_distance(normalization_factor_list[0])
        self.coverage_distance_computer.normalize_distance(normalization_factor_list[1])
        self.pacbio_distance_computer.normalize_distance(normalization_factor_list[2]) 

        # recompute the average distance

        #get item_name_list : 
        item_list1 = list(self.coverage_distance_computer.distance_table.keys())

        # initialize an entry for item_list 1 and get the list of their corresponding entries
        for item1 in item_list1 :

            self.average_distance_table[item1] = {}
            item_list2 = list(self.coverage_distance_computer.distance_table[item1].keys())
          
            for item2 in item_list2 :

                self.average_distance_table[item1][item2] = float(self.ponderation_list[0])*self.tetranu_distance_computer.distance_table[item1][item2] + float(self.ponderation_list[1])*self.coverage_distance_computer.distance_table[item1][item2] + float(self.ponderation_list[2])*self.pacbio_distance_computer.distance_table[item1][item2] + float(self.ponderation_list[3])*self.busco_distance_computer.distance_table[item1][item2]

        return self.average_distance_table



if __name__ == '__main__':
    pass






