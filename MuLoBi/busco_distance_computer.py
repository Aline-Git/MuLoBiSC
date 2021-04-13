#!/usr/bin/env python
# -*- coding: utf-8 -*-

" module containing buscoDistanceComputer, daugther of distanceComputer"


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


class buscoDistanceComputer(distanceComputer):

    """
class buscoDistanceComputer: this object compute the distance between two sets (lists) of DNA sequences
    """

    # ------------------------------------------------------------------ #
    # Constructors/Destructors                                           #
    # ------------------------------------------------------------------ #

    def __init__(self):
        """__init__: the initialization of coverageDistanceComputer is the same as the mother distanceComputer """

        distanceComputer.__init__(self)
    
        #self.busco_info_stored = False

        # dictionnary to store busco information
        self.input_busco = {}

      



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

        # check that the contigs in the lists have a non empty busco_dict.
        for item in list(set(input_list1 + input_list2)) :
            if isinstance(item, Contig):
                item_name = item.record.id 

                # if the item dont have a busco entry associated yet, fill it.
                if item.busco_dict == {}:
                    item.busco_dict = self.input_busco[item_name] # check if the name correspond exaclty, else correct 
        
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
                
                # compute busco distance between the two items
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
                        
                        intersection = 0
                        union = 0

                    # change busco to numerical
                        for busco in item1.busco_dict :
                            if item1.busco_dict[busco] == 'Complete' or item1.busco_dict[busco] == 'Duplicated' or item1.busco_dict[busco] == 'Undefined':
                                busco_value_item1 = 1
                            elif item1.busco_dict[busco] == 'Fragmented' :
                                busco_value_item1 = 0.5
                            else :
                                busco_value_item1 = 0	    

                            if item2.busco_dict[busco] == 'Complete' or item2.busco_dict[busco] == 'Duplicated' or item2.busco_dict[busco] == 'Undefined':
                                busco_value_item2 = 1
                            elif item2.busco_dict[busco] == 'Fragmented' :
                                busco_value_item2 = 0.5
                            else :
                                busco_value_item2 = 0	    
                   
                            somme = busco_value_item1 + busco_value_item2

                            if somme >= 2 :
                                union += 1
                                intersection += 1
                            elif somme == 1.5 :
                                union += 1
                                intersection += 0.5
                            elif somme == 1 :
                                union +=1

                        
                            # compute the intersection/union
                            try : 
                                self.distance_table[item_name1][item_name2] = float(intersection)/float(union)
                                #print('intersection :', intersection, 'union : ', union)

                            except : 
                                self.distance_table[item_name1][item_name2] = 0
                                #print('intersection :', intersection, 'union : ', union)     

        return self.distance_table



    def get_feature(self,filename, contig_collection):
        """ get_feature : store the coverage from the file given in input """

        print('get busco feature form ', filename)

        # this function reads a coverage table provided in input and stores it in a double dictionnary self.input_pacbio
        table_reader = tableReaderPrinter()
        self.input_busco = table_reader.read_file(filename,'\t')

        # attribute the input_busco to every contig 
        for contig_name in contig_collection.dict : 
            try :
                contig_collection.dict[contig_name].busco_dict = self.input_busco[contig_name] 
            except : 
                pass

        
        return 0
     

    def complete_dictionary(self,contig_collection):
        """ complete the missing entries in the busco dictionary """

        compte_missing = 0
        # compare the key of the input table with the contig collection, add an entry if it is missing
        for contig_name in contig_collection.dict :
            if contig_name not in self.input_busco:
                self.input_busco[contig_name] = {}
                compte_missing += 1
                print('missing busco in : ', contig_name)

                # complete the new entry with 'Missing' values
                for busco in self.input_busco[list(self.input_busco.keys())[0]]:
                    self.input_busco[contig_name][busco] = 'Missing'
     
        print('****\n',compte_missing, ' missing busco entries where completed')
        return 0
        

    def add_busco_tag(self,contig_collection):
        """ complete the missing entries in the busco dictionary """

        is_busco_empty = True

        # add a busco tag to the contig
        for contig_name in contig_collection.dict :
            for busco in self.input_busco[list(self.input_busco.keys())[0]]:
                if self.input_busco[contig_name][busco] == 'Complete' or self.input_busco[contig_name][busco] == 'Duplicated' : 
                    is_busco_empty = False
            if is_busco_empty :
                contig_collection.dict[contig_name].tags.append('no_busco')

        return 0



##########################################################################


if __name__ == '__main__':
    pass

