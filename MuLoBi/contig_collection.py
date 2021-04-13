#!/usr/bin/env python
# -*- coding: utf-8 -*-

" This module contains the class ContigCollection"

##########################################################################
import operator
##########################################################################


class ContigCollection(object):

    """
class ContigCollection : objects of type ContigCollection contains a collection of contigs (type Contig) stored in a dictionnary with the record.id as key and the contigs as value. The names of the contigs are stored in a list according in increasing length of contig.
    """

    # ------------------------------------------------------------------ #
    # Constructors/Destructors                                           #
    # ------------------------------------------------------------------ #

    def __init__(self):
        """__init__: a contig collection is a list ordonnee of contig names and a dictionnary of contigs (of type record) 
        """

        # Members ---------------------- #

        # list
        self.sorted_name_list = []
	
	# the contig of type Contig is stored in a dictionnary 
        self.dict = {}



    def __del__(self):
        """__del__: not implemented """
        pass



    def addContig(self,new_contig):

        """ addContig : add a new contig to the contig_collection """

        self.dict[new_contig.record.id] = new_contig
            
	# add the name according to the length of the sequence

        return 0


    def getContig(self,contig_id):

        """ getContig : returns the contig corresponding to the contig_id given in argument  """

        return self.dict[contig_id]



    def sortContigsBySize(self):

        """ sortContigsBySize : put the contig names of the contig_collection in order, beginning by the biggest """
        # temportary dict with contig_id as key and length as value
        temp_dict = {}

        for contig_id in self.dict :
            temp_dict[contig_id] = len(self.dict[contig_id].record.seq)

        # sort the dictionnary by value
        sorted_dict = dict(sorted(temp_dict.items(), key=operator.itemgetter(1),reverse=True))

        # uncomment the next line to check the sorting:
        #print('dictionnaire dans l ordre :',sorted_dict)

        # construct the sorted list
        for item in sorted_dict :
            self.sorted_name_list.append(item)

        return 0




##########################################################################


if __name__ == '__main__':
    test = ClusterCollection()
