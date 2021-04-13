#!/usr/bin/env python
# -*- coding: utf-8 -*-

" this module contains the class distance computer mother of tetranuDistanceComputer, ..."


##########################################################################

##########################################################################


class distanceComputer(object):

    """
class distanceComputer: general class of distance computer that will compute and store distances in a dictionnary as a matrix.
    """

    # ------------------------------------------------------------------ #
    # Constructors/Destructors                                           #
    # ------------------------------------------------------------------ #

    def __init__(self):
        """__init__: initialize the distance computer with an empty dictionnary"""

        # Members ---------------------- #

	# the distance table of the distance computer will be stored in a dictionnary
        self.distance_table = {}



    def __del__(self):
        """__del__: not implemented """
        pass

    # ------------------------------------------------------------------ #
    # Methods                                                            #
    # ------------------------------------------------------------------ #

    # public:

    def compute_distance(self,input_list1,input_list2,filename = None,ponderation_list=None):
        """compute_distance : interface function """
        pass

    def read(self):
        """read: interface function """
        pass



    def normalize_distance(self,normalization_factor):
        """normalize_disance : normalizes the distance table of the computer with the given factors """
        for item1 in self.distance_table :
            for item2 in self.distance_table[item1]:
                self.distance_table[item1][item2] = self.distance_table[item1][item2]/float(normalization_factor)

        return


##########################################################################


if __name__ == '__main__':
    test = distanceComputer()
