#!/usr/bin/env python
# -*- coding: utf-8 -*-

" module containing class Contig"


class Contig(object):

    """
class Contig: an object of type SeqRecord plus aditionnal features. It is defined by a number and initialized with a contig.
    """

    # ------------------------------------------------------------------ #
    # Constructors/Destructors                                           #
    # ------------------------------------------------------------------ #

    def __init__(self,record):
        """__init__: a contig is initialized with a record. Divers information will be stored in the object Contig """

        # Members ---------------------- #

        # SeqRecord object
        self.record = record

        # typical 4mer frequency vector of the contig
        self.d_4mers_median = {}

        # 4mer frequency table on the sliding window
        self.d_4mers_matrix = {}      

        # coverage vector stores the coverage of the illumina extractions A and B 
        self.coverage_dict = {}

        # coverage vector stores the coverage of the pacbio extractions
        self.pacbio_dict = {}

        # busco_dictionnary stores the busco features with the busco name as key
        self.busco_dict = {}

        # a list to stores the different tags of the contigs
        self.tags = []




    def __del__(self):
        """__del__: not implemented """
        pass






##########################################################################


if __name__ == '__main__':
    test = Contig()
