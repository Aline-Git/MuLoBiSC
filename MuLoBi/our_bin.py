#!/usr/bin/env python
# -*- coding: utf-8 -*-

" module containing class Bin"

##########################################################################
import operator
##########################################################################

class Bin(object):

    """
class Bin: a bin is a list of contigs (objects of type SeqRecord). It is defined by a number and initialized with a contig.
    """

    # ------------------------------------------------------------------ #
    # Constructors/Destructors                                           #
    # ------------------------------------------------------------------ #

    def __init__(self,number,contig):
        """__init__: a bin is initialized with a first contig and the number of the bin. It also contains dictionnaries to store the informations about the 4_mers statistics, coverage illumina, coverage pacbio and busco. """

        # Members ---------------------- #

        # number of the bin
        self.number = number

        # number of the bin
        self.length = 0

        # list of the contig of the bins
        self.contig_list = []

	# attribute the contig to the bin
        self.contig_list.append(contig)

        # typical 4mer frequency vector of the bin (is reinitialized each time a contig is added)
        self.d_4mers_median = {}

        # 4mer frequency table on the sliding window over the bin
        self.d_4mers_matrix = {}  

        # coverage vector stores the coverage of the illumina extractions A and B  
        self.coverage_dict = {}

        # coverage vector stores the coverage of the pacbio reads
        self.pacbio_dict = {}

        # busco_dictionnary stores the busco features with the busco name as key
        self.busco_dict = {}
    


    def __del__(self):
        """__del__: not implemented """
        pass



    def evaluate_N50_contigs(self):
        """ evaluate_N50_contigs :  tags the contigs of the bin with 'bigger_N50' if their size is bigger than the N50 of the bin at the end of the binning 
        """     
  
        # store the contig names and their size in a temporary dictionary temp_dict
        temp_dict = {}
        for contig in self.contig_list :
            temp_dict[contig.record.id] = len(contig.record.seq)

        # sort the temp dictionnary by value
        sorted_dict = dict(sorted(temp_dict.items(), key=operator.itemgetter(1),reverse=True))
        # uncomment the next line to check the sorting
        #print('sorted_dict',sorted_dict)

        # set the sum of the len of contigs to zero and N50_depasse to False
        sum_len_contigs = 0
        N50_depasse = False

        # go through the contig names in the dictionnary sorted by size
        for item in sorted_dict :
            sum_len_contigs += temp_dict[item]

            # if N50_depasse is False, search for the contig corresponding to item and tag it with 'bigger_N50'
            if not N50_depasse :
                for contig in self.contig_list :
                    if contig.record.id == item :
                        contig.tags.append('bigger_N50')
                        break

            # if the sum is bigger than half of the length of the contig set N50_depasse to True
            if sum_len_contigs > 0.5 * self.length :
                N50_depasse = True

        return 0





##########################################################################


if __name__ == '__main__':
    test = Sequence()
