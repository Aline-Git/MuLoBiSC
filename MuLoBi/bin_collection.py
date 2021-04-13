#!/usr/bin/env python
# -*- coding: utf-8 -*-

" This module contains the class BinCollection "

##########################################################################
from our_bin import Bin
from table_reader_printer import tableReaderPrinter
##########################################################################


class BinCollection(object):

    """
class BinCollection: the bin collection is a dictionnary of Bins, with a number as key. It will evolved during the binning.
    """

    # ------------------------------------------------------------------ #
    # Constructors/Destructors                                           #
    # ------------------------------------------------------------------ #

    def __init__(self):
        """__init__: a bin collection is initialized as an empty dictionnary """

        # Members ---------------------- #

        # dictionary of bins
        self.dict = {}
        
	# current_number  
        self.current_number = 0
        self.current_number = int(self.current_number)


    def __del__(self):
        """__del__: not implemented"""
        pass


    def create_new_bin(self,contig):

        """ __create_new_bin__ : add a bin to the bin collection with a contig and a number as a name (the number is also contained in the bin. """

	# create a new bin having as number the current number and add it to the present collection 
        self.dict[self.current_number]=Bin(self.current_number,contig)

        # at the begining, the attributes of the bin are those of the only contig member
        self.dict[self.current_number].d_4mers_median = contig.d_4mers_median
        self.dict[self.current_number].d_4mers_matrix = contig.d_4mers_matrix
        self.dict[self.current_number].coverage_dict = contig.coverage_dict
        self.dict[self.current_number].pacbio_dict = contig.pacbio_dict
        self.dict[self.current_number].busco_dict = contig.busco_dict

        # with one contig in the bin, the length is the one of hte contig
        self.dict[self.current_number].length = len(contig.record.seq)

        # let this at the end of the function !!!
        self.current_number += 1
        print('new bin init with contig : ',contig.record.id)

        return 0   



    def add_contig_to_bin(self, contig, bin_number):

        """__add_contig_to_bin__: add the contig contig_name to the bin bin_name and actualize all the features of the bin """

        # append the contig to the contig_list of the bin number bin_number
        bin_number = int(bin_number)
        self.dict[bin_number].contig_list.append(contig)

        # the 4-mers matrix is actualized by combining the contig matrix to the bin matrix

        self.addTwoFrequencyMatrices(contig,bin_number)


        # the 4-mer median frequency vector is reinitialized
        self.dict[bin_number].d_4mers_median = {}

        # actualise the bin 'bin_number' illumina coverage_dict
        for feature in  self.dict[bin_number].coverage_dict:
            self.dict[bin_number].coverage_dict[feature] = (float(contig.coverage_dict[feature]) * len(contig.record.seq) + float(self.dict[bin_number].coverage_dict[feature]) * self.dict[bin_number].length) / (len(contig.record.seq) + self.dict[bin_number].length)


        # actualise the bin 'bin_number' pacbio_dict
        for feature in  self.dict[bin_number].pacbio_dict:
            self.dict[bin_number].pacbio_dict[feature] = (float(contig.pacbio_dict[feature]) * len(contig.record.seq) + float(self.dict[bin_number].pacbio_dict[feature]) * self.dict[bin_number].length) / (len(contig.record.seq) + self.dict[bin_number].length)


        # actualise the bin 'bin_number' busco_dict
        for busco in self.dict[bin_number].busco_dict:
            self.dict[bin_number].busco_dict[busco] = self.merge_busco(contig.busco_dict[busco], self.dict[bin_number].busco_dict[busco])

        # actualise bin length
        self.dict[bin_number].length = self.dict[bin_number].length + len(contig.record.seq) 


        return 0



    def addTwoFrequencyMatrices(self,contig,bin_number):

        #contig.d_4mers_matrix
        print('add two matrices')
        print('contig : ',contig.record.id, '   bin : ', bin_number)

        # find the maximum value in self.dict[bin_number].contig.d_4mers_matrix
        maximum_window = max(list(self.dict[int(bin_number)].d_4mers_matrix['TTTC'].keys()))
        #print('maximum_window',maximum_window)
        # add this value to the second key in a copy named temp_dict
        temp_dict = {}
        for oligonu in contig.d_4mers_matrix :
            temp_dict[oligonu] = {}

            # find the resolution (if the contig is as long as the sliding window, it will create an error: 
            try :
                resolution = int(list(contig.d_4mers_matrix[oligonu])[1])- int(list(contig.d_4mers_matrix[oligonu])[0])
            except : 
                resolution = 200

            for start_pos in contig.d_4mers_matrix[oligonu]:
                temp_dict[oligonu][start_pos + int(maximum_window)+resolution] = contig.d_4mers_matrix[oligonu][start_pos]
   
            # add the values of the contig dictionnary in the bin dictionnary
            self.dict[int(bin_number)].d_4mers_matrix[oligonu].update(temp_dict[oligonu])

        return 0



    def merge_busco(self,busco1,busco2):
        busco = ''
        if  busco1 =='Duplicated' or busco2 =='Duplicated':
            busco ='Duplicated'

        elif busco1 == 'Missing':
            busco = busco2

        elif busco2 == 'Missing' :
            busco = busco1

        elif busco1 == 'Fragmented':
            if busco2 == 'Fragmented' or busco2 == 'Complete' : 
                busco = 'Undefined'
            elif busco2 == 'Undefined':
                busco = 'Duplicated'
            else :
                print('exception1')

        elif busco2 == 'Fragmented':
            if busco1 == 'Complete' : 
                busco = 'Undefined'
            elif busco1 == 'Undefined':
                busco = 'Duplicated'
            else :
                print('exception2')
                print('busco1 : ', busco1 , ' busco2 : ',busco2)

        elif busco1 == 'Undefined':
            if busco2 =='Undefined' or busco2 == 'Complete' : 
                busco = 'Duplicated'
            else :
                print('exception3')

        elif busco1 == 'Complete':
            if busco2 =='Undefined' or busco2 == 'Complete' : 
                busco = 'Duplicated'
		    
            else :
                print('exception4')

        else :
            print('problem to determinate busco ')
            print('busco1 : ', busco1 , ' busco2 : ',busco2)
            busco = 'NA'

        return busco



          

##########################################################################


if __name__ == '__main__':
   # test = Cluster('test_cluster_head_id', 'ATTCAGTAAT')
   pass

