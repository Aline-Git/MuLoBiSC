#!/usr/bin/env python
# -*- coding: utf-8 -*-

" This module contains the class tablePrinter"


##########################################################################

##########################################################################


class tableReaderPrinter(object):

    """
class tableReaderPrinter: this type of object is implemented to transfer a table from a file to a double dictionnary and conversly
    """

    # ------------------------------------------------------------------ #
    # Constructors/Destructors                                           #
    # ------------------------------------------------------------------ #

    def __init__(self):
        """__init__: initialize the distance computer given a table which can be empty """

        # Members ---------------------- #

        self.filename = ''


    def __del__(self):
        """__del__: not implemented """
        print('\nDeleting table printer')

    # ------------------------------------------------------------------ #
    # Methods                                                            #
    # ------------------------------------------------------------------ #

    # public:
    def read_file(self,filename, separator):
        """read_file: read a table from a file and return it as a dictionnary of dictionnary """

        print('read from file : ',filename)

        # open the file given in input
        input_file = open(filename,'r')
        table = {}
       
        # store the entries of the header
        header = input_file.readline().replace('\n','').split(separator)

        # go through the entries of the input file
        line = input_file.readline()
        while line :
            split_line = line.replace('\n','').split(separator)
            #print('split line : ', split_line)

            # the name of the contig is put as a key of the table
            item_name = split_line[0]
            table[item_name] = {}
         
            # the corresponding values fro this contig are stored in the table
            for i in range(1,len(split_line)):  
                table[item_name][header[i]]= split_line[i]

            line = input_file.readline()
        input_file.close()

        return table



    def output_in_file(self,filename, table):
        """output_in_file: print a table which is given as a double dictionary, with the first keys as header and the second keys in the first column """
        
        output_file = open(filename,'w')
        
        # go through the keys of the table and print them as header
        header = ' '
        nb_item2 = 0

        for item2 in table[list(table.keys())[0]] :
            header += '\t' + str(item2)
            nb_item2 += 1
        output_file.write(header + '\n')

        # print the elements of the table :
        for item1 in table :
            output_line = item1
            print('output line for item1 : ',item1)
            # add missing values so that the matrix is not half full
            for item in table :
                print('item : ',item)
                if item not in table[item1] and item != item1:
                    output_line += '\t' + str(table[item][item1])

            for item2 in table[item1] :
                print('item2 : ',item2)

                output_line += '\t' + str(table[item1][item2])
            output_file.write(output_line + '\n')
                
        output_file.close()


        return 0



##########################################################################


if __name__ == '__main__':
    test = tablePrinter()
