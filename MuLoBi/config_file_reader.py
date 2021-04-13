#!/usr/bin/env python
# -*- coding: utf-8 -*-

" module config_file_reader contains the definition of the class ConfigFileReader"

##########################################################################
from dotdict import Dotdict
##########################################################################


class ConfigFileReader(object):

    """
class ConfigFileReader: can read and parse options given in a configuration file
    """

    # ------------------------------------------------------------------ #
    # Constructors/Destructors                                           #
    # ------------------------------------------------------------------ #

    def __init__(self,filename):
        """__init__: the config file reader is initialized with 4 members : a filename and the 2 dictionnaries where the information will be stored plus the dictionnary of arguments to return """

        # Members ---------------------- #

	# string filename
        self.filename = filename

	# dictionnary dict
	# self.dict[option_name] = [option_value]
        self.dict = {}
	
	# dictionnary conversion_dict
        self.conversion_dict = {}

	# dictionnary arguments
        self.args = {}

	
    def __del__(self):
        """__del__: not implemented """
        pass

    # ------------------------------------------------------------------ #
    # Methods                                                            #
    # ------------------------------------------------------------------ #

    # public:

    def add_argument(self, arg_short, arg_long, help, required):
        """ add the arguments that can be passed to the program by the configuration file. """
	
	# the argument names and the usage is stored in the self.dict dictionnary
        self.dict[arg_long.replace('--','')] = [help, required]
	# the correspondance between the short and the long version of the argument name is stored in the self.conversion_dict
        self.conversion_dict[arg_long.replace('--','')] = arg_short.replace('-','')

        return 0 



    def read(self):
        """read: this function read the information in the config file and store it in the self.dict. In the future the error message can be enhanced by providing the list of valid options for example et eventuellement sortir du programme"""
        
        config_file = open(self.filename,'r')
	
	# read the configuration file line by line from the start
        config_line = config_file.readline()
	
        while config_line:
	    # cont will switch to false if an option is not valid, then the option is not stored
            cont = True

	    # short and long argument names are initialized
            arg_long = ''
            arg_short = ''

	    # the argument name is extracted from the current line of the configuration file
            argument = config_line.split('\t')[0]

	    # if it is the long name 
            if argument[:2]=='--': 
                arg_long = argument.replace('--','')

		# check if the option name is valid

		# if not print an error message
                if arg_long not in self.dict : 
                    print('error : ', arg_long, ' is not a valid option.')

		# if yes the short option name is extracted from the conversion dictionnary
                else : 
                    arg_short = self.conversion_dict[arg_long]

	    # same if the option name is the short version
            elif argument[:1] == '-':
                arg_short = argument.replace('-','')
                if arg_short not in self.conversion_dict.values():
                    print('error : ', arg_short, ' is not a valid option.')
                    print('check there are only tabs in config files and not some spaces')
                else :
                    for key in self.conversion_dict:
                        if self.conversion_dict[key] == arg_short:
                            arg_long = key

	    # if it is not valid, print an error message and set cont to False
            else :
                print('error: ', argument, ' is not a valid option.')
                cont = False

	    # if the option is valid
            if cont :
		# get the option value
                arg_value = config_line.split('\t')[1].replace('\n','')
		# store it in the argument dictionnary
                self.args[arg_long] = arg_value


            config_line = config_file.readline()
	   
        config_file.close() 

        return 0



    def parse_args(self):
        """return a dotdictonnary with arguments values : dict.option_name = option_value"""

	# allow acces to the values with dot (to be compatible with the ArgumentParser)
        self.args = Dotdict(self.args)


        return self.args





##########################################################################


if __name__ == '__main__':
    test = ConfigFileReader('/home/aline/spe_repository_aa/project/class_diagram/input_files/config_clustering_cd_hit.txt')

