#!/usr/bin/env python
# -*- coding: utf-8 -*-

" module containing class dotdict "

##########################################################################

##########################################################################


class Dotdict(dict):

    """
dot.notation access to dictionary attributes, found on internet :  

http://stackoverflow.com/questions/2352181/how-to-use-a-dot-to-access-members-of-dictionary

    """



    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    



