#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 09:45:29 2020

@author: adamo010
"""

import numpy as np
import csv
import subprocess
import os
import pandas as pd
import fnmatch
import shutil
import glob
import re
import copy


origlist = pd.read_csv("need_institution_assignments.csv")

origlist2 = origlist.applymap(str)
#don't know if this does anything

origlist2["affiliation string"] = origlist2["affiliation string"].str.split(",", expand=False)      
#this splits the "affiliation string" column in orglist2 into a list of strings, where elements are 
#defined by where the commas are in the original string. Expand is set to False so a list of strings
#will be returned instead of a dataframe

#now, need to find some way to create a new column with only the last element in each list in the column

end_aff_string = []

for value in origlist2["affiliation string"]:
    print(type(value))
    end_aff_string.append(value[-1])

origlist2["end_aff_string"] = end_aff_string    

###########VICTORY!!!! Thanks Sambhawa for the idea to just pull the last element from the string

origlist2.to_csv("need_inst_assignment_last_string.csv")


##############Junk
#for row in origlist['affiliation string']:
    #print(row)
    #row.str.split(",")
#error: AttributeError: 'str' object has no attribute 'str'