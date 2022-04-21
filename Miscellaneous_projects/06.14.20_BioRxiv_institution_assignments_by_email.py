#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 12:41:59 2020

@author: adamo010
"""

#poach from 02.14.20_BioRxiv_institution_assignments
#doc source is /Users/adamo010/OneDrive/Documents/Blekhman_lab/Rich BioRxiv paper 2020

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


origlist = pd.read_csv("unknown_email_domains_BA2.csv")
#use version 2 b/c it is clean and does not contain special characters, which don't work with CSVs

#copy the domain column to preserve the original
origlist["domain_original"] = origlist["domain"]
origlist["domain"] = origlist["domain"].str.split(".", expand=False)      
#this splits the "domain" column (a list of email addresses) in orglist2 into a list of strings, where elements are 
#defined by where the periods are in the original string. Expand is set to False so a list of strings
#will be returned instead of a dataframe

#I'm not sure if this is the best option, but here we are.

#the last element of the list contains information on country; need to use a dictionary of key-value pairs with known values and
#their respective countries to scan through the dataframe and auto-input this information.

#now, need to find some way to create a new column with only the last element in each list in the column

#Maybe I should get a list of all endings, and figure out which country codes I need to look up

#only keep rows where 'domain' is not NAN (I think this is all but the first row)
origlist2 = copy.deepcopy(origlist)
origlist2 = origlist2[origlist2['domain'].notna()]

country_email = []

for row in origlist2["domain"]:   
    #print(type(value))
    print(row[-1])
    country_email.append(row[-1])
del row

unique_country_emails = []
for elem in country_email:
    if elem not in unique_country_emails:
        unique_country_emails.append(elem)
del elem        
        
with open('email_ends.txt', 'w') as f:
    for item in unique_country_emails:
        f.write("%s\n" % item)    
del f
del item
        
#I cheated and used vlookup in Excel to assign IANA (Internet Assigned Numbers Authority) to the endings of email addresses
email_iana_codes = pd.read_csv("iana_email_assignments.csv")    

#create a dictionary from email_iana_codes where keys are the country codes from emails, and values are from the country names
with open('iana_email_assignments.csv', mode='r') as infile:
    reader = csv.reader(infile)
    with open('new_iana_email_assignments.csv', mode='w') as outfile:
        writer = csv.writer(outfile)
        domain_dict = {rows[0]:rows[1] for rows in reader}
del infile
del outfile
del reader
del writer

#create a new column with last element from 'domain' column lists (which is the country code)
origlist2["country_domain"] = country_email

#create a new column, country_from_domain, that contains the keys from new_dict whose values correspond to country_domain column
origlist2['country_from_domain'] = origlist2['country_domain'].map(domain_dict)

#dynamite. very handy.

origlist2.to_csv("updated_unknown_email_domains_BA.csv")

########################################################I HATE THIS FUCK IT
#use https://github.com/hipo/university-domains-list
#this has a list of all domains and names of universities and countries in the world, apparently
#it's a .json file and I don't know what that is or how to handle it.

import json
#read file
with open('university-domains-list-master/world_universities_and_domains.json', 'r') as myfile:
    master_domain_database=myfile.read()
#parse file
master_domain_database_obj = json.loads(master_domain_database)

#well, this is really a challenge. Each line seems to be a dictionary
for line in master_domain_database_obj:
    print(line)

#let's try it as a pandas series
master_domain_df = pd.read_json(r'university-domains-list-master/world_universities_and_domains.json')    
#hell. yes. Screw .json files  

master_domain_df.columns   
#columns are 'web_pages', 'name', 'alpha_two_code', 'state-province', 'domains', 'country'    

master_domain_df_slimmed = master_domain_df.drop(['web_pages', 'alpha_two_code', 'state-province'], axis = 1)    
    
#worth noting that domains is a list, b/c multiple domains may be associated with a single university
#yeah okay, this is going to be a problem
#need to duplicate rows where domains have multiple entries and have a separate row for each domain

master_domain_df_slim2= master_domain_df_slimmed.explode('domains')
#huh.
#only added 80ish new rows, which is good. 

#master_domain_df_slimmed.domains.apply(pd.Series).merge(master_domain_df_slimmed, left_index = True, right_index=True)

for elem in master_domain_df_slim2['domains']:
    print(type(elem))
    
#so yesterday I did a merge on master_domain_df_slim2 and origlist2: this was frustrating and stuff happened that I didn't understand.
#now, I'm going to make some dictionaries and apply them
    
#columns in origlist2:
origlist2.columns    
#domain', 'count', 'institution', 'country', 'notes', 'domain_original', 'country_domain', 'country_from_domain'
master_domain_df_slim2.columns
#'name', 'domains', 'country'

#domain is shared. Need to make a couple of dictionaries:
#create a new column, country_from_domain, that contains the keys from new_dict whose values correspond to country_domain column
#key is the domain, value is the country OR the institution name
github_country_dict = pd.Series(master_domain_df_slim2.country.values,index=master_domain_df_slim2.domains).to_dict()
github_institution_dict = pd.Series(master_domain_df_slim2.name.values,index=master_domain_df_slim2.domains).to_dict()

origlist2['country_from_github'] = origlist2['domain_original'].map(github_country_dict)
origlist2['institution_from_github'] = origlist2['domain_original'].map(github_institution_dict)

#okay, now I need to merge info from the two countries columns: country_from_domain and country_from_github. In cases where they disagree,
#I'd llike to keep the value from country_from_domain

origlist3 = copy.deepcopy(origlist2)

# replacing na values in college with No_country because NA VALUES FUCK EVERYTHING UP
origlist3["country_from_domain"].fillna("no_country", inplace = True) 
origlist3["country_from_github"].fillna("no_country", inplace = True) 

#maybe just... add them both together in a list then only keep the first element?
origlist3['country_from_both'] = origlist3[["country_from_domain", "country_from_github"]].values.tolist()

combined_country = []

for row in origlist3['country_from_both']:
    if row[0] != "no_country":
        combined_country.append(row[0])  
    elif row[0] == "no_country" and row[1] != "no_country":
        combined_country.append(row[1])
    elif row[0] == "no_country" and row[1] == "no_country":    
        combined_country.append('none')

origlist3 ['final_country']= combined_country            

final_list = origlist3[['domain_original', 'count', 'institution_from_github','final_country', 'notes']].copy()
final_list = final_list.rename({'domain_original': 'domain', 'institution_from_github': 'institution', 'final_country': "country"}, axis =1)

final_list.to_csv("unknown_email_domains_autofilled.csv", encoding='utf-8-sig')

final_list.to_csv("unknown_email_domains_autofilled.txt", header=True, index=False, sep='\t', mode='a')
a.to_csv('xgboost.txt', header=True, index=False, sep='\t', mode='a')



#############################JUNK##############################JUNK##############################JUNK##############################    

#now need a clean, assignment-free version of the BioRxiv email lists
domains_to_assign = pd.read_csv("unknown_email_domains_BA2.csv")
#use version 2 b/c it is clean and does not contain special characters, which don't work with CSVs
domains_to_assign.columns   
#columns are 'domain', 'count', 'institution', 'country', 'notes'

#I would like to look up the domain from domains_to_assign in master_domain_df_slimmed and return the name (to institution)
#and country (to country) from master_domain_df_slimmed to domains_to_assign

master_domain_df_slim2= master_domain_df_slim2.rename(columns={"domains": "domain"})

domains_to_assign['institution_assigned'] = domains_to_assign.domain.map(master_domain_df_slim2.set_index('domain').institution)    

results = domains_to_assign.merge(master_domain_df_slim2, on='domain', how="right")    
#that... kind of worked. only kept ~450 rows???
results2 = master_domain_df_slim2.merge(domains_to_assign, on='domain', how="right")    
#same thing. Clearly, this is just keeping the values that are present in both lists
#adding the 'how' statement... is not that helpful.
#it makes results2 1620 rows, while domains_to_assign is 1611.
#I do not understand. 

results.to_csv("github_domain_based_ids.csv")
#this identified 449 of 1620 unknown domains. Just going by country abbreviations got 1005. This 1005 does not include any US institutions.


    