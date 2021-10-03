'''
## NAME:
  Entrez.py
  
## LANGUAGE & VERSION:
  python 3.8.5
  
## AUTHOR:
  Phabel Antonio Lopez Delgado <phabel2001@gmail.com>
  
## DATE:
  October 02nd, 2021.
  
## DESCRIPTION & LOGIC:
  This script solves two tasks. The first one is to access
  two specific fiels and subfields from the database "protein"
  and print their description. The second one automatizes an
  esearch with a string of key terms, retrieves files' IDs and
  writes a new file with them.
  
## USAGE:
  Entrez.py python 3.8.5


## INPUT - OUTPUT:

    INPUT: 
        Reference email
        search_terms
        

    OUTPUT: A file with all IDs
        according to search_terms


## SOFTWARE REQUIREMENTS:
    python3
    BioPython  -->  Library Entrez

## EXTRA COMMENTS:
   Remeber to always add the reference email.
    
    
## LAST MODIFICATION:
  Phabel Antonio Lopez Delgado: October 02nd, 2021  [Creation]

## SOURCE:
  GitHub: https://github.com/phabel-LD/python_classII/blob/master/Tasks/Entrez.py
  
'''

# Library
from Bio import Entrez


# Main Code

# Task 1: Print fileds' description
# email reference
Entrez.email = "phabel@lcg.unam.mx"
# Create handle and record with info from 'protein' DB
handle = Entrez.einfo(db = "protein")
record = Entrez.read(handle)
handle.close()
print(record["DbInfo"]["FieldList"]["ECNO"])
print(record["DbInfo"]["LinkList"]["protein_protein_small_genome"])


# Task 2: Use esearch to write IDs file
# Create term
search_terms = "(Amaranta Manrique[Author]) AND (alacranes[Title] OR etica[Title])"
# Access info with Entrez and esearch
handle = Entrez.esearch(term = search_terms)
record = Entrez.read(handle)
handle.close()
IDs = record["IdList"]
print(IDs)
# Create file with all ID's
file_name = "IDs.txt"
with open(file_name, "w") as file:
    for id_name in IDs:
        file_name.write(id_name)

