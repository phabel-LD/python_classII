'''
## NAME:
  Entrez.py
  
## LANGUAGE & VERSION:
  python 3.8.5
  
## AUTHOR:
  Phabel Antonio Lopez Delgado <phabel2001@gmail.com>
  
## DATE:
  October 10th, 2021.
  
## DESCRIPTION & LOGIC:
  This script solves two tasks. The first one is to access
  two specific fiels and subfields from the database "protein"
  and print their description. The second one automatizes an
  esearch with a string of key terms, retrieves files' IDs and
  writes a new file with them. Then uses the id's in the file
  to access the abstracts and citations from PubMed and write
  a second file.
  
## USAGE:
  Entrez.py python 3.8.5


## INPUT - OUTPUT:

    INPUT: 
        Reference email
        search_terms
        

    OUTPUT: A file with all IDs
        according to search_terms.
        A file with all abstracts
        and citations per ID


## SOFTWARE REQUIREMENTS:
    python3
    BioPython  -->  Library Entrez

## EXTRA COMMENTS:
   Remeber to always add the reference email.
    
    
## LAST MODIFICATION:
  Phabel Antonio López Delgado: October 02nd, 2021  [Creation]
  Phabel Antonio López Delgado: October 10th, 2021  [Upgrade to Part II]

## SOURCE:
  GitHub: https://github.com/phabel-LD/python_classII/blob/master/Tasks/Entrez.py
  
'''
# Check for permissions
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

# Library
from Bio import Entrez, SeqIO


# Main Code ############

# Part I ###############

# Task 1: Print fileds' description
# email reference
Entrez.email = "phabel@lcg.unam.mx"
# Create handle and record with info from 'protein' DB
handle = Entrez.einfo(db = "protein")
record = Entrez.read(handle)
handle.close()

# Print description for FieldList "ECNO"
for field in record["DbInfo"]["FieldList"]:
    if field["Name"]  == "ECNO":
        print(field["Description"])

# Print description for LinkList "protein_protein_small_genome"
for field in record["DbInfo"]["LinkList"]:
    if field["Name"]  == "protein_protein_small_genome":
        print(field["Description"])
        

# Task 2: Use esearch to write IDs file
# Create term
author = "Jie Su"
title_word1 = "microRNA"
title_word2 = "mice"
search_terms = f"{author} AND ({title_word1} OR {title_word2})"
# Access info with Entrez and esearch
count = int(record["DbInfo"]["Count"])
handle = Entrez.esearch(db = "pubmed", term = search_terms, retmax = count)
record = Entrez.read(handle)
handle.close()
IDs = record["IdList"]
print(f"IDs({len(IDs)}): {IDs}")
# Create file with all ID's
with open("IDs.txt", "w") as file_ids:
    for id_name in IDs:
        file_ids.write(f"{id_name}\n")

# Part II ###############
# Open abst_file in "w" mode, and IDs.txt just with "r+" mode
with open("Abstracts.txt", "w") as abst_file:
    with open("IDs.txt", "r+") as file_ids:
        # Get the info (abstract and citations) for every id in the file
        for id_name in file_ids:
            # Use efetch to retrieve abstracts 
            fetch_handle = Entrez.efetch(db="pubmed", id = id_name,
                                rettype="abstract", retmode="text")
            data = fetch_handle.read()
            fetch_handle.close()
            abst_file.write(data)
            
            # Get citations
            citations_results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc",
                                                 LinkName = "pubmed_pmc_refs", from_uid = id_name))
            # Write every citations in abst_file
            for citation in citations_results[0]["LinkSetDb"]:
                abst_file.write(str(citations_results[0]["LinkSetDb"]))
