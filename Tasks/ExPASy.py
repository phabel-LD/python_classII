'''
## NAME:
  ExPASy.py
  
## LANGUAGE & VERSION:
  python 3.8.5
  
## AUTHOR:
  Phabel Antonio Lopez Delgado <phabel2001@gmail.com>
  
## DATE:
  October, 2021.
  
## DESCRIPTION & LOGIC:
  This script uses BioPython tools ExPASy, SwissProt and
  Entrez in order to analyse SwissProt IDs. The main aim is
  to receive an ID list and GO terms list, retrieve
  SwissProt files for every ID, and write an analysis file
  if the GO terms correspond to the ID.

## USAGE:
  ExPASy.py python 3.8.5
  
## ARGUMENTS:
  The script receives no args from Terminal.

## INPUT - OUTPUT:
   Input: List of IDs and list of GO terms
   Output: Analysis files .txt for every valid ID.
       
## EXAMPLES:
   Input(Inscript):
       id_list = ["O87765", "O23729", "O23435", "000000"]
        go_list = ["C:cytosol", "F:pyroglutamyl-peptidase activity", "P:protein folding"]
       
   Output:
        -->Starting ExPASy analysis for A0A0K2RVI7_9BETC...
        A0A0K2RVI7_9BETC_GO.txt ready for A0A0K2RVI7_9BETC
        ->ID: A0A0K2RVI7_9BETC ExPASy analysis completed!

        (...)

        -->Starting ExPASy analysis for B8I7R6_RUMCH...
        No GOs from the list were found in B8I7R6_RUMCH.
        ->ID: B8I7R6_RUMCH ExPASy analysis completed!

        ->-> Full ExPASy analysis completed! <-<-
   
## SOFTWARE REQUIREMENTS:
    python3
    
## FUNCTIONS: There are six functions, five are
subfunctions of the sixth main one. Every function
is documented in the Functions section.
  
## EXTRA COMMENTS:
    This scripts needs BioPython installed.
    The script is robust to the next Errors:
    UnboundLocalError, ValueError. Look for details
    in Functions section.
    
## LAST MODIFICATION:
  Phabel Antonio Lopez Delgado: October, 2021. [Creation]

## SOURCE:
  GitHub: https://github.com/phabel-LD/python_classII/Tasks
'''

# Libraries #######################
from Bio import ExPASy
from Bio import SwissProt
from Bio import Entrez
from Bio.ExPASy import Prosite
import ssl


# Entrez email ####################
Entrez.email = "phabel@lcg.unam.mx"


# Ensure perimissions #############
ssl._create_default_https_context = ssl._create_unverified_context


# Functions #######################

# Get GOs
def get_go(cross_refs):
    '''This function gets the GO terms in the ID analyzed.
    
    Parameters:
        cross_refs: A list derived from a Bio.SwissProt.Record object.
                    Contains valuable info, including the ID's GO terms.               
    Return:
        go_terms: A list with the ID's GO terms to be validated afterwards.'''
    go_terms = []
    # Extract every GO
    for term in cross_refs:
        if "GO" in term:
            go_terms.append(term)
    return(go_terms)


# Validate GOs
def go_val(go_terms, go_list):
    '''This function validates the GO terms found by 'get_go()',
    by comparing them with the query GO list.
        
        Parameters:
            go_terms: A list with ID's GO terms to be validated.
            go_list: A list with query GO terms to search and used to
                    validate 'go_terms'.
        Return:
            go_found: A list with the validated GO terms from 'go_terms'
                    found in 'go_list'.'''
    go_found = []
    for go in go_list:
        # Assure GO query terms are in ID's GO list
        for term in go_terms:
            if go in term:
                go_found.append(term)
    return(go_found)


# Get one abstract
def get_abstract(record, ID):
    '''This function retrieves an abstract from the PubMed references
        for every ID with validated GO terms.
        
        Parameters:
            record: A Bio.SwissProt.Record object from the ID to be
                analyzed.
            ID: ID: A string with the ID from ExPASy to analyze.
    '''
    # Get one id number from Pubmed. Watch out for errors.
    try:
        # Get one PubMed reference id
        for refer in record.references:
            for db in refer.references:
                if "PubMed" in db:
                    id_num = db[1]
                    break
        # Get one abstract
        handle = Entrez.efetch(db = "pubmed", id = id_num, rettype = "abstract",
                          retmode = "text")
        abstract = handle.read()
        handle.close()
        return(abstract)
    except(UnboundLocalError):
        # No abstract
        print(f'No PubMed abstract available for {ID}.')
        return(0)
    pass

# Get PDOC
def get_pdoc(record, ID):
    '''This function gets de Prosite documentation (pdoc)
    for each validated ID.
        
        Parameters:
            record: A Bio.SwissProt.Record object from the ID
                to be analyzed.
            ID: A string with the ID from ExPASy to analyze.
        Return: A string with the Prosite documentation (pdoc).
    '''
    # Get Prosite ID
    for refer in record.cross_references:
        if "PROSITE" in refer:
            prosite_id = refer[1]
            break
    try: 
        # Access UniProt
        uniprot_handle = ExPASy.get_prosite_raw(prosite_id)
        uniprot = Prosite.read(uniprot_handle)
        # Get UniProt Documentation-pdoc
        pdoc = uniprot.pdoc
        return(pdoc)
    # No PDOC available
    except(UnboundLocalError):
        return(0)


# Write GO file
def go_file(ID, go_found, record):
    '''This function writes the ExPASy analysis .txt files for each ID whose
        GOs have been validated. Calls another function: 'get_abstract()'.
        
        Parameters: 
            ID: A string with the ID from ExPASy to analyze.
            go_found: A list derived from 'go_val' with the
                validated GO terms for the ID.
            record: A Bio.SwissProt.Record object from the ID
                to be analyzed.
        Return: No return value.
    '''
    with open(f'{record.entry_name}_GO.txt', 'w') as file:
        # ID and protein name
        file.write(f'->ID: {ID}.\n->PROTEIN NAME: {record.description}')
        # GOs and their definitions 
        for go in go_found:
            file.write(f'\n->GOs: {str(go)}')
        # Organism
        file.write(f'\n->ORGANISM: {record.organism}')
        # SubCell localization
        for comment in record.comments:
            if "SUBCELLULAR LOCATION" in comment:
                file.write(f'\n->{comment}')
        # Abstract
        abstract = get_abstract(record, ID)
        if not abstract:
            file.write("\n->No PubMed abstract available")
        else:
            file.write("\n->ABSTRACT:\n")
            file.write(abstract)
        # PDOC
        pdoc = get_pdoc(record, ID)
        if not pdoc:
            file.write("\n->No PDOC available.")
        else:
            file.write(f'\n->UniProt PDOC: {pdoc}')
        # Update user
        print(f'{record.entry_name}_GO.txt ready for {ID}')
    pass


# Main ExPASy analysis
def expasy_analysis(id_list, go_list):
    '''This main function directs the whole ExPASy anaylisis.
    Calls other four functions.
    
    Parameters:
        id_list: A list with all IDs query to analyze.
        go_list: A list with all GO terms to validate for every ID.
        
    Return: No return value.'''
    # Apply to every ID
    for ID in id_list:
        print(f'\n-->Starting ExPASy analysis for {ID}...')
        try:
            # Get record
            handle = ExPASy.get_sprot_raw(ID)
            record = SwissProt.read(handle)
            handle.close()
            # Get GO terms
            go_terms = get_go(record.cross_references)
            # Compare with go list to get GOs from the query list found in terms
            go_found = go_val(go_terms, go_list)
            # Write file if GOs were found
            if go_found:
                go_file(ID, go_found, record)
            else:
                print(f'No GOs from the list were found in {ID}.')
        except(ValueError):
            print(f'ID: {ID} not available.')
        finally:
            print(f'->ID: {ID} ExPASy analysis completed!')
    # Finish process
    print("\n->-> Full ExPASy analysis completed! <-<-")
    pass


# Main Code ######################

go_list = ["GO:0046755", "GO:0046761",
              "GO:0046760", "GO:0039702",
              "GO:0046765", "GO:0046762"]
id_list = ["A0A0K2RVI7_9BETC", "A8R4D4_9BETC",
                 "POLG_YEFV1", "POLG_DEN1W",
                 "Q6W352_CVEN9", "D9SV67_CLOC7",
                 "A9KSF7_LACP7", "B8I7R6_RUMCH"]


expasy_analysis(id_list, go_list)

