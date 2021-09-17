'''
## NAME:
  GenBankSummary.py
  
## LANGUAGE & VERSION:
  python 3.8.5
  
## AUTHOR:
  Phabel Antonio Lopez Delgado <phabel2001@gmail.com>
  
## DATE:
  September 16th, 2021.
  
## DESCRIPTION & LOGIC:
  This script uses BioPython's libbs: Bio & BioSeq to analize
  a GenBank file and provide a summary from its metadata.

## USAGE:
  GenBankSummary.py python 3.8.5
  
## ARGUMENTS:
  The script receives no args from Terminal

## INPUT - OUTPUT:
   Input: List of genes from script
   Output: GenBank file summary
       
## EXAMPLES:
   Input:
       genes = ['L', 'N', 'P']
       
   Output:
        Organism: Isfahan virus
        Date: 13-AUG-2018
        Sample country: ['Iran:Isfahan province']
        Isolate: ['Phlebotomus papatasi']
        Gene name: ['N']
        DNA: ATGACTTCTGTAGTA...
        RNA: AUGACUUCUGUAGUA...
        Protein: MTSVV...
        Gene name: ['P']
        DNA: ATGTCTCGACTCAAC...
        RNA: AUGUCUCGACUCAAC...
        Protein: MSRLN...
        Gene name: ['L']
        DNA: ATGGATGAGTACTCT...
        RNA: AUGGAUGAGUACUCU...
        Protein: MDEYS...
   
## SOFTWARE REQUIREMENTS:
    python3
    
## FUNCTIONS:
  Methods:
        GenBankSummary(path, genes):
            Receives a path (GenBank file) and parses it
            using SeqIO.parse(). Receives a list of genes
            to search in the file. Then gives information
            about de files metadata: organism, date, country,
            isolate and then gives a 15 NTs DNA, RNA and
            5 aa protein from de CDS for each gene.

## EXTRA COMMENTS:
    This scripts needs BioPython installed.
    
## LAST MODIFICATION:
  Phabel Antonio Lopez Delgado: September 16th, 2021  [Creation]

## SOURCE:
  GitHub: https://github.com/phabel-LD/python_classII/Tasks
'''

##############################################################################
## Libraries
##############################################################################
from Bio.Seq import Seq
from Bio import SeqIO

##############################################################################
## Function
##############################################################################
def GenBankSummary(path, genes):
    for gb_record in SeqIO.parse(path, "genbank"):
        ## Organism
        organism = gb_record.annotations['organism']
        print(f'Organism: {organism}')
        ## Date
        date = gb_record.annotations['date']
        print(f'Date: {date}')
        ## Sample country
        country = gb_record.features[0].qualifiers['country']
        print(f'Sample country: {country}')
        ## Isolate number
        isolate = gb_record.features[0].qualifiers['isolation_source']
        print(f'Isolate: {isolate}')
        ## Find CDS indexes in features for every gene
        indexes = []
        index = 2
        while (index < len(gb_record.features)):
            for gene in genes:
                if (gb_record.features[index].qualifiers['gene'][0] == gene):
                    indexes.append(index)
                    break
            ## Just take CDS features, ignore genes features
            index += 2
        ## Sequence analysis
        for index in indexes:
            ## Gene name
            gene_name = gb_record.features[index].qualifiers['gene']
            print(f'Gene name: {gene_name}')            
            ## 15 NTs DNA --> RNA --> Prot
            start = gb_record.features[index].location.nofuzzy_start
            end = gb_record.features[index].location.nofuzzy_end
            CDS = gb_record.seq[start:end]
            DNA = CDS[:15]
            print(f'DNA: {DNA}...')
            RNA = DNA.transcribe()
            print(f'RNA: {RNA}...')
            prot = DNA.translate()
            print(f'Protein: {prot}...')
    pass

##############################################################################
## Main Code
##############################################################################
path = "virus.gb"
genes = ['L', 'N', 'P']
GenBankSummary(path, genes)
