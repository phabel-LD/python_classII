'''
## NAME:
  ArchivosPBD.py
  
## LANGUAGE & VERSION:
  python 3.8.5
  
## AUTHOR:
  Phabel Antonio Lopez Delgado <phabel2001@gmail.com>
  
## DATE:
  September 22nd, 2021.
  
## DESCRIPTION & LOGIC:
  This script includes a function that receives a PBD file,
  a chain and an aminoacid to look for. Returns a list with
  all the residues of a particular amioacid in tha chain.

  If no chain 'chain' is in the file, a KeyError is caught.

## USAGE:
  ArchivosPBD.py python 3.8.5
  
## ARGUMENTS:
  No arguments are taken.

## INPUT - OUTPUT:

    INPUT: 3 parameters
        PBD file path
        aminoacid_list
        chains

    OUTPUT: A list with all aminoacid residues in the chain and a summary
    
## EXAMPLES:

  --> Script
    file = '1kcw.pdb'
    aminoacid_list = ['CYS', 'PHE']
    chains = ['A', 'B']

  --> Terminal
    For aminoacid CYS in chain A in file 1kcw.pdb...
    There are 14 CYS residues: 
    [<Residue CYS het=  resseq=155 icode= >, ..., <Residue CYS het=  resseq=1021 icode= >]

    For aminoacid PHE in chain A in file 1kcw.pdb...
    There are 49 PHE residues: 
     [<Residue PHE het=  resseq=61 icode= >, ..., <Residue PHE het=  resseq=1010 icode= >]

    No chain B in file 1kcw.pdb. Try other combination.
    There are no CYS in B
    No chain B in file 1kcw.pdb. Try other combination.
    There are no PHE in B

## SOFTWARE REQUIREMENTS:
    python3
    BioPython  -->  Library Bio
    
## FUNCTIONS:

  get_residue(path, chain, residue):
      Takes the path to parse the PBD file,
      access to chain 'chain' and adds to a
      list all 'residue'aminoacids. If a
      KeyError is caught, other combination
      must be tried.

## EXTRA COMMENTS:
   The function is robust to KeyError only,
   thus the path must be right.
    
    
## LAST MODIFICATION:
  Phabel Antonio Lopez Delgado: September 22nd, 2021  [Creation]

## SOURCE:
  GitHub: https://github.com/phabel-LD/python_classII/blob/master/Tasks/ArchivosPBD.py
  
'''

# Library
from Bio import PDB


# Function
def get_residue(path, chain, residue):
    """
    Takes the path to parse the PBD file, access to chain
    'chain' and adds to a list all 'residue'aminoacids.
    If a KeyError is caught, either the residue or chain
    do not exist, thus other combination must be tried.
    """
    try:
        parser = PDB.PDBParser(QUIET = True)
        pdb_structure = parser.get_structure('pdb_file', path)
        model = pdb_structure[0]
        aminoacids = []
        for aminoacid in model[chain]:
            if (aminoacid.get_resname() == residue):
                aminoacids.append(aminoacid)
        return (aminoacids)
    except(KeyError):
        print(f'No chain {chain} in file {path}. Try other chain.')
        pass


# Main Code
file = '1kcw.pdb'
aminoacid_list = ['CYS', 'PHE']
chains = ['A', 'B']
# Iterate for all chains
for chain in chains:
    # Iterate for all aminoacids
    for aminoacid in aminoacid_list:
        aa_list = get_residue(file, chain, aminoacid)
        if (aa_list):
            print(f'For aminoacid {aminoacid} in chain {chain} in file {file}...')
            print(f'There are {len(aa_list)} {aminoacid} residues: \n {aa_list}\n')
        # No aminoacid
        else:
            print(f'There are no {aminoacid} in {chain}')

