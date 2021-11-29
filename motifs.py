'''
NAME
    motifs.py
VERSION
    2.1.2
AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>
DESCRIPTION
    This module is part of prot_motifs package which obtains the number of proteins motifs in a protein and
    the sequences themselves. The module allows to search disulfide bonds, alpha helices, beta sheets and
    specific motifs given by the user.
CATEGORY
     PDB files analysis
USAGE
    .\ds_bonds.py -i [str str ...] -d [float] -a [str] -b [str] -m [str str str]
ARGUMENTS
    -i, --input (strings separated by whitespaces)
            PDB files paths
    -d, --disulfide (float)
            This option indicates to find disulfide bonds and receives the distance between sulfur atoms
            implicated in the bond, in Armstrongs, or uses 8 as default by passing -1
    -a, --alpha (str in double quotes)
            This option indicates to find alpha helices in the protein and receives a search pattern or
            uses the default option by passing 'default'
    -b, --beta (str in double quotes)
            This option indicates to find beta sheets in the protein and receives a search pattern or
            uses the default option by passing 'default'
    -m, --motif (pattern str in double quotes, int as minimal length, and the motif name as str, separated
                by whitespaces)
            This option indicates to find new motifs in the protein and receives the search pattern,
            the minimal length of the motif in residues and the motif name.
INPUT
    Requires the relative or absolute paths to files with .pdb extension and format.
    Depending on the arguments passed, the additional information provided.
OUTPUT
    The script execution  prints the output dictionaries generated for each file protein
    containing the protein PDB ID as name, the number of motifs found and the motifs. And returns
    each motif sequences number found in each protein.
EXAMPLES
    Example 1.1: receives in terminal:
         .\motifs.py -i 1kcw.pdb 1fat.pdb -d 2 -a  default -b "/BT.{1,3}H/B" -m "V.{4}" 3 beta_turns
    and prints:
        {'name': '1kcw', 'num_bonds': 1, 'di_bonds': [[155, 181, 1.9980831, <Model id=0>, <Chain id=A>]]}
        {'name': '1kcw', 'num_helix': 4, 'helix_seqs': [[['RIYHSHIDAPKD',..., 'KVNKDDEEFIE'],
                                                                <Model id=0>, <Chain id=A>]]}
        {'name': '1kcw', 'num_sheets': 0, 'sheets_seqs': [[[], <Model id=0>, <Chain id=A>]]}
        {'name': '1kcw', 'num_beta_turns': 50, 'beta_turns_seqs': [[['VDTEH',..., 'VLQN-'], <Model id=0>,
                                                                                        <Chain id=A>]]}
        {'name': '1kcw', 'num_bonds': 1, 'num_helix': 4, 'num_sheets': 0, 'num_beta_turns': 50}
        {'name': '1fat', 'num_bonds': 0, 'di_bonds': []}
        {'name': '1fat', 'num_helix': 4, 'helix_seqs': [[['KTTRWDFVNGE'], <Model id=0>, <Chain id=A>],
                                                        [['KTTRWDFVNGE'], <Model id=0>, <Chain id=B>],
                                                        [['KTTRWDFVNGE'], <Model id=0>, <Chain id=C>],
                                                        [['KTTRWDFVNGE'], <Model id=0>, <Chain id=D>]]}
        {'name': '1fat', 'num_sheets': 0, 'sheets_seqs': [[[], <Model id=0>, <Chain id=A>], [[], <Model id=0>,
                            <Chain id=B>], [[], <Model id=0>, <Chain id=C>], [[], <Model id=0>, <Chain id=D>]]}
        {'name': '1fat', 'num_beta_turns': 60, 'beta_turns_seqs': [[['VSSSG',...,'VLSWS'], <Model id=0>, <Chain id=A>],
                                                                   [['VSSSG',..,'VLSWS'], <Model id=0>, <Chain id=B>],
                                                                   [['VSSSG',...,'VLSWS'], <Model id=0>, <Chain id=C>],
                                                                   [['VSSSG',...,'VLSWS'], <Model id=0>, <Chain id=D>]]}
        {'name': '1fat', 'num_bonds': 0, 'num_helix': 4, 'num_sheets': 0, 'num_beta_turns': 60}
SOURCE
    https://github.com/daianna21/python_class/blob/master/tareas/Entrez_search.py
'''
from Bio import PDB
import argparse
import re

'''
# Create the argument parser
my_parser = argparse.ArgumentParser(description="Program that obtains protein motifs")
# Add arguments: choose which motifs to search
# Request input files
my_parser.add_argument("-i", "--input",  nargs='+',
                    type=str,
                    help="List of the proteins files paths, separated by whitespace and comma",
                    required=True)
# Give distance between S atoms
my_parser.add_argument("-d", "--disulfide",
                    type=float,
                    help="Distance between S-S atoms. Write -1 to use the default value: 8",
                    required=False)
# Give a pattern to search helices
my_parser.add_argument("-a", "--alpha",
                    type=str,
                    help="Sequence pattern to search alpha helices. Use default to search standard pattern",
                    required=False)
# Give a pattern to search sheets
my_parser.add_argument("-b", "--beta",
                    type=str,
                    help="Sequence pattern to search beta sheets. Use default to search standard pattern",
                    required=False)
# Give a new motif to search, its minimal length and name
my_parser.add_argument("-m", "--motif", nargs='+',
                    type=str,
                    help="Motif to search, minimal length of the motif sequence and its name",
                    required=False)

# Execute the parse_args() method
args = my_parser.parse_args()
'''

# Dictionary of single letter code of aa
aa_code = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
 

def ds_bond(path, prot_name, distance=8):
    """
    This function obtains potential disulfide bonds between sulfur atoms in a protein sequence from a pdb file.
            Parameters:
                    path (str): relative or absolute path to .pdb file.
                    prot_name (str): the file protein name as a PDB ID
                    distance (float): the maximal distance between sulfur atoms that form disulfide bonds
                                    (in Armstrongs). 8 is taken as default.
            Prints:
                    prot_dic (dict): dictionary with the protein name, number of disulfide bonds in
                                    the whole protein, and the bonds themselves as lists with the
                                    implicated Cys residues IDs, the distance between their sulfur
                                    atoms, and the chain and model they belong to.
            Returns:
                    num_bonds (str): string indicating the number of bonds
                    c (int): number of bonds found
    """
    # Lists of bonds
    bonds=[]
    # Number of bonds
    c=0
    # Create the PDB parser and ignore warnings
    parser = PDB.PDBParser(QUIET=True)
    # Get the structure from the file
    struct=parser.get_structure('protein',path)
    # For each model, get each chain
    for model in struct:
        for chain in model:
            cys_chain=[]
            # Search Cys residues in the chain and store them
            for residue in chain:
                if residue.get_resname() == 'CYS':
                    cys_chain.append(residue)
            pairs = []
            # For each unique pair of Cys residues get its ID's and its sulfur atoms
            for cys_1 in cys_chain:
                for cys_2 in cys_chain:
                    if not cys_1==cys_2:
                        # Distance between sulfur atoms less than 8 A for unique pairs of Cys
                        if (cys_1['SG']-cys_2['SG'])<distance and [cys_1.get_id()[1], cys_2.get_id()[1]] \
                                not in pairs and [cys_2.get_id()[1], cys_1.get_id()[1]] not in pairs:
                            pairs.append([cys_1.get_id()[1], cys_2.get_id()[1]])
                            bonds.append([cys_1.get_id()[1], cys_2.get_id()[1], cys_1['SG']-cys_2['SG'],
                                          model, chain])
                            c=c+1
    # Output dictionary
    prot_dic={'name':prot_name, 'num_bonds': c, 'di_bonds':bonds}
    print(prot_dic)
    return('num_bonds', c)



def al_helix(path, prot_name, reg_exp="\B[HKR].{0,4}[HKR]{1}[^HRKPSG]{1,3}[DE]{1}.{0,4}[DE]\B"):
    """
    This function obtains potential alpha helices in a protein from a pdb file according to a search pattern
            Parameters:
                    path (str): relative or absolute path to .pdb file
                    prot_name (str): the file protein name as a PDB ID
                    reg_exp (str): regular expression in double quotes that describes the sequence pattern to
                                   search in the protein sequence. By default it uses "\B[HKR].{0,4}[HKR]{1}
                                   [^HRKPSG]{1,3} [DE]{1}.{0,4}[DE]\B" searching sequences that start with
                                   positive aa, end with negative ones and have positive aa 1 or 3 residues
                                   away from negatives.
            Prints:
                    prot_dic (dict): dictionary with the protein name, number of alpha helices found under
                                    the given patter, and the helices themselves as lists of the
                                    alpha helices sequences, and the model and chain they belong to.
            Returns:
                    num_alpha (str): string indicating the number of alpha helices
                    c (int): number of alpha helices found
    """
    # Count and store helices
    c=0
    motifs = []
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure('protein', path)
    # Iterate each model, chain, and residue
    for model in struct:
        for chain in model:
            aa_seq = []
            helix_seqs = []
            # Obtain the aa sequence of the protein
            for residue in chain:
                if residue.resname in aa_code.keys():
                    aa_seq.append(aa_code[residue.resname])
                # If element is not an aa
                else:
                    aa_seq.append('-')
            aa_seq = ''.join(aa_seq)
            aa_seq = str(aa_seq)
            # Search the pattern in the protein sequence
            matches = re.findall(reg_exp, aa_seq)
            # Save only helices of at least 10 residues
            for helix in matches:
                if len(helix) >= 10:
                    helix_seqs.append(helix)
                    c=c+1
            # Helices list
            motifs.append([helix_seqs, model, chain])
    prot_dic = {'name': prot_name, 'num_helix': c, 'helix_seqs': motifs}
    print(prot_dic)
    return('num_helix',c)


def b_sheets(path, prot_name, reg_exp="\B.{1,4}[V,I,T,F,Y,W]+.{1,4}[V,I,T,F,Y,W]+.{1,4}"):
    """
    This function obtains potential beta sheets in a protein from a pdb file according to a search pattern
            Parameters:
                    path (str): relative or absolute path to .pdb file
                    prot_name (str): the file protein name as a PDB ID
                    reg_exp (str): regular expression in double that describes the sequence pattern to
                                   search in the protein sequence. By default it uses "\B.{1,4}[V,I,T,F,Y,W]
                                   +.{1,4} [V,I,T,F,Y,W]+.{1,4}", searching sequences with aromatic aa and
                                   beta branched aa in the middle of the sequence.
            Prints:
                    prot_dic (dict): dictionary with the protein name, number of beta sheets found under
                                    the given patter, and the sheets themselves as lists of the
                                    beta sheets sequences, and the model and chain they belong to.
            Returns:
                    num_beta (str): string indicating the number of beta sheets
                    c (int): number of beta sheets found
    """
    # Count and store beta sheets
    c=0
    motifs = []
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure('protein', path)
    # Get the protein sequence
    for model in struct:
        for chain in model:
            aa_seq = []
            for residue in chain:
                if residue.resname in aa_code.keys():
                    aa_seq.append(aa_code[residue.resname])
                else:
                    aa_seq.append('-')
            aa_seq = ''.join(aa_seq)
            aa_seq = str(aa_seq)
            # Search the pattern in the sequence
            matches = re.findall(reg_exp, aa_seq)
            sheets_seqs = []
            # Save only sequences larger than 2 residues
            for sheet in matches:
                if len(sheet) >= 3:
                    sheets_seqs.append(sheet)
                    c=c+1
            motifs.append([sheets_seqs, model, chain])
    prot_dic = {'name': prot_name, 'num_sheets': c, 'sheets_seqs': motifs}
    print(prot_dic)
    return('num_sheets',c)


def others(path, prot_name, size, reg_exp, motif_name):
    """
    This function obtains potential sequences of a certain motif in a protein from a pdb file
            Parameters:
                    path (str): relative or absolute path to .pdb file
                    prot_name (str): the file protein name as a PDB ID
                    size (int): minimal size of the motif sequences (in residues)
                    reg_exp (str): regular expression in double that describes the motif pattern to
                                   search in the protein sequence. No default.
                    motif_name(str): motif name
            Prints:
                    prot_dic (dict): dictionary with the protein name, number of motif sequences found
                                    under the given patter, and the sequences themselves as lists of the
                                    motif sequences, and the model and chain they belong to.
            Returns:
                    num_[motif name] (str): string indicating the number of the motif sequences
                    c (int): number of those motifs found
    """
    # Count and store the motifs
    c=0
    motifs = []
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure('protein', path)
    # Obtain the protein sequence
    for model in struct:
        for chain in model:
            aa_seq = []
            motif_seqs = []
            for residue in chain:
                if residue.resname in aa_code.keys():
                    aa_seq.append(aa_code[residue.resname])
                else:
                    aa_seq.append('-')
            aa_seq = ''.join(aa_seq)
            aa_seq = str(aa_seq)

            # Search the particular motif pattern in the sequence
            matches = re.findall(reg_exp, aa_seq)
            for sheet in matches:
                if len(sheet) >= size:
                    motif_seqs.append(sheet)
                    c=c+1
            motifs.append([motif_seqs, model, chain])
    prot_dic = {'name': prot_name, f'num_{motif_name}': c, f'{motif_name}_seqs': motifs}
    print(prot_dic)
    return(f'num_{motif_name}',c)


def unique_dict(path,prot_name, args_disulfide, args_alpha, args_beta, args_motif):
    """
    This function makes a count of the number of the motifs given found in the protein.
            Parameters:
                    path (str): relative or absolute path to .pdb file
                    prot_name (str): the file protein name as a PDB ID
                    args_disulfide (class str): -d arguments given by the user. 'None' if they don't exist
                    args_alpha (class str): -a arguments given by the user. 'None' if they don't exist
                    args_beta (class str): -b arguments given by the user. 'None' if they don't exist
                    args_motif (class str): -m arguments given by the user. 'None' if they don't exist
            Prints & Returns:
                    dic (dict): dictionary with the protein name and the number of each of the given motifs
                                of the protein
    """
    # Create dictionary
    dic={'name': prot_name,}
    # Check if -d option is given
    if not args_disulfide == None:
        # Default distance value
        if args_disulfide == -1:
            count=ds_bond(path, prot_name)
        # Use distance value given
        else:
            count=ds_bond(path, prot_name, args_disulfide)
        # Add key and value
        dic[count[0]]=count[1]


    # Check if -a argument is given
    if not args_alpha == None:
        # Use default pattern
        if args_alpha == 'default':
            count=al_helix(path, prot_name)
        # Use pattern given
        else:
            count = al_helix(path, prot_name, args_alpha)
        # Add key and value
        dic[count[0]]=count[1]

    # Check if -b argument is given
    if not args_beta == None:
        # Use default pattern
        if args_beta == 'default':
           count= b_sheets(path, prot_name)
        # Use pattern given
        else:
            count=b_sheets(path, prot_name, args_beta)
        # Add key and value
        dic[count[0]]=count[1]

    # Check if -m argument is given
    if not args_motif == None:
        count = others(path, prot_name, int(args_motif[1]), args_motif[0], args_motif[2])
        # Add key and value
        dic[count[0]]=count[1]

    print(dic)
    return (dic)

'''
# Check there are arguments for at least one motif
if not args.alpha and not args.disulfide and not args.beta and not args.motif:
    print('At least one type of motif is required')

# For each path, check if they exist
for path in args.input:
    try:
        # Check they are .pdb files
        if path.endswith('.pdb'):
            # Parse protein name from the path file
            if not '/' in path:
                prot_name=str(path).split('.')[0]
            else:
                prot_name = str(path).split('/')[-1]
                prot_name = str(prot_name).split('.')[0]

            # Obtain motifs for each protein
            unique_dict(path,prot_name, args.disulfide, args.alpha, args.beta, args.motif)
        else:
            print(path,' file must have .pdb format')
    except FileNotFoundError as ex:
        print(path,' : File not found')
'''
