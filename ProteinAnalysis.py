'''
## NAME:
  ProteinAnalysis.py
  
## LANGUAGE & VERSION:
  python 3.8.5
  
## AUTHORS:
  Phabel Antonio Lopez Delgado <phabel2001@lcg.unam.mx>
  Daianna Gonzales Padilla <daianna#lcg.unam.mx>
  
## DATE:
  November, 2021.
  
## DESCRIPTION & LOGIC:
  This script uses BioPython tools numpy, pandas, seaborn,
  matplotlib and argparse to run a functional protein analysis
  between a protein query and series of proteins. In order to
  asses functional relatedness to a query protein from the analyzed.

## USAGE:
  ProteinAnalysis.py python 3.8.5
  
## ARGUMENTS & HELP:

Protein functional analysis and comparison

options:
  -h, --help            show this help message and exit
  -int_matx_df          Prints intersection-matrix dataframe
  -means_df             Prints means vector as dataframe.
  -best                 Prints best match protein.
  -heatmap              Prints intersection-matrix heatmap.
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        List of the proteins files paths, separated by
                        whitespace and comma
  -d DISULFIDE, --disulfide DISULFIDE
                        Distance between S-S atoms. Write -1 to use the
                        default value: 8
  -a ALPHA, --alpha ALPHA
                        Sequence pattern to search alpha helices. Use default
                        to search standard pattern
  -b BETA, --beta BETA  Sequence pattern to search beta sheets. Use default to
                        search standard pattern
  -m MOTIF [MOTIF ...], --motif MOTIF [MOTIF ...]
                        Motif to search and minimal length of the motif
                        sequence

## INPUT - OUTPUT:

   Input: List of protein files to be analyzed (-i,), distance between S-S atoms (-d),
           sequence pattern to search alpha helices (-a), sequence pattern to search
           beta sheets (-b), motif to search and minimal length of the motif (-m)

   
   Output: Prints motifs found.Prints functional analysis stadistics, depending on the analysis argument
           given (see ARGUMENTS & HELP). Saves my_protein object as an instance of ProtAnalysis class.

       
## EXAMPLES:

   Input: (From Terminal)
       python3 ProteinAnalysis.py -i 1kcw.pdb 1fat.pdb 3jbz.pdb 1kbe.pdb 4g68.pdb 1hp8.pdb
       -d 2 -a default -b default -int_matx_df -means_df -best -heatmap

       Query protein name as key (str): 1kbe

   Output: (Std output)
        {'name': '1kcw', 'num_bonds': 1, 'di_bonds': [[155, 181, 1.9980831, <Model id=0>, <Chain id=A>]]}
        {'name': '1kcw', 'num_helix': 4, 'helix_seqs': [[['RIYHSHIDAPKD', 'KEKEKHIDRE', 'KVDKDNEDFQE', 'KVNKDDEEFIE'], <Model id=0>, <Chain id=A>]]}
        (...)
        {'name': '1hp8', 'num_bonds': 0, 'num_helix': 0, 'num_sheets': 1}
        Proteins to be analyzed: [{'name': '1kcw', 'num_bonds': 1, 'num_helix': 4, 'num_sheets': 61}, {'name': '1fat', 'num_bonds': 0, 'num_helix': 4, 'num_sheets': 60}, {'name': '3jbz', 'num_bonds': 0, 'num_helix': 5, 'num_sheets': 58}, {'name': '1kbe', 'num_bonds': 0, 'num_helix': 0, 'num_sheets': 3}, {'name': '4g68', 'num_bonds': 0, 'num_helix': 3, 'num_sheets': 81}, {'name': '1hp8', 'num_bonds': 0, 'num_helix': 0, 'num_sheets': 1}]
        Valid protein names are: ['1kcw', '1fat', '3jbz', '1kbe', '4g68', '1hp8']
        Query protein name as key (str): 1kbe

        -------- Starting analysis --------

        Intersection matrix as np.array:

        [[0.         0.04918033]
         [0.         0.05      ]
         [0.         0.05172414]
         [0.         0.        ]
         [0.         0.03703704]
         [0.         0.33333333]] 


        Intersection matrix's row means as np.array:

        [0.02459016 0.025      0.02586207 0.         0.01851852 0.16666667] 


        Intersection matrix as dataframe:

              num_helix  num_sheets
        1kcw        0.0    0.049180
        1fat        0.0    0.050000
        3jbz        0.0    0.051724
        1kbe        0.0    0.000000
        4g68        0.0    0.037037
        1hp8        0.0    0.333333 

        Intersection matrix means vector as dataframe:

          1kbe
        1kcw  0.024590
        1fat  0.025000
        3jbz  0.025862
        1kbe  0.000000
        4g68  0.018519
        1hp8  0.166667 

        Printing best match...
        Best match           1kbe
        1hp8  0.166667

        Printing heatmap...

        -------- Analysis compleyed -------- 
   

## SOFTWARE REQUIREMENTS:
    python3
    argparse
    numpy
    pandas
    seaborn
    matplotlib.pyplot
    motifs.py
    
## FUNCTIONS: There are 8 functions, each one is necesary
    for a step in the analysis. Their documentation is in
    their docstrings. Many other functions are imported in
    motifs from <motifs.py>
  
## EXTRA COMMENTS:
   This script imports the module motifs.py

## LAST MODIFICATION:
  Phabel Antonio Lopez Delgado & Daianna Gonzáles Padilla: November, 2021. [Creation]

## SOURCE:
  GitHub: https://github.com/phabel-LD/python_classII/
  
'''


##################################################################
# Libraries
##################################################################

import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import motifs as mo

##################################################################
# Functions
##################################################################

# Dictionary of single letter code of aa necesaary for further analysis
aa_code = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def comparison_intersection(protein_A, protein_B):
    '''
    This function makes the intersection comparsion between two
    proteins.
    Parameters:
        protein_A (dict): Dictionary with info from the first protein,
                              including motif frequences.
        protein_B (dict): Dictionary with info from the first protein,
                              including motif frequences.
    Return:
        intersections (np.array): Array with intersection matrix.
    '''
    motifs_A = np.array(protein_A.keys())
    values_A = np.array(list(protein_A.values())[2:])
    motifs_B = np.array(protein_B.keys())
    values_B = np.array(list(protein_B.values())[2:])
    # Start Vector pairwise comparison
    intersections = []
    for index in range(0, len(values_A)):
        # Avoid division by 0.
        if max(values_A[index], values_B[index]) == 0:
            intersection = 0
        else:
            intersection = min(values_A[index], values_B[index])/ max(values_A[index], values_B[index])
        intersections.append(intersection)
    # Intersection vector is ready
    intersections = np.array(intersections)
    return(intersections)


def comparison_mean(protein_A, protein_B):
    '''
    This function makes the comparison mean from the instersection analysis
    of two proteins
    Parameters:
        protein_A (dict): Dictionary with info from the first protein,
                              including motif frequences.
        protein_B (dict): Dictionary with info from the first protein,
                              including motif frequences.
    Return:
        mean (np.array): Mean intersection value.
    '''
    intersections = comparison_intersection(protein_A, protein_B)
    # Mean intersection value as numpy array for further analysis
    mean = np.mean(intersections)
    return(mean)

def intersections_matrix(proteins, query_prot):
    '''
    This function makes the intersection matrix (generalized comparisson)
    between a query protein and a list of proteins.
    Parameters:
        proteins (list): A list of dictionaries, each with info. of a protein.
        query_prot (dict): A dictionary with the info. of protein to be analyzed.
    Return:
        comparison_matrix (np.array): Comparison matrix with intersection values.
    '''
    # Initialize matrix
    comparison_matrix = []
    # Analyze every query protein
    for index in range(0, len(proteins)):
        if (proteins[index] == query_prot):
            intersection = np.zeros(len(comparison_intersection(proteins[index], query_prot)))
        else:
            intersection = comparison_intersection(proteins[index], query_prot)
        comparison_matrix.append(intersection)
    # Create numpy array for further analysis
    comparison_matrix = np.array(comparison_matrix)
    return(comparison_matrix)

def intersections_means(proteins, query_prot):
    '''
    This function takes the row means from an intersection matrix.
    Parameters:
        proteins (list): A list of dictionaries, each with info. of a protein.
        query_prot (dict): A dictionary with the info. of protein to be analyzed.
    Return:
        means (np.array): Array of 1D with the row means.
    '''
    # Get matrix
    matrix = intersections_matrix(proteins, query_prot)
    # Get rows means
    means = np.mean(matrix, axis = 1)
    return(means)

def plot_heatmap(proteins, query_prot):
    '''
    This function plots a heatmap from the info of an intersection matrix, onced
    applied.
    Parameters:
        proteins (list): A list of dictionaries, each with info. of a protein.
        query_prot (dict): A dictionary with the info. of protein to be analyzed.
    Return: None
    '''
    # Get matrix
    matrix = intersections_matrix(proteins, query_prot)
    # Get x labels: motifs
    x = [ motif for motif in list(query_prot.keys())[2:] ]
    # Get y labels: proteins
    y = [prot["name"] for prot in proteins]
    # Get name of analyzed protein
    query_name = query_prot["name"]
    # Create Heatmap & plot
    heatmap_matrix = sns.heatmap(matrix, cmap="YlGnBu",
                             xticklabels = x, yticklabels = y)
    plt.title(f"{query_name} functional analysis")
    plt.xlabel("Motifs")
    plt.ylabel("Proteins")
    plt.show()
    return

def get_means_df(proteins, query_prot):
    '''
    This function takes the mean array from a intersection matrix an creates
    a dataframe
    Parameters:
        proteins (list): A list of dictionaries, each with info. of a protein.
        query_prot (dict): A dictionary with the info. of protein to be analyzed.
    Return:
        means_df (pd.dataframe): Dataframe with the row means.
    '''
    # Get means
    means = intersections_means(proteins, query_prot)
    # Create dataframe
    means_df = pd.DataFrame(means.T, columns = [query_prot["name"]])
    # Upgrae header and indexes
    old_index = [x for x in means_df.index]
    new_index = [prot["name"] for prot in proteins]
    means_df = means_df.rename(index=dict(zip(old_index, new_index)))
    return(means_df)

def intersect_matrix_df(proteins, query_prot):
    '''
    This function makes the intersection matrix (generalized comparisson)
    between a query protein and a list of proteins and create a dataframe
    from it.
    Parameters:
        proteins (list): A list of dictionaries, each with info. of a protein.
        query_prot (dict): A dictionary with the info. of protein to be analyzed.
    Return:
        intersection_df (pd.dataframe): Dataframe of the comparison matrix.
    '''
    # Get intersection matrix
    matrix = intersections_matrix(proteins, query_prot)
    # Rename columns according to motif names
    x = [ motif for motif in list(query_prot.keys())[2:] ]
    # Create dataframe
    intersection_df = pd.DataFrame(matrix, columns = x )
    # Upgrade indexes
    old_index = [x for x in intersection_df.index]
    new_index = [prot["name"] for prot in proteins]
    intersection_df = intersection_df.rename(index=dict(zip(old_index, new_index)))
    return(intersection_df)

def print_best_match(proteins, query_prot):
    '''
    This function determines the analyzed protein's best match as the protein with the
    highest intersection mean value. It prints the result.
    Parameters:
        proteins (list): A list of dictionaries, each with info. of a protein.
        query_prot (dict): A dictionary with the info. of protein to be analyzed.
    Return: None
    '''
    # Get mean dataframe
    means_df = get_means_df(proteins, query_prot)
    # Get best match motif represented as an index
    best_match = means_df.idxmax()
    # Show result
    print(f"Best match {means_df.loc[best_match] }")
    return


##################################################################
# Arguments & Parser
##################################################################

# Parser
parser = argparse.ArgumentParser(description = "Protein functional analysis and comparison")

# Arguments
# Print matrix dataframe
parser.add_argument("-int_matx_df",
                    action = "store_true",
                    help = "Prints intersection-matrix dataframe",
                    required = False)
# Print means dataframe
parser.add_argument("-means_df",
                    action = "store_true",
                    help = "Prints means vector as dataframe.",
                    required = False)
# Print best match protein for query
parser.add_argument("-best",
                    action = "store_true",
                    help = "Prints best match protein.",
                    required = False)
# Print heatmap
parser.add_argument("-heatmap",
                    action = "store_true",
                    help = "Prints intersection-matrix heatmap.",
                    required = False)
# Request input files
parser.add_argument("-i", "--input",  nargs='+',
                    type=str,
                    help="List of the proteins files paths, separated by whitespace and comma",
                    required=True)
# Give distance between S atoms
parser.add_argument("-d", "--disulfide",
                    type=float,
                    help="Distance between S-S atoms. Write -1 to use the default value: 8",
                    required=False)
# Give a pattern to search helices
parser.add_argument("-a", "--alpha",
                    type=str,
                    help="Sequence pattern to search alpha helices. Use default to search standard pattern",
                    required=False)
# Give a pattern to search sheets
parser.add_argument("-b", "--beta",
                    type=str,
                    help="Sequence pattern to search beta sheets. Use default to search standard pattern",
                    required=False)
# Give a new motif to search and its minimal length
parser.add_argument("-m", "--motif", nargs='+',
                    type=str,
                    help="Motif to search and minimal length of the motif sequence",
                    required=False)


# Parse arguments
args = parser.parse_args()

# Assign argument values
matrix_df_action = args.int_matx_df
means_df_action = args.means_df
best_action = args.best
heatmap_action = args.heatmap

##################################################################
# ProtAnalysis class
################################################################## 
class ProtAnalysis():
    ''' This class stores data from the protein analysis. '''
    def __init__(self, query_name, query_prot, proteins, matrix, means):
        self.query_name = query_name
        self.query_prot = query_prot
        self.proteins = proteins
        self.matrix = matrix
        self.means = means


##################################################################
# Main Code
##################################################################

# Get protein dictionaries
if not args.alpha and not args.disulfide and not args.beta and not args.motif:
    print('At least one type of motif is required')

# For each path, check if they exist
print("\n-------- Motifs -------- \n")
proteins = []
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
            protein = mo.unique_dict(path, prot_name, args.disulfide, args.alpha, args.beta, args.motif)
            proteins.append(protein)
        else:
            print(f"\n{path} file must have .pdb format\n")
    except FileNotFoundError as ex:
        print(path,' : File not found')

# Get protein names
protein_names = []
for protein in proteins:
    protein_names.append(protein["name"])

# Ask for query protein name
print(f"Proteins to be analyzed: {proteins}")
print(f"Valid protein names are: {protein_names}")
query_name = input("Query protein name as key (str): ")
query_valid = False

# Validate query protein
if query_name in protein_names:
    for protein in proteins:
        if protein["name"] == query_name:
            query_prot = protein
            query_valid = True
# Query protein is not found
if not query_valid:
    print(f"No query protein named {query_name} in searched protein names. Valid names are: {protein_names}")
    # Exit script with error
    exit(1)

# Start Analysis
print("\n-------- Starting analysis --------")

# Default analysis
# Print intersection matrix as np.array
matrix = intersections_matrix(proteins, query_prot)
print("\nIntersection matrix as np.array:\n")
print(matrix, "\n\n")

# Print intersection matrix's row means as np.array
means = intersections_means(proteins, query_prot)
print("Intersection matrix's row means as np.array:\n")
print(means, "\n\n")

# Create object ProtAnalysis
my_protein = ProtAnalysis(query_name, query_prot, proteins, matrix, means)


# Analysis with arguments
if matrix_df_action:
    intersection_df = intersect_matrix_df(my_protein.proteins, my_protein.query_prot)
    print("\nIntersection matrix as dataframe:\n")
    print(intersection_df, "\n\n")

if means_df_action:
    means_df = get_means_df(my_protein.proteins, my_protein.query_prot)
    print("\nIntersection matrix means vector as dataframe:\n")
    print(means_df, "\n\n")

if best_action:
    print("Printing best match...")
    print_best_match(my_protein.proteins, my_protein.query_prot)
    print("\n\n")

if heatmap_action:
    print("Printing heatmap...")
    plot_heatmap(my_protein.proteins, my_protein.query_prot)
    print("\n\n")

print("-------- Analysis compleyed -------- \n\n")

