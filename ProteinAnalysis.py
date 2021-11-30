'''
## NAME:
  ProteinAnalysis.py
  
## LANGUAGE & VERSION:
  python 3.8.5
  
## AUTHOR:
  Phabel Antonio Lopez Delgado <phabel2001@gmail.com>
  
## DATE:
  November, 2021.
  
## DESCRIPTION & LOGIC:
  This script uses BioPython tools numpy, pandas, seaborn,
  matplotlib and argparse to run a functional protein analysis
  between a protein query and series of proteins. In order to
  asses funciontal relatedness.

## USAGE:
  ProteinAnalysis.py python 3.8.5
  
## ARGUMENTS:

  -h, --help
          Show this help message and exit
  -int_matx_df
          Prints intersection-matrix dataframe
  -means_df
          Prints means vector as dataframe.
  -best
          Prints best match protein.
  -heatmap
          Prints intersection-matrix heatmap.

## INPUT - OUTPUT:

   Input: 

   
   Output: Prints functional analysis stadistics.

       
## EXAMPLES:

   Input:



   Output:


## SOFTWARE REQUIREMENTS:
    python3
    argparse
    numpy
    pandas
    seaborn
    matplotlib.pyplot
    
## FUNCTIONS: There are 8 functions, each one is necesary
    for a step in the analysis. Their documentation is in
    their docstrings.
  
## EXTRA COMMENTS:


    
    
## LAST MODIFICATION:
  Phabel Antonio Lopez Delgado: November, 2021. [Creation]

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

##################################################################
# Functions
##################################################################

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
parser.add_argument("-int_matx_df",
                    action = "store_true",
                    help = "Prints intersection-matrix dataframe",
                    required = False)

parser.add_argument("-means_df",
                    action = "store_true",
                    help = "Prints means vector as dataframe.",
                    required = False)

parser.add_argument("-best",
                    action = "store_true",
                    help = "Prints best match protein.",
                    required = False)

parser.add_argument("-heatmap",
                    action = "store_true",
                    help = "Prints intersection-matrix heatmap.",
                    required = False)

# Parse arguments
args = parser.parse_args()

# Assign argument values
matrix_df_action = args.int_matx_df
means_df_action = args.means_df
best_action = args.best
heatmap_action = args.heatmap


##################################################################
# Main Code
##################################################################

# Protein motif cuantification
prot1 = {"name": "A", "bonds": 5, "motif1": 4, "motif2": 5, "motif3": 0, "motif4": 12, "motif5": 1}
prot2 = {"name": "B", "bonds": 2, "motif1": 6, "motif2": 5, "motif3": 0, "motif4": 10, "motif5": 8}
prot3 = {"name": "C", "bonds": 6, "motif1": 0, "motif2": 1, "motif3": 3, "motif4": 12, "motif5": 18}
prot4 = {"name": "D", "bonds": 8, "motif1": 4, "motif2": 7, "motif3": 9, "motif4": 13, "motif5": 0}


#class ProteinData:
    


# Proteins to analyze as a list of dicts
proteins = [prot1, prot2, prot3, prot4]

# Set query protein to analyze
query_prot = prot3
# Set query protein's name as string
query_name = query_prot["name"]

# Start nnalysis
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

# Analysis with arguments
if matrix_df_action:
    intersection_df = intersect_matrix_df(proteins, query_prot)
    print("\nIntersection matrix as dataframe:\n")
    print(intersection_df, "\n\n")

if means_df_action:
    means_df = get_means_df(proteins, query_prot)
    print("\nIntersection matrix means vector as dataframe:\n")
    print(means_df, "\n\n")

if best_action:
    print("Printing best match...")
    print_best_match(proteins, query_prot)
    print("\n\n")

if heatmap_action:
    print("Printing heatmap...")
    plot_heatmap(proteins, query_prot)
    print("\n\n")

print("-------- Analysis completed -------- \n\n")

