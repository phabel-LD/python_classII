'''
## NAME:
  dtype_array.py
  
## LANGUAGE & VERSION:
  python 3.8.5
  
## AUTHOR:
  Phabel Antonio Lopez Delgado <phabel2001@gmail.com>
  
## DATE:
  October, 2021.
  
## DESCRIPTION & LOGIC:
  This script uses Numpy to customize arrays with dtype,
  in order to get arryas with differente data types.

## USAGE:
  dtype_array.py python 3.8.5
  
## ARGUMENTS:
  The script receives no args from Terminal.
       
## EXAMPLES:
    sort_35C_gain = np.sort(gain_dtype, order = '35C_gain')
    print(sort_35C_gain)

        [('Gen_2', 0.45454545, 0.71428571)
        ('Gen_4', 2.15      , 0.71666667)
        ('Gen_3', 1.75      , 0.77777778)
        ('Gen_1', 0.7       , 1.16666667)]


## SOFTWARE REQUIREMENTS:
    python3

# Libraries
    numpy

## LAST MODIFICATION:
  Phabel Antonio Lopez Delgado: October, 2021. [Creation]

## SOURCE:
  GitHub: https://github.com/phabel-LD/python_classII/Tasks/dtype_array.py
'''
# Libraries #####################
import numpy as np


# Main Code #####################

# Original production matrix
production = np.array([[5,3], [11,7], [4,9], [2,6]])
print(production)

# Dtype production matrix. Use str and float64
production_dtype = np.array( [('Gen_1', 5.0, 3.0), ('Gen_2', 11.0, 7.0),
                       ('Gen_3', 4.0, 9.0), ('Gen_4', 2.0, 6.0)],
                      dtype = [('Gen', (np.str_, 5)),
                               ('30C', np.float64),
                               ('35C', np.float64) ] )
print(production_dtype)

# Original costs matrix
costs = np.array([3.5, 5, 7, 4.3])
print(costs)

# Dtype costs matrix. Use str and float64
costs_dtype = np.array( [('Gen_1', 3.5), ('Gen_2', 5),
                       ('Gen_3', 7), ('Gen_4', 4.3)],
                      dtype = [('Gen', (np.str_, 5)),
                               ('Cost', np.float64)] )
print(costs_dtype)

# Get gain matrix with original matrixes
gain = (costs / production.T).T
print(gain)

# Adjust dtype gain matrix. Use str and float64
gain_dtype = np.array( [('Gen_1', 0.7, 1.16666667), ('Gen_2', 0.45454545, 0.71428571),
                       ('Gen_3', 1.75, 0.77777778), ('Gen_4', 2.15, 0.71666667)],
                      dtype = [('Gen', (np.str_, 5)),
                               ('30C_gain', np.float64),
                               ('35C_gain', np.float64) ] )
print(gain_dtype)

# Order dtype gain matrix...
# by Gene str
sort_gene = np.sort(gain_dtype, order = 'Gen')
print(sort_gene)
# by 30C float64
sort_30C_gain = np.sort(gain_dtype, order = '30C_gain')
print(sort_30C_gain)
# by 35C float64
sort_35C_gain = np.sort(gain_dtype, order = '35C_gain')
print(sort_35C_gain)

