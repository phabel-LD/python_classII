'''
## NAME:
  ClassesOOP.py

## LANGUAGE & VERSION:
  python 3.8.5

## AUTHOR:
  Phabel Antonio Lopez Delgado <phabel2001@gmail.com>

## DATE:
   September 03rd, 2021.

## DESCRIPTION & LOGIC:
   This script exemplifies the properties of OOP:
   heredability, overriding & polimorphisms

## USAGE:
  ClassesOOP.py python 3.8.5

## EXAMPLES:
   No exmaples of input, just output from print:
        {'ID': 1234, 'length_kb': 3.5, 'GC_content': 0.65, 'Tm': 24, 'species': 'E. coli'}
        {'ID': 1234, 'length_kb': 3.5, 'GC_content': 0.65, 'Tm': 24, 'species': 'E. coli', 'DNA': True, 'RNA': False}
        {'ID': 2341, 'length_kb': 3.4, 'GC_content': 0.63, 'Tm': 28, 'species': 'Coronavirus', 'RNA': True, 'DNA': False, 'strand': 'Single(+)'}
        {'ID': 65432, 'length_kb': 3.5, 'GC_content': 0.6, 'Tm': 30, 'species': 'Coronavirus', 'strand_type': 'RNA'}
        RNA

## SOFTWARE REQUIREMENTS:
    python3

## SOURCE:
  GitHub: https://github.com/phabel-LD/python_classII

'''

## Tarea 1: Crear clases propias con métodos. Explicar uso.

## 1) Stablish my class, which can be used in order to characterize a
##    DNA or RNA sequence from any known species. It is meant to receive a 
##    full characterization with several parameters, so it can be used to fill
##    in databases.
class Dna_sequence():
  RNA = False
  DNA = False
  
  ## Constructor
  def __init__(self, ID, length_kb, GC_content, Tm, species):
    self.ID = ID
    self.length_kb = length_kb
    self.GC_content = GC_content
    self.Tm = Tm
    self.species = species
    
  ## Assign DNA type
  def DNA(self):
    self.DNA = True
    self.RNA = False
    
  ## Assign RNA type
  def RNA(self):
    self.RNA = True
    self.DNA = False
    
  ## Access DNA/RNA type
  def strand(self):
    print(self.DNA)
    print(self.RNA)

seq1 = Dna_sequence(1234, 3.5, .65, 24, 'E. coli')
print(seq1.__dict__)
seq1.DNA()
print(seq1.__dict__)

## 2) Apply inheritance to characterize a virus strand and biol. cycle. Like the
##    previous class, it is for a full virus (like SarsCov2) characterization for
##    further experiments and databases annotation.
class Virus(Dna_sequence):
  ## Strand attribute: single/doble
  strand_type = ''
  
  ## Assign cycle type
  def lysogenic(self):
    self.cycle = 'Lysogenic'
    
  def lytic(self):
    self.cycle = 'Lytic'

SarsCov2 = Virus(2341, 3.4, .63, 28, 'Coronavirus')
SarsCov2.RNA()
SarsCov2.strand = 'Single(+)'
print(SarsCov2.__dict__)

## 3) Apply overriding to access a virus strand type DNA or RNA. Herency from Virus()
## 4) Create a polymorphism by changing the strand() method to override the use
##    of strand_type
##    This last class can be implemented with a new approach for genomic analysis
##    of the virus
class Virus_strand(Virus):
  
  def strand(self):
    print(self.strand_type)

SarsCov2 = Virus_strand(65432, 3.5, 0.6, 30, 'Coronavirus')
SarsCov2.strand_type = 'RNA'
print(SarsCov2.__dict__)
SarsCov2.strand()


