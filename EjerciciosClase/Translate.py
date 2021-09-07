from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio import SeqUtils
from Bio.SeqUtils import nt_search, GC, molecular_weight

sequence = "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"
inicio = Seq("ATG")
position = nt_search(str(sequence), inicio)
for ORF in range(1, len(position)):  ## Position 0 is patron
    ## Find every sequence per starting position
    sec_prot = sequence[ORF:]
    protein = sec_prot.translate(to_stop=True)
    print(protein)
