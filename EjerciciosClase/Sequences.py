## Sequences

from Bio.Seq import Seq
from Bio.Seq import MutableSeq

seqobj = Bio.Seq.Seq("ATCGTG")
print(str(seqobj))

mutable = MutableSeq(seqobj)
seqobj[0] = "T"
