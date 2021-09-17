from Bio import SeqIO

filename = ".../seq.nt.fa"
## Check every seqrecord from file out
for seq_record in SeqIO.parse(filename, "fasta"):
    print("ID {}".format(seq_record.id))
    print("len{}".format(len(seq_record)))
    print("Traduccion {}".format(seq_record.seq.translate(to_stop=False)))

## Read file in dict
id_dict = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
print(id_dict["seq4"].seq.transcribe())
for i in id_dict:
    print(">{}".format(i))
    ## Sequence to modify
    sec = id_dict[i].seq
    for codon in re.findall(r".{3}", str(sec[1:len(sec)])):
        print(codon, end = "\t")
    print("\n")
