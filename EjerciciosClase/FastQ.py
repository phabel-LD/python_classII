#############  FastQ
n = 0
for record in SeqIO.parse("files/sample.fastq", "fastq"):
    if n < 2:
        print("%s %s" % (record.seq))
        n += 1
    else: break
print(record.letter_annotation["phred_quality"])


for record in SeqIO.parse("files/sample.fastq", "fastq"):
    mean = sum(record.letter_annotations["phred quality"]) / len(record.letter_annotations["phred quality"])
    if (promedio < umbral):
        temp = (promedio, record.seq)
        mala_calidad.append(temp)

print("{} secuencias con promedio menor a umbral: {}\n".format(len(mala_calidad, umbral)))

