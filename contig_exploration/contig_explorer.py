from Bio import SeqIO

# print lengths of all contigs
sequences = open("scobliq.fasta", "rU")
for record in SeqIO.parse(sequences, "fasta"):
    print "\n"    
    print record.id   
    print len(record.seq)
sequences.close()



