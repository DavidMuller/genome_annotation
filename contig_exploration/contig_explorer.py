from Bio import SeqIO

# print lengths of all contigs
sequences = open("scobliq.fasta", "rU")
for record in SeqIO.parse(sequences, "fasta"):
    print "\n"    
    print record.id    
    print len(record.seq)
sequences.close()


"""
# print contigs longer than 125000...
sequences = open("scobliq.fasta", "rU")
total = 0
for record in SeqIO.parse(sequences, "fasta"):
    if len(record.seq) > 125000:
        print "\n"    
        print record.id    
        print len(record.seq)
sequences.close()
"""
