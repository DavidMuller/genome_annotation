# create a fasta file with only our contigs

from Bio import SeqIO

our_contigs = {}
our_contigs["Contig243"] = 'Y'
our_contigs["Contig133"] = 'Y'
our_contigs["Contig170"] = 'Y'
our_contigs["Contig236"] = 'Y'
our_contigs["Contig323"] = 'Y'
our_contigs["Contig410"] = 'Y'
our_contigs["Contig180"] = 'Y'
our_contigs["Contig181"] = 'Y'
our_contigs["Contig179"] = 'Y'
our_contigs["Contig178"] = 'Y'
our_contigs["Contig175"] = 'Y'
our_contigs["Contig172"] = 'Y'
our_contigs["Contig163"] = 'Y'
our_contigs["Contig155"] = 'Y'
our_contigs["Contig154"] = 'Y'
our_contigs["Contig127"] = 'Y'
our_contigs["Contig132"] = 'Y'
our_contigs["Contig131"] = 'Y'
our_contigs["Contig128"] = 'Y'
our_contigs["Contig129"] = 'Y'
our_contigs["Contig130"] = 'Y'
our_contigs["Contig157"] = 'Y'
our_contigs["Contig171"] = 'Y'
our_contigs["Contig173"] = 'Y'
our_contigs["Contig174"] = 'Y'
our_contigs["Contig177"] = 'Y'
our_contigs["Contig237"] = 'Y'
our_contigs["Contig238"] = 'Y'
our_contigs["Contig239"] = 'Y'
our_contigs["Contig240"] = 'Y'
our_contigs["Contig241"] = 'Y'
our_contigs["Contig242"] = 'Y'
our_contigs["Contig411"] = 'Y'
our_contigs["Contig412"] = 'Y'
our_contigs["Contig413"] = 'Y'
our_contigs["Contig324"] = 'Y'
our_contigs["Contig325"] = 'Y'
our_contigs["Contig326"] = 'Y'
our_contigs["Contig158"] = 'Y'
our_contigs["Contig159"] = 'Y'
our_contigs["Contig160"] = 'Y'
our_contigs["Contig161"] = 'Y'
our_contigs["Contig162"] = 'Y'


# append only our contigs to "our_contig_list"
our_contig_list = []
sequences = open("scobliq.fasta", "rU")
for record in SeqIO.parse(sequences, "fasta"):
    if our_contigs.has_key(record.id):    
        our_contig_list.append(record)
sequences.close()

# create a fasta file with these contigs
SeqIO.write(our_contig_list, "our_contigs.fasta", "fasta")
