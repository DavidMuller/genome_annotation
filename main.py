import argparse
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline


################
# "Main":
################

# build the command-line parser
parser = argparse.ArgumentParser(description='CSE 182: Genome annotation project')
parser.add_argument('contig_data_base', help='fasta data base with all contigs', type=argparse.FileType('r'))
parser.add_argument('protein_data_base', help='Protein seq ', type=argparse.FileType('r'))

# parse command line arguments
args = parser.parse_args()  


for record in SeqIO.parse(args.contig_data_base, "fasta"):
    if record.id == "Contig130":
        blastx_cline = NcbiblastxCommandline(query=record.seq, db=args.protein_data_base)
        stdout, stderr = blastx_cline()
    



