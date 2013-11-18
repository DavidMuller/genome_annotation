from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline
import argparse
import os
import glob

# Make a directory if it doesn't already exist. 
# Parameters: the name of the directory.  
def make_directory(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name) 


# Given a file, or a the full path to the file,
# return only the file name with no file extension.
def return_only_name(full_path):
    path_and_name = os.path.split(full_path)
    name_and_extension = os.path.splitext(path_and_name[1])
    return name_and_extension[0]


# Create individual fasta files from a multi-fasta file.
# Parameters: a multi-fasta file, name for a directory
# to put all the new individual fasta files in.
def make_individual_fatsa_files(multi_fasta, fasta_dir):
    make_directory(fasta_dir)
    # Create an individual fasta file for each of our contigs
    for record in SeqIO.parse(multi_fasta,"fasta"):
        file_name = record.id    
        SeqIO.write(record, fasta_dir + "/" + file_name +".fasta", "fasta")


# Run blastx on a directory containg fasta files. Put the
# output in a directory specified by blastx_output_dir.
# Parameters: name of a directory with all the fasta files,
# the protein data-base to use, and the name for an output directory 
# for blastx.
def blastx_fasta_files(fasta_dir, protein_db, blastx_output_dir):
    make_directory(blastx_output_dir)
    all_files = glob.glob(fasta_dir + "/*")
    for fasta_file in all_files:
        contig_name = return_only_name(fasta_file)
        saved_output = blastx_output_dir + "/" + contig_name + ".xml"
        print "Running blastx on " + contig_name + "..."
        blastx_cline = NcbiblastxCommandline(cmd='blastx', query=fasta_file, db=protein_db, outfmt=5, out = saved_output)
        stdout, stderr = blastx_cline()



################
# "Main":
################

# build the command-line parser
parser = argparse.ArgumentParser(description='CSE 182: Genome annotation project')
parser.add_argument('contig_data_base', help='fasta data base with all contigs')
parser.add_argument('protein_data_base', help='Protein db ')

# parse command line arguments
args = parser.parse_args()  

individual_fasta_dir = "individual_contigs"
blastx_output = "blastx_output"

print "\nCreating Individual Fasta files for all our contigs:"
print "Saving them in folder " + "'" + individual_fasta_dir + "'"

# break up all the contigs into individual contig files 
make_individual_fatsa_files(args.contig_data_base, individual_fasta_dir)

print "\nRunning Blastx on each of our contigs:"
print "Saving output in folder " + "'" + blastx_output + "'"

# run blastx on all contigs
blastx_fasta_files(individual_fasta_dir, args.protein_data_base, blastx_output)





