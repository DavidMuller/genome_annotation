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


# Given an id for a protein record from our Chlamy 
# fasta database, return the id of that protein.  
# (Eg: jgi|Chlre4|28893|estExt_fgenesh1_pg.C_80374
# returns 28893)
def return_id_number(protein_id):
    id_components = protein_id.split("|")
    return id_components[2]


# Given blast id for a protein record, return the id
# of that protein.  
# (Eg: gnl|BL_ORD_ID|10482 jgi|Chlre4|153454|Chlre2_kg.scaffold_58000016
# returns 153454)
def return_id_number_blast(protein_id):
    id_components = protein_id.split("|")
    return id_components[4]



# Create individual fasta files from a multi-fasta file.
# Parameters: a multi-fasta file, name for a directory
# to put all the new individual fasta files in.
def make_individual_fatsa_files(multi_fasta, fasta_dir):
    make_directory(fasta_dir)
    # Create an individual fasta file for each of our contigs
    for record in SeqIO.parse(multi_fasta,"fasta"):
        file_name = record.id
        SeqIO.write(record, fasta_dir + "/" + file_name +".fasta", "fasta")


# Create individual protein files from a multi-fasta file.
# Parameters: a multi-fasta protein file, name for a directory
# to put all the new protein fasta files in.
def make_individual_protein_files(protein_fasta, fasta_dir):
    make_directory(fasta_dir)
    # Create an individual fasta file for each of our contigs
    for record in SeqIO.parse(protein_fasta,"fasta"):
        file_name = return_id_number(record.id)
        SeqIO.write(record, fasta_dir + "/" + file_name +".fasta", "fasta")



# Run blastx on a directory contaning fasta files. Put the
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


        if contig_name == "Contig325": #our shortest contig

            print "Running blastx on " + contig_name + "..."
            blastx_cline = NcbiblastxCommandline(cmd='blastx', query=fasta_file, db=protein_db, outfmt=5, out = saved_output)
            print blastx_cline            
            stdout, stderr = blastx_cline()



# Given a blastx alignment, determine if it
# "matched well." 
# Parameters: an e value.
# Return: True, if it passes our thresholding 
# strategy.  False otherwise.   
def blastx_filter(e):
    if e < .001:       
        return True
    else:
        return False
    

# Given hsps, determine the genomic end points
# of the hit.
# Parameters: an hsp object, or a list of them.
# Return: The genomic endpoints of the hit:
# start, end.
def find_end_points(hsps):
    begin = []
    end = []
    # each hsp_data is one hit for the protein
    for hsp_data in hsps:
        begin.append(hsp_data.query_start)
        end.append(hsp_data.query_end)
    lowest_start = min(begin)
    highest_end = max(end)
    return lowest_start, highest_end



# Filter blastx output..run Exonerate on FULL LENGTH contigs
# and whatver proteins they have significant matches with.
# Parameters: blastx_output_dir--the blastx output, exonerate_out_dir--
# to put the exonerate output, contig_dir--where the individual contig
# files are, protein_dir--where the individiual protein files are.
def exonerate(blastx_output_dir, exonerate_out_dir, contig_dir, protein_dir):
    # make sure we have a directory for exonerate output
    make_directory(exonerate_out_dir)
    
    # grab all blastx output files
    all_files = glob.glob(blastx_output_dir + "/*")

    # look at the xml output files for all contigs    
    for xml_output in all_files:
        contig_number = return_only_name(xml_output)        
        opened_file = open(xml_output)
        blast_record = NCBIXML.read(opened_file)
    
        # look at all protein hits
        descriptions_alignments = zip(blast_record.descriptions, blast_record.alignments)
        for element in descriptions_alignments:
            description = element[0]            
            alignment = element[1]
            hsps = alignment.hsps
            protein_length = alignment.length
            significant = blastx_filter(description.e)
        
            # only continue analysis with proteins that pass our filter        
            if significant:
                protein_of_interest = description.title
                protein_file_id = return_id_number_blast(protein_of_interest)
                protein_file = protein_dir + "/" + protein_file_id + ".fasta"
                contig_file = contig_dir + "/" + contig_number + ".fasta"

                #save output as Contig#_ProteinID.gff (eg: Contig325_188153.gff)
                output_file = exonerate_out_dir + "/" + contig_number + "_" + protein_file_id + ".gff"

                print "Running Exonerate on",contig_number,"and",protein_of_interest
                exonerate_argument = "exonerate --model protein2genome -q " + protein_file + " -t " + contig_file + " --showtargetgff no --showsugar no --showalignment no --showvulgar yes --verbose 0 > " + output_file  
                os.system(exonerate_argument)

                    

################
# "Main":
################

# build the command-line parser
parser = argparse.ArgumentParser(description='CSE 182: Genome annotation project.')
parser.add_argument('-out', required=True, help ='g for GFF output, a for a multi_fasta file of predicted protein, c for cDNA')
parser.add_argument('-in', help='Fasta data base with all contigs')
parser.add_argument('-p', required=True, help='Protein ortholog data base file')
parser.add_argument('-m', required=False, help='Optional cDNA file')

# parse command line arguments
args = parser.parse_args()  

individual_contig_dir = "individual_contigs"
individual_protein_dir = "individual_proteins"
blastx_output = "blastx_output"
exonerate_output = "exonerate_output"


print "\nCreating individual fasta files for all our contigs:"
print "Saving them in folder " + "'" + individual_contig_dir + "' ..."

# break up all the contigs into individual contig files 
make_individual_fatsa_files(args.genomeseq_file, individual_contig_dir)



print "\nCreating individual fasta files for proteins in",args.p
print "Saving them in folder " + "'" + individual_protein_dir + "' ..."

# break up all the proteins into individual protein files 
make_individual_protein_files(args.p, individual_protein_dir)



print "\nRunning blastx on each of our contigs:"
print "Saving output in folder " + "'" + blastx_output + "'"

# run blastx on all contigs
blastx_fasta_files(individual_fasta_dir, args.p, blastx_output)



print "\nRunning Exonerate on each contig and protein sequence that matched well"
print "Saving output in folder " + "'" + exonerate_output + "'"

# filter blastx output, and run exonerate appropriately
exonerate(blastx_output, exonerate_output, individual_contig_dir, individual_protein_dir)




