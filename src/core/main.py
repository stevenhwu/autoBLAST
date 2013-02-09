"""
Created on Jan 17, 2013

@author: Steven Wu


at ./src folder, run
python core/main.py file_name outfile
e.g.
python core/main.py /home/sw167/Postdoc/Project_Fungi/data/test.fasta /home/sw167/Postdoc/Project_Fungi/data/test_out.txt

"""

from core.run_BLAST import runBLAST
import argparse

import sys
import pickle


def main(args):

    print args

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("file_name", help="FASTA file for all seqs")
    parser.add_argument("outfile", help="Output file")
    parser.add_argument("-d", "--debug", action="count", default=0,
                        help="increase debugging level")

#    file_name = "/home/sw167/Postdoc/Project_Fungi/autoBLAST/data/test.fasta"
#    file_name = "/home/sw167/Postdoc/Project_Fungi/autoBLAST/data/test_blastx.fasta"
    outfile = "/home/sw167/Postdoc/Project_Fungi/data/test_out.txt"
    #    args = parser.parse_args(["-h"])

    file_name = "/home/sw167/Postdoc/Project_Fungi/autoBLAST/data/Euro_Blast_All.fas"
    outfile = "/home/sw167/Postdoc/Project_Fungi/autoBLAST/data/Euro_Blast_All_130202.out"

    args = parser.parse_args([file_name, outfile])
#    args = parser.parse_args()


#    main(args)


    BLAST = runBLAST(args.file_name, args.outfile)
#    BLAST.run()
#    BLAST.custom_Update()
#    BLAST.check2()
    taxa_list_file = "/home/sw167/Postdoc/Project_Fungi/autoBLAST/data/Euro_Blast_All_130202.out.Taxa_67Genes"
    outfile = "/home/sw167/Postdoc/Project_Fungi/autoBLAST/data/Euro_Blast_All_130202.out.Seqs_67genes"

    taxa_file = open(taxa_list_file, "r")
    efetch_taxa_list = pickle.load(taxa_file)
    exist_ref_seq_list = open("/home/sw167/Postdoc/Project_Fungi/autoBLAST/data/Eurotio_alltaxa.list", "r")
#        infile = open(file_name, "r")

    efetch_taxa_set = set()
    for taxa in efetch_taxa_list:
        taxa = taxa.lower()
        efetch_taxa_set.add(taxa)

    exist_taxa = set()
    for line in exist_ref_seq_list:
        line = line.lower().strip()
        if line.startswith("rock_"):
            line = line[5:len(line)]
        line = line.replace("_", " ")
        exist_taxa.add(line)
    print len(exist_taxa)
    master_taxa_list = set(efetch_taxa_set.difference(exist_taxa))
#    print len(exist_taxa), exist_taxa
#    print len(efetch_taxa_set), efetch_taxa_set
#    print len(master_taxa_list), master_taxa_list

#    master_taxa_list = set(["cladosporium cladosporioides"])
#    BLAST.get_seqs_from_taxa_list(master_taxa_list, outfile)




    pass
# try:
#   with open('filename') as f: pass
# except IOError as e:
#   print 'Oh dear.'
#
# if not os.path.exists(filename):
#    file(filename, 'w').close()
# Alternatively:
#
# file(filename, 'w+').close()
