"""
Created on Jan 17, 2013

@author: Steven Wu
"""

from Bio import Entrez, SeqIO, Blast
from Bio.Blast import NCBIWWW, NCBIXML
from Bio._py3k import _as_bytes
from core.run_BLAST import runBLAST
import argparse
import re
import sys
import urllib
import warnings

def main(args):

    print args
    BLAST = runBLAST(args.infile, args.outfile)

    BLAST.run()
    dict_id = {
               "S2": [u'333384728', u'380350107', u'380350105', u'299835223'],
               "S3": [u'384474978', u'333384728', u'333384729', u'333384727', u'333384725']

               }
#    BLAST.get_organism(dict_id)

#    list_of_organisms = ['Ascomycota sp. H31', 'Exophiala lecanii-corni', 'Chaetothyriales sp. M-Cre1-2', 'Fungal endophyte sp. H114', 'fungal sp. GMG_C6', 'Fungal endophyte sp. H125', 'uncultured Exophiala', 'uncultured Ascomycota', 'Exophiala sp. NH1238', 'Exophiala pisciphila', 'Capronia kleinmondensis', 'Exophiala castellanii', 'Exophiala mesophila', 'uncultured soil fungus', 'Exophiala sp. NH512', 'uncultured fungus', 'Rhinocladiella atrovirens', 'uncultured Herpotrichiellaceae', 'fungal sp. H8501']

#    b = ['Exophiala+lecanii-corni', 'uncultured+Ascomycota', 'Capronia+kleinmondensis', 'uncultured+soil+fungus', 'Ascomycota+sp.+H31', 'Chaetothyriales+sp.+M-Cre1-2', 'fungal+sp.+GMG_C6', 'uncultured+Exophiala', 'Fungal+endophyte+sp.+H114', 'fungal+sp.+H8501', 'Homo+sapiens', 'Pongo+abelii', 'uncultured+Herpotrichiellaceae', 'Fungal+endophyte+sp.+H125', 'Pan+troglodytes', 'Exophiala+sp.+NH1238', 'Exophiala+sp.+NH512', 'uncultured+fungus', 'Exophiala+mesophila', 'Exophiala+pisciphila', 'Exophiala+castellanii', 'Rhinocladiella+atrovirens']
#    BLAST.get_seqs_from_names(list_of_organisms)
    # # TDO: filter
    # # del no_of_seq == 1
    # Taxa:Exophiala castellanii - ITS:9, 16S:10, 18S:2 + G_A + G_B
#




    sys.exit(-1)
#    BLAST.get_organism(["407312812", "407312813"])
#    print "p2"
#    BLAST.get_organism()


if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="FASTA file for all seqs")
    parser.add_argument("outfile", help="Output file")
    parser.add_argument("-d", "--debug", action="count", default=0,
                        help="increase debugging level")


    infile = "/home/sw167/Postdoc/Project_Fungi/data/test.fasta"
    outfile = "/home/sw167/Postdoc/Project_Fungi/data/test_out.txt"
    #    args = parser.parse_args(["-h"])
    args = parser.parse_args([infile, outfile])
#    args = parser.parse_args()

    main(args)
    pass
