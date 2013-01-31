"""
Created on Jan 17, 2013

@author: Steven Wu


at ./src folder, run
python core/main.py infile outfile
e.g.
python core/main.py /home/sw167/Postdoc/Project_Fungi/data/test.fasta /home/sw167/Postdoc/Project_Fungi/data/test_out.txt

"""

from core.run_BLAST import runBLAST
import argparse

import sys


def main(args):

    print args
    BLAST = runBLAST(args.infile, args.outfile)
    BLAST.run()


if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="FASTA file for all seqs")
    parser.add_argument("outfile", help="Output file")
    parser.add_argument("-d", "--debug", action="count", default=0,
                        help="increase debugging level")


#    infile = "/home/sw167/Postdoc/Project_Fungi/data/test.fasta"
#    outfile = "/home/sw167/Postdoc/Project_Fungi/data/test_out.txt"
    #    args = parser.parse_args(["-h"])
#    args = parser.parse_args([infile, outfile])
    args = parser.parse_args()

    main(args)
    pass
