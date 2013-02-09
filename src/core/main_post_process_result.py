"""
Created on Feb 1, 2013

@author: sw167
"""
import pickle
keep_key = set(["18S", "28S", "18S_mt", "RPB1", "RPB2", "MCM7" , "ITS"  ])  # ITS

keep_key_2 = ["ZZ"]  # ="rrna", "rna", "ribosomal"]
# 26S like.... "28S":[["28s"], ["large subunit ribosomal"], ["large subunit rrna"]],
#

remove_keyword_2 = ["mrna", "whole genome", "calmodulin", "cmd", "cytb",
                     "from patent", "cytochrome", "intron",
                     "microsatellite", "histone", "recombinant",
                     "polyketide synthase", "pks", "elongase", "ef1", "elo1",
                     "enzyme mixture", "nonsense", "mating", "biogenesis", "chaperonin",
                     "amds", "acetamidase",
                     "heat shock", "glyceraldehyde-3-phosphate",
                     "sterol c14 reductase (erg3)",
                     "faca", "rodlet", "roda",
                     "tryptophan synthase", "trpc", "sugar transporter",
                     "afls", "aflh", "aflo", "aflp", "aflr", "mfs1", "alfe", "alfg", "aflw",

                      "method of detecting of fungus",
                     "aspergillus oryzae gene",
                     "laundry methods",
                     "genes from a gene cluster",
                     "polypeptides having beta-glucosidase activity and polynucleotides encoding same",
                     "a system for detection and identification",
                     "promoter variants for expressing genes in a fungal cell",
                     "methods for increasing expression",
                     "detergent compositions",
                     "pcr system for detecting specific mold",
                     "material and method for detecting fungi"
                     "inflatum transposon restless deletion derivative",
                     "antifungal protein and usage thereof",
                     ]



EQUAL_SEP = "==============================\n"

def post_proc_1(file_name):
    infile = open(file_name, "r")
    outfile = open(file_name + ".1", "w")

    master_taxa_list = list()
    output, hit_output, taxa_line = [""] * 3
    keep_count, hit_count, remove2_count, Keep_2_count = [0] * 4
    taxa_count = 0
    master_taxa_count = [0] * 9
    print master_taxa_count
    for i, line in enumerate(infile):
#        line = line.strip()
#        if i == 100000:
#            s9s.exit(-1)
#        print line,
        if line == EQUAL_SEP:

#            for j in range(9):
# #                if keep_count == j :
#                if keep_count >= j:
#                    master_taxa_count[j] += 1

            if  keep_count > 5:  # hit_count > 1:
                taxa_count += 1
                master_taxa_list.append(taxa_line)
                if remove2_count > 0:
                    output += "Remove2:%d\n" % remove2_count

                if Keep_2_count > 0:
                    output += "Keep_key_2:%d\n" % Keep_2_count

                output += "%s\n%s\n" % (hit_output, EQUAL_SEP)
                outfile.write(output)

            output, hit_output, taxa_line = [""] * 3
            keep_count, hit_count, remove2_count, Keep_2_count = [0] * 4
        elif line.find("Count:") is not -1 and line.find("Taxa:") is not -1:
            output += line
            index = line.find("Taxa:") + 5
#            print line
            taxa_line = line[index:len(line)].strip()


        elif line.startswith("[") and line.endswith("]\n"):
            output += line
        else:
            index = line.find(":")
            header = line[0:index]
            if header in keep_key:
                output += line
                keep_count += 1
            elif header == "Remove":
                output += line
            elif header == "Hit":
                line_l = line.lower()
                keep_hit = True
                for key in remove_keyword_2:

                    if line_l.find(key) is not -1:
                        remove2_count += 1
                        keep_hit = False
                        break

                for key in keep_key_2:
                    if line_l.find(key) is not -1:
                        Keep_2_count += 1
                        keep_hit = False
                        break


                if keep_hit:
                    hit_output += line
                    hit_count += 1
    print "T_C: ", taxa_count
    print master_taxa_count
    for taxa in sorted(master_taxa_list):
        print taxa
    tempfile = open("/home/sw167/Postdoc/Project_Fungi/autoBLAST/data/Euro_Blast_All_130202.out.Taxa_67Genes", "w")
    pickle.dump(sorted(master_taxa_list), tempfile)
    outfile.write(str(taxa_count))
    outfile.close()


if __name__ == '__main__':
    file_name = "/home/sw167/Postdoc/Project_Fungi/autoBLAST/data/Euro_Blast_All_130202.out"

    post_proc_1(file_name)
