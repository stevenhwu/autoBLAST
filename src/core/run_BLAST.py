"""
Created on Jan 30, 2013

@author: sw167


"""
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from urllib2 import HTTPError
import math
import pickle
import re
import sys
import time
import urllib
import warnings
import os.path



MAX_BLAST = 500
RE_GI_PATTERN = re.compile("gi\|(\d+)\|\w+\|.*\|")
E_VALUE_THRESH = 1e-5
Entrez.email = "Ather@hotmail.com"
EQUAL_SEP = "==============================\n"

keep_keys = {

    "18S":[["18s"], ["small subunit ribosomal"], ["small subunit rrna"]],
    "28S":[["28s"], ["large subunit ribosomal"], ["large subunit rrna"]],
    "RPB1":[["rpb1"]],  # blastx
    "RPB2":[["rpb2"], ["rpbii"]],  # blastx
    "MCM7":[["mcm7"]],  # blastx
    "18S_mt":[["18s", "mitochondrial"],
              ["small subunit ribosomal", "mitochondrial"],
              ["small subunit rrna", "mitochondrial"] ],
    "ITS":[["internal transcribed spacer"], ["its"]]
         }
ALL_GENE_SET = set(keep_keys.keys())

program_seq = {
    "blastn" : ["SSU", "LSU", "ITS"],
    "blastx" : ["RPB", "MCM7"]
    }
remove_keys = [
    "elongation",
    "tubulin",
    "actin",
    "cytb",
    "whole genome shotgun",
    "rpb3",
    "tsr1",
    "cct8",

    "mrna", "whole genome", "calmodulin", "cmd", "cytb",
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


def get_gi_id(string):

#    string = "gi|434860974|gb|KC254027.1|"
    title_match = re.search(RE_GI_PATTERN, string)
    gi_id = title_match.group(1)
    return gi_id


def _increment_count(keyword_count, key):

    try:
        keyword_count[key] += 1
    except KeyError:
        keyword_count[key] = 1

RE_GI_ID_FROM_SEQIDS = re.compile("gi\|(\d+)")
# ['emb|HE605217.1|', 'gi|380350107']
# ['emb|HE605215.1|', 'gi|380350105']
# ['lcl|comp9774-28S', 'gb|AF050275.1|AF050275', 'gi|4206351']


def get_gi_id_from_other_seqids(id_list):

#    if len(id_list) == 2:
#        match = re.search(RE_GI_ID_FROM_SEQIDS, id_list[1])
#        gi_id = match.group(1)
#    else:
    for gid in id_list:
        match = re.search(RE_GI_ID_FROM_SEQIDS, gid)
        if match is not None:
            gi_id = match.group(1)
    try:
        return gi_id
    except:
        print id_list


def _append_to_dict_list(dicts, key, value):

    try:
        dicts[key].append(value)
    except:
        dicts[key] = [value]




class runBLAST(object):

    def __init__(self, file_name, outfile):
        self.outfile_name_prefix = outfile

        self.record_index = SeqIO.index(file_name, "fasta")
#        self.outfile = open(outfile, "w")

    def run(self):

        self.outfile = open(self.outfile_name_prefix, "w")
        self.data_to_gid = self.blast_seqs();
        tempfile = open(self.outfile_name_prefix + "_data_to_gid", "w")
        pickle.dump(self.data_to_gid, tempfile)
        tempfile.close()
#        tempfile = open(self.outfile_name_prefix + "_data_to_gid", "r")
#        self.data_to_gid = pickle.load(tempfile)


        self.gid_to_organism = self.get_organism(self.data_to_gid)
        tempfile = open(self.outfile_name_prefix + "_gid_to_organism", "w")
        pickle.dump(self.gid_to_organism, tempfile)
        tempfile.close()

        self.all_organisms = self.merge_to_unique_organism(self.gid_to_organism)
        self.taxa_to_data = self.map_organism_to_data(self.all_organisms,
                                                       self.gid_to_organism, self.data_to_gid)
        tempfile = open(self.outfile_name_prefix + "_taxa_to_data", "w")
        pickle.dump(self.taxa_to_data, tempfile)
        tempfile.close()

        self.get_seqs_from_names(self.all_organisms, self.taxa_to_data)
        self.outfile.write("\n==end==\n")
        self.outfile.close()

    def custom_Update(self):
        tempfile = open(self.outfile_name_prefix + "_data_to_gid", "r")
        self.data_to_gid = pickle.load(tempfile)
        tempfile.close()

        tempfile = open(self.outfile_name_prefix + "_gid_to_organism", "r")
        self.gid_to_organism = pickle.load(tempfile)
        tempfile.close()

        tempfile = open(self.outfile_name_prefix + "_taxa_to_data", "r")
        self.taxa_to_data = pickle.load(tempfile)
        tempfile.close()

        self.all_organisms = self.merge_to_unique_organism(self.gid_to_organism)

        self.outfile = open(self.outfile_name_prefix, "r+")
        set_all_organisms = set(self.all_organisms)
        print 'checking taxa name', len(set_all_organisms)
        taxa_count = 0
        for line in self.outfile:
            index = line.find("Taxa:")
            if index is not -1:
                taxa_name = line[index + 5:len(line)].strip()
#                print taxa_name
                set_all_organisms.remove(taxa_name)
                taxa_count += 1
        print "Removed ", taxa_count, "\tend: ", len(set_all_organisms)
        other_taxa = list(set_all_organisms)

        self.get_seqs_from_names(other_taxa, self.taxa_to_data)
        self.outfile.write("\n==end==\n")
        self.outfile.close()




    def check2(self):
        tempfile = open("/home/sw167/Postdoc/Project_Fungi/autoBLAST/data/Euro_BLAST_130201/Euro_Blast_All.out_data_to_gid", "r")
        data_to_gid = pickle.load(tempfile)
#        for key, value in data_to_gid.iteritems():
#            print len(value), "  \t", key

        tempfile = open("/home/sw167/Postdoc/Project_Fungi/autoBLAST/data/Euro_BLAST_130201/Euro_Blast_All.out_gid_to_organism", "r")
        gid_to_taxa = pickle.load(tempfile)
        for key, value in gid_to_taxa.iteritems():
            print len(value), "  \t", key, value



    def blast2(self):

 #        File = open("output"+x+".txt","w")
        fasta_string = open(self.infile).read()  # or make the names fasta1.fasta and just do open(i).read
        print(fasta_string)
        database = "nr"
        program = "blastn"
        parameters = [
         ('DATABASE', database),

          ('PROGRAM', program),
          # ('PSSM',pssm), - It is possible to use PSI-BLAST via this API?
          ('QUERY', fasta_string),

        ('CMD', 'Put'),
          ]
        query = [x for x in parameters if x[1] is not None]
        message = (urllib.urlencode(query))
        print (query)
        print(message)
        result_handle = NCBIWWW.qblast("blastn", "nr", fasta_string, hitlist_size=10)

        blast_records = NCBIXML.parse(result_handle)
        # or blast_record = NCBIXML.read(result_handle) if you only have one seq in file
        E_VALUE_THRESH = 0.001
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        print "alignment:", alignment.title
                        print "e-value:", hsp.expect


    def blast_seqs(self):
        """
            # similarity
            # 90 -> 88
            # coverage
            # 85 - >83 (query_end - query_start+1)/seq_length
        """
        print "BLAST seqs..."
        gi_id = 0
        data_to_gid = dict()
        print len(self.record_index)
        for key in self.record_index:

            program, db = "blastn", "nt"

            if key.find("RPB") is not -1 or key.find("MCM7") is not -1:
                program, db = "blastx", "nr"
            print "%s using %s..." % (key, program),
            seq = self.record_index[key].seq
            seq_len = float(len(seq))

            for i in range(50):
                try:
                    result_handle = NCBIWWW.qblast(program, db, self.record_index[key].seq,
                                           hitlist_size=250, megablast=True)
                    print "done! ",
                    break
                except:
                    print "sleep for 3s...",
                    time.sleep(3)
                    continue
            print "parsing records..."
            blast_record = NCBIXML.read(result_handle)
#            or, if you have lots of results (i.e. multiple query sequences):
#            blast_records = NCBIXML.parse(result_handle)
            tempfile = open(self.outfile_name_prefix + ".blastResult_" + key, "w")
            pickle.dump(blast_record, tempfile)
            tempfile.close()
            list_gi_id = list()
            for alignment in blast_record.alignments:

                for hsp in alignment.hsps:
                    coverage = (hsp.query_end - hsp.query_start + 1) / seq_len
                    identity = hsp.identities / float(hsp.align_length)

#                    if program == "blastx":
#                        print hsp.query
#                        print hsp.match
#                        print hsp.sbjct
#                        print "++++++++"
#                        print hsp.match.count("|"), len(hsp.match)
#                        print hsp.query_start, hsp.query_end
#                        print hsp.sbjct_start, hsp.sbjct_end


                    if identity > 0.88 and coverage > 0.88 and hsp.expect < E_VALUE_THRESH:
#                        print "========="
#                        print identity, coverage
#                        print alignment
#
#                        print hsp, "\n", hsp.expect, hsp.identities, hsp.positives, hsp.gaps, hsp.align_length
#                        print hsp.gaps, hsp.align_length, 1 - (hsp.gaps / float(hsp.align_length))
                        if program == "blastx":
                            print hsp.query
                            print hsp.match
                            print hsp.sbjct
                            print "++++++++"
                            print hsp.match.count("|"), len(hsp.match)
                            print hsp.query_start, hsp.query_end
                            print hsp.sbjct_start, hsp.sbjct_end

                        gi_id = get_gi_id(alignment.hit_id)
                        if gi_id is not None:
                            list_gi_id.append(gi_id)
            data_to_gid[key] = list_gi_id
#            break
        return data_to_gid


    def get_organism(self, data_to_gid):
        print "Pool all organisms... ",
        list_id = [item for sublist in data_to_gid.values() for item in sublist]

        list_id = list(set(list_id))
        print len(list_id)




        list_all_id = []
        max_count = math.ceil(len(list_id) / float(MAX_BLAST))
        count = 0
        while count < max_count:
            temp = list_id[(count * MAX_BLAST):((count + 1) * MAX_BLAST)]

            list_all_id.append(temp)
            count += 1

        gid_to_organism = dict()
        for ids in list_all_id:

            handle = Entrez.efetch(db="nucleotide", id=ids, rettype="gb",
                                   retmode="xml", retmax=100000)
            data = Entrez.parse(handle)

            for i, t_dict in enumerate(data):
                gid = get_gi_id_from_other_seqids(t_dict["GBSeq_other-seqids"])
    #            print(t_dict["GBSeq_other-seqids"], gid)
#                print "i:\t" + str(i)
                gid_to_organism[gid] = t_dict["GBSeq_organism"]

        return gid_to_organism

    def merge_to_unique_organism(self, gid_to_organism):
        all_organisms = list(set(gid_to_organism.values()))
        # http://stackoverflow.com/questions/1207406/remove-items-from-a-list-while-iterating-in-python
        all_organisms[:] = [temp for temp in all_organisms if temp.find("uncultured") == -1]

        return all_organisms

    def get_seqs_from_names(self, list_of_organisms, taxa_to_raw_data):
        print "Get seq def..."
        result = dict()
#        output = ""
#        print "List_of_organisms:\t", list_of_organisms
        for taxa in list_of_organisms:
            qterm = "%s[Organism]" % taxa

            handle = Entrez.esearch(db="nucleotide", term=qterm, rettype="gb",
                                    retmode="xml", retmax=100000)
            data = Entrez.read(handle)
            count = int(data["Count"])
#            print (data["Count"], type(data["Count"]), str(data["Count"]), int(data["Count"]))
            if count > 1:

                all_id = data["IdList"]
                print count, taxa  # , len(all_id)  # ,all_id_list
                try:
                    self.outfile.write("Count: %d\tTaxa: %s\n" % (count, taxa))
                    self.outfile.write("%s\n" % str(taxa_to_raw_data[taxa]))
                except:
                    print "Error: check\t", taxa


                keyword_count = dict()

                list_all_id = []
                max_count = math.ceil(len(all_id) / float(MAX_BLAST))
#                print len(all_id), max
                count = 0
                while count < max_count:
                    temp = all_id[(count * MAX_BLAST):((count + 1) * MAX_BLAST)]

                    list_all_id.append(temp)
                    count += 1

                for ids in list_all_id:
                    self._get_id_to_name_parser(keyword_count, ids)

                for key, value in keyword_count.iteritems():
                    self.outfile.write("%s:%s\n" % (key, value))
                self.outfile.write("\n==============================\n")
        # TODO: need fix/avoid this
        # Traceback (most recent call last):
        #  File "/home/sw167/Postdoc/Project_Fungi/autoBLAST/src/core/main.py", line 52, in <module>
        #    main(args)
        #  File "/home/sw167/Postdoc/Project_Fungi/autoBLAST/src/core/main.py", line 26, in main
        #    BLAST.custom_Update()
        #  File "/home/sw167/Postdoc/Project_Fungi/autoBLAST/src/core/run_BLAST.py", line 187, in custom_Update
        #    self.get_seqs_from_names(other_taxa, self.taxa_to_data)
        #  File "/home/sw167/Postdoc/Project_Fungi/autoBLAST/src/core/run_BLAST.py", line 397, in get_seqs_from_names
        #    self._get_id_to_name_parser(keyword_count, ids)
        #  File "/home/sw167/Postdoc/Project_Fungi/autoBLAST/src/core/run_BLAST.py", line 412, in _get_id_to_name_parser
        #    for j, seq in enumerate(seqs):
        #  File "/usr/lib64/python2.7/site-packages/Bio/Entrez/Parser.py", line 197, in parse
        #    text = handle.read(BLOCK)
        #  File "/usr/lib64/python2.7/socket.py", line 380, in read
        #    data = self._sock.recv(left)
        #  File "/usr/lib64/python2.7/httplib.py", line 541, in read
        #    return self._read_chunked(amt)
        #  File "/usr/lib64/python2.7/httplib.py", line 592, in _read_chunked
        #    value.append(self._safe_read(amt))
        #  File "/usr/lib64/python2.7/httplib.py", line 649, in _safe_read
        #    raise IncompleteRead(''.join(s), amt)
        # httplib.IncompleteRead: IncompleteRead(114 bytes read, 910 more expected)




    def _get_id_to_name_parser(self, keyword_count, id_list):


        output = ""
        handle = Entrez.efetch(db="nucleotide", id=id_list,
                               rettype="gb", retmode="xml", retmax=100000)
        seqs = Entrez.parse(handle)

        for j, seq in enumerate(seqs):
            seq_def = seq["GBSeq_definition"].lower()
#            print "j\t%d" % j
            is_print = True

            keep = dict()
            for key, terms in keep_keys.iteritems():

                match = [[seq_def.find(x) > -1 for x in term ] for term in terms]
                keep[key] = any([ all(x) for x in match])

            if keep["ITS"]:
                is_print = False
                _increment_count(keyword_count, "ITS")

            elif keep["18S_mt"]:
                is_print = False
                _increment_count(keyword_count, "18S_mt")
            else:
                for key in keep:
                    if keep[key]:
                        is_print = False
                        _increment_count(keyword_count, key)

            if is_print:
                for word in remove_keys:
                    index = seq_def.find(word)
                    if index is not -1:
                        is_print = False
                        _increment_count(keyword_count, "Remove")
                        break

            if is_print:
                result = "Hit:%s\n" % (seq_def)
                output += result

        self.outfile.write(output)
        return keyword_count

    def map_organism_to_data(self, all_organisms, gid_to_organism, data_to_gid):

#        all_organisms = ['Exophiala castellanii', 'Exophiala pisciphila', 'Exophiala mesophila']
#        all_organism = {'333384777': 'Exophiala mesophila', '333384776': 'Exophiala mesophila', '380350134': 'uncultured Herpotrichiellaceae', '332903373': 'uncultured fungus', '295883724': 'uncultured soil fungus', '333384768': 'Exophiala mesophila', '359372668': 'Exophiala mesophila', '324034947': 'fungal sp. H8501', '193297429': 'Fungal endophyte sp. H114', '333384799': 'Exophiala castellanii', '305682536': 'Exophiala sp. NH1238', '407312813': 'Exophiala mesophila', '407312812': 'Exophiala mesophila', '407312811': 'Exophiala mesophila', '407312810': 'Exophiala mesophila', '305854671': 'uncultured fungus', '333384771': 'Exophiala mesophila', '34740141': 'Rhinocladiella atrovirens', '384474978': 'Exophiala castellanii', '333384724': 'Exophiala castellanii', '333384725': 'Exophiala castellanii', '333384726': 'Exophiala castellanii', '333384727': 'Exophiala castellanii', '387583275': 'uncultured soil fungus', '327387952': 'Exophiala mesophila', '380350107': 'Exophiala castellanii', '333384729': 'Exophiala castellanii', '333384766': 'Exophiala mesophila', '401836251': 'uncultured Exophiala', '299777240': 'uncultured fungus', '27464851': 'Exophiala mesophila', '116042296': 'Exophiala mesophila', '380350145': 'uncultured Herpotrichiellaceae', '333384723': 'Exophiala castellanii', '380350129': 'uncultured Herpotrichiellaceae', '90102407': 'uncultured soil fungus', '116042294': 'Exophiala mesophila', '116042295': 'Exophiala mesophila', '323404636': 'Chaetothyriales sp. M-Cre1-2', '299835224': 'Exophiala mesophila', '193297430': 'Fungal endophyte sp. H125', '189909435': 'Capronia kleinmondensis', '434860974': 'Exophiala mesophila', '333384728': 'Exophiala castellanii', '4206349': 'Exophiala pisciphila', '380350105': 'Exophiala castellanii', '379990855': 'Exophiala sp. NH512', '154082167': 'uncultured fungus', '333384767': 'Exophiala mesophila', '299835223': 'Exophiala castellanii'}
#        all_gi_id = {'S2': [u'380350107', u'380350105', u'299835223', u'384474978', u'333384728', u'333384729', u'333384727', u'333384725', u'333384726', u'299835224', u'333384724', u'333384799', u'333384723', u'434860974', u'407312810', u'327387952', u'407312812', u'407312813', u'333384767', u'333384768', u'333384766', u'333384776', u'407312811', u'116042295', u'333384771', u'333384777', u'116042296', u'116042294', u'359372668', u'4206349', u'27464851', u'380350145', u'380350129', u'401836251', u'305854671', u'379990855', u'295883724', u'34740141', u'305682536', u'323404636', u'380350134', u'189909435', u'154082167', u'387583275', u'299777240', u'193297429', u'193297430', u'90102407', u'332903373', u'324034947']}
        print "maping organism to data...",
        match_taxa_to_gid = dict()
        for taxa in all_organisms:
            for key, value in gid_to_organism.iteritems():  # all_orgamism
                if value == taxa:
                    _append_to_dict_list(match_taxa_to_gid, taxa, key)

        taxa_to_data = dict()
        for taxa in all_organisms:
            for ids in match_taxa_to_gid[taxa]:
                for key, value in data_to_gid.iteritems():
                    if ids in value:
                        _append_to_dict_list(taxa_to_data, taxa, key)
#        print(self.taxa_to_data)
        for key in taxa_to_data:
            taxa_to_data[key] = list(set(taxa_to_data[key]))
        print "Done!"
        return taxa_to_data


    def get_seqs_from_taxa_list(self, master_taxa_list, outfile_name):


        self.outfile = open(outfile_name, "w")
        self.out_temp = open("/home/sw167/Postdoc/Project_Fungi/autoBLAST/data/aa_Euro_Blast_All_130202.out__possibleStrainName", "w")
        output = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("taxa", "strain", "gene_key", "gid", "gene_length", "gene_seq", "gene_def")
        self.outfile.write(output)
#        output = ""
#        print "List_of_organisms:\t", list_of_organisms
        for taxa in master_taxa_list:
            print "taxa: %s. " % taxa,
            file_path = "/home/sw167/Postdoc/Project_Fungi/autoBLAST/data/Euro_Blast_All_130202.out.BLAST_seqs." + taxa
            if os.path.exists(file_path):
                print "load file....",
                tempfile = open(file_path, "r")
                keep_seqs_list = pickle.load(tempfile)
                tempfile.close()
            else:
                qterm = "%s[Organism]" % taxa
                handle = Entrez.esearch(db="nucleotide", term=qterm, rettype="gb",
                                        retmode="xml", retmax=100000)
                data = Entrez.read(handle)

                all_id = data["IdList"]
                print "efetch.. %d" % len(all_id)
                list_all_id = []
                max_count = math.ceil(len(all_id) / float(MAX_BLAST))

                count = 0
                while count < max_count:
                    temp = all_id[(count * MAX_BLAST):((count + 1) * MAX_BLAST)]
                    list_all_id.append(temp)
                    count += 1

                keep_seqs_list = dict()
                for ids in list_all_id:
                    keep_seqs_list = self.get_seqs_from_ids(taxa, keep_seqs_list, ids)
#            print keep_seqs_list
                tempfile = open(file_path, "w")
                pickle.dump(keep_seqs_list, tempfile)
                tempfile.close()


            keep_key = self.seq_selector(keep_seqs_list)
            for key in keep_key:

                v = keep_seqs_list[key]
#                print v
                output = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (v["taxa"], v["strain"], v["gene_key"], v["gid"], v["gene_length"], v["gene_seq"], v["gene_def"])
                self.outfile.write(output)
            self.outfile.flush()
            print "==========\n"
#            print "FINAL_KEY", keep_key, "\n\n"
#        sys.exit(-1)
        self.out_temp.close()
        self.outfile.close()


    def seq_selector(self, seqs_list):

        # # select up to 3 lagrest non-overlapping subset

        RE_KEY_TERMS = [
                        re.compile("(?P<g1>cbs)\s*[:#]?\s*(?P<g2>\d+\.?\d+)", re.IGNORECASE),
                       re.compile("(?P<g1>aftol)\s*(-id)?\s*(?P<g2>\d+)", re.IGNORECASE) ,
                       re.compile("(?P<g1>atcc)\s*[:#]?\s*(?P<g2>[\d]+)", re.IGNORECASE) ,
                       re.compile("(?P<g1>atcc mya)\-?\s*[:#]?\s*(?P<g2>[\d]+)", re.IGNORECASE) ,
                       re.compile("(?P<g1>athum)\s*[:#]?\s*(?P<g2>\d+)", re.IGNORECASE) ,
                       re.compile("(?P<g1>daom)\s*[:#]?\s*(?P<g2>\d+)", re.IGNORECASE) ,
                       re.compile("(?P<g1>dto)\s*[:#]?\s*(?P<g2>\d+)", re.IGNORECASE) ,
                       re.compile("(?P<g1>nrrl\s*:?_?a?-?)\s*[_:#-]?\s*(?P<g2>[\d]+)", re.IGNORECASE) ,

                       ]


        strain_count2 = dict()
        all_group_index_set = set()
        all_possible_genes = set()
        all_strain_genes = set()

        for key, seq_info in seqs_list.iteritems():
#            print key, seq_info["gene_key"], seq_info["gene_def"]
            gene_def = seq_info["gene_def"]
            gene_key = seq_info["gene_key"]
            all_possible_genes.update([gene_key])

            is_match = False
            for pattern in RE_KEY_TERMS:

                match = pattern.search(gene_def)
                if match:
                    is_match = True
                    all_strain_genes.update([gene_key])
#                    print (match.group(1), match.group(2))
                    group_index = "%s %s" % (match.group("g1"), match.group("g2"))
#                    print group_index
                    all_group_index_set.update([group_index])
                    seqs_list[key]["strain"] = group_index
                    try:
                        strain_count2[group_index]["count"] += 1
                    except:
                        strain_count2[group_index] = dict()
                        strain_count2[group_index]["count"] = 1
                        strain_count2[group_index]["genes"] = list()
                        strain_count2[group_index]["gid"] = list()
                    finally:
                        strain_count2[group_index]["genes"].append(gene_key)
                        strain_count2[group_index]["gid"].append(key)
                    break

            if not is_match:
                check_list = ["cbs", "aftol", "atcc", "athum", "daom", "dto", "nrrl"]
                check_result = [True for check in check_list if gene_def.find(check) is not -1]
                if check_result:
                    warnings.warn("unmatched seq:\n %s" % gene_def)
                    print check_result
                    sys.exit(-1)

            if gene_def.find("strain") is not -1:
                index = gene_def.find("strain")
                output = "possible strain key: %s...%s\n" % (gene_def[0:30], gene_def[index:(index + 30)])
                self.out_temp.write(output)


        keep_gid = list()
        target_set = set()
        if strain_count2 != dict():
            max_count = max([v["count"] for v in strain_count2.values()])

            hit_max = dict()
            for k, v in strain_count2.iteritems():
                if v["count"] is max_count:
                    hit_max[k] = strain_count2[k]

            to_add_gene = hit_max.popitem()
            keep_gid.extend(to_add_gene[1]["gid"])
            target_set.update(to_add_gene[1]["genes"])
            print keep_gid, target_set,
            while target_set != all_strain_genes:

                missing_set = all_strain_genes.difference(target_set)

                temp_dict = self._get_max_intersection(missing_set, strain_count2)
                add_dict = self._get_max_intersection(target_set, temp_dict)
    #            print "========"
    #            for k, v in temp_dict.iteritems():
    #                print k, v
    #            print "=="
    #            for k, v in add_dict.iteritems():
    #                print k, v

                to_add_gene = add_dict.popitem()
                keep_gid.extend(to_add_gene[1]["gid"])
                target_set.update(to_add_gene[1]["genes"])

    #            keep_gid =
    #        target_set = set(to_add_gene[1]["genes"])

            print " to ", keep_gid, target_set

        if all_possible_genes != all_strain_genes:
#            while target_set != all_possible_genes:

            missing_set = all_possible_genes.difference(target_set)
            print "single:", missing_set
            for key, seq_info in seqs_list.iteritems():
                if seq_info["gene_key"] in missing_set:
                    missing_set.remove(seq_info["gene_key"])
                    keep_gid.append(seq_info["gid"])
#                    print seq_info
#                    print missing_set

        return keep_gid

# CBS ###.###
# CBS ####
#AFTOL ###
#AFTOL-id ###

    def _get_max_intersection(self, target_set, temp_dict):
        max_count = 0
        target_dict = dict()
        for k, v in temp_dict.iteritems():
            counter = len(target_set.intersection(v["genes"]))
            if counter > max_count:
                target_dict = dict()
                target_dict[k] = v
                max_count = counter
            elif counter == max_count:
                target_dict[k] = v
        return target_dict


    def get_seqs_from_ids(self, taxa, keep_seqs_list, id_list):


#        print "efetch taxa: %s\tlength:%d" % (taxa, len(id_list))

        for _ in range(10):
            try:
                handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="gb", retmode="xml", retmax=100000)
                seqs = Entrez.parse(handle)


                break
            except HTTPError as e:
                print "sleep for 5s.....", e
                time.sleep(5)
                continue

        for j, seq in enumerate(seqs):
            is_keep = False
            seq_def = seq["GBSeq_definition"].lower()
            gene_key = ""

            keep = dict()
            for key, terms in keep_keys.iteritems():

                match = [[seq_def.find(x) > -1 for x in term ] for term in terms]
                keep[key] = any([ all(x) for x in match])

            if keep["ITS"]:
                is_keep = True
                gene_key = "ITS"

            elif keep["18S_mt"]:
                is_keep = True
                gene_key = "18S_mt"
            else:
                for key in keep:
                    if keep[key]:
                        is_keep = True
                        gene_key = key
#                        _increment_count(keyword_count, key)
            if is_keep:
                gene_seq = seq["GBSeq_sequence"].upper()
                gid = get_gi_id_from_other_seqids(seq["GBSeq_other-seqids"])
                keep_seqs_list[gid] = {"taxa": taxa,
                                       "strain": "",
                                       "gene_key": gene_key,
                                       "gene_seq": gene_seq,
                                       "gid": gid,
                                       "gene_def": seq_def,
                                       "gene_length": len(gene_seq)
                                       }

#                print keep_seqs_list[gid]["gene_key"], keep_seqs_list[gid]["gene_def"]
        return keep_seqs_list



