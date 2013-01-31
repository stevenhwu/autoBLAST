"""
Created on Jan 30, 2013

@author: sw167
"""
import re
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import urllib

RE_GI_PATTERN = re.compile("gi\|(\d+)\|\w+\|.*\|")
E_VALUE_THRESH = 1e-5
Entrez.email = "Ather@hotmail.com"
keep_keys = {
#    "5.8S":[["5.8S"]],
    "18S":[["18S"], ["small subunit ribosomal"], ["small subunit rRNA"]],
    "28S":[["28S"], ["large subunit ribosomal"], ["large subunit rRNA"]],
    "RPB1":[["RPB1"]],
    "RPB2":[["RPB2"], ["RPBII"]],
    "MCM7":[["MCM7"]],
    "18S_mt":[["18S", "mitochondrial"],
              ["small subunit ribosomal", "mitochondrial"],
              ["small subunit rRNA", "mitochondrial"] ],
#    "ITS":[["internal transcribed spacer", "5.8S"],
#           ["internal transcribed spacer", "28S"]],
    "ITS":[["internal transcribed spacer"], ["ITS"]]
         }

remove_keys = [
    "elongation",
    "tubulin",
    "actin",
    "cytb",
    "Whole genome shotgun",
    "rpb3",
    "tsr1",
    "cct8"]


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
def get_gi_id_from_other_seqids(id_list):

    if len(id_list) == 2:
        match = re.search(RE_GI_ID_FROM_SEQIDS, id_list[1])
        gi_id = match.group(1)
    else:
        for id in id_list:
            match = re.search(RE_GI_PATTERN, id)
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

    def __init__(self, infile, outfile):
        self.infile = infile
        self.record_index = SeqIO.index(infile, "fasta")
        self.outfile = open(outfile, "w")

    def run(self):

        self.data_to_gid = self.blast_seqs();
        self.gid_to_organism = self.get_organism(self.data_to_gid)
        self.all_organisms = self.merge_to_unique_organism(self.gid_to_organism)
        self.taxa_to_data = self.map_organism_to_data(self.all_organisms,
                                                       self.gid_to_organism, self.data_to_gid)


        self.get_seqs_from_names(self.outfile, self.all_organisms, self.taxa_to_data)
        print "==end==\n\n"

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
        for key in self.record_index:

            print key
            seq = self.record_index[key].seq
            seq_len = float(len(seq))

            result_handle = NCBIWWW.qblast("blastn", "nt", self.record_index[key].seq,
                                           hitlist_size=250, megablast=True)
            blast_record = NCBIXML.read(result_handle)
#            or, if you have lots of results (i.e. multiple query sequences):
#            blast_records = NCBIXML.parse(result_handle)

            list_gi_id = list()
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    coverage = (hsp.query_end - hsp.query_start + 1) / seq_len
                    identity = hsp.identities / float(hsp.align_length)

                    if identity > 0.88 and coverage > 0.83 and hsp.expect < E_VALUE_THRESH:
#                        print "========="
#                        print identity, coverage
#                        print alignment
#
#                        print hsp, "\n", hsp.expect, hsp.identities, hsp.positives, hsp.gaps, hsp.align_length
#                        print hsp.gaps, hsp.align_length, 1 - (hsp.gaps / float(hsp.align_length))
#
#                        print hsp.query
#                        print hsp.match
#                        print hsp.sbjct
#                        print "++++++++"
#                        print hsp.match.count("|"), len(hsp.match)
#                        print hsp.query_start, hsp.query_end
#                        print hsp.sbjct_start, hsp.sbjct_end

                        gi_id = get_gi_id(alignment.hit_id)
                        if gi_id is not None:
                            list_gi_id.append(gi_id)
            data_to_gid[key] = list_gi_id
#            break
        return data_to_gid


    def get_organism(self, data_to_gid):
        print "Pool all organisms..."
        list_id = [item for sublist in data_to_gid.values() for item in sublist]

        list_id = list(set(list_id))

        handle = Entrez.efetch(db="nucleotide", id=list_id, rettype="gb", retmode="xml")
        data = Entrez.parse(handle)

        gid_to_organism = dict()
        for i, t_dict in enumerate(data):

            gid = get_gi_id_from_other_seqids(t_dict["GBSeq_other-seqids"])
#            print(t_dict["GBSeq_other-seqids"], gid)
            gid_to_organism[gid] = t_dict["GBSeq_organism"]

        return gid_to_organism

    def merge_to_unique_organism(self, gid_to_organism):
        all_organisms = list(set(gid_to_organism.values()))
        # http://stackoverflow.com/questions/1207406/remove-items-from-a-list-while-iterating-in-python
        all_organisms[:] = [temp for temp in all_organisms if temp.find("uncultured") == -1]

        return all_organisms

    def get_seqs_from_names(self, outfile, list_of_organisms, taxa_to_raw_data):
        print "Get seq def..."
        result = dict()
#        print "List_of_organisms:\t", list_of_organisms
        for taxa in list_of_organisms:
            qterm = "%s[Organism]" % taxa

            handle = Entrez.esearch(db="nucleotide", term=qterm, rettype="gb", retmode="xml", retmax=10000)
            data = Entrez.read(handle)
            count = int(data["Count"])
#            print (data["Count"], type(data["Count"]), str(data["Count"]), int(data["Count"]))
            if count > 1:

# # FIXME: something wrong with data["count"], count == 46 but 47 seqs
                print count, taxa  # , len(data["IdList"])  , data["IdList"]
                outfile.write("Count: %d\tTaxa: %s\n" % (count, taxa))
                handle = Entrez.efetch(db="nucleotide", id=data["IdList"], rettype="gb", retmode="xml")
                seqs = Entrez.parse(handle)
                keyword_count = dict()
                for j, seq in enumerate(seqs):
                    seq_def = seq["GBSeq_definition"]

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

                    if is_print:
                        result = "Hit:%s\n" % (seq_def)
                        outfile.write(result)


                outfile.write(str(taxa_to_raw_data[taxa]))
                outfile.write("\n")


                for key, value in keyword_count.iteritems():
                    outfile.write("%s:%s\n" % (key, value))
                outfile.write("\n==============================\n\n")


    def map_organism_to_data(self, all_organisms, gid_to_organism, data_to_gid):

#        all_organisms = ['Exophiala castellanii', 'Exophiala pisciphila', 'Exophiala mesophila']
#        all_organism = {'333384777': 'Exophiala mesophila', '333384776': 'Exophiala mesophila', '380350134': 'uncultured Herpotrichiellaceae', '332903373': 'uncultured fungus', '295883724': 'uncultured soil fungus', '333384768': 'Exophiala mesophila', '359372668': 'Exophiala mesophila', '324034947': 'fungal sp. H8501', '193297429': 'Fungal endophyte sp. H114', '333384799': 'Exophiala castellanii', '305682536': 'Exophiala sp. NH1238', '407312813': 'Exophiala mesophila', '407312812': 'Exophiala mesophila', '407312811': 'Exophiala mesophila', '407312810': 'Exophiala mesophila', '305854671': 'uncultured fungus', '333384771': 'Exophiala mesophila', '34740141': 'Rhinocladiella atrovirens', '384474978': 'Exophiala castellanii', '333384724': 'Exophiala castellanii', '333384725': 'Exophiala castellanii', '333384726': 'Exophiala castellanii', '333384727': 'Exophiala castellanii', '387583275': 'uncultured soil fungus', '327387952': 'Exophiala mesophila', '380350107': 'Exophiala castellanii', '333384729': 'Exophiala castellanii', '333384766': 'Exophiala mesophila', '401836251': 'uncultured Exophiala', '299777240': 'uncultured fungus', '27464851': 'Exophiala mesophila', '116042296': 'Exophiala mesophila', '380350145': 'uncultured Herpotrichiellaceae', '333384723': 'Exophiala castellanii', '380350129': 'uncultured Herpotrichiellaceae', '90102407': 'uncultured soil fungus', '116042294': 'Exophiala mesophila', '116042295': 'Exophiala mesophila', '323404636': 'Chaetothyriales sp. M-Cre1-2', '299835224': 'Exophiala mesophila', '193297430': 'Fungal endophyte sp. H125', '189909435': 'Capronia kleinmondensis', '434860974': 'Exophiala mesophila', '333384728': 'Exophiala castellanii', '4206349': 'Exophiala pisciphila', '380350105': 'Exophiala castellanii', '379990855': 'Exophiala sp. NH512', '154082167': 'uncultured fungus', '333384767': 'Exophiala mesophila', '299835223': 'Exophiala castellanii'}
#        all_gi_id = {'S2': [u'380350107', u'380350105', u'299835223', u'384474978', u'333384728', u'333384729', u'333384727', u'333384725', u'333384726', u'299835224', u'333384724', u'333384799', u'333384723', u'434860974', u'407312810', u'327387952', u'407312812', u'407312813', u'333384767', u'333384768', u'333384766', u'333384776', u'407312811', u'116042295', u'333384771', u'333384777', u'116042296', u'116042294', u'359372668', u'4206349', u'27464851', u'380350145', u'380350129', u'401836251', u'305854671', u'379990855', u'295883724', u'34740141', u'305682536', u'323404636', u'380350134', u'189909435', u'154082167', u'387583275', u'299777240', u'193297429', u'193297430', u'90102407', u'332903373', u'324034947']}

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
        return taxa_to_data
#        print (self.taxa_to_data)
#                sys.exit(-1)
#            print("\n")

