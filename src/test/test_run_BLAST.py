"""
Created on Jan 30, 2013

@author: sw167
"""
import unittest
from core.run_BLAST import runBLAST
from unittest.suite import TestSuite
from core.utils.path_utils import get_data_dir
import math

class TestRunBLAST(unittest.TestCase):

    def setUp(self):
        datadir = get_data_dir()

        self.infile = datadir + "test.fasta"
        self.outfile = datadir + "test_temp.txt"
    #    args = parser.parse_args(["-h"])

        self.BLAST = runBLAST(self.infile, self.outfile)


    def tearDown(self):
        pass


    def test_get_organism(self):

        data_to_gid = {"EC": ['333384728', '380350107', ],
                       "EM": ["407312812", "333384767"],
                       }
        expected = {'333384728':"Exophiala castellanii",
                    '380350107':"Exophiala castellanii",

                    "407312812":"Exophiala mesophila",
                    "333384767":"Exophiala mesophila",
                    }
        gid_to_taxa = self.BLAST.get_organism(data_to_gid)
        self.assertEqual(expected, gid_to_taxa)


        data_to_gid = {"EC": ['333384728', '380350107', '380350105'],
                       "EM": ["407312812", "333384767"],
                       "EL": ["406857581", "407312809"]
                       }
        expected = {'333384728':"Exophiala castellanii",
                    '380350107':"Exophiala castellanii",
                    '380350105':"Exophiala castellanii",

                    "407312812":"Exophiala mesophila",
                    "333384767":"Exophiala mesophila",
                    "406857581":"Exophiala lecanii-corni",
                    "407312809":"Exophiala lecanii-corni"
                    }
        gid_to_taxa = self.BLAST.get_organism(data_to_gid)
        self.assertEqual(expected, gid_to_taxa)

    def test_create_list_organism(self):
        gid_to_taxa = {'333384728':"Exophiala castellanii",
                       '380350107':"Exophiala castellanii",
                       '380350105':"Exophiala castellanii",

                       "407312812":"Exophiala mesophila",
                       "333384767":"Exophiala mesophila",
                       "406857581":"Exophiala lecanii-corni",
                       "407312809":"Exophiala lecanii-corni"
                }
        list_organism = self.BLAST.merge_to_unique_organism(gid_to_taxa)
        expected = ["Exophiala castellanii", "Exophiala mesophila", "Exophiala lecanii-corni"]
        self.assertEqual(expected, list_organism)

        gid_to_taxa = {"1":"1", "2":"2", "3":"3", "4":"1", "5":"1", "6":"a", "7":"aa",
                "10":"uncultured 1",
                "11":"uncultured a",
                }
        list_organism = self.BLAST.merge_to_unique_organism(gid_to_taxa)
        expected = set(["1", "2", "3", "a", "aa"])
        self.assertEqual(expected, set(list_organism))
        for t in expected:
            print t

    def test_map_organism_to_data(self):

#        self.taxa_to_raw_data = self.map_organism_to_search(self.all_organisms,
#                                                       self.gid_to_organism, self.data_to_gid)
        all_taxa = ["A", "B", "C"]
        gid_to_organism = {"1":"A", "2":"A", "3":"A",
                           "11":"B", "111":"B",
                           "9":"C", "91":"C", "92":"C", "93":"C"
                           }
        data_to_gid = {"D1":["1", "11", "111"],
                      "D9":["9", "91", "92", "93"],
                      "DX":["1", "11", "2", "91"]
                      }
        result = self.BLAST.map_organism_to_data(all_taxa, gid_to_organism, data_to_gid)
        expected = {"A":set(["D1", "DX"]),
                    "B":set(["D1", "DX"]),
                    "C":set(["D9", "DX"])
                    }

        for key, value in result.iteritems():
            self.assertEqual(expected[key], set(value))


    def test_get_seq_from_name(self):
        taxa = ["fungal endophyte"]
        self.BLAST.get_seqs_from_names(taxa, None)

if __name__ == "__main__":
    print "main"

    suite = TestSuite()
    suite.addTest(TestRunBLAST("test_get_seq_from_name"))
#    suite.addTest(TestGoConnector("test_GoConnector_long"))
    unittest.TextTestRunner().run(suite)
#
#    all = range(31)
#    new_list = []
#    max = math.ceil(len(all) / 10.0)
#    print all, max
#    count = 0
#    while count < max:
#        t = all[(count * 10):((count + 1) * 10) ]
#        print t
#        new_list.append(t)
#        count += 1
#    print "=="
#    print new_list
#
#
#    data = xrange(100)
#    data2 = range(100)
#    data3 = set(range(100))
#    stmt = """
#        for i in data:
#            i in data
#       """
#
#
#    import timeit
# #    a = timeit.Timer(stmt, setup="data = xrange(100)")
# # #    print a.timeit()
# #    print a.repeat()
# #    print timeit.repeat(stmt, setup="data = range(100)")
# #    print timeit.repeat(stmt, setup="data  = set(range(100)) ")
#
#
#    print timeit.timeit('char in text', setup='text = "sample string"; char = "g"')
#    print timeit.timeit("char in text", setup="text = 'sample string'; char = 'g'")
#
# # if __name__ == '__main__':
