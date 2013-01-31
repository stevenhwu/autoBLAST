'''


#fully auto method
cd ./src
python -m unittest discover -v

 #no ".py", auto call testRunner
 #create proper test suite, which allow the testRunner to call it automatically
 #need add from * import *
cd ./src
python -m unittest -v test.main
'''
from test import *
from test.utils import *
import unittest


class TestAll(unittest.TestCase):

    CWD = os.getcwd()

if __name__ == '__main__':
    unittest.main(verbosity=2)
#    suite = unittest.TestLoader().loadTestsFromTestCase(TestAll)
#    suite = unittest.TestSuite()
#    suite.addTests(unittest.TestLoader().discover(os.getcwd()))
#    suite.addTests(unittest.TestLoader().loadTestsFromModule(TestAll))
#    unittest.TextTestRunner(verbosity=2).run(suite)
