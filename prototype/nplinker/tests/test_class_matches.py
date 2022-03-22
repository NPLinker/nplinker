import unittest

import sys
import os
base_path = os.path.join('..', '..', '..', 'prototype')
sys.path.append(base_path)


from nplinker.class_info.class_matches import ClassMatches


class TestClassMatches(unittest.TestCase):
    def test_running(self):
        mibig_class_file = '../data/MIBiG2.0_compounds_with_AS_BGC_CF_NPC_classes.txt'
        class_matches = ClassMatches(mibig_class_file)
        for i, elem in enumerate([class_matches.class_matches,
                                  class_matches.class_matches_counts,
                                  class_matches.bgc_class_names,
                                  class_matches.chem_class_names]):
            self.assertTrue(len(elem) != 0, f"Element {i} failed to load")


if __name__ == '__main__':
    unittest.main()
