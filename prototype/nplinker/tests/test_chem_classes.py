import unittest

import sys
import os

# todo: change imports when making correct package
base_path = os.path.join('..', '..', '..', 'prototype')
sys.path.append(base_path)
from nplinker.class_info.chem_classes import CanopusResults, \
    MolNetEnhancerResults


class TestCanopusResults(unittest.TestCase):
    def setUp(self):
        canopus_test_dir = "test_files"
        can_res = CanopusResults(canopus_test_dir, 'gnps_test_data')
        self._cr = can_res

    def test_running(self):
        for i, elem in enumerate([
                self._cr.spectra_classes, self._cr.spectra_classes_names,
                self._cr.spectra_classes_names_inds, self._cr.molfam_classes,
                self._cr.molfam_classes_names,
                self._cr.molfam_classes_names_inds]):
            self.assertTrue(len(elem) != 0, f"Element {i} failed to load")

    def test_can_treemap(self):
        # todo
        pass


class TestMolNetEnhancerResults(unittest.TestCase):
    def setUp(self):
        mne_test_dir = "test_files"
        mne_res = MolNetEnhancerResults(mne_test_dir)
        self._mr = mne_res

    def test_running(self):
        for i, elem in enumerate([
                self._mr.spectra2molfam, self._mr.spectra_classes_names,
                self._mr.spectra_classes_names_inds, self._mr.molfam_classes]):
            self.assertTrue(len(elem) != 0, f"Element {i} failed to load")

        test_spec_classes = self._mr.spectra_classes('130522')
        self.assertTrue(len(test_spec_classes) != 0,
                        f"func spectra_classes() not working")
        self.assertTrue(self._mr.spectra_classes('31') is None,
                        "handling 'no matches' goes wrong, MFs/spectra without"
                        " classes end up in mne results")



if __name__ == '__main__':
    unittest.main()
