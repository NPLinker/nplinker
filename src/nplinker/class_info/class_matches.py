# Copyright 2021 The NPLinker Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import logging
import os
from collections import Counter
from collections import defaultdict
import pandas as pd


logger = logging.getLogger(__name__)


class ClassMatches:
    """Holds all info concerning class matches (based on known bgc-structure links in MIBiG)."""

    def __init__(self, mibig_classes_file):
        """Read mibig classes file and convert into count and scoring tables.

        Args:
            mibig_classes_file: str, filepath of the MIBiG data as made by
            https://github.com/louwenjjr/mibig_classifications
        """
        self._mibig_classes_file = mibig_classes_file
        (
            self._class_matches,
            self._class_matches_counts,
            self._bgc_class_names,
            self._chem_class_names,
        ) = ({}, {}, [], [])

        if os.path.isfile(self._mibig_classes_file):
            self._read_mibig_classes()
            self._get_class_counts()
            self._get_scoring_tables()
            logger.info("Loaded MIBiG classes, and class matching tables")
        else:
            logger.warn("MIBiG classes not found, class matches not loaded!")

        pd.options.display.float_format = "{:,.3f}".format  # adjust pd formatting

        self._bigscape_mibig_conversion = {
            "PKSI": "Polyketide",
            "PKSother": "Polyketide",
            "NRPS": "NRP",
            "RiPPs": "RiPP",
            "Saccharides": "Saccharide",
            "Others": "Other",
            "Terpene": "Terpene",
            "PKS-NRP_Hybrids": "PKS-NRP_Hybrids",
        }

        self._as_conversion = {
            "NAGGN": "other",
            "NAPAA": "other",
            "RRE-containing": "bacteriocin",
            "RiPP-like": "bacteriocin",
            "cf_fatty_acid": "fatty_acid",
            "cf_putative": "other",
            "cf_saccharide": "saccharide",
            "guanidinotides": "fused",
            "lanthipeptide-class-i": "lanthipeptide",
            "lanthipeptide-class-ii": "lanthipeptide",
            "lanthipeptide-class-iii": "lanthipeptide",
            "lanthipeptide-class-iv": "lanthipeptide",
            "lanthipeptide-class-v": "lanthipeptide",
            "lantipeptide": "lanthipeptide",
            "linaridin": "lanthipeptide",
            "lipolanthine": "lanthipeptide",
            "nrps": "NRPS",
            "otherks": "hglE-KS",
            "prodigiosin": "other",
            "pyrrolidine": "other",
            "ranthipeptide": "bacteriocin",
            "redox-cofactor": "other",
            "t1pks": "T1PKS",
            "t2pks": "T2PKS",
            "t3pks": "T3PKS",
            "thioamide-NRP": "other",
            "thioamitides": "bacteriocin",
            "transatpks": "transAT-PKS",
        }

    def get_gcf_as_classes(self, gcf, cutoff=0.5):
        """Get antismash classes for a gcf if antismash class occurs in more than <cutoff> of gcf.

        Args:
            - gcf: GCF NPLinker object
            - cutoff: float - fraction of the GCF that needs to contain the class for it to be counted
        Returns:
            List of str, the antismash classes present in the gcf
        """
        # todo: move to GCF object?
        gcf_size = len(gcf.bgcs)
        unlist_all_products = [
            product for bgc in gcf.bgcs for product in bgc.product_prediction.split(".")
        ]
        sorted_as_classes = Counter(unlist_all_products).most_common()
        # keep if in more than half of bgcs?
        cutoff = 0.5
        size_cutoff = gcf_size * cutoff
        filtered_as_classes = []
        for product in sorted_as_classes:
            if product[1] >= size_cutoff:
                filtered_as_classes.append(product[0])
        return filtered_as_classes

    def convert_as_classes(self, init_as_classes: list):
        """Convert AS classes to class names that are in scoring table.

        Args:
            - init_as_classes: list of str, the initial antismash class names
        Returns:
            List of str: converted antismash classes with _as_conversion_table
        """
        as_classes = []
        for as_class in init_as_classes:
            as_conversion = self.as_conversion.get(as_class)
            if as_conversion:
                as_classes.append(as_conversion)
            else:
                as_classes.append(as_class)
        return as_classes

    def _read_mibig_classes(self):
        """Read mibig file to dict of list {chem_id: [bgc_classes, chem_classes]}.

        Returns:
            dict(str, list(list(str))) - {chem_id: [bgc_classes, chem_classes]}
        """
        classes_dict = {}
        with open(self._mibig_classes_file) as inf:
            header = inf.readline()
            for line in inf:
                elems = line.strip().split("\t")
                chem_id = elems.pop(0)
                class_base = elems.pop(0).split(",")
                classes = [cls.partition(":")[0] for cls in class_base]
                sub_classes = [cls for cls in class_base if cls.split(":")[1]]
                as_classes = elems.pop(0).split(",")

                bgc_classes = [classes, sub_classes, as_classes]
                chem_classes = [chem_cls.split("; ") for chem_cls in elems[2:]]
                classes_dict[chem_id] = [bgc_classes, chem_classes]
        self._mibig_classes = classes_dict
        # add header info
        s_h = header.strip().split("\t")

        self._bgc_class_names = ["mibig_classes"] + s_h[1:3]
        self._chem_class_names = s_h[5:]

        return self._mibig_classes

    def _get_class_counts(self):
        """Aggregate pairwise class matrices for all compounds.

        Returns:
            Recurring defaultdict of {bgc_cat: chem_cat: bgc_c: chem_c: int}
        """

        def _rec_dd():
            """Initialises a recurring defaultdict."""
            return defaultdict(_rec_dd)

        result = _rec_dd()
        # loop through each mibig compound
        for mibig_chem_id, (bgc_classes, chem_classes) in self._mibig_classes.items():
            # get all combinations of classes for this compound
            for i, bgc_cat in enumerate(self.bgc_class_names):
                init_bgc_class = bgc_classes[i]
                if not init_bgc_class or init_bgc_class == [""]:
                    continue

                bgc_class = init_bgc_class[:]  # if no exceptions, just assign classes

                # do some cleanup for mibig classes
                if bgc_cat == "mibig_classes":
                    # group pks-nrp hybrids for MIBiG classes
                    hyb_count = len(
                        [
                            1
                            for init_bgc_c in init_bgc_class
                            if any(
                                [
                                    test in init_bgc_c.lower()
                                    for test in ["nrp", "pks", "polyketide"]
                                ]
                            )
                        ]
                    )
                    if hyb_count >= 2:
                        # if hybrid, reconstruct the bgc_class
                        bgc_class = ["PKS-NRP_Hybrids"]
                        # append other classes if there are more
                        for init_bgc_c in init_bgc_class:
                            if not any(
                                [
                                    test in init_bgc_c.lower()
                                    for test in ["nrp", "pks", "polyketide"]
                                ]
                            ):
                                bgc_class.append(init_bgc_c)

                    # replace Alkaloid with Other in bgc_class
                    bgc_class = ["Other" if bgc_c == "Alkaloid" else bgc_c for bgc_c in bgc_class]

                for j, chem_cat in enumerate(self.chem_class_names):
                    chem_class = chem_classes[j]
                    if not chem_class or chem_class == [""]:
                        continue

                    for bgc_c in bgc_class:
                        for chem_c in chem_class:
                            try:
                                result[bgc_cat][chem_cat][bgc_c][chem_c] += 1
                            except TypeError:
                                result[bgc_cat][chem_cat][bgc_c][chem_c] = 1
        self._class_count_dict = result
        return result

    def _get_scoring_tables(self):
        """Makes dicts linked to pd.DataFrame that stores counts/scores.

        The resulting dataframes (tables) are stored in _class_matches and
        _class_matches_counts.

        They should be accessed from column to row so that if you want a
        count/score from BGC class -> chem class you should query the
        _class_matches as _class_matches[bgc_class_category]\
        [chem_class_category][bgc_class_name][chem_class_name]
        """
        class_matching_tables = {}
        class_matching_counts = {}  # store the counts in df/get rid of defaultdicts
        for bgc_key, bgc_chem_counts in self._class_count_dict.items():
            for chem_key, counts in bgc_chem_counts.items():
                # init entries in dict
                if bgc_key not in class_matching_tables:
                    class_matching_tables[bgc_key] = {}
                    class_matching_counts[bgc_key] = {}
                if chem_key not in class_matching_tables:
                    class_matching_tables[chem_key] = {}
                    class_matching_counts[chem_key] = {}
                # add matching tables as DataFrames
                counts_df = pd.DataFrame.from_dict(counts)
                class_matching_tables[bgc_key][chem_key] = (
                    counts_df / counts_df.sum(axis=0)
                ).fillna(0)
                class_matching_counts[bgc_key][chem_key] = counts_df.fillna(0)
                class_matching_tables[chem_key][bgc_key] = (
                    counts_df.T / counts_df.sum(axis=1)
                ).fillna(0)
                class_matching_counts[chem_key][bgc_key] = counts_df.T.fillna(0)
        self._class_matches = class_matching_tables
        self._class_matches_counts = class_matching_counts
        return class_matching_tables

    @property
    def class_matches(self):
        return self._class_matches

    @property
    def class_matches_counts(self):
        return self._class_matches_counts

    @property
    def bgc_class_names(self):
        return self._bgc_class_names

    @property
    def chem_class_names(self):
        return self._chem_class_names

    @property
    def bigscape_mibig_conversion(self):
        return self._bigscape_mibig_conversion

    @property
    def as_conversion(self):
        return self._as_conversion
