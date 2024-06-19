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
import glob
import logging
import os
from collections import Counter
from canopus import Canopus
from canopus.classifications_to_gnps import analyse_canopus


logger = logging.getLogger(__name__)


# load Chem_class_predictions (canopus, molnetenhancer are loaded)
# for canopus, check if results can be converted with canopus_treemap
# otherwise use the pre-existing output of canopus
class ChemClassPredictions:
    """Class for storing results for chemical class predictions of spectra.

    Currently, CANOPUS and MolNetEnhancer results are loaded
    """

    def __init__(self, canopus_dir, mne_dir, gnps_dir):
        """Load classes with CanopusResults, MolNetEnhancerResults.

        Args:
            canopus_dir: str, canopus_dir found in root_dir of nplinker project
            mne_dir: str, mne_dir found in root_dir of nplinker project
            gnps_dir: str, root dir where all gnps info is found
        """
        self._canopus = CanopusResults(canopus_dir, gnps_dir)
        self._molnetenhancer = MolNetEnhancerResults(mne_dir)

        class_predict_options = []
        if self._canopus.spectra_classes:
            class_predict_options.append("canopus")
        if self._molnetenhancer.spectra2mf:
            class_predict_options.append("molnetenhancer")
        if class_predict_options:
            class_predict_options = ["mix", "main"] + class_predict_options
        self._class_predict_options = class_predict_options

    @property
    def canopus(self):
        return self._canopus

    @property
    def molnetenhancer(self):
        return self._molnetenhancer

    @property
    def class_predict_options(self):
        """The available class predictions."""
        return self._class_predict_options


class CanopusResults:
    """Class for storing canopus results.

    The results from the canopus dir are read and combined with the MN from GNPS
    using canopus_treemap: github.com/louwenjjr/canopus_treemap/tree/master/canopus
    This creates the two files that are read for the spectra and mfs:
        -cluster_index_classifications.txt
        -component_index_classifications.txt

    If canopus_treemap somehow fails (like for too old gnps results), the canopus
    output is read as is from canopus_summary.tsv
    """

    def __init__(self, canopus_dir, gnps_dir):
        """Read the class info from root_dir/canopus.

        Args:
            canopus_dir: str, canopus_dir found in root_dir of nplinker project
            gnps_dir: str, root dir where all gnps info is found

        If possible, convert canopus output with canopus_treemap. Otherwise
        canopus output is used directly. If converting fails, probably the MN
        version is too old, and canopus output is read directly
        """
        self._canopus_dir = canopus_dir
        self._gnps_dir = gnps_dir
        self._mf_classes, self._mf_classes_names, self._mf_classes_names_inds = (
            None,
            None,
            None,
        )
        self._spectra_classes, self._spectra_classes_names, self._spectra_classes_names_inds = (
            None,
            None,
            None,
        )
        if os.path.isdir(self._canopus_dir):
            self._read_all_classes()
        else:
            logger.info(
                f"No CANOPUS results present at {self._canopus_dir}. "
                f"(set run_canopus=true in the .toml to run CANOPUS)"
            )

    def _read_all_classes(self):
        """Wrapper to read all canopus output and store them in the object."""
        ci_file = False  # find cluster_index_classifications in canopus_dir
        compi_file = False  # same for component_index_classifications
        canopus_files = glob.glob(os.path.join(self._canopus_dir, "*"))
        possible_ci_file = [
            canopus_file
            for canopus_file in canopus_files
            if canopus_file.endswith("cluster_index_classifications.txt")
        ]
        possible_compi_file = [
            canopus_file
            for canopus_file in canopus_files
            if canopus_file.endswith("component_index_classifications.txt")
        ]
        if possible_ci_file:
            ci_file = possible_ci_file[0]
        if possible_compi_file:
            compi_file = possible_compi_file[0]

        if not os.path.isfile(ci_file):
            logger.info("Converting canopus output using canopus_treemap")
            try:
                analyse_canopus(self._canopus_dir, self._gnps_dir, self._canopus_dir)
            except FileNotFoundError as er1:
                logger.warning(f"{er1}. CANOPUS output is missing in {self._canopus_dir}")
            except Exception as er2:
                logger.warning(
                    f"canopus_treemap failed with: {er2}. Probably the MN "
                    f"version from GNPS is too old. Will attempt to read "
                    f"directly from canopus_dir"
                )
            else:
                # find processed canopus files again
                canopus_files = glob.glob(os.path.join(self._canopus_dir, "*"))
                possible_ci_file = [
                    canopus_file
                    for canopus_file in canopus_files
                    if canopus_file.endswith("cluster_index_classifications.txt")
                ]
                possible_compi_file = [
                    canopus_file
                    for canopus_file in canopus_files
                    if canopus_file.endswith("component_index_classifications.txt")
                ]
                if possible_ci_file:
                    ci_file = possible_ci_file[0]
                if possible_compi_file:
                    compi_file = possible_compi_file[0]

        if os.path.isfile(ci_file):
            spectra_classes_names, spectra_classes = self._read_spectra_classes(ci_file)

            if os.path.isfile(compi_file):
                mf_classes_names, mf_classes = self._read_mf_classes(compi_file)
                self._mf_classes = mf_classes
                self._mf_classes_names = mf_classes_names
                self._mf_classes_names_inds = {elem: i for i, elem in enumerate(mf_classes_names)}
        else:
            # use canopus output correctly (only for spectra)
            logger.info(
                "Attempting to read spectra classes directly from "
                "canopus_dir (canopus_summary.tsv)"
            )
            spectra_classes_names, spectra_classes = self._read_spectra_classes_directly()
            # mfs have to be added later with info about mf <- spectra
            # this happens with transfer_spec_classes_to_mfs() in loader.py

        self._spectra_classes = spectra_classes
        self._spectra_classes_names = spectra_classes_names
        self._spectra_classes_names_inds = {elem: i for i, elem in enumerate(spectra_classes_names)}

    def _read_spectra_classes(self, input_file):
        """Read canopus classes for spectra, return classes_names, classes.

        Args:
            input_file: str, cluster_index_classifications.txt
        Returns:
            Tuple of:
            - ci_classes_names: list of str - the names of each different level
            - ci_classes: dict of {str, lists of tuple(str, float)} - per spectrum index (key) the classes for each level
                where each level is a list of (class_name, score) sorted on best choice so index 0 is the best
                class prediction for a level. When no class is present, instead of Tuple it will be None for that level.
        """
        ci_classes = {}  # make a dict {ci: [[(class,score)]]}
        ci_classes_header = None
        ci_classes_names = []

        if os.path.isfile(input_file):
            logger.info(f"reading canopus results for spectra from {input_file}")
            with open(input_file) as inf:
                ci_classes_header = inf.readline().strip().split("\t")
                for line in inf:
                    line = line.strip("\n").split("\t")
                    classes_list = []
                    for lvl in line[3:]:
                        lvl_list = []
                        for l_class in lvl.split("; "):
                            if l_class:
                                l_class = l_class.split(":")
                                c_tup = tuple([l_class[0], float(l_class[1])])
                            else:
                                c_tup = None  # default value for class value
                            lvl_list.append(c_tup)
                        classes_list.append(lvl_list)
                    ci_classes[line[1]] = classes_list
        else:
            logger.warn(
                "could not load cluster_index_classifications.txt; missing from canopus_dir"
            )

        if ci_classes_header:
            #  todo: rename the output from the canopus script directly
            ci_classes_names = [f"cf_{elem}" for elem in ci_classes_header[3:-3]] + [
                f"npc_{elem}" for elem in ci_classes_header[-3:]
            ]
        return ci_classes_names, ci_classes

    def _read_spectra_classes_directly(self):
        """Read canopus classes directly from canopus output and save output.

        Returns:
            Tuple of:
            - can_classes_names: list of str - the names of each different level
            - can_classes: dict of {str, lists of tuple(str, float/None)} - per spectrum index (key) the classes for each level
                where each level is a list of (class_name, score) sorted on best choice so index 0 is the best
                class prediction for a level. When no class is present, instead of Tuple it will be None for that level.
        """
        can_classes = {}  # make a dict {spec: [[(class,score)]]}
        can_classes_header = None
        can_classes_names = []

        input_file = os.path.join(self._canopus_dir, "canopus_summary.tsv")
        sep = "\t"
        if os.path.isfile(input_file):
            logger.info(f"reading canopus results for spectra directly from " f"{input_file}")
            with open(input_file) as inf:
                can_classes_header = inf.readline().strip().split(sep)
                # get [4:8] - level 5 subclass        class   superclass
                starti = 4
                stopi = 8
                for line in inf:
                    line = line.strip("\n").split(sep)
                    # make tuples with scores as 1
                    classes_list = [[(cls, 1)] if cls else [None] for cls in line[starti:stopi]]
                    # assuming that all MGFs have ScanNumber headers..
                    spec_index = line[0].split("ScanNumber")[-1]
                    classes_list.reverse()  # reverse so superclass is first
                    can_classes[spec_index] = classes_list

            # use canopus_treemap to produce NPClassifier classes
            # TODO: probably change when sirius v5 comes out
            logger.info("Using canopus_treemap to get NPC classes")
            canopus_workspace = Canopus(sirius=self._canopus_dir)
            npc_file = os.path.join(self._canopus_dir, "npc_summary.tsv")
            canopus_workspace.npcSummary().to_csv(npc_file, sep=sep)
            with open(npc_file) as inf:
                npc_can_classes_header = inf.readline().strip().split(sep)
                # get [2:7] - name,directoryName,pathway,pathwayProbability,superclass,superclassProbability,class,classProbability,ClassyFirePrediction,ClassyFirePredictionProbability
                npc_starti = 2
                npc_stopi = 7
                for line in inf:
                    line = line.strip("\n").split(sep)
                    # make tuples with scores
                    npc_classes_list = [
                        [(cls, float(c_score))] if cls != "N/A" else [None]
                        for cls, c_score in zip(
                            line[npc_starti:npc_stopi:2], line[npc_starti + 1 : npc_stopi + 1 : 2]
                        )
                    ]
                    # assuming that all MGFs have ScanNumber headers..
                    spec_index = line[1].split("ScanNumber")[-1]
                    can_classes[spec_index].extend(npc_classes_list)
        else:
            logger.warn(
                f"could not load canopus results for spectra; "
                f"{input_file} missing from canopus_dir"
            )

        if can_classes_header:
            # important that the names here match with MIBiG
            can_classes_names = [f"cf_{elem}" for elem in can_classes_header[starti:stopi]]
            can_classes_names.reverse()  # reverse so superclass is first
            npc_can_classes_names = [
                f"npc_{elem}" for elem in npc_can_classes_header[npc_starti:npc_stopi:2]
            ]
            can_classes_names.extend(npc_can_classes_names)

            # make sure that non-existing NPC predictions are added as [None]
            # to all spectra - some spectra where not included in npc_summary
            len_classes_names = len(can_classes_names)
            for spec_index, classes_list in can_classes.items():
                diff = len_classes_names - len(classes_list)
                if diff > 0:
                    for _ in range(diff):
                        classes_list.append([None])

            # save output, next time clusterindex_classifications.txt is read
            # save just the info that we need in the format that is parsed in
            # read_spectra_classes()
            output_file = os.path.join(self._canopus_dir, "cluster_index_classifications.txt")
            with open(output_file, "w") as outf:
                header = ["componentindex", "cluster index", "formula"] + [
                    name.partition("_")[-1] for name in can_classes_names
                ]
                outf.write("\t".join(header) + "\n")
                for spec, classes in can_classes.items():
                    formatted_classes = []
                    for lvl in classes:
                        lvl_classes = []
                        for tup in lvl:
                            if not tup:
                                lvl_classes.append("")
                            else:
                                lvl_classes.append(f"{tup[0]}:{tup[1]:.3f}")
                        formatted_classes.append("; ".join(lvl_classes))
                    output_l = ["-", spec, "-"] + formatted_classes
                    outf.write("\t".join(output_l) + "\n")
        return can_classes_names, can_classes

    def _read_mf_classes(self, input_file):
        """Read canopus classes for mfs, return classes_names, classes.

        Args:
            input_file: str, component_index_classifications.txt
        Returns:
            Tuple of:
            - compi_classes_names: list of str - the names of each different level
            - compi_classes: dict of {str: lists of tuple(str, float)} - per mf index (key) the classes for each level
                where each level is a list of (class_name, fraction) sorted on best choice so index 0 is the best
                class prediction for a level. When no class is present, instead of Tuple it will be None for that level.
        """
        compi_classes = {}  # make a dict {compi: [[(class,score)]]}
        compi_classes_header = None
        compi_classes_names = []

        if os.path.isfile(input_file):
            with open(input_file) as inf:
                compi_classes_header = inf.readline().strip().split("\t")
                for line in inf:
                    line = line.strip("\n").split("\t")
                    classes_list = []
                    for lvl in line[2:]:
                        lvl_list = []
                        for l_class in lvl.split("; "):
                            if l_class:
                                l_class = l_class.split(":")
                                c_tup = tuple([l_class[0], float(l_class[1])])
                            else:
                                c_tup = None  # default value for class value
                            lvl_list.append(c_tup)
                        classes_list.append(lvl_list)
                    compi_classes[line[0]] = classes_list
        else:
            logger.warn(
                "could not load component_index_classifications.txt; missing from canopus_dir"
            )

        if compi_classes_header:
            #  todo: rename the output from the canopus script directly
            compi_classes_names = [f"cf_{elem}" for elem in compi_classes_header[2:-3]] + [
                f"npc_{elem}" for elem in compi_classes_header[-3:]
            ]
        return compi_classes_names, compi_classes

    def transfer_spec_classes_to_mfs(self, mfs, fraction_cutoff=0.0):
        """Set _mf_classes(_names) from spectra_classes and return classes.

        This can be used in the _loader to get mf classes when the GNPS MN
        version is too old and canopus_treemap fails to work directly.

        Args:
            mfs: list of MolecularFamily from the NPLinker space
            fraction_cutoff: float, cut-off for the fraction of class terms
                needed to be included in the mf
        Returns:
            dict of {str: lists of tuple(str, float)} - per mf (key) the classes for each level
                where each level is a list of (class_name, fraction) sorted on best choice so index 0 is the best
                class prediction for a level. When no class is present, instead of Tuple it will be None for that level.
        """
        self._mf_classes_names = self._spectra_classes_names
        self._mf_classes_names_inds = self._spectra_classes_names_inds
        mf_classes = {}

        for mf in mfs:
            fid = mf.id  # the key
            spectra = mf.spectra
            # if singleton family, format like 'fid_spectrum-id'
            if fid.startswith("singleton-"):
                spec_id = spectra[0].id
                fid += f"_{spec_id}"
            len_mf = len(spectra)

            classes_per_spectra = []
            for spec in spectra:
                spec_classes = self.spectra_classes.get(spec.id)
                if spec_classes:  # account for spectra without prediction
                    classes_per_spectra.append(spec_classes)

            if not classes_per_spectra:
                continue  # no spectra with classes for this mf

            sorted_classes = []
            for i, class_level in enumerate(self._mf_classes_names):
                # 1. aggregate classes from all spectra for this class level
                classes_cur_level = []
                for spec_classes in classes_per_spectra:
                    try:
                        for class_tup in spec_classes[i]:
                            if class_tup:
                                classes_cur_level.append(class_tup[0])
                    except IndexError:
                        print(self._mf_classes_names)
                        print(i, class_level)
                        print(classes_per_spectra)
                        print(spec_classes)
                        print(spectra)
                # 2. count the instances of each class in the MF per class lvl
                counts_cur_level = Counter(classes_cur_level)
                # 3. calculate fraction and sort high to low, filter out Nones
                fraction_tups = sorted(
                    (
                        (cls, count / len_mf)
                        for cls, count in counts_cur_level.most_common()
                        if count / len_mf >= fraction_cutoff
                    ),
                    key=lambda x: x[1],
                    reverse=True,
                )
                if not fraction_tups:
                    fraction_tups = [None]
                sorted_classes.append(fraction_tups)
            mf_classes[fid] = sorted_classes

        self._mf_classes = mf_classes
        return mf_classes

    def show(self, objects):
        """Show a table of predicted chemical compound classes for spectrum/MF.

        Args:
              objects: list of Spectrum or MolecularFamily objects
        Returns:
            pandas.DataFrame of objects (rows) vs classes (columns)
        """
        pass

    @property
    def spectra_classes(self):
        return self._spectra_classes

    @property
    def spectra_classes_names(self):
        return self._spectra_classes_names

    @property
    def spectra_classes_names_inds(self):
        return self._spectra_classes_names_inds

    @property
    def mf_classes(self):
        return self._mf_classes

    @property
    def mf_classes_names(self):
        return self._mf_classes_names

    @property
    def mf_classes_names_inds(self):
        return self._mf_classes_names_inds


class MolNetEnhancerResults:
    """Class for storing MolNetEnhancer results.

    The input file for ClassyFire results is read from the molnetenhancer directory:
        - ClassyFireResults_Network.txt
    """

    def __init__(self, mne_dir):
        """Read the class info from file in root_dir/molnetenhancer/.

        Args:
            mne_dir: str, mne_dir found in root_dir of nplinker project
        """
        cf_classes_names, mf_classes, spectra2mf = self._read_cf_classes(mne_dir)
        self._spectra2mf = spectra2mf
        self._mf_classes = mf_classes
        self._spectra_classes_names = cf_classes_names  # if NPC gets implemented, add here
        self._spectra_classes_names_inds = {elem: i for i, elem in enumerate(cf_classes_names)}

    def _read_cf_classes(self, mne_dir):
        r"""Read ClassyFireResults_Network.txt in molnetenhancer dir.

        Args:
            mne_dir: str, mne_dir found in root_dir of nplinker project
        Returns:
            tuple of:
            -list of str - names of the classes in order
            -dict of {str: [(str, float)]} - linking mfs to (classes, scores) in order of names,
                singleton families are denoted with S[\d]+
            -dict of {str:str} - linking spectra to mfs
        """
        columns = []
        mne_component_dict = {}
        mne_cluster2component = {}
        # flexible finding of CF results
        input_file = "not_found"
        wanted_file = "ClassyFireResults_Network.txt"
        possible_files = glob.glob(os.path.join(mne_dir, "*")) + glob.glob(
            os.path.join(mne_dir, "*", "*")
        )
        try:
            input_file = [
                pos_file for pos_file in possible_files if pos_file.endswith(wanted_file)
            ][0]
        except IndexError:
            pass

        if not os.path.isfile(input_file):
            logger.info(
                f"No MolNetEnhancer result present at {mne_dir}. "
                f"(run it on GNPS and download it here if you want to "
                f"use it)"
            )
            return columns, mne_component_dict, mne_cluster2component

        with open(input_file) as inf:
            logger.info(f"reading molnetenhancer results from {mne_dir}")
            header = inf.readline().strip().split("\t")
            # get the columns that are of interest to us
            columns = [
                "cf_direct_parent" if col == "CF_Dparent" else col.lower()
                for i, col in enumerate(header[3:])
                if i % 2 == 0
            ]
            for line in inf:
                line = line.strip("\n").split("\t")
                cluster = line.pop(0)
                component = line.pop(0)
                nr_nodes = line.pop(0)
                # todo: make it easier to query classes of singleton families
                # if singleton family, format like '-1_spectrum-id' like canopus results
                # Note that the singleton families id is "singleton-" + spectrum-id.
                if nr_nodes == "1":
                    component = f"-1_{cluster}"
                class_info = []
                # get class names and scores in (class, score)
                for i in range(0, len(line), 2):
                    cur_class = line[i]
                    if not cur_class or cur_class == "no matches":
                        # add None instead of tuple when no match at this lvl
                        class_tup = None
                        # catch when there is 'no matches' at superclass lvl ->
                        # no prediction and dont add to dict (same as canopus)
                        if columns[int(i / 2)] == "cf_kingdom":
                            break
                    else:
                        class_tup = (line[i], float(line[i + 1]))

                    class_info.append(class_tup)
                if component not in mne_component_dict and class_info:
                    mne_component_dict[component] = class_info
                mne_cluster2component[cluster] = component

        return columns, mne_component_dict, mne_cluster2component

    def spectra_classes(self, spectrum_id):
        """Return classes by relating spectrum_id in the mf_classes.

        Args:
            spectrum_id: int/str, spectrum_id - ints will be converted to str
        """
        classes = []
        if isinstance(spectrum_id, int):
            spectrum_id = str(spectrum_id)
        mf_id = self.spectra2mf.get(spectrum_id)
        if mf_id:
            classes = self.mf_classes.get(mf_id)
        return classes

    @property
    def spectra2mf(self):
        return self._spectra2mf

    @property
    def spectra_classes_names(self):
        return self._spectra_classes_names

    @property
    def spectra_classes_names_inds(self):
        return self._spectra_classes_names_inds

    @property
    def mf_classes(self):
        return self._mf_classes
