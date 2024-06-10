import logging
import time
from nplinker.genomics import BGC
from nplinker.genomics import GCF
from nplinker.metabolomics import Spectrum
from nplinker.strain import StrainCollection
from .abc import ScoringBase
from .link_graph import LinkGraph
from .score import Score
from .scoring_method import ScoringMethod


logger = logging.getLogger(__name__)


class NPClassScoring(ScoringBase):
    name = ScoringMethod.NPLCLASS.value

    def __init__(self, npl):
        super().__init__(npl)
        self.cutoff = 0.25
        self.method_options = npl.chem_classes.class_predict_options
        self.method = self.method_options[0]
        self.num_results = 1  # how many scores do you want for each link
        # take care what targets are as this method can link both bgc/gcf to
        # both spectrum/gcf
        self.equal_targets = False  # if obj GCF, target is spec not MF
        self.both_targets = False  # if obj GCF, target both spec and MF
        # filter out spectra without a score due to missing spectrum classes
        self.filter_missing_scores = False
        self._target_no_scores = set()

    def _npclass_score(self, obj, target, method="mix", obj_classes=None, target_classes=None):
        """Return sorted class link scores for scoring obj and target.

        The input objects can be any mix of the following NPLinker types:
            - BGC
            - GCF
            - Spectrum
            - MolecularFamily

        Classes will be retrieved from NPLinker object with _get_gen_classes
        and _get_met_classes or the classes can be predetermined and given as
        input with the function, this should correspond to obj and
        target (if obj is BGC then  obj_classes should be classes for the BGC).
        Can be:
            - BGC/GCF classes: {'as_classes': list(str), 'bigscape_class': str}
            - Spectrum/mf classes: tuple of (classes, class_names_indices),
            where classes is a list of list of tuples/None, where each
            tuple is a class and a score (str, float), and class_names_indices
            is a list of ints that relate to the name of a class ontology lvl

        Args:
            obj: one of the possible input objects
            target: one of the possible input objects
            method: str, which classification method to use for spectra. default is 'main'
                options:
                    -'main': use the main method - currently canopus
                    -'mix': use main method first, when no main classes present, use the others
                    if present
                    -'canopus': use only canopus class predictions
                    -'molnetenhancer': use only molnetenhancer class predictions
            obj_classes: default - None, or classes for one of the input types
            target_classes: default - None, or classes for one of the input types

        Returns:
            List of tuple of
            (score, obj_class_lvl, target_class_lvl, obj_class, target_class)
            (float, str, str, str, str)
            List will be empty if either BGC or spectrum classes are missing
        """
        # todo: make subroutines
        # assess what is obj and target
        spec_like = obj
        bgc_like = target
        spec_like_classes_tup = obj_classes
        bgc_like_classes_dict = target_classes
        bgc_to_spec = False
        if isinstance(obj, BGC) or isinstance(obj, GCF):
            bgc_like = obj
            spec_like = target
            bgc_like_classes_dict = obj_classes
            spec_like_classes_tup = target_classes
            bgc_to_spec = True

        # assess method - move to get_links?
        assert method in self.method_options, (
            f"NPClass method should be one of method options: {self.method_options}, if your method is not "
            + "in the options check if the class predictions (canopus, etc.) are loaded correctly"
        )

        # gather correct classes if not provided, dict for bgcs, tup for spec
        if not bgc_like_classes_dict:
            bgc_like_classes_dict = self._get_gen_classes(bgc_like)
        if not spec_like_classes_tup:
            spec_like_classes_tup = self._get_met_classes(spec_like, method)
        # unpack spec_like classes - both are lists
        spec_like_classes, spec_like_classes_names_inds = spec_like_classes_tup

        scores = []  # this will be returned if one of the class sides is absent
        std_score = 0  # if link not recorded in scores (mibig) return this score
        # loop through classes that are possible to link (names in class_match object)
        for bgc_class_name in self.npl.class_matches.bgc_class_names:
            if bgc_class_name == "mibig_classes":
                # treat specially as bigscape class needs to be translated to mibig class
                bigscape_class = bgc_like_classes_dict["bigscape_class"]
                # convert bigscape class to mibig class
                bgc_like_classes = [
                    self.npl.class_matches.bigscape_mibig_conversion.get(bigscape_class)
                ]
            else:
                bgc_like_classes = bgc_like_classes_dict.get(bgc_class_name)
            if bgc_like_classes and spec_like_classes:  # check for classes from both sides
                for bgc_class in bgc_like_classes:
                    for chem_class_name in self.npl.class_matches.chem_class_names:
                        # does info exist for this spectrum class level, return index for class level
                        spec_class_level_i = spec_like_classes_names_inds.get(chem_class_name)
                        if spec_class_level_i:
                            spec_class_tup = spec_like_classes[spec_class_level_i]
                            if spec_class_tup:  # if there is a class at this lvl
                                # is a tuple of (name, score) so take [0]
                                spec_class = spec_class_tup[0]
                                # determine direction of scoring: BGC -> spectrum
                                if bgc_to_spec:
                                    score = (
                                        self.npl.class_matches.class_matches[bgc_class_name][
                                            chem_class_name
                                        ]
                                        .get(bgc_class, {})
                                        .get(spec_class, std_score)
                                    )
                                    result_tuple = (
                                        score,
                                        bgc_class_name,
                                        chem_class_name,
                                        bgc_class,
                                        spec_class,
                                    )
                                else:  # spectrum -> BGC
                                    score = (
                                        self.npl.class_matches.class_matches[chem_class_name][
                                            bgc_class_name
                                        ]
                                        .get(spec_class, {})
                                        .get(bgc_class, std_score)
                                    )
                                    result_tuple = (
                                        score,
                                        chem_class_name,
                                        bgc_class_name,
                                        spec_class,
                                        bgc_class,
                                    )
                                scores.append(result_tuple)
        return sorted(scores, reverse=True)

    def _get_targets(self, test_id):
        """Get the targets based upon instance of test_id, returns list of targets.

        Args:
            test_id: one of the NPLinker objects: BGC, GCF, Spectrum, mf
        Returns:
            List of one or more of one of the NPLinker objects
        """
        if isinstance(test_id, BGC) or isinstance(test_id, GCF):  # genome side
            return self._get_targets_genomics(test_id)
        else:  # metabolome side
            return self._get_targets_metabolomics(test_id)

    def _get_targets_metabolomics(self, test_id):
        if self.both_targets:
            targets = self.npl.bgcs + self.npl.gcfs
        elif isinstance(test_id, Spectrum):
            if self.equal_targets:
                targets = self.npl.bgcs
            else:
                targets = self.npl.gcfs
        else:  # obj are mf
            if self.equal_targets:
                targets = self.npl.gcfs
            else:
                targets = self.npl.bgcs
        return targets

    def _get_targets_genomics(self, test_id):
        if self.both_targets:  # no matter BGC or GCF take both spec and MF
            targets = self.npl.spectra + self.npl.mfs
        elif isinstance(test_id, BGC):  # obj are BGC
            if self.equal_targets:  # take
                targets = self.npl.spectra
            else:
                targets = self.npl.mfs
        else:  # obj are GCF
            if self.equal_targets:
                targets = self.npl.mfs
            else:
                targets = self.npl.spectra
        return targets

    def _get_gen_classes(self, bgc_like, gcf_as_cutoff=0.5):
        """Get classes for genomics objects.

        Args:
            bgc_like: BGC or GCF object from NPLinker input objects
            gcf_as_cutoff: float - if GCF, get antismash classes with
                class_matches.get_gcf_as_classes if the class occurs in more
                than this fraction of the GCF, default = 0.5
        Returns:
            Genome-based classes as a dict with antismash and bigscpae classes
            {'as_classes': list(str), 'bigscape_class': str}
        """
        # assess if bgc or gcf
        is_bgc = isinstance(bgc_like, BGC)
        if is_bgc:
            # get parent gcf for bgc
            bgc_like_gcf = [
                gcf for gcf in self.npl.gcfs if bgc_like.id in [b.id for b in gcf.bgcs]
            ][0]
            # gather AS classes and convert to names in scoring dict
            as_classes = self.npl.class_matches.convert_as_classes(
                bgc_like.product_prediction.split(".")
            )
            bgc_like_classes_dict = {
                "bigscape_class": bgc_like_gcf.bigscape_class,
                # str - always one bigscape class right?
                "as_classes": as_classes,
            }  # list(str)
        else:
            as_classes = self.npl.class_matches.convert_as_classes(
                self.npl.class_matches.get_gcf_as_classes(bgc_like, gcf_as_cutoff)
            )
            bgc_like_classes_dict = {
                "bigscape_class": bgc_like.bigscape_class,
                # str - always one bigscape class right?
                "as_classes": as_classes,
            }  # list(str)
        return bgc_like_classes_dict

    def _get_met_classes(self, spec_like, method="mix"):
        """Get chemical classes for a Spectrum or mf based on method.

        Args:
            spec_like: Spectrum or mf, one of the NPLinker input types
            method: str, one of the appropriate methods for chemical class
                predictions (mix, canopus...), default='mix'
        Returns:
            tuple of (classes, class_names_indices),
            where classes is a list of list of tuples/None, where each
            tuple is a class and a score (str, float), and class_names_indices
            is a list of ints that relate to the name of a class ontology lvl
        """
        # assess if spectrum or mf
        is_spectrum = isinstance(spec_like, Spectrum)

        # gather classes for spectra, using right method
        # choose the main method here by including it as 'main' in the method parameter
        use_canopus = (
            "main" in method or "canopus" in method or "mix" in method
        ) and "canopus" in self.method_options
        use_mne = (
            "molnetenhancer" in method or "mix" in method
        ) and "molnetenhancer" in self.method_options
        spec_like_classes, spec_like_classes_names, spec_like_classes_names_inds = (
            None,
            None,
            None,
        )
        # the order in which the classes are read, determines the priority (now: first canopus, then mne)
        if use_canopus and not spec_like_classes:
            if is_spectrum:
                # list of list of tuples/None - todo: add to spectrum object?
                # take only 'best' (first) classification per ontology level
                all_classes = self.npl.chem_classes.canopus.spectra_classes.get(spec_like.id)
                if all_classes:
                    spec_like_classes = [
                        cls_per_lvl
                        for lvl in all_classes
                        for i, cls_per_lvl in enumerate(lvl)
                        if i == 0
                    ]
                spec_like_classes_names_inds = (
                    self.npl.chem_classes.canopus.spectra_classes_names_inds
                )
            else:  # mf
                fam_id = spec_like.family.id
                if fam_id.startswith("singleton-"):  # account for singleton families
                    fam_id += f"_{spec_like.spectra[0].id}"
                all_classes = self.npl.chem_classes.canopus.mf_classes.get(fam_id)
                if all_classes:
                    spec_like_classes = [
                        cls_per_lvl
                        for lvl in all_classes
                        for i, cls_per_lvl in enumerate(lvl)
                        if i == 0
                    ]
                spec_like_classes_names_inds = self.npl.chem_classes.canopus.mf_classes_names_inds
        if use_mne and not spec_like_classes:
            # if mne or when main/canopus does not get classes
            if is_spectrum:
                spec_like_classes = self.npl.chem_classes.molnetenhancer.spectra_classes(
                    spec_like.id
                )
            else:  # mf
                fam_id = spec_like.family.id
                if fam_id.startswith("singleton"):  # account for singleton families
                    fam_id += f"_{spec_like.spectra[0].id}"
                spec_like_classes = self.npl.chem_classes.molnetenhancer.mf_classes.get(fam_id)
            # classes are same for mf and spectrum so names are irrespective of is_spectrum
            spec_like_classes_names_inds = (
                self.npl.chem_classes.molnetenhancer.spectra_classes_names_inds
            )
        return spec_like_classes, spec_like_classes_names_inds

    @classmethod
    def setup(cls, npl):
        """Perform any one-off initialisation required (will only be called once)."""
        logger.info("Set up NPClassScore scoring")
        met_options = npl.chem_classes.class_predict_options
        logger.info(f"Please choose one of the methods from {met_options}")
        if not met_options:
            logger.warn(
                "There are no methods available! This is probably "
                "because no class predictions (canopus, etc.) were "
                "present or they are not loaded correctly"
            )
        else:
            logger.info(f"Currently the method '{met_options[0]}' is selected")
        # todo: give info about parameters

    def get_links(self, *objects, **parameters):
        # TODO: replace some attributes with parameters
        """Given a set of objects, return link information."""
        # todo: pickle results
        logger.info("Running NPClassScore...")
        begin = time.time()
        first_obj = objects[0]
        targets = self._get_targets(first_obj)
        obj_is_gen = isinstance(first_obj, BGC) or isinstance(first_obj, GCF)

        # only get target classes once for each target here
        if obj_is_gen:  # obj is genome so get metabolome classes for target
            targets_classes = [self._get_met_classes(target, self.method) for target in targets]
        else:
            targets_classes = [self._get_gen_classes(target) for target in targets]

        # TODO: implement the computation of common strains between objects and targets
        common_strains = StrainCollection()

        logger.info(
            f"Calculating NPClassScore for {len(objects)} objects to "
            f"{len(targets)} targets ({len(common_strains)} pairwise "
            f"interactions that share at least 1 strain). This might "
            f"take a while."
        )

        lg = LinkGraph()
        results = {}
        for obj in objects:
            results[obj] = {}
            # get obj class
            if obj_is_gen:
                obj_classes = self._get_gen_classes(obj)
            else:
                obj_classes = self._get_met_classes(obj, self.method)

            for target, target_classes in zip(targets, targets_classes):
                self._create_object_link(
                    obj_is_gen,
                    common_strains,
                    lg,
                    obj,
                    obj_classes,
                    target,
                    target_classes,
                    parameters,
                )

        # info about spectra/MFs with missing scoring
        len_missing = len(self._target_no_scores)
        if len_missing > 0:
            filter_msg = "kept"
            if self.filter_missing_scores:
                filter_msg = "filtered out"
            logger.warning(
                f"{len_missing} targets have no NPClassScore "
                f"prediction due to missing class predictions and are "
                f"{filter_msg} by default. Adjust .filter_missing_scores "
                f"to change."
            )

        logger.info(f"NPClassScore completed in {time.time() - begin:.1f}s")
        return lg

    def _create_object_link(
        self, obj_is_gen, common_strains, lg, obj, obj_classes, target, target_classes, parameters
    ):
        # only consider targets that have shared strains
        common_tup = (obj, target)
        if obj_is_gen:
            common_tup = (target, obj)
        shared_strains = common_strains.get(common_tup)
        if shared_strains is not None:  # todo: fix for bgcs
            full_score = self._npclass_score(obj, target, self.method, obj_classes, target_classes)[
                : self.num_results
            ]
            try:
                npclassscore = full_score[0][0]
            except IndexError:
                # no score is found due to missing classes for spectra
                self._target_no_scores.add(target)  # keep track
                if not self.filter_missing_scores:
                    lg.add_link(obj, target, Score(self.name, full_score, parameters))
            else:
                if npclassscore > self.cutoff:
                    lg.add_link(obj, target, Score(self.name, full_score, parameters))

    def format_data(self, data):
        """Given whatever output data the method produces, return a readable string version."""
        # data or full_score is a list of tuples, here return just NPClassScore
        formatted_data = None  # default when there is no score (missing class)
        if data:
            # there is a score
            formatted_data = f"{data[0][0]:.3f}"
        return formatted_data

    def sort(self, objects, reverse=True):
        """Given a list of objects, return them sorted by link score."""
        return sorted(objects, key=lambda objlink: objlink[self], reverse=reverse)
