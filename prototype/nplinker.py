import sys
import logging
import pickle

import numpy as np

from metabolomics import Spectrum
from metabolomics import MolecularFamily

from genomics import GCF

from data_linking import DataLinks
from data_linking import RandomisedDataLinks
from data_linking import LinkLikelihood
from data_linking import LinkFinder
from data_linking import SCORING_METHODS

from config import Config, Args
from loader import DatasetLoader

from logconfig import LogConfig
logger = LogConfig.getLogger(__file__)

class NPLinker(object):

    OBJ_CLASSES = [Spectrum, MolecularFamily, GCF]
    # enumeration for accessing results of LinkFinder.get_links, which are (3, num_links) arrays:
    # - R_SRC_ID: the ID of an object that was supplied as input to get_links
    # - R_DST_ID: the ID of an object that was discovered to have a link to an input object
    # - R_SCORE: the score for the link between a pair of objects
    R_SRC_ID, R_DST_ID, R_SCORE = range(3)

    def __init__(self, userconfig=None):
        """Initialise an NPLinker instance, automatically loading a dataset and generating scores.

        NPLinker instances can be configured in multiple ways, in ascending order of priority:
            1. A global user-level default configuration file in TOML format, found in the directory:
                    $XDG_CONFIG_HOME/nplinker/nplinker.toml
            2. A local TOML configuration file
            3. Command-line arguments / supplying a manually constructed dict instance

        The user-level configuration file will be created when you first create an NPLinker instance. 
        It contains sensible default values for each setting and is intended to be copied and edited to
        produce dataset-specific configuration files, which will then override any parameters shared
        with the user-level file. To load such a file, simply set the "userconfig" parameter to a 
        string containing the filename. 

        It's also possible to selectively override configuration file parameters by 
        supplying command-line arguments (if running nplinker.py as a script), or by passing
        a dict with a structure corresponding to the configuration file format to this method.

        Some examples may make the various possible combinations a bit clearer:
            # load the default/user-level configuration file and nothing else
            npl = NPLinker()

            # override the default file with a different one
            npl = NPLinker('myconfig.toml')

            # the same thing but running NPLinker as a script
            > python nplinker.py --config "myconfig.toml"

            # use the defaults from the user-level config while modifying the root path
            # to load the dataset from (this is the minimum you would need to change in the
            # default config file)
            npl = NPLinker({'dataset': {'root': '/path/to/dataset'}})

            # the same thing running NPLinker as a script
            > python nplinker.py --dataset.root /path/to/dataset
    
        Args:
            userconfig: supplies user-defined configuration data. May take one of 3 types:
                - None: just load the user-level default configuration file
                - str: treat as filename of a local configuration file to load 
                        (overriding the defaults)
                - dict: contents will be used to override values in the dict generated 
                        from loading the configuration file(s)
        """

        # if userconfig is None => create Config() from empty dict
        # if userconfig is a string => create a dict with 'config' key and string as filename
        # if userconfig is a dict => pass it to Config() directly
        if userconfig is None:
            userconfig = {}
        elif isinstance(userconfig, str):
            userconfig = {'config': userconfig}
        elif not isinstance(userconfig, dict):
            raise Exception('Invalid type for userconfig (should be None/str/dict, found "{}")'.format(type(userconfig)))

        self._config = Config(userconfig)

        # configure logging based on the supplied config params
        LogConfig.setLogLevelStr(self._config.config['loglevel'])
        logfile = self._config.config['logfile']
        if len(logfile) > 0:
            logfile_dest = logging.FileHandler(logfile)
            LogConfig.setLogDestination(logfile_dest)

        # this object takes care of figuring out the locations of all the relevant files/folders
        # and will show error/warning if any missing (depends if optional or not)
        self._loader = DatasetLoader(self._config.config)
    
        self._spectra = []
        self._bgcs = []
        self._gcfs = []
        self._strains = []
        self._metadata = {}
        self._families = []
        self._mibig_bgc_dict = {}

        self._bgc_lookup = {}
        self._spec_lookup = {}

        self._datalinks = None
        self._linkfinder = None

        self._random_scores = {}
        self._links = {}

        self._repro_data = {}
        repro_file = self._config.config['repro_file']
        if len(repro_file) > 0:
            self._save_repro_data(repro_file)

    def _collect_repro_data(self):
        """Creates a dict containing data to aid reproducible runs of nplinker.

        This method creates a dict containing various bits of information about
        the current execution of nplinker. This data will typically be saved to
        a file in order to aid reproducibility using :func:`_save_repro_data`. 

        TODO describe contents

        Returns:
            A dict containing the information described above
        """

        self._repro_data = {}
        # TODO best way to embed latest git commit hash? probably with a packaging script...
        # TODO versions of all Python dependencies used (just rely on Pipfile.lock here?)
        
        # insert command line arguments
        self._repro_data['args'] = {}
        for i, arg in enumerate(sys.argv):
            self._repro_data['args'][i] = arg

        # TODO anything else to include here?

        return self._repro_data

    def _save_repro_data(self, filename):
        self._collect_repro_data()
        with open(filename, 'wb') as repro_file:
            # TODO is pickle the best format to use?
            pickle.dump(self._repro_data, repro_file)
            logger.info('Saving reproducibility data to {}'.format(filename))

    def load_data(self, new_bigscape_cutoff=None):
        """Loads the basic components of a dataset.

        This method is responsible for loading the various pieces of the supplied dataset into
        memory and doing any initial parsing/object manipulation required. After it completes,
        applications can access the lists of GCFs, Spectra, MolecularFamilies and strains 
        using the corresponding properties of the NPLinker class.

        Returns:
            True if successful, False otherwise
        """

        # clear any existing link/score data
        self.clear_links()

        # typical case where load_data is being called with no params
        if new_bigscape_cutoff is None:
            logger.debug('load_data (normal case, full load)')
            if not self._loader.load():
                return False
        else:
            logger.debug('load_data with new cutoff = {}'.format(new_bigscape_cutoff))
            self._loader._bigscape_cutoff = new_bigscape_cutoff
            self._loader._load_genomics()

        self._spectra = self._loader.spectra
        self._molfams = self._loader.molfams
        self._bgcs = self._loader.bgcs
        self._gcfs = self._loader.gcfs
        self._mibig_bgc_dict = self._loader.mibig_bgc_dict
        self._strains = self._loader.strains

        logger.debug('Generating lookup tables: genomics')
        self._bgc_lookup = {}
        for i, bgc in enumerate(self._bgcs):
            self._bgc_lookup[bgc.name] = i
        
        self._gcf_lookup = {}
        for i, gcf in enumerate(self._gcfs):
            self._gcf_lookup[gcf.id] = i

        # don't need to do these two if cutoff changed (indicating genomics data
        # was reloaded but not metabolomics)
        if new_bigscape_cutoff is None:
            logger.debug('Generating lookup tables: metabolomics')
            self._spec_lookup = {}
            for i, spec in enumerate(self._spectra):
                self._spec_lookup[spec.spectrum_id] = i

            self._molfam_lookup = {}
            for i, molfam in enumerate(self._molfams):
                self._molfam_lookup[molfam.id] = i

        logger.debug('load_data: completed')
        return True


    def _generate_scores(self, dl, lf):
        # if metcalf scoring enabled
        if self.scoring.metcalf.enabled:
            lf.metcalf_scoring(dl, type='spec-gcf')
            lf.metcalf_scoring(dl, type='fam-gcf')

        # if hg scoring enabled
        if self.scoring.hg.enabled:
            lf.hg_scoring(dl, type='spec-gcf')
            lf.hg_scoring(dl, type='fam-gcf')

        # if likescore scoring enabled
        if self.scoring.likescore.enabled:
            ll = LinkLikelihood()
            ll.calculate_likelihoods(dl, type='spec-gcf')
            ll.calculate_likelihoods(dl, type='fam-gcf')
            lf.likelihood_scoring(dl, ll, type='spec-gcf')
            lf.likelihood_scoring(dl, ll, type='fam-gcf')

    def process_dataset(self, random_count=None):
        """Construct the DataLinks and LinkFinder objects from the loaded dataset.

        Deals with initialising all the objects used for scoring/linking, and also
        currently manages the creation of randomised scoring matrices for the different
        object types. 

        Args:
            random_count: number of randomised instances to create. The default "None"
                effectively means "use the previously configured value" (e.g. from a
                configuration file). Otherwise the value must be >0, and will override
                and previously configured value.
        """

        if len(self._spectra) == 0 or len(self._gcfs) == 0 or len(self._strains) == 0:
            logger.info('process_dataset: calling load_data')
            self.load_data()

        self._datalinks = DataLinks()
        self._datalinks.load_data(self._spectra, self._gcfs, self._strains)
        self._datalinks.find_correlations()
        self._linkfinder = LinkFinder()

        # generate the actual scores from the standard DataLinks object
        logger.debug('Generating scores, enabled methods={}'.format(self.scoring.enabled()))
        self._generate_scores(self._datalinks, self._linkfinder)

        # create the randomised scoring matrices, and stack them together for later use
        if random_count is not None and random_count < 1:
            raise Exception('random_count must be None or >0 (value={})'.format(random_count))

        if random_count is None:
            random_count = self.scoring.random_count

        logger.debug('Generating {} randomised instances for scoring'.format(random_count))
        rdatalinks = [RandomisedDataLinks.from_datalinks(self._datalinks, True) for x in range(random_count)]
        rlinkfinders = [LinkFinder() for x in range(random_count)]
        logger.debug('Generating randomised scores, enabled methods={}'.format(self.scoring.enabled()))
        
        for i in range(random_count):
            self._generate_scores(rdatalinks[i], rlinkfinders[i])
        
        logger.debug('Finished generating randomised scores')
        # TODO this can be a bit slow?
        # create a set of dicts indexed by class to store the matrices
        self._random_scores = {}
        for m in self.scoring.enabled():
            rand_scores = \
            {
                Spectrum:           np.zeros((len(self._spectra), 0)),
                MolecularFamily:    np.zeros((len(self._families), 0)),
                GCF:                {
                                        Spectrum:           np.zeros((0, len(self._gcfs))),
                                        MolecularFamily:    np.zeros((0, len(self._gcfs))),
                                    }
            }

            for i in range(random_count):
                rand_scores[GCF][Spectrum] = np.r_[rand_scores[GCF][Spectrum], rlinkfinders[i].get_scores(m.name, 'spec-gcf')]
                rand_scores[GCF][MolecularFamily] = np.r_[rand_scores[GCF][MolecularFamily], rlinkfinders[i].get_scores(m.name, 'fam-gcf')]
                rand_scores[Spectrum] = np.c_[rand_scores[Spectrum], rlinkfinders[i].get_scores(m.name, 'spec-gcf')]
                if len(self._families) > 0:
                    rand_scores[MolecularFamily] = np.c_[rand_scores[MolecularFamily], rlinkfinders[i].get_scores(m.name, 'fam-gcf')]

            self._random_scores[m.name] = rand_scores
        logger.debug('Scoring matrices created')

        self.clear_links()
        return True

    def clear_links(self):
        """Clear any previously retrieved links"""
        self._links = {}

    def get_links(self, objects, scoring_method, datalinks=None):
        """Collects sets of links between a given set of objects and their counterparts.

        This method allows an application to pass in one or more objects of type GCF, Spectrum,
        or MolecularFamily, and retrieve the set of links to their counterpart objects (e.g.
        GCF for Spectrum input) that meet the scoring criteria for the selected scoring
        method (see ScoringConfig class / self.scoring)
        
        Returned scores can then be accessed using NPLinker.all_links or NPLinker.links_for_obj.

        Args:
            objects: An instance of one of the above types, or a list containing any mix of them.
            scoring_method: selects the scoring method to be used (e.g. metcalf, likescore). The value
                of this parameter should be one of the ScoringConfig.<method> members, which also contain the
                configurable parameters for each method.
            datalinks: The DataLinks instance used for scoring. Defaults to None, in which case the
                standard internal DataLinks instance is used. 

        Returns:
            If either a single instance or a list with uniform type is used as input, the return value
            is a list containing the objects for which links were found to satisfy the scoring criteria.
            In the case where a mixed list of objects is used, the return value will instead be a 3 element
            list, containing the individual results for the classes GCF, Spectrum, and MolecularFamily 
            (in that order). 
        """
        # single object, pass it straight through
        if not isinstance(objects, list):
            logger.debug('get_links: single object')
            return self._get_links([objects], datalinks, scoring_method)

        # check if the list contains one or multiple types of object
        type_counts = {NPLinker.OBJ_CLASSES[t]: 0 for t in range(len(NPLinker.OBJ_CLASSES))}
        for o in objects:
            try:
                type_counts[type(o)] += 1
            except KeyError:
                raise Exception('Unsupported type passed to get_links ({})'.format(type(o)))

        if max(type_counts.values()) == len(objects):
            # if the list seems to have a uniform type, again pass straight through
            logger.debug('get_links: uniformly typed list ({})'.format(type(objects[0])))
            return self._get_links(objects, datalinks, scoring_method)

        # finally, if it's a mix of objects, extract them into separate lists and call 
        # the internal method multiple times instead
        obj_lists = {t: [] for t in NPLinker.OBJ_CLASSES}
        for o in objects:
            obj_lists[type(o)].append(o)

        # keep consistent ordering here: GCF, Spectrum, MolecularFamily
        obj_lists = [obj_lists[GCF], obj_lists[Spectrum], obj_lists[MolecularFamily]]

        results = [[] for t in NPLinker.OBJ_CLASSES]
        logger.debug('get_links: multiple lists')
        for i in range(len(NPLinker.OBJ_CLASSES)):
            if len(obj_lists[i]) > 0:
                results[i].append(self._get_links(obj_lists[i], datalinks, scoring_method))

        return results

    def _store_object_links(self, src, dst, score, scoring_method):
        # structure of self._links is:
        # src objects as keys, values are scoring methods
        #   scoring methods as keys, values are dst objects
        #       dst objects as keys, values are link scores
        if src not in self._links:
            self._links[src] = {m: {} for m in SCORING_METHODS}

        if scoring_method.name not in self._links[src]: 
            raise Exception('Should never happen!')

        self._links[src][scoring_method.name][dst] = score

    def _get_links(self, objects, datalinks, scoring_method):
        if len(objects) == 0:
            return []

        if self._linkfinder is None:
            raise Exception('Need to call process_dataset first')

        if scoring_method.name not in SCORING_METHODS:
            raise Exception('unknown scoring method "{}"'.format(scoring_method.name))

        datalinks = datalinks or self._datalinks
        input_type = type(objects[0])

        # if using the randomised scores to generate a significance threshold for the links, 
        # then just call get_links without a cutoff value so it will return all scores
        # for the given objects, which we can then apply the perecentile threshold to (this
        # seemed easier than modifying the LinkFinder implementation)

        # TODO this block can probably be simplified 
        if scoring_method.name == 'metcalf':
            if scoring_method.sig_percentile < 0 or scoring_method.sig_percentile > 100:
                raise Exception('sig_percentile invalid! Expected 0-100, got {}'.format(scoring_method.sig_percentile))
            logger.debug('_get_links for metcalf scoring')
            results = self._linkfinder.get_links(datalinks, objects, scoring_method.name, None)
            # if using percentile option, need to filter the links based on that before continuing...
            results = self._get_links_percentile(input_type, objects, results, scoring_method.name, scoring_method.sig_percentile)
        elif scoring_method.name == 'likescore':
            results = self._linkfinder.get_links(datalinks, objects, scoring_method.name, scoring_method.cutoff)
        elif scoring_method.name == 'hg':
            results = self._linkfinder.get_links(datalinks, objects, scoring_method.name, scoring_method.prob)
        else:
            raise Exception('Not handled yet')

        scores_found = set()

        if input_type == GCF:
            logger.debug('get_links: result type is GCF, inputs={}, results={}'.format(len(objects), len(results)))
            # for GCF input, results contains two arrays of shape (3, x), 
            # which contain spec-gcf and fam-gcf links respectively 
            result_gcf_spec, result_gcf_fam = results[0], results[1]

            for res, type_ in [(result_gcf_spec, Spectrum), (result_gcf_fam, MolecularFamily)]:
                if res.shape[1] == 0:
                    logger.debug('Found no links for {}, [{}]'.format(objects, type_))
                    continue # no results


                # TODO: self._families is never populated atm, so this breaks with MolecularFamily results!
                if type_ == MolecularFamily:
                    logger.warning('FIXME: MolecularFamily support')
                    continue
                logger.debug('Found {} links for {}, [{}]'.format(res.shape[1], objects, type_))

                lookup = {gcf.id: objects.index(gcf) for gcf in objects}

                # for each entry in the results (each Spectrum or MolecularFamily)
                for j in range(res.shape[1]):
                    # extract the ID of the object and get the object itself
                    obj_id = int(res[NPLinker.R_DST_ID, j])
                    obj = self._spectra[obj_id] if type_ == Spectrum else self._families[obj_id]

                    # retrieve the GCF object too
                    index = lookup[res[NPLinker.R_SRC_ID, j]]
                    gcf = objects[index]

                    # finally, create the entry in self._links and record that this GCF
                    # has at least one link associated with it
                    self._store_object_links(gcf, obj, res[NPLinker.R_SCORE, j], scoring_method)
                    scores_found.add(gcf)

        else:
            logger.debug('get_links: result type is Spectrum/MolFam, inputs={}, results={}'.format(len(objects), len(results)))
            # for non-GCF input, result is a list containing a single array, shape (3, x)
            # where x is the total number of links found
            results = results[0]
            if results.shape[1] == 0:
                logger.debug('Found no links for {}'.format(objects))
                return []

            logger.debug('Found {} links for {}'.format(results.shape[1], objects))

            lookup = {obj.id: objects.index(obj) for obj in objects}

            # for each entry in the results (each GCF)
            for j in range(results.shape[1]):
                # extract the ID of the GCF and use that to get the object
                gcf_id = int(results[NPLinker.R_DST_ID, j])
                gcf = self._gcfs[gcf_id]

                # retrieve the source object similarly
                index = lookup[results[NPLinker.R_SRC_ID, j]]
                obj = objects[index]

                # finally, create the entry in self._links and record that this Spectrum or
                # MolecularFamily has at least one link associated with it
                self._store_object_links(obj, gcf, results[NPLinker.R_SCORE, j], scoring_method)
                scores_found.add(obj)

        return list(scores_found)

    def _get_links_percentile(self, input_type, objects, results, scoring_method_name, sig_percentile):
        # self._random_scores contains the randomised and stacked scores for each given type of 
        # object, so pull those out and apply np.percentile with the given threshold to the correct 
        # matrices here. the results object from get_links will contain *ALL* potential candiate
        # links at this stage - e.g. for a Spectrum input, it will return all GCFs

        # retrieve the IDS for all the supplied objects
        obj_ids = [obj.id for obj in objects]

        if input_type == GCF:
            # 2-element list returned with spec-gcf, fam-gcf results
            results_spec, results_fam = results

            # lookup the appropriate random scores
            rand_spec_scores = self._random_scores[scoring_method_name][input_type][Spectrum]
            rand_fam_scores = self._random_scores[scoring_method_name][input_type][MolecularFamily]

            # calculate the threshold values for each input object
            spec_perc_thresholds = [np.percentile(rand_spec_scores[:, id], sig_percentile) for id in obj_ids]
            fam_perc_thresholds = [np.percentile(rand_fam_scores[:, id], sig_percentile) for id in obj_ids]

            perc_results = [np.zeros((3, 0)), np.zeros((3, 0))]

            data = [(results_spec, spec_perc_thresholds, len(self._spectra))]
            data.append((results_fam, fam_perc_thresholds, len(self._families)))

            for d, (full_results, thresholds, num_objs) in enumerate(data):
                for i in range(len(objects)):
                    # want to get the indices of full_results where a) the source ID matches the
                    # current object and b) the score exceeds the threshold for that object
                    obj_indices = np.intersect1d(np.where(full_results[NPLinker.R_SCORE, :] >= thresholds[i]),
                                                 np.where(full_results[NPLinker.R_SRC_ID, :] == objects[i].id))

                    # concat the selected columns to the final result
                    perc_results[d] = np.c_[perc_results[d], full_results[:, obj_indices]]

            return [perc_results[0], perc_results[1]]

        elif input_type == Spectrum or input_type == MolecularFamily:
            # here results will be a list which in turn contains a single array, so just extract that
            results = results[0]

            # lookup the appropriate random scores
            rand_scores = self._random_scores[scoring_method_name][input_type]

            # calculate the threshold values for each input object
            perc_thresholds = [np.percentile(rand_scores[id, :], sig_percentile) for id in obj_ids]

            perc_results = np.zeros((3, 0))

            for i in range(len(objects)):
                # want to get the indices of full_results where a) the source ID matches the
                # current object and b) the score exceeds the threshold for that object
                obj_indices = np.intersect1d(np.where(results[NPLinker.R_SCORE, :] >= perc_thresholds[i]),
                                             np.where(results[NPLinker.R_SRC_ID, :] == obj_ids[i]))

                # build up a new matrix containing only the sufficiently high scoring columns
                perc_results = np.c_[perc_results, results[:, obj_indices]]

            return [perc_results]

        logger.error('Bad input type in _get_links_percentile??')
        return None

    def get_common_strains(self, objects_a, objects_b, filter_no_shared=True):
        # this is a dict with structure:
        #   (Spectrum/MolecularFamily, GCF) => list of strain indices
        common_strains = self._datalinks.common_strains(objects_a, objects_b, filter_no_shared)

        # replace the lists of strain indices with actual strain objects
        for objpair in common_strains.keys():
            common_strains[objpair] = [self._strains[x] for x in common_strains[objpair]]
    
        return common_strains

    def has_bgc(self, name):
        return name in self._bgc_lookup

    def lookup_bgc(self, name):
        if name not in self._bgc_lookup:
            return None

        return self._bgcs[self._bgc_lookup[name]]

    def lookup_spectrum(self, name):
        if name not in self._spec_lookup:
            return None

        return self._spectra[self._spec_lookup[name]]

    @property
    def strains(self):
        return self._strains

    @property
    def bgcs(self):
        return self._bgcs

    @property
    def gcfs(self):
        return self._gcfs

    @property
    def spectra(self):
        return self._spectra

    @property
    def molfams(self):
        return self._molfams

    @property
    def metadata(self):
        return self._metadata

    @property
    def mibig_bgc_dict(self):
        return self._mibig_bgc_dict

    @property
    def datalinks(self):
        """
        Returns the main DataLinks object used for scoring etc
        """
        return self._datalinks

    @property
    def rdatalinks(self):
        """
        Returns the list of RandomisedDataLinks objects
        """
        return self._rdatalinks

    @property
    def scoring(self):
        """
        Returns the ScoringConfig object used to manipulate the scoring process
        """
        return self._config.scoring

    @property
    def repro_data(self):
        """
        Returns the dict of reproducibility data
        """
        return self._repro_data

    @property
    def all_links(self):
        """
        Returns the links for all objects and scoring methods
        """
        return self._links

    def links_for_obj(self, obj, scoring_method, sort=True, type_=None):
        # the type_ part is only really useful for GCFs, to pick between MolFam/Spectrum results
        if obj not in self._links:
            logger.warning('No links found for supplied object. Either get_links was not called first or no results were found')
            return None
        links = [(obj, score) for obj, score in self._links[obj][scoring_method.name].items() if (type_ is None or isinstance(obj, type_))]
        if sort:
            links.sort(key=lambda x: x[1], reverse=True)
        return links

if __name__ == "__main__":
    # can set default logging configuration this way...
    LogConfig.setLogLevel(logging.DEBUG)

    # initialise NPLinker from the command-line arguments
    npl = NPLinker(Args().get_args())

    # load the dataset
    if not npl.load_data():
        print('Failed to load the dataset!')
        sys.exit(-1)

    # make sure metcalf scoring is enabled 
    npl.scoring.metcalf.enabled = True

    # generate all datalinks objects etc
    if not npl.process_dataset(random_count=1):
        print('Failed to process dataset')
        sys.exit(-1)

    # pick a GCF to get links for
    test_gcf = npl.gcfs[8]

    # set metcalf scoring significance percentile to 99%
    npl.scoring.metcalf.sig_percentile = 99

    # call get_links and set the scoring method to metcalf by referencing the 
    # scoring.metcalf object
    result = npl.get_links(test_gcf, scoring_method=npl.scoring.metcalf)

    # check if any links were found
    if test_gcf not in result:
        print('No links found!')
        sys.exit(0)

    # retrieve the specific links for this object (get_links() supports passing
    # multiple objects as input and so doesn't return them directly). since multiple
    # scoring methods may be in use, need to supply the desired scoring object as a
    # parameter. By default this will return objects of any available type (so Spectrum
    # or MolecularFamily for GCF input), but the type_ parameter can be used to select
    # a single class if needed.
    test_gcf_links = npl.links_for_obj(test_gcf, npl.scoring.metcalf, type_=Spectrum)

    # print the objects (spectra in this example) and their scores
    for obj, score in test_gcf_links:
        print('{} : {}'.format(obj, score))
        # also retrieve common strains. the return value here is a dict indexed
        # by object tuples (either (Spectrum, GCF) or (MolecularFamily, GCF) depending
        # on the inputs), with the values being lists of shared strains (strings)
        common_strains = npl.get_common_strains(test_gcf, obj)
        strain_names = list(common_strains.values())[0]
        if len(strain_names) == 0:
            print('   No shared strains!')
        else:
            print('   {} shared strains: {}'.format(len(strain_names), strain_names))

    print('{} total links found'.format(len(test_gcf_links)))
