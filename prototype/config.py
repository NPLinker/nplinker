import argparse
import os
import glob
from types import SimpleNamespace
from shutil import copyfile
from collections import Mapping

import toml
from xdg import XDG_CONFIG_HOME

from data_linking import SCORING_METHODS

from logconfig import LogConfig
logger = LogConfig.getLogger(__file__)

class Args(object):

    def __init__(self):
        self.parser = argparse.ArgumentParser(description='nplinker arguments', epilog='Note: command-line arguments will override '
                                              'arguments from configuration files')
        self.parser.add_argument('-c', '--config', help='Path to a .toml configuration file', metavar='path')
        self.parser.add_argument('-d', '--dataset', help='Root path for the dataset to be loaded', metavar='path')
        self.parser.add_argument('-l', '--loglevel', help='Logging verbosity level: DEBUG, INFO, WARNING, ERROR', metavar='loglevel')
        self.parser.add_argument('-f', '--logfile', help='Redirect logging from stdout to this file', metavar='logfile')

        self.parser.add_argument('--antismash-format', help='Antismash file format (currently must be "flat"!)', metavar='format')
        self.parser.add_argument('--repro-file', help='Filename to store reproducibility data in', metavar='filename')

        self.parser.add_argument('-r', '--scoring.random', help='Number of randomized instances to create during scoring', metavar='num')
        # TODO better just leaving these in config file?
        self.parser.add_argument('--scoring.metcalf.enabled', type=bool, help='Metcalf scoring enabled/disabled', metavar='true|false')
        self.parser.add_argument('--scoring.metcalf.sig_percentile', type=int, help='Metcalf scoring percentile threshold value (0-100)', metavar='val')
        self.parser.add_argument('--scoring.hg.enabled', type=bool, help='Hypergeometric scoring enabled/disabled', metavar='true|false')

        self.parser.add_argument('--scoring.hg.prob', type=float, help='Hypergeometric scoring threshold (0-1.0)', metavar='val')
        self.parser.add_argument('--scoring.likescore.enabled', type=bool, help='Likescore scoring enabled/disabled', metavar='true|false')
        self.parser.add_argument('--scoring.likescore.cutoff', type=int, help='Likescoring cutoff threshold value', metavar='val')
        self.args = self.parser.parse_args()

    def get_args(self):
        # restructure the arguments into a dict with the same nested structure as Config class expects
        orig = vars(self.args)
        args = {}
        for k, v in orig.items():
            # ignore any params with no value given
            if v is None:
                continue

            # values with non-dotted names can get inserted directly
            if k.find('.') == -1:
                args[k] = v
            else:
                # otherwise add a nested dict for each dotted part, then
                # insert the actual value on the innermost level
                parts = k.split('.')
                root = args
                for p in parts[:-1]:
                    if p not in root:
                        root[p] = {}
                    root = root[p]
                root[parts[-1]] = v
        return args

class DatasetLoader(object):

    def __init__(self, root_dir, overrides):
        """Given a root directory for a dataset, attempts to construct the filenames
            that are necessary to load everything into nplinker.
        """

        self.root = root_dir
        spec_dir = os.path.join(root_dir, 'spectra')

        # just <root_dir>/mibig_json
        if 'mibig_json_dir' in overrides:
            self.mibig_json_dir = overrides['mibig_json_dir']
        else:
            self.mibig_json_dir = os.path.join(root_dir, 'mibig_json')

        if 'mgf_file' in overrides:
            self.mgf_file = overrides['mgf_file']
        else:
            # should have an MGF file named <root_dir>/spectra/<something>.mgf
            self.mgf_file = self._find_via_glob(os.path.join(spec_dir, '*.mgf'), 'MGF file')

        # annotation file(s) may not exist
        if 'annotation_files' in overrides:
            self.annotation_files = overrides['annotation_files']
        else:
            self.annotation_files = glob.glob(os.path.join(spec_dir, '*.annotations.*'))

        # both of these should exist

        if 'edges_file' in overrides:
            self.edges_file = overrides['edges_file']
        else:
            self.edges_file = self._find_via_glob(os.path.join(spec_dir, "*.pairsinfo"), 'edges file')

        if 'nodes_file' in overrides:
            self.nodes_file = overrides['nodes_file']
        else:
            self.nodes_file = self._find_via_glob(os.path.join(spec_dir, '*.out'), 'nodes file')

        if 'bigscape_dir' in overrides:
            self.bigscape_dir = overrides['bigscape_dir']
        else:
            self.bigscape_dir = os.path.join(root_dir, 'bigscape')

        if 'antismash_dir' in overrides:
            self.antismash_dir = overrides['antismash_dir']
        else:
            self.antismash_dir = os.path.join(root_dir, 'antismash')

    def _find_via_glob(self, path, file_type):
        try:
            filename = glob.glob(path)[0]
            return filename
        except (OSError, IndexError) as e:
            # "from None" suppresses the traceback for the original exception, which isn't really needed
            raise Exception('ERROR: unable to find {} in path "{}"'.format(file_type, path)) from None

    def key_paths(self):
        return [self.mgf_file, self.edges_file, self.nodes_file, self.bigscape_dir, self.antismash_dir]

    def __repr__(self):
        return 'Root={}\n   MGF={}\n   EDGES={}\n   NODES={}\n   BIGSCAPE={}\n   MIBIG_JSON={}\n   ANTISMASH={}\n   ANNOTATIONS={}'.format(
                self.root, self.mgf_file, self.edges_file, self.nodes_file, self.bigscape_dir, self.mibig_json_dir, self.antismash_dir, self.annotation_files)


class Config(object):
    """Wrapper for all NPLinker configuration options"""

    DEFAULT_CONFIG = 'nplinker.toml'

    def __init__(self, config_dict):
        self.default_config_path = os.path.join(XDG_CONFIG_HOME, 'nplinker', Config.DEFAULT_CONFIG)
        if not os.path.exists(self.default_config_path):
            logger.debug('Creating default config file')
            os.makedirs(os.path.join(XDG_CONFIG_HOME, 'nplinker'), exist_ok=True)
            copyfile(os.path.join(os.path.dirname(__file__), Config.DEFAULT_CONFIG), self.default_config_path)

        # load the default per-user config file, then check for one provided as an argument
        # and if present use it to override the defaults
        logger.debug('Parsing default config file: {}'.format(self.default_config_path))
        config = toml.load(open(self.default_config_path, 'r'))
        if 'config' in config_dict:
            logger.debug('Loading user config {}'.format(config_dict['config']))
            user_config = toml.load(open(config_dict['config'], 'r'))
            config.update(user_config)
            del config_dict['config']

        # remaining values in the dict should override the existing ones from config files
        # however if running non-interactively, argparse will set values of all non-specified
        # options to None and don't want to wipe out existing settings, so do things this way
        def update(d, u):
            for k, v in u.items():
                if isinstance(v, Mapping):
                    d[k] = update(d.get(k, {}), v)
                elif v is not None:
                    d[k] = v
            return d
        config = update(config, config_dict)

        # store the scoring configuration in a dedicated object so it can be more
        # easily manipulated while the app is running, then remove from the overall config
        self.scoring = ScoringConfig(config['scoring'])
        del config['scoring']
        self.config = config

class ScoringConfig(object):

    def __init__(self, scoringdict):
        self._methods = {}
        for m in SCORING_METHODS:
            scoringdict[m]['name'] = m
            self._methods[m] = SimpleNamespace(**scoringdict[m])
            setattr(self, m, self._methods[m])
        self.random_count = scoringdict['random_count']

    def get(self, m):
        if m in self._methods:
            return self._methods[m]

        raise Exception('Unknown scoring method "{}"'.format(m))

    def enabled(self):
        """
        Returns a list of the enabled scoring methods
        """
        return [self._methods[x] for x in SCORING_METHODS if self._methods[x].enabled]

    def all(self):
        """
        Returns a list of all available scoring methods
        """
        return list(self._methods.keys())

