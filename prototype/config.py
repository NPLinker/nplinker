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
        def bool_checker(x):
            return str(x).lower() == 'true'

        self.parser = argparse.ArgumentParser(description='nplinker arguments', epilog='Note: command-line arguments will override '
                                              'arguments from configuration files')
        self.parser.add_argument('-c', '--config', help='Path to a .toml configuration file', metavar='path')
        self.parser.add_argument('-d', '--dataset.root', help='Root path for the dataset to be loaded', metavar='path')
        self.parser.add_argument('-l', '--loglevel', help='Logging verbosity level: DEBUG, INFO, WARNING, ERROR', metavar='loglevel')
        self.parser.add_argument('-f', '--logfile', help='Redirect logging from stdout to this file', metavar='logfile')

        self.parser.add_argument('--antismash-format', help='Antismash file format (currently must be "flat"!)', metavar='format')
        self.parser.add_argument('--repro-file', help='Filename to store reproducibility data in', metavar='filename')

        self.parser.add_argument('-r', '--scoring.random', help='Number of randomized instances to create during scoring', metavar='num')
        # TODO better just leaving these in config file?
        self.parser.add_argument('--scoring.metcalf.enabled', type=bool_checker, help='Metcalf scoring enabled/disabled', metavar='true|false')
        self.parser.add_argument('--scoring.metcalf.sig_percentile', type=int, help='Metcalf scoring percentile threshold value (0-100)', metavar='val')
        self.parser.add_argument('--scoring.hg.enabled', type=bool_checker, help='Hypergeometric scoring enabled/disabled', metavar='true|false')

        self.parser.add_argument('--scoring.hg.prob', type=float, help='Hypergeometric scoring threshold (0-1.0)', metavar='val')
        self.parser.add_argument('--scoring.likescore.enabled', type=bool_checker, help='Likescore scoring enabled/disabled', metavar='true|false')
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
            if config_dict['config'] is not None:
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

        if 'dataset' not in config:
            raise Exception('No dataset defined in configuration!')

        if 'dataset.root' in config:
            root = config['dataset.root']
            logger.debug('Dataset root is being set to "{}"'.format(root))
            config['dataset']['root'] = root
            del config['dataset.root']
            if root is None or not os.path.exists(root):
                raise Exception('Dataset path "{}" not found or not accessible'.format(root))

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
        Returns a list of the currently enabled scoring methods
        """
        return [self._methods[x] for x in SCORING_METHODS if self._methods[x].enabled]

    def enabled_names(self):
        """
        Returns a list of the names of all enabled scoring methods
        """
        return [m.name for m in self.enabled()]

    def all(self):
        """
        Returns a list of all supported scoring methods
        """
        return list(self._methods.values())

    def all_names(self):
        """
        Returns a list of the names of all supported scoring methods
        """
        return list(self._methods.keys())

    def __repr__(self):
        r = ''
        for n in self._methods.keys():
            r += '[{}]\n'.format(n)
            vals = self._methods[n]
            for vn, vv in vals.__dict__.items():
                if vn == 'name':
                    continue

                r += '   {} = {}\n'.format(vn, vv)
            r += '------------------\n'

        return r


