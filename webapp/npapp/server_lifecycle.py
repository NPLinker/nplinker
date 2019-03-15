import sys
import os
import csv

import bokeh.palettes as bkp

# add path to nplinker/prototype 
# TODO probably will need changed at some point!
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../prototype'))

from nplinker import NPLinker
from logconfig import LogConfig

class NPLinkerHelper(object):
    """
    This is just a simple class to wrap up the various objects involved in 
    loading and plotting a dataset via nplinker/bokeh
    """
    def __init__(self):
        LogConfig.setLogLevelStr('DEBUG')

    def load(self):
        # initialise nplinker and load the dataset, using a config file in the webapp dir
        self.nplinker = NPLinker(os.path.join(os.path.dirname(__file__), 'nplinker_webapp.toml'))
        if not self.nplinker.load_data():
            raise Exception('Failed to load data')

        if not self.nplinker.process_dataset():
            raise Exception('Failed to process dataset')

        # load the BGC TSNE csv file from the webapp's /data dir
        self.bgc_data = {'x': [], 'y': [], 'strain': [], 'name': [], 'gcf': []}
        uniq_gcfs = set()
        gcf_lookup = {}
        with open(os.path.join(os.path.dirname(__file__), 'data/crusemann-bgc-tsne.csv'), 'r') as csvfile:
            csvr = csv.reader(csvfile)
            for l in csvr:
                (name, x, y) = l
                # NOTE cast to float or nothing will appear!
                if self.nplinker.has_bgc(name):
                    self.bgc_data['x'].append(float(x))
                    self.bgc_data['y'].append(float(y))
                    self.bgc_data['name'].append(name)
                    self.bgc_data['strain'].append(name) # TODO
                    # self.bgc_data['radius'].append(0.45) 
                    gcf = self.nplinker.lookup_bgc(name).parent.id
                    if gcf not in uniq_gcfs:
                        gcf_lookup[gcf] = len(uniq_gcfs)
                        uniq_gcfs.add(gcf)
                    self.bgc_data['gcf'].append(gcf)

                # else:
                #     print('Missing: {}'.format(name))


        # load the spectra TSNE csv file from the webapp's /data dir
        self.spec_data = {'x': [], 'y': [], 'name': [], 'family': []}
        uniq_fams = set()
        fam_lookup = {}
        with open(os.path.join(os.path.dirname(__file__), 'data/crusemann-spectra-tsne.csv')) as csvfile:
            csvr = csv.reader(csvfile)
            for l in csvr:
                (name, x, y) = l
                # self.spec_data['radius'].append(0.4)
                spec = self.nplinker.lookup_spectrum(name)
                if spec is None:
                    print('*** LOOKUP FAILED: spec name={}'.format(name))
                else:
                    # family = spec.family.family_id
                    family = spec.family
                    self.spec_data['name'].append(name)
                    self.spec_data['x'].append(float(x))
                    self.spec_data['y'].append(float(y))
                    self.spec_data['family'].append(family)
                    if family not in uniq_fams:
                        fam_lookup[family] = len(uniq_fams)
                        uniq_fams.add(family)

        print('Unique families: {}'.format(len(uniq_fams)))
        total = len(uniq_fams)
        cmap = []
        while total > 0:
            c = bkp.inferno(min(total, 256))
            total -= len(c)
            cmap.extend(c)

        self.spec_data['fill'] = []
        for i in range(len(self.spec_data['name'])):
            # if self.spec_data['family'][i] == '-1':
            #     self.spec_data['fill'].append('#ffeeeeee')
            # else:
            self.spec_data['fill'].append(cmap[fam_lookup[self.spec_data['family'][i]]])


        total = len(uniq_gcfs)
        cmap = []
        print('Unique GCFs: {}'.format(len(uniq_gcfs)))
        while total > 0:
            c = bkp.viridis(min(total, 256))
            total -= len(c)
            cmap.extend(c)

        self.bgc_data['fill'] = []
        for i in range(len(self.bgc_data['name'])):
            self.bgc_data['fill'].append(cmap[gcf_lookup[self.bgc_data['gcf'][i]]])
            


        self.bgc_indices = {}
        self.spec_indices = {}
        for i, bgc in enumerate(self.bgc_data['name']):
            self.bgc_indices[bgc] = i
        for i, spec in enumerate(self.spec_data['name']):
            self.spec_indices[spec] = i

nh = NPLinkerHelper()

def on_server_loaded(server_context):
    """
    This is called when the server is initially loaded, so use it to load the dataset,
    generate scores and store the results in an instance of NPLinkerHelper so that
    it only needs to be done once.
    """
    print('on_server_loaded')
    nh.load()

def on_session_created(session_context):
    """
    This is called for every new session that is created. To make the already-loaded
    objects visible in the context of the bokeh document, access it through the 
    context object and add the NPLinkerHelper instance as an attribute.
    """
    print('on_session_created')
    setattr(session_context._document, 'nh', nh)
