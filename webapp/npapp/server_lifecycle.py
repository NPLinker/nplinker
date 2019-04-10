import sys
import os
import csv

import bokeh.palettes as bkp

# add path to nplinker/prototype 
# TODO probably will need changed at some point!
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../prototype'))

from nplinker import NPLinker
from genomics import load_mibig_map, MiBIGBGC
from logconfig import LogConfig

class NPLinkerHelper(object):
    """
    This is just a simple class to wrap up the various objects involved in 
    loading and plotting a dataset via nplinker/bokeh
    """
    def __init__(self):
        LogConfig.setLogLevelStr('DEBUG')

    def load(self):
        data_dir = os.path.join(os.path.dirname(__file__), 'data')

        # initialise nplinker and load the dataset, using a config file in the webapp dir
        self.nplinker = NPLinker(os.path.join(os.path.dirname(__file__), 'nplinker_webapp.toml'))
        if not self.nplinker.load_data():
            raise Exception('Failed to load data')

        # TODO this should be handled better/elsewhere?
        # hacky
        mibig_file = os.path.join(data_dir, 'mibig_gnps_links_q1.csv')
        mibig_map = load_mibig_map(mibig_file)
        # mibig_bgcs = [MiBIGBGC(name, product) for (name, product) in mibig_map.items()]
        # c = len(self.nplinker.bgcs)
        # self.nplinker._bgcs.extend(mibig_bgcs)
        # for i, mb in enumerate(mibig_bgcs):
        #     self.nplinker._bgc_lookup[mb.name] = i + c

        if not self.nplinker.process_dataset():
            raise Exception('Failed to process dataset')

        # load each of the BGC TSNE csv files from the /data dir
        # these are all named crusemann-bgc-tsne-<name>.csv
        self.bgc_data = {}
        self.available_gcfs = {}

        for f in os.listdir(data_dir):
            if not f.startswith('crusemann-bgc-tsne-'):
                continue

            fname = os.path.join(data_dir, f)
            print('Loading BGC TSNE from {}'.format(fname))
            fid = f.replace('crusemann-bgc-tsne-', '')[:-4]

            bgc_data = {'x': [], 'y': [], 'strain': [], 'name': [], 'gcf': []}
            uniq_gcfs = set()
            gcf_lookup = {}
            self.available_gcfs[fid] = set()
            
            with open(fname, 'r') as csvfile:
                csvr = csv.reader(csvfile)
                headers = next(csvr)
                missing = 0
                for l in csvr:
                    if len(l) == 4:
                        # ignore bgc type col in some files
                        (name, x, y, _) = l
                    else:
                        (name, x, y) = l

                    # NOTE cast to float or nothing will appear!
                    if self.nplinker.has_bgc(name):
                        bgc_data['x'].append(float(x))
                        bgc_data['y'].append(float(y))
                        bgc_data['name'].append(name)
                        bgc_data['strain'].append(name) # TODO
                        gcf_obj = self.nplinker.lookup_bgc(name).parent
                        self.available_gcfs[fid].add(gcf_obj)
                        gcf_id = gcf_obj.id
                        if gcf_id not in uniq_gcfs:
                            gcf_lookup[gcf_id] = len(uniq_gcfs)
                            uniq_gcfs.add(gcf_id)
                        bgc_data['gcf'].append(gcf_id)
                    else:
                        missing += 1
                        # print('Missing: {}'.format(name))

                total = len(uniq_gcfs)
                cmap = []
                while total > 0:
                    c = bkp.d3['Category20c'][20]
                    total -= len(c)
                    cmap.extend(c)

                bgc_data['fill'] = []
                for i in range(len(bgc_data['name'])):
                    bgc_data['fill'].append(cmap[gcf_lookup[bgc_data['gcf'][i]]])
            
                self.bgc_data[fid] = bgc_data

                if missing > 0:
                    print('{} had {} missing BGCs!'.format(fid, missing))

        # load the spectra TSNE csv file from the webapp's /data dir
        self.spec_data = {'x': [], 'y': [], 'name': [], 'family': []}
        uniq_fams = set()
        fam_lookup = {}
        with open(os.path.join(os.path.dirname(__file__), 'data/crusemann-spectra-tsne.csv')) as csvfile:
            csvr = csv.reader(csvfile)
            headers = next(csvr)
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
            c = bkp.d3['Category20'][20]
            total -= len(c)
            cmap.extend(c)

        self.spec_data['fill'] = []
        for i in range(len(self.spec_data['name'])):
            # fix singletons to a single obvious colour
            if self.spec_data['family'][i] == '-1':
                self.spec_data['fill'].append('#000000')
            else:
                self.spec_data['fill'].append(cmap[fam_lookup[self.spec_data['family'][i]]])

        self.bgc_indices = {}
        self.spec_indices = {}
        for i, spec in enumerate(self.spec_data['name']):
            self.spec_indices[spec] = i

        for fid in self.bgc_data.keys():
            self.bgc_indices[fid] = {}
            for i, bgc_name in enumerate(self.bgc_data[fid]['name']):
                self.bgc_indices[fid][bgc_name] = i

        # provide a way to quickly look up the list of GCFs containing a particular BGC
        self.bgc_gcf_lookup = {}
        for gcf in self.nplinker.gcfs:
            for bgc in gcf.bgc_list:
                if bgc in self.bgc_gcf_lookup:
                    self.bgc_gcf_lookup[bgc].add(gcf)
                else:
                    self.bgc_gcf_lookup[bgc] = set([gcf])

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
