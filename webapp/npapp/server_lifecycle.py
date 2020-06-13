import sys
import os
import csv
import time

import bokeh.palettes as bkp
from bokeh.models import ColumnDataSource

# add path to nplinker/prototype (for local dev)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../prototype'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../nplinker/prototype'))

from nplinker.nplinker import NPLinker
from nplinker.genomics import load_mibig_map, MiBIGBGC
from nplinker.metabolomics import SingletonFamily
from nplinker.logconfig import LogConfig
from nplinker.layout import create_genomics_graph, create_metabolomics_graph

from tables_init import TableData, TableSessionData

class NPLinkerHelper(object):
    """
    This is just a simple class to wrap up the various objects involved in 
    loading and plotting a dataset via nplinker/bokeh
    """
    def __init__(self):
        LogConfig.setLogLevelStr('DEBUG')
        self.nx_scale = 11
        self.data_dir = os.path.join(os.path.dirname(__file__), 'data')

    def _construct_bgc_data(self, fid, nodes_xy, edges_start, edges_end):
        # self.bgc_data in the TSNE code is a dict of dicts, indexed by TSNE
        # filename at the top level. each entry within that is a dict with 
        # columndatasource structure: a list for each of x,y,strain,name,gcf 
        # columns, which can be passed directly to ColumnDataSource constructor
        # here can simplify a bit as only one "layout" to deal with
        self.bgc_data = {}
        # this is usually a dict with one entry per TSNE filename, and contains
        # a set of available GCF objects within each TSNE projection
        self.available_gcfs = {fid: set()}
        bgc_data = {'index': [], 'x': [], 'y': [], 'strain': [], 'name': [], 'gcf': [], 'gcf_name': [], 'mibig': [], 'prodtype': []}
        uniq_gcfs = set()
        gcf_lookup = {}
        bgcindex = 0
        for bgc in self.nplinker.bgcs:
            if bgc.id not in nodes_xy:
                print('??? {}'.format(bgc))
                continue

            bgc_xy = nodes_xy[bgc.id]
            bgc_data['index'].append(bgcindex)
            bgc_data['x'].append(bgc_xy[0])
            bgc_data['y'].append(bgc_xy[1])
            
            bgc_data['name'].append(bgc.name)
            # TODO is this right thing to do here?
            bgc_data['strain'].append(bgc.strain.id if bgc.strain is not None else 'MiBIGBGC') 
            bgc_data['mibig'].append(isinstance(bgc, MiBIGBGC))
            bgcindex += 1
            gcf_obj = bgc.parent
            self.available_gcfs[fid].add(gcf_obj)
            gcf_id = gcf_obj.id
            if gcf_id not in uniq_gcfs:
                gcf_lookup[gcf_id] = len(uniq_gcfs)
                uniq_gcfs.add(gcf_id)
            bgc_data['gcf'].append(gcf_id)
            bgc_data['gcf_name'].append(gcf_obj.gcf_id)
            bgc_data['prodtype'].append(gcf_obj.product_type)

        total = len(uniq_gcfs)
        cmap = []
        # colour by GCF
        # while total > 0:
        #     c = bkp.d3['Category20c'][20]
        #     total -= len(c)
        #     cmap.extend(c)
        #bgc_data['fill'] = []
        #for i in range(len(bgc_data['name'])):
        #    if bgc_data['mibig'][i] == True:
        #        bgc_data['fill'].append('#e2e0c0')
        #    else:
        #        bgc_data['fill'].append(cmap[gcf_lookup[bgc_data['gcf'][i]]])

        # colour by product type
        bgc_data['fill'] = []
        prodtypes = self.nplinker.product_types
        cmap = bkp.d3['Category10'][len(prodtypes)+1]
        for i in range(len(bgc_data['name'])):
            if bgc_data['mibig'][i] == True:
                bgc_data['fill'].append(cmap[-1])
            else:
                bgc_data['fill'].append(cmap[prodtypes.index(bgc_data['prodtype'][i])])
            
        self.bgc_cmap = {pt: cmap[prodtypes.index(pt)] for pt in prodtypes}
        self.bgc_data[fid] = bgc_data

        bgc_edge_data = {'start': [], 'end': [], 'mibig': [], 'xs': [], 'ys': []}
        bgc_edge_lookup = {}
        for i in range(len(edges_start)):
            start_index = edges_start[i]
            end_index = edges_end[i]

            li = len(bgc_edge_data['xs'])
            for index in [start_index, end_index]:
                if index not in bgc_edge_lookup:
                    bgc_edge_lookup[index] = [(li, start_index if index == end_index else start_index)]
                else:
                    bgc_edge_lookup[index].append((li, start_index if index == end_index else start_index))

            bgc_edge_data['mibig'].append(bgc_data['mibig'][start_index] or bgc_data['mibig'][end_index])

            bgc_edge_data['start'].append(start_index)
            bgc_edge_data['end'].append(end_index)
            
            bgc_edge_data['xs'].append([bgc_data['x'][start_index], bgc_data['x'][end_index]])
            bgc_edge_data['ys'].append([bgc_data['y'][start_index], bgc_data['y'][end_index]])

        self.bgc_edge_data = {}
        self.bgc_edge_data[fid] = bgc_edge_data
        self.bgc_positions = nodes_xy
        self.bgc_edge_lookup = bgc_edge_lookup

    def _construct_spec_data(self, fid, nodes_xy, edges_start, edges_end):
        # self.spec_data is the same as self.bgc_data above, but holding spectrum
        # info rather than BGC info...
        self.spec_data = {'index': [], 'x': [], 'y': [], 'name': [], 'family': [], 'singleton': [], 'parent_mass': []}
        uniq_fams = set()
        fam_lookup = {}

        specindex = 0
        for spec in self.nplinker.spectra:
            family = spec.family.family_id
            if spec.id not in nodes_xy:
                continue
            spec_xy = nodes_xy[spec.id]
            self.spec_data['index'].append(specindex)
            self.spec_data['x'].append(spec_xy[0])
            self.spec_data['y'].append(spec_xy[1])
            self.spec_data['name'].append(spec.spectrum_id)
            self.spec_data['family'].append(family)
            self.spec_data['singleton'].append(isinstance(spec.family, SingletonFamily))
            self.spec_data['parent_mass'].append(spec.parent_mz)
            specindex += 1
            # precache JCAMP data
            spec.to_jcamp_str()
            if family not in uniq_fams:
                fam_lookup[family] = len(uniq_fams)
                uniq_fams.add(family)

        print('Unique families: {}'.format(len(uniq_fams)))
        total = len(uniq_fams)

        # colouring by MolFam
        # cmap = []
        # while total > 0:
        #     c = bkp.d3['Category20'][20]
        #     total -= len(c)
        #     cmap.extend(c)

        # self.spec_data['fill'] = []
        # for i in range(len(self.spec_data['name'])):
        #     # fix singletons to a single obvious colour
        #     if self.spec_data['family'][i] == '-1':
        #         self.spec_data['fill'].append('#000000')
        #     else:
        #         self.spec_data['fill'].append(cmap[fam_lookup[self.spec_data['family'][i]]])

        # colouring by parent mass
        parent_mass_categories = [(None, 150), (150, 300), (300, 500), (500, 700), (700, 900), (900, 1100), (1100, None)]
        def mass_to_cmap_index(pm):
            for i in range(len(parent_mass_categories)):
                pmin, pmax = parent_mass_categories[i]
                if pmax is None:
                    return i
                if pm >= pmax:
                    continue
                if (pmin is not None and pm > pmin) or pmin is None:
                    return i

        cmap = bkp.d3['Category10'][7]
        self.spec_data['fill'] = []
        for i in range(len(self.spec_data['name'])):
            self.spec_data['fill'].append(cmap[mass_to_cmap_index(self.spec_data['parent_mass'][i])])
        self.spec_cmap = []
        for i, limits in enumerate(parent_mass_categories):
            pmin, pmax = limits
            if pmin is not None and pmax is not None:
                self.spec_cmap.append(('{:.0f}--{:.0f}'.format(pmin, pmax), cmap[i]))
            elif pmin is None:
                self.spec_cmap.append(('< {:.0f}'.format(pmax), cmap[i]))
            elif pmax is None:
                self.spec_cmap.append(('> {:.0f}'.format(pmin), cmap[i]))

        spec_edge_data = {'start': [], 'end': [], 'singleton': [], 'xs': [], 'ys': []}
        spec_edge_lookup = {}
        for i in range(len(edges_start)):
            start_index = edges_start[i]
            end_index = edges_end[i]

            li = len(spec_edge_data['xs'])
            for index in [start_index, end_index]:
                if index not in spec_edge_lookup:
                    spec_edge_lookup[index] = [(li, start_index if index == end_index else start_index)]
                else:
                    spec_edge_lookup[index].append((li, start_index if index == end_index else start_index))

            spec_edge_data['singleton'].append(self.spec_data['singleton'][start_index] or self.spec_data['singleton'][end_index])

            spec_edge_data['start'].append(start_index)
            spec_edge_data['end'].append(end_index)

            spec_edge_data['xs'].append([self.spec_data['x'][start_index], self.spec_data['x'][end_index]])
            spec_edge_data['ys'].append([self.spec_data['y'][start_index], self.spec_data['y'][end_index]])

        print('# edges = {}'.format(len(spec_edge_data['start'])))
        self.spec_edge_data = spec_edge_data
        self.spec_positions = nodes_xy
        self.spec_edge_lookup = spec_edge_lookup

    def get_bgc_edges(self, selected_node_indices):
        """
        Given a list of selected BGC indices, return a list of edge indices
        so that only the appropriate subset of edges are shown as part of the 
        selection. 
        """
        edgeset = set()
        for n_index in selected_node_indices:
            if n_index in self.bgc_edge_lookup:
                # value is a list of (edge index, node index) tuples
                node_edge_info = self.bgc_edge_lookup[n_index]
                # to avoid "dangling" edges, also filter out any edges between 
                # nodes that are not both in the selected set
                edgeset.update([ni[0] for ni in node_edge_info if ni[1] in selected_node_indices and ni[1] != n_index])

        return list(edgeset)

    def get_spec_edges(self, selected_node_indices):
        """
        Given a list of selected spectrum indices, return a list of edge indices
        so that only the appropriate subset of edges are shown as part of the 
        selection. 
        """
        edgeset = set()
        for n_index in selected_node_indices:
            if n_index in self.spec_edge_lookup:
                # value is a list of (edge index, node index) tuples
                node_edge_info = self.spec_edge_lookup[n_index]
                # to avoid "dangling" edges, also filter out any edges between 
                # nodes that are not both in the selected set
                edgeset.update([ni[0] for ni in node_edge_info if ni[1] in selected_node_indices and ni[1] != n_index])

        return list(edgeset)

    def plot_x_range(self, extended=False, padding=1.2):
        return (-self.nx_scale * padding, self.nx_scale) 

    def plot_y_range(self, extended=False, padding=1.2):
        if not extended:
            return (-self.nx_scale * padding, self.nx_scale)

        return (-self.nx_scale * 3 * padding, self.nx_scale)

    def load(self):
        self.load_nplinker()

        self.load_genomics()
        self.load_metabolomics()

        self.load_tables()

        self.load_rosetta()

    def load_rosetta(self):
        # trigger rosetta scoring setup during initial load as the first run 
        # can take a significant amount of time. the results are cached and
        # pickled for subsequent use.
        print('Preprocessing for rosetta scoring...')
        rs_obj = self.nplinker.scoring_method('rosetta')
        print('Finished preprocessing for rosetta scoring')

    def load_nplinker(self):
        # initialise nplinker and load the dataset, using a config file in the webapp dir
        datapath = os.getenv('NPLINKER_CONFIG', os.path.join(os.path.join(os.path.dirname(__file__), 'nplinker_webapp.toml')))
        print('DATAPATH: {}'.format(datapath))
        self.nplinker = NPLinker(datapath)
        if not self.nplinker.load_data():
            raise Exception('Failed to load data')

    def load_metabolomics(self):
        print('Creating metabolomics graph')
        met_points, met_edges = create_metabolomics_graph(self.nplinker.spectra, width=self.nx_scale * 2, singletons=True, split_singletons=True)
        met_edge_start = met_edges['start']
        met_edge_end = met_edges['end']

        fid = '<networkx>' # to replace previous BGC TSNE file selection for now
        self._construct_spec_data(fid, met_points, met_edge_start, met_edge_end)

        self.spec_indices = {}
        for i, spec in enumerate(self.spec_data['name']):
            self.spec_indices[spec] = i

    def load_genomics(self):
        print('Creating genomics graph')
        gen_points, gen_edges = create_genomics_graph(self.nplinker.bgcs, width=self.nx_scale * 2, mibig=True, split_mibig=False)
        gen_edge_start = gen_edges['start']
        gen_edge_end = gen_edges['end']

        fid = '<networkx>' # to replace previous BGC TSNE file selection for now
        self._construct_bgc_data(fid, gen_points, gen_edge_start, gen_edge_end)

        self.bgc_indices = {}
        for fid in self.bgc_data.keys():
            self.bgc_indices[fid] = {}
            for i, bgc_name in enumerate(self.bgc_data[fid]['name']):
                self.bgc_indices[fid][bgc_name] = i

        # provide a way to quickly look up the list of GCFs containing a particular BGC
        self.bgc_gcf_lookup = {}
        for gcf in self.nplinker.gcfs:
            for bgc in gcf.bgcs:
                if bgc in self.bgc_gcf_lookup:
                    self.bgc_gcf_lookup[bgc].add(gcf)
                else:
                    self.bgc_gcf_lookup[bgc] = set([gcf])

    def load_tables(self):
        print('Loading tables data')
        self.table_data = TableData(self)
        self.table_data.setup()

    def load_tsne(self):
        # initialise nplinker and load the dataset, using a config file in the webapp dir
        self.nplinker = NPLinker(os.path.join(os.path.dirname(__file__), 'nplinker_webapp.toml'))
        if not self.nplinker.load_data():
            raise Exception('Failed to load data')

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
                        bgc = self.nplinker.lookup_bgc(name)
                        # TODO is this right thing to do here?
                        bgc_data['strain'].append(bgc.strain if bgc.strain is not None else 'MiBIGBGC') 

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
                    # precache JCAMP data
                    spec.to_jcamp_str()
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
                self.bgc_indices[fid][self.nplinker.lookup_bgc(bgc_name)] = i

        # provide a way to quickly look up the list of GCFs containing a particular BGC
        self.bgc_gcf_lookup = {}
        for gcf in self.nplinker.gcfs:
            for bgc in gcf.bgcs:
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

    print('==================================')
    print('NPLinker server loading completed!')
    print('==================================')

def on_session_created(session_context):
    """
    This is called for every new session that is created. To make the already-loaded
    objects visible in the context of the bokeh document, access it through the 
    context object and add the NPLinkerHelper instance as an attribute.
    """
    print('on_session_created')
    args = session_context.request.arguments
    # TODO this should ideally not change the cutoff for every session?
    if 'cutoff' in args:
        val = int(args['cutoff'][0])
        if val != nh.nplinker.bigscape_cutoff:
            print('*** Updating cutoff value to {} and reloading data'.format(val))
            nh.nplinker.load_data(new_bigscape_cutoff=val)
            nh.nplinker.process_dataset()
            nh.load_genomics()
    setattr(session_context._document, 'nh', nh)
    setattr(session_context._document, 'table_session_data', TableSessionData(nh.table_data))

def on_session_destroyed(session_context):
    print('on_session_destroyed')
    session_data = getattr(session_context._document, 'table_session_data')
    # closes the database connection used by the Linker class
    session_data.linker.close()
