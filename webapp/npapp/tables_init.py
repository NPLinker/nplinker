import os
import pickle

import pandas as pd

from bokeh.models import Button, CustomJS, ColumnDataSource
from bokeh.models.widgets import DataTable, TableColumn, Div

from linker import Linker
from tables_functions import create_links, get_table_info, NA_ID, NA_TEXT

class TableData(object):

    def __init__(self, nh):
        self.nh = nh

    def setup(self):
        """
        This method sets up all the data structures needed to manage the tables
        interface to the webapp
        """

        self.local_cache = os.path.join(os.getenv('HOME'), 'nplinker_data', 'tables_links', self.nh.nplinker.dataset_id)
        os.makedirs(self.local_cache, exist_ok=True)

        pickled_scores_path = os.path.join(self.local_cache, 'scores.pckl')
        last_cutoff_path = os.path.join(self.local_cache, 'cutoff.pckl')
        pickled_links_path = os.path.join(self.local_cache, 'links.pckl')
        database_path = os.path.join(self.local_cache, 'linker.sqlite3')

        # need to handle the user changing the metcalf cutoff for the tables 
        # in the config file between runs here, so the last used value is stored
        # and compared against the current one. if it changes, have to regenerate
        # all the data structs below
        last_cutoff_val = None
        current_cutoff_val = self.nh.nplinker._loader.webapp_scoring_cutoff()
        if os.path.exists(last_cutoff_path):
            last_cutoff_val = pickle.load(open(last_cutoff_path, 'rb'))

        if last_cutoff_val != current_cutoff_val:
            print('*** Metcalf cutoff for tables has been changed (old={}, new={}), will regenerate data structs'.format(last_cutoff_val, current_cutoff_val))
            for path in [pickled_scores_path, pickled_links_path, database_path]:
                if os.path.exists(path):
                    os.unlink(path)
        pickle.dump(current_cutoff_val, open(last_cutoff_path, 'wb'))

        if not os.path.exists(pickled_scores_path):
            # get the metcalf scoring links between spectra:gcf and gcf:spectra pairs
            metcalf_scoring = self.nh.nplinker.scoring_method('metcalf')
            metcalf_scoring.cutoff = self.nh.nplinker._loader.webapp_scoring_cutoff()
            print('Metcalf tables cutoff is {}'.format(metcalf_scoring.cutoff))
            spec_links = self.nh.nplinker.get_links(self.nh.nplinker.spectra, metcalf_scoring)
            pickle.dump(spec_links, open(pickled_scores_path, 'wb'))
        else:
            spec_links = pickle.load(open(pickled_scores_path, 'rb'))


        spectra = list(spec_links.links.keys())
        gcfs = spec_links.get_all_targets()
        bgcs = []
        molfams = []

        print('len(spec_links) = {}, gcfs={}'.format(len(spec_links), len(gcfs)))

        _spec_indices, _gcf_indices, _bgc_indices, _molfam_indices = {NA_ID: 0}, {NA_ID: 0}, {NA_ID: 0}, {NA_ID: 0}

        # construct pandas dataframes and bokeh datatables for each object class
        # (this is probably fast enough to not be worth pickling like the links)
        self.bgc_data = dict(bgc_pk=[NA_ID], name=[NA_TEXT], product_type=[NA_TEXT])
        self.gcf_data = dict(gcf_pk=[NA_ID], gcf_id=[NA_TEXT], product_type=[NA_TEXT], nbgcs=[NA_TEXT], metcalf_score=[NA_TEXT])
        self.spec_data = dict(spec_pk=[NA_ID], spectrum_id=[NA_TEXT], family=[NA_TEXT], metcalf_score=[NA_TEXT])
        self.molfam_data = dict(molfam_pk=[NA_ID], family_id=[NA_TEXT], nspectra=[NA_TEXT])

        gcfs = sorted(gcfs, key=lambda gcf: gcf.id)
        bgcs_seen = set()
        for gcf in gcfs:
            self.gcf_data['gcf_pk'].append(gcf.id)
            self.gcf_data['gcf_id'].append(gcf.gcf_id)
            self.gcf_data['product_type'].append(gcf.product_type)
            self.gcf_data['nbgcs'].append(len(gcf.bgcs))
            self.gcf_data['metcalf_score'].append(NA_TEXT)
            _gcf_indices[gcf] = len(_gcf_indices)
            for bgc in gcf.bgcs:
                if bgc not in bgcs_seen:
                    bgcs_seen.add(bgc)
                    _bgc_indices[bgc] = len(bgcs_seen)
                    bgcs.append(bgc)

        bgcs = sorted(bgcs, key=lambda bgc: bgc.id)
        for bgc in bgcs:
            self.bgc_data['bgc_pk'].append(bgc.id)
            self.bgc_data['name'].append(bgc.name)
            self.bgc_data['product_type'].append(bgc.product_prediction)

        molfams_seen = set()
        for spec in spectra:
            self.spec_data['spec_pk'].append(spec.id)
            self.spec_data['spectrum_id'].append(spec.spectrum_id)
            self.spec_data['family'].append(spec.family.family_id)
            self.spec_data['metcalf_score'].append(NA_TEXT)
            _spec_indices[spec] = len(_spec_indices)
            if spec.family not in molfams_seen:
                self.molfam_data['molfam_pk'].append(spec.family.id)
                self.molfam_data['family_id'].append(spec.family.family_id)
                self.molfam_data['nspectra'].append(len(spec.family.spectra))
                molfams_seen.add(spec.family)
                _molfam_indices[spec.family] = len(molfams)
                molfams.append(spec.family)


        self.bgc_df = pd.DataFrame(self.bgc_data)
        self.gcf_df = pd.DataFrame(self.gcf_data)
        self.spec_df = pd.DataFrame(self.spec_data)
        self.molfam_df = pd.DataFrame(self.molfam_data)

        bgcs = list(bgcs)
        molfams = list(molfams)

        print('DataFrames created')

        # create links between GCF and BGC objects (or load pickled version)
        # note the +1s are because the tables code uses 0 for NA
        if not os.path.exists(pickled_links_path):
            index_mappings_1, index_mappings_2 = {}, {}

            print('Constructing link data structures')

            bgc_to_gcf_indices = []
            for gcf in gcfs:
                for bgc in gcf.bgcs:
                    if bgc not in bgcs_seen:
                        continue
                    # bgc_to_gcf_indices.append((bgc.id + 1, gcf.id + 1))
                    bgc_to_gcf_indices.append((_bgc_indices[bgc] + 0, _gcf_indices[gcf] + 0))
                    index_mappings_1[_bgc_indices[bgc] + 0] = bgc.id
                    index_mappings_2[_gcf_indices[gcf] + 0] = gcf.id

            self.bgc_gcf = create_links(self.bgc_df, self.gcf_df, 'bgc_pk', 'gcf_pk', bgc_to_gcf_indices, index_mappings_1, index_mappings_2)
            print(' + bgc_gcf')

            # links between Spectrum and MolFam objects
            # note the +1s are because the tables code uses 0 for NA
            index_mappings_1, index_mappings_2 = {}, {}
            molfam_to_spec_indices = []
            tmp = set(spectra)
            # for molfam in self.nh.nplinker.molfams:
            for molfam in molfams:
                for spec in molfam.spectra:
                    # molfam_to_spec_indices.append((molfam.id + 0, spec.id + 0))
                    if spec not in tmp:
                        continue
                    molfam_to_spec_indices.append((_molfam_indices[molfam] + 0, _spec_indices[spec] + 0))
                    index_mappings_1[_molfam_indices[molfam] + 0] = molfam.id
                    index_mappings_2[_spec_indices[spec] + 0] = spec.id
            self.mf_spectra = create_links(self.molfam_df, self.spec_df, 'molfam_pk', 'spec_pk', molfam_to_spec_indices, index_mappings_1, index_mappings_2)
            print(' + mf_spectra')

            # links between GCF<==>Spectrum objects based on metcalf scores, via BGCs
            index_mappings_1, index_mappings_2 = {}, {}
            spec_to_bgc_indices = []
            tmp = set()
            tmpbgcs = {}
            for spec, result in spec_links.links.items():
                for gcf, data in result.items():
                    #print('spec {} <==> gcf {}'.format(spectrum.id, gcf.id))
                    for bgc in gcf.bgcs:
                        # spec_to_bgc_indices.append((spec.id + 0, bgc.id + 0))
                        # avoid dup entries
                        if (_spec_indices[spec] + 0, _bgc_indices[bgc] + 0) in tmp:
                            continue
                        spec_to_bgc_indices.append((_spec_indices[spec] + 0, _bgc_indices[bgc] + 0))
                        tmp.add((_spec_indices[spec] + 0, _bgc_indices[bgc] + 0))
                        index_mappings_1[_spec_indices[spec] + 0] = spec.id
                        index_mappings_2[_bgc_indices[bgc] + 0] = bgc.id
            self.spectra_bgc = create_links(self.spec_df, self.bgc_df, 'spec_pk', 'bgc_pk', spec_to_bgc_indices, index_mappings_1, index_mappings_2)
            print(' + spectra_bgc')

            # pickle the data structs to avoid doing the above again
            pickle.dump((self.bgc_gcf, self.mf_spectra, self.spectra_bgc), open(pickled_links_path, 'wb'))
        else:
            self.bgc_gcf, self.mf_spectra, self.spectra_bgc = pickle.load(open(pickled_links_path, 'rb'))

        # combine the data and link informations into a list of dicts in the format the
        # linker class requires 
        self.table_info = get_table_info(self.molfam_df, self.mf_spectra, self.spec_df, self.spectra_bgc, self.bgc_df, self.bgc_gcf, self.gcf_df)
        Linker(self.table_info, self.nh.nplinker.dataset_id, database_path, do_init=True)

        print('TableData setup completed')

class TableSessionData(object):

    def __init__(self, tabledata):
        self.data = tabledata
        self.widgets = []

    def setup(self):
        # currently setting fixed widths here, not ideal, but using auto-sizing seems to
        # place limits on how far you can resize the columns for some reason...
        self.bgc_cols = [TableColumn(field='bgc_pk', title='ID', width=15), 
                         TableColumn(field='name', title='BGC name', width=220), 
                         TableColumn(field='product_type', title='Product type', width=75)]

        self.gcf_cols = [TableColumn(field='gcf_pk', title='ID', width=15), 
                         TableColumn(field='gcf_id', title='GCF ID', width=75), 
                         TableColumn(field='product_type', title='Product type', width=75),
                         TableColumn(field='nbgcs', title='#bgcs', width=75), 
                         TableColumn(field='metcalf_score', title='Metcalf score', width=75)]

        self.spec_cols = [TableColumn(field='spec_pk', title='ID', width=15), 
                          TableColumn(field='spectrum_id', title='Spectrum ID', width=75), 
                          TableColumn(field='family', title='Family ID', width=75),
                          TableColumn(field='metcalf_score', title='Metcalf score', width=75)]

        self.molfam_cols = [TableColumn(field='molfam_pk', title='ID', width=15), 
                            TableColumn(field='family_id', title='MolFam ID', width=75), 
                            TableColumn(field='nspectra', title='#spectra', width=75)]

        self.bgc_ds = ColumnDataSource(self.data.bgc_df)
        self.gcf_ds = ColumnDataSource(self.data.gcf_df)
        self.spec_ds = ColumnDataSource(self.data.spec_df)
        self.molfam_ds = ColumnDataSource(self.data.molfam_df)

        self.bgc_dt = DataTable(source=self.bgc_ds, columns=self.bgc_cols, name='table_bgcs', sizing_mode='scale_width', width=300, fit_columns=False, index_width=15)
        self.gcf_dt = DataTable(source=self.gcf_ds, columns=self.gcf_cols, name='table_gcfs', sizing_mode='scale_width', width=300, fit_columns=False, index_width=15)
        self.spec_dt = DataTable(source=self.spec_ds, columns=self.spec_cols, name='table_spectra', sizing_mode='scale_width', width=300, fit_columns=False, index_width=15)
        self.molfam_dt = DataTable(source=self.molfam_ds, columns=self.molfam_cols, name='table_molfams', sizing_mode='scale_width', width=300, fit_columns=False, index_width=15)

        # start building a list of widgets that the main module will add to the
        # bokeh document object
        self.widgets = []
        self.widgets.append(self.bgc_dt)
        self.widgets.append(self.gcf_dt)
        self.widgets.append(self.spec_dt)
        self.widgets.append(self.molfam_dt)

        data_sources = {
            'molfam_table': self.molfam_ds,
            'spec_table': self.spec_ds,
            'bgc_table': self.bgc_ds,
            'gcf_table': self.gcf_ds
        }
        self.data_sources = data_sources

        self.data_tables = {
            'molfam_table': self.molfam_dt, 
            'spec_table': self.spec_dt,
            'bgc_table': self.bgc_dt, 
            'gcf_table': self.gcf_dt
        }

        # https://stackoverflow.com/questions/31824124/is-there-a-way-to-save-bokeh-data-table-content
        download_button_code = """
            function table_to_csv(source) {
                const columns = Object.keys(source.data)
                const nrows = source.get_length()
                const lines = [columns.join(',')]

                for (let i = 0; i < nrows; i++) {
                    let row = [];
                    for (let j = 0; j < columns.length; j++) {
                        const column = columns[j]
                        row.push(source.data[column][i].toString())
                    }
                    lines.push(row.join(','))
                }
                return lines.join('\\n').concat('\\n')
            }

            const filename = name + '_data.csv'
            filetext = table_to_csv(source)
            const blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' })

            //addresses IE
            if (navigator.msSaveBlob) {
                navigator.msSaveBlob(blob, filename)
            } else {
                const link = document.createElement('a')
                link.href = URL.createObjectURL(blob)
                link.download = filename
                link.target = '_blank'
                link.style.visibility = 'hidden'
                link.dispatchEvent(new MouseEvent('click'))
            }
        """

        self.molfam_dl_button = Button(label='Download as CSV', name='molfam_dl_button')
        self.spec_dl_button = Button(label='Download as CSV', name='spec_dl_button')
        self.bgc_dl_button = Button(label='Download as CSV', name='bgc_dl_button')
        self.gcf_dl_button = Button(label='Download as CSV', name='gcf_dl_button')

        self.molfam_dl_button.callback = CustomJS(args=dict(name='molfam', source=self.molfam_ds), code=download_button_code)
        self.spec_dl_button.callback = CustomJS(args=dict(name='spectrum', source=self.spec_ds), code=download_button_code)
        self.bgc_dl_button.callback = CustomJS(args=dict(name='bgc', source=self.bgc_ds), code=download_button_code)
        self.gcf_dl_button.callback = CustomJS(args=dict(name='gcf', source=self.gcf_ds), code=download_button_code)

        self.widgets.append(self.molfam_dl_button)
        self.widgets.append(self.spec_dl_button)
        self.widgets.append(self.bgc_dl_button)
        self.widgets.append(self.gcf_dl_button)

        self.resetting = False

        def table_callback(name):
            cb_obj = self.data_sources[name]
            selected_indices = cb_obj.selected.indices
            print('table {}, selected = {}'.format(name, selected_indices))
            if self.resetting:
                print('skipping remainder of callback, resetting is True')
                return

            self.linker.removeConstraints(name)

            if name in ['bgc_table', 'gcf_table']:
                self.molfam_dt.selectable = len(selected_indices) == 0
                self.spec_dt.selectable = len(selected_indices) == 0
                self.bgc_dt.selectable = True
                self.spec_dt.selectable = True
            else:
                self.molfam_dt.selectable = True
                self.spec_dt.selectable = True
                self.gcf_dt.selectable = len(selected_indices) == 0
                self.bgc_dt.selectable = len(selected_indices) == 0

            for idx in selected_indices:
                self.linker.addConstraint(name, idx, self.data_sources[name])

            print('table_callback: query')
            self.linker.query()
            print('table_callback: updating data sources')
            self.linker.updateDataSources(self.data_sources)
            print('table_callback {} done'.format(name))

        def reset_tables_callback(unused):
            print('reset callback')
            # this is used to prevent the table_callback being triggered during
            # this callback, which doesn't break anything but does repeat some
            # operations and waste time
            self.resetting = True

            for table_name, table_ds in self.data_sources.items():
                table_ds.selected.indices = []
                self.linker.removeConstraints(table_name)

            self.linker.query()
            self.linker.updateDataSources(self.data_sources)
            print('reset callback done')

            self.molfam_dt.selectable = True
            self.spec_dt.selectable = True
            self.gcf_dt.selectable = True
            self.bgc_dt.selectable = True
            self.tables_score_met.disabled = False
            self.tables_score_gen.disabled = False

            self.resetting = False

        self.molfam_ds.selected.on_change('indices', lambda a, b, c: table_callback('molfam_table'))
        self.spec_ds.selected.on_change('indices', lambda a, b, c: table_callback('spec_table'))
        self.bgc_ds.selected.on_change('indices', lambda a, b, c: table_callback('bgc_table'))
        self.gcf_ds.selected.on_change('indices', lambda a, b, c: table_callback('gcf_table'))

        # the buttons to trigger scoring
        self.tables_score_met = Button(label='Show scores for selected spectra', name='tables_score_met', button_type='success')
        self.tables_score_gen = Button(label='Show scores for selected BGCs', name='tables_score_gen', button_type='success')
        self.widgets.append(self.tables_score_met)
        self.widgets.append(self.tables_score_gen)

        selection_change_callback_code = """
            var value = '50%';
            if(source.selected.indices.length == 0)
                value = '';

            if (name === 'bgc_table' || name === 'gcf_table') {
                $('#spec_table').css('opacity', value);
                $('#molfam_table').css('opacity', value);
                btn_met.disabled = true;
                btn_gen.disabled = false;
            } else {
                $('#bgc_table').css('opacity', value);
                $('#gcf_table').css('opacity', value);
                btn_met.disabled = false;
                btn_gen.disabled = true;
            }
        """
        self.molfam_ds.selected.js_on_change('indices', CustomJS(code=selection_change_callback_code, args=dict(name='molfam_table', btn_met=self.tables_score_met, btn_gen=self.tables_score_gen, source=self.molfam_ds)))
        self.spec_ds.selected.js_on_change('indices', CustomJS(code=selection_change_callback_code, args=dict(name='spec_table', btn_met=self.tables_score_met, btn_gen=self.tables_score_gen, source=self.spec_ds)))
        self.bgc_ds.selected.js_on_change('indices', CustomJS(code=selection_change_callback_code, args=dict(name='bgc_table', btn_met=self.tables_score_met, btn_gen=self.tables_score_gen, source=self.bgc_ds)))
        self.gcf_ds.selected.js_on_change('indices', CustomJS(code=selection_change_callback_code, args=dict(name='gcf_table', btn_met=self.tables_score_met, btn_gen=self.tables_score_gen, source=self.gcf_ds)))

        # CustomJS callback code for resetting the selection (and opacity) states
        # of all the tables 
        reset_selection_code = """    
            // remove selections from all tables
            var tablenames = ['molfam_table', 'spec_table', 'bgc_table', 'gcf_table'];
            for(let t=0;t<tablenames.length;t++) {
                const table = tablenames[t];
                // restore opacity
                $('#' + table).css('opacity', '');
            }
        """
        self.tables_reset = Button(label='Clear selections', name='tables_reset', button_type='danger')
        self.tables_reset.js_on_click(CustomJS(args=dict(), code=reset_selection_code))
        self.tables_reset.on_click(reset_tables_callback)
        self.widgets.append(self.tables_reset)

        self.linker = Linker(self.data.table_info, self.data.nh.nplinker.dataset_id, os.path.join(self.data.local_cache, 'linker.sqlite3'))
