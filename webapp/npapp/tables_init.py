import os
import pickle

import pandas as pd

from bokeh.models import Button, CustomJS, ColumnDataSource
from bokeh.models.widgets import DataTable, TableColumn, Div

from linker import Linker
from tables_functions import create_links, get_table_info, NA

class TableData(object):

    def __init__(self, nh):
        self.nh = nh

    def setup(self):
        """
        This method sets up all the data structures needed to manage the tables
        interface to the webapp
        """

        # first step: get the metcalf scoring links between spectra:gcf and gcf:spectra pairs
        self.spec_links = self.nh.nplinker.get_links(self.nh.nplinker.spectra, self.nh.nplinker.scoring_method('metcalf'))
        self.gcf_links = self.nh.nplinker.get_links(self.nh.nplinker.gcfs, self.nh.nplinker.scoring_method('metcalf'))

        # construct pandas dataframes and bokeh datatables for each object class
        # (this is probably fast enough to not be worth pickling like the links)
        self.bgc_data = dict(bgc_pk=[NA], name=[NA], product_type=[NA])
        self.gcf_data = dict(gcf_pk=[NA], gcf_id=[NA], product_type=[NA], nbgcs=[NA])
        self.spec_data = dict(spec_pk=[NA], spectrum_id=[NA], family=[NA])
        self.molfam_data = dict(molfam_pk=[NA], family_id=[NA], nspectra=[NA])

        for bgc in self.nh.nplinker.bgcs:
            self.bgc_data['bgc_pk'].append(bgc.id)
            self.bgc_data['name'].append(bgc.name)
            self.bgc_data['product_type'].append(bgc.product_prediction)

        for gcf in self.nh.nplinker.gcfs:
            self.gcf_data['gcf_pk'].append(gcf.id)
            self.gcf_data['gcf_id'].append(gcf.gcf_id)
            self.gcf_data['product_type'].append(gcf.product_type)
            self.gcf_data['nbgcs'].append(len(gcf.bgcs))

        for spec in self.nh.nplinker.spectra:
            self.spec_data['spec_pk'].append(spec.id)
            self.spec_data['spectrum_id'].append(spec.spectrum_id)
            self.spec_data['family'].append(spec.family.family_id)

        for molfam in self.nh.nplinker.molfams:
            self.molfam_data['molfam_pk'].append(molfam.id)
            self.molfam_data['family_id'].append(molfam.family_id)
            self.molfam_data['nspectra'].append(len(molfam.spectra))

        self.bgc_df = pd.DataFrame(self.bgc_data)
        self.gcf_df = pd.DataFrame(self.gcf_data)
        self.spec_df = pd.DataFrame(self.spec_data)
        self.molfam_df = pd.DataFrame(self.molfam_data)

        # currently setting fixed widths here, not ideal, but using auto-sizing seems to
        # place limits on how far you can resize the columns for some reason...
        self.bgc_cols = [TableColumn(field='bgc_pk', title='ID', width=15), 
                         TableColumn(field='name', title='BGC name', width=220), 
                         TableColumn(field='product_type', title='Product type', width=75)]

        self.gcf_cols = [TableColumn(field='gcf_pk', title='ID', width=15), 
                         TableColumn(field='gcf_id', title='GCF ID', width=75), 
                         TableColumn(field='product_type', title='Product type', width=75)]

        self.spec_cols = [TableColumn(field='spec_pk', title='ID', width=15), 
                          TableColumn(field='spectrum_id', title='Spectrum ID', width=76), 
                          TableColumn(field='family', title='Family ID', width=75)]

        self.molfam_cols = [TableColumn(field='molfam_pk', title='ID', width=15), 
                            TableColumn(field='family_id', title='MolFam ID', width=75), 
                            TableColumn(field='nspectra', title='#spectra', width=75)]

        self.bgc_ds = ColumnDataSource(self.bgc_df)
        self.gcf_ds = ColumnDataSource(self.gcf_df)
        self.spec_ds = ColumnDataSource(self.spec_df)
        self.molfam_ds = ColumnDataSource(self.molfam_df)

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

        self.molfam_dl_button = Button(label='Download MolFam data', name='molfam_dl_button')
        self.spec_dl_button = Button(label='Download Spectrum data', name='spec_dl_button')
        self.bgc_dl_button = Button(label='Download BGC data', name='bgc_dl_button')
        self.gcf_dl_button = Button(label='Download GCF data', name='gcf_dl_button')

        self.molfam_dl_button.callback = CustomJS(args=dict(name='molfam', source=self.molfam_ds), code=download_button_code)
        self.spec_dl_button.callback = CustomJS(args=dict(name='spectrum', source=self.spec_ds), code=download_button_code)
        self.bgc_dl_button.callback = CustomJS(args=dict(name='bgc', source=self.bgc_ds), code=download_button_code)
        self.gcf_dl_button.callback = CustomJS(args=dict(name='gcf', source=self.gcf_ds), code=download_button_code)

        self.widgets.append(self.molfam_dl_button)
        self.widgets.append(self.spec_dl_button)
        self.widgets.append(self.bgc_dl_button)
        self.widgets.append(self.gcf_dl_button)

        def table_callback(name):
            cb_obj = self.data_sources[name]
            selected_indices = cb_obj.selected.indices
            print('table {}, selected = {}'.format(name, selected_indices))
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

            self.linker.query()
            self.linker.updateDataSources(self.data_sources)

        def reset_tables_callback(unused):
            for table_name, table_ds in self.data_sources.items():
                table_ds.selected.indices = []
                self.linker.removeConstraints(table_name)

            self.linker.query()
            self.linker.updateDataSources(self.data_sources)

            self.molfam_dt.selectable = True
            self.spec_dt.selectable = True
            self.gcf_dt.selectable = True
            self.bgc_dt.selectable = True


        self.molfam_ds.selected.on_change('indices', lambda a, b, c: table_callback('molfam_table'))
        self.spec_ds.selected.on_change('indices', lambda a, b, c: table_callback('spec_table'))
        self.bgc_ds.selected.on_change('indices', lambda a, b, c: table_callback('bgc_table'))
        self.gcf_ds.selected.on_change('indices', lambda a, b, c: table_callback('gcf_table'))

        selection_change_callback_code = """
            var value = '50%';
            if(source.selected.indices.length == 0)
                value = '';

            if (name === 'bgc_table' || name === 'gcf_table') {
                $('#spec_table').css('opacity', value);
                $('#molfam_table').css('opacity', value);
            } else {
                $('#bgc_table').css('opacity', value);
                $('#gcf_table').css('opacity', value);
            }
        """
        self.molfam_ds.selected.js_on_change('indices', CustomJS(code=selection_change_callback_code, args=dict(name='molfam_table', source=self.molfam_ds)))
        self.spec_ds.selected.js_on_change('indices', CustomJS(code=selection_change_callback_code, args=dict(name='spec_table', source=self.spec_ds)))
        self.bgc_ds.selected.js_on_change('indices', CustomJS(code=selection_change_callback_code, args=dict(name='bgc_table', source=self.bgc_ds)))
        self.gcf_ds.selected.js_on_change('indices', CustomJS(code=selection_change_callback_code, args=dict(name='gcf_table', source=self.gcf_ds)))

        # pickle the links structs as they don't change for a given dataset and can take
        # some time to generate when the webapp is instantiated
        # TODO proper location for the pickled data file
        pickled_links = 'table_links_{}.pckl'.format(self.nh.nplinker.dataset_id)
        if not os.path.exists(pickled_links):
            # 1. links between GCF and BGC objects
            # note the +1s are because the tables code uses 0 for NA
            bgc_to_gcf_indices = []
            for gcf in self.nh.nplinker.gcfs:
                for bgc in gcf.bgcs:
                    bgc_to_gcf_indices.append((bgc.id + 1, gcf.id + 1))
            self.bgc_gcf = create_links(self.bgc_df, self.gcf_df, 'bgc_pk', 'gcf_pk', bgc_to_gcf_indices)

            # 2. links between Spectrum and MolFam objects
            # note the +1s are because the tables code uses 0 for NA
            molfam_to_spec_indices = []
            for molfam in self.nh.nplinker.molfams:
                for spec in molfam.spectra:
                    molfam_to_spec_indices.append((molfam.id + 1, spec.id + 1))
            self.mf_spectra = create_links(self.molfam_df, self.spec_df, 'molfam_pk', 'spec_pk', molfam_to_spec_indices)

            # 3. links between GCF<==>Spectrum objects based on metcalf scores, via BGCs
            spec_to_bgc_indices = []
            for spectrum in self.spec_links:
                # lookup the GCF that this spectrum has a link to 
                linked_gcfs = self.nh.nplinker.links_for_obj(spectrum, self.nh.nplinker.scoring.metcalf)
                for gcf, score in linked_gcfs:
                    #print('spec {} <==> gcf {}'.format(spectrum.id, gcf.id))
                    for bgc in gcf.bgcs:
                        spec_to_bgc_indices.append((spectrum.id + 1, bgc.id + 1))
            self.spectra_bgc = create_links(self.spec_df, self.bgc_df, 'spec_pk', 'bgc_pk', spec_to_bgc_indices)

            pickle.dump((self.bgc_gcf, self.mf_spectra, self.spectra_bgc), open(pickled_links, 'wb'))
        else:
            self.bgc_gcf, self.mf_spectra, self.spectra_bgc = pickle.load(open(pickled_links, 'rb'))

        # add the buttons to trigger scoring
        self.tables_score_met = Button(label='Show scores for selected spectra', name='tables_score_met', button_type='success')
        self.tables_score_gen = Button(label='Show scores for selected BGCs', name='tables_score_gen', button_type='success')
        self.widgets.append(self.tables_score_met)
        self.widgets.append(self.tables_score_gen)

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
        self.tables_reset.js_on_click(CustomJS(args={}, code=reset_selection_code))
        self.tables_reset.on_click(reset_tables_callback)
        self.widgets.append(self.tables_reset)

        # combine the data and link informations into a list of dicts in the format the
        # linker expects. this object is passed in to the linker by a callback after 
        # the webapp loads (see below & main.py)
        self.table_info = get_table_info(self.molfam_df, self.mf_spectra, self.spec_df, self.spectra_bgc, self.bgc_df, self.bgc_gcf, self.gcf_df)
        self.linker = Linker(self.table_info, self.nh.nplinker.dataset_id)
