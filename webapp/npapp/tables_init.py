import os
import pickle

import pandas as pd

from bokeh.models import Button, CustomJS, ColumnDataSource
from bokeh.models.widgets import DataTable, TableColumn

from tables_functions import create_links, get_table_info, NA

class TableData(object):

    def __init__(self, nh):
        self.nh = nh

    def setup(self):
        # get the metcalf scoring links between spectra:gcf and gcf:spectra pairs
        self.spec_links = self.nh.nplinker.get_links(self.nh.nplinker.spectra, scoring_method=self.nh.nplinker.scoring.metcalf)
        self.gcf_links = self.nh.nplinker.get_links(self.nh.nplinker.gcfs, scoring_method=self.nh.nplinker.scoring.metcalf)

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

        self.bgc_cols = [TableColumn(field='bgc_pk', title='ID'), TableColumn(field='name', title='BGC name'), TableColumn(field='product_type', title='Product type')]
        self.gcf_cols = [TableColumn(field='gcf_pk', title='ID'), TableColumn(field='gcf_id', title='GCF ID'), TableColumn(field='product_type', title='Product type')]
        self.spec_cols = [TableColumn(field='spec_pk', title='ID'), TableColumn(field='spectrum_id', title='Spectrum ID'), TableColumn(field='family', title='Family ID')]
        self.molfam_cols = [TableColumn(field='molfam_pk', title='ID'), TableColumn(field='family_id', title='MolFam ID'), TableColumn(field='nspectra', title='#spectra')]

        self.bgc_ds = ColumnDataSource(self.bgc_df)
        self.gcf_ds = ColumnDataSource(self.gcf_df)
        self.spec_ds = ColumnDataSource(self.spec_df)
        self.molfam_ds = ColumnDataSource(self.molfam_df)

        self.bgc_dt = DataTable(source=self.bgc_ds, columns=self.bgc_cols, width=300, name='table_bgcs', visible=False)
        self.gcf_dt = DataTable(source=self.gcf_ds, columns=self.gcf_cols, width=300, name='table_gcfs', visible=False)
        self.spec_dt = DataTable(source=self.spec_ds, columns=self.spec_cols, width=300, name='table_spectra', visible=False)
        self.molfam_dt = DataTable(source=self.molfam_ds, columns=self.molfam_cols, width=300, name='table_molfams', visible=False)

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

        # custom JS callback codes to call the linker when a data table is clicked
        code = """    
            // get linker object    
            const linker = window.shared['linker'];
            if (!linker) {
                alert('Please click Load Data first!');
                cb_obj.indices = [];
                return;
            }
            
            // set selections       
            const selected_indices = cb_obj.indices;
            console.log(table_name, selected_indices);
            linker.removeConstraints(table_name);
            for(let i = 0; i < selected_indices.length; i++) {
                const idx = selected_indices[i];
                console.log("%o %o %o", table_name, idx, data_sources[table_name]);
                linker.addConstraint(table_name, idx, data_sources[table_name]);        
            }
            
            // query the linker and update data sources with the query result
            linker.query();
            linker.updateDataSources(data_sources);    
        """
        self.molfam_ds.selected.js_on_change('indices', CustomJS(args={'table_name': 'molfam_table', 'data_sources': data_sources}, code=code))
        self.spec_ds.selected.js_on_change('indices', CustomJS(args={'table_name': 'spec_table', 'data_sources': data_sources}, code=code))
        self.bgc_ds.selected.js_on_change('indices', CustomJS(args={'table_name': 'bgc_table', 'data_sources': data_sources}, code=code))
        self.gcf_ds.selected.js_on_change('indices', CustomJS(args={'table_name': 'gcf_table', 'data_sources': data_sources}, code=code))

        # pickle the links as they don't change for a given dataset and can take
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

        self.tables_score_met = Button(label='Show scores for selected spectra', name='tables_score_met', visible=False)
        self.tables_score_gen = Button(label='Show scores for selected GCFs', name='tables_score_gen', visible=False)
        self.widgets.append(self.tables_score_met)
        self.widgets.append(self.tables_score_gen)

        # combine the data and link informations into a table info
        # this will be passed to the linker object when the Load button is clicked
        self.table_info = get_table_info(self.molfam_df, self.mf_spectra, self.spec_df, self.spectra_bgc, self.bgc_df, self.bgc_gcf, self.gcf_df)
        table_widgets = [self.bgc_dt, self.gcf_dt, self.spec_dt, self.molfam_dt]
        button_widgets = [self.tables_score_met, self.tables_score_gen]
        self.load_tables = Button(label='Load table data', name='load_tables', sizing_mode='stretch_both')
        self.load_tables.callback = CustomJS(args=dict(tableInfo=self.table_info, tableWidgets=table_widgets, buttonWidgets=button_widgets), code="""
            // disable and hide the button
            cb_obj.disabled = true;
            cb_obj.visible = false;

            // create linker object and store in window.shared 
            const linker = new Linker(tableInfo);
            window.shared['linker'] = linker;

            // make the tables and buttons visible
            var i = 0;
            for(; i<tableWidgets.length;i++) {
                var table = tableWidgets[i];
                table.visible = true;
            }
            i = 0;
            for(; i<buttonWidgets.length;i++) {
                var button = buttonWidgets[i];
                button.visible = true;
            }
        """)
        self.widgets.append(self.load_tables)
