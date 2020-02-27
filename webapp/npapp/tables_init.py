import os
import pickle

import pandas as pd

from bokeh.models import Button, CustomJS, ColumnDataSource
from bokeh.models.widgets import DataTable, TableColumn, Div

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

        self.data_tables = {
            'molfam_table': self.molfam_dt, 
            'spec_table': self.spec_dt,
            'bgc_table': self.bgc_dt, 
            'gcf_table': self.gcf_dt
        }

        # custom JS callback to call the linker when a data table is clicked
        # (note only triggered for the table that's actually clicked, updates
        # to the others don't trigger this)
        code = """
            // get linker object
            const linker = window.shared['linker'];
            if (!linker) {
                alert('Please click Load Data first!');
                cb_obj.indices = [];
                return;
            }

            console.log("js_on_change for %o", table_name);

            // this is a hacky way to show that the tables on the opposite side
            // from the clicked table are now "disabled", by setting the opacity
            // to 50% via CSS (but only if something was actually selected)
            if(cb_obj.indices.length > 0) {
                if(table_name == "molfam_table" || table_name == "spec_table") {
                    var tables = $(".gen-table");
                    for(var i=0;i<tables.length;i++) {
                        $('#' + tables[i].id).css('opacity', '50%');
                    }
                } else {
                    var tables = $(".met-table");
                    for(var i=0;i<tables.length;i++) {
                        $('#' + tables[i].id).css('opacity', '50%');
                    }
                }
            }
            
            // propagate selections to the linker
            const selected_indices = cb_obj.indices;
            linker.removeConstraints(table_name);
            for(let i = 0; i < selected_indices.length; i++) {
                const idx = selected_indices[i];
                linker.addConstraint(table_name, idx, data_sources[table_name]);
            }

            // query the linker and update data sources with the query result
            linker.query();
            linker.updateDataSources(data_sources);
        """

        # make the genomics tables non-selectable when a metabolomics table is clicked
        def met_tables_clicked(attr, old, new):
            self.gcf_dt.selectable = False
            self.bgc_dt.selectable = False

        # make the metabolomics tables non-selectable when a genomics table is clicked
        def gen_tables_clicked(attr, old, new):
            self.molfam_dt.selectable = False
            self.spec_dt.selectable = False

        # NOTE: passing in DataTable widgets here in the same way as data_sources causes bokeh JS to throw weird exceptions... 
        self.molfam_ds.selected.js_on_change('indices', CustomJS(args={'table_name': 'molfam_table', 'data_sources': data_sources}, code=code))
        self.spec_ds.selected.js_on_change('indices', CustomJS(args={'table_name': 'spec_table', 'data_sources': data_sources}, code=code))
        self.bgc_ds.selected.js_on_change('indices', CustomJS(args={'table_name': 'bgc_table', 'data_sources': data_sources}, code=code))
        self.gcf_ds.selected.js_on_change('indices', CustomJS(args={'table_name': 'gcf_table', 'data_sources': data_sources}, code=code))

        self.molfam_ds.selected.on_change('indices', met_tables_clicked)
        self.spec_ds.selected.on_change('indices', met_tables_clicked)
        self.bgc_ds.selected.on_change('indices', gen_tables_clicked)
        self.gcf_ds.selected.on_change('indices', gen_tables_clicked)

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
        self.tables_score_met = Button(label='Show scores for selected spectra', name='tables_score_met')
        self.tables_score_gen = Button(label='Show scores for selected BGCs', name='tables_score_gen')
        self.widgets.append(self.tables_score_met)
        self.widgets.append(self.tables_score_gen)

        # CustomJS callback code for resetting the selection (and opacity) states
        # of all the tables 
        reset_selection_code = """    
            // get linker object    
            const linker = window.shared['linker'];
            if (!linker) {
                alert('Please click Load Data first!');
                cb_obj.indices = [];
                return;
            }

            
            // remove selections from all tables
            var tablenames = ['molfam_table', 'spec_table', 'bgc_table', 'gcf_table'];
            for(let t=0;t<tablenames.length;t++) {
                const table = tablenames[t];
                // TODO is this the best way to do this??
                //linker.removeConstraints(table);
                data_sources[table].selected.indices = [];
                
                // restore opacity
                $('#' + table).css('opacity', '');
            }
            linker.query();
            linker.updateDataSources(data_sources);

        """
        self.tables_reset = Button(label='Clear selections', name='tables_reset')
        self.tables_reset.js_on_click(CustomJS(args={'data_sources': data_sources, 'data_tables': self.data_tables}, code=reset_selection_code))
        self.widgets.append(self.tables_reset)

        # combine the data and link informations into a list of dicts in the format the
        # linker expects. this object is passed in to the linker by a callback after 
        # the webapp loads (see below & main.py)
        self.table_info = get_table_info(self.molfam_df, self.mf_spectra, self.spec_df, self.spectra_bgc, self.bgc_df, self.bgc_gcf, self.gcf_df)

        # hack alert: bokeh doesn't provide any obvious way to trigger a JavaScript callback
        # from Python code other than attaching a CustomJS object to a property of some 
        # object, and then changing that property. This has limitations, such as not 
        # working during the loading process. That's very annoying, because it makes it
        # difficult to load the data into the tables (i.e. transferring data from 
        # Python -> JS) without the original solution of asking the user to click a button
        # to trigger the CustomJS callback necessary. 
        #
        # I finally found a very ugly but usable workaround thanks to this issue:
        #    https://github.com/bokeh/bokeh/issues/8728
        # which can be summed up as:
        # 1. create an empty/hidden div on the page
        # 2. attach the callback to be executed to its 'text' property
        # 3. add a callback on the document object 
        # 4. use *that* callback to modify the div text...
        # 5. ... which then triggers the CustomJS callback... 
        # 6. ... and finally remove the document callback 
        # 
        # Not very nice but it does work

        # this is the empty div
        self.dummydiv = Div(text='', name='dummydiv')
        self.widgets.append(self.dummydiv)

        # this callback is used to initialise the code behind the data tables, which 
        # is all Javascript. It must be passed the set of tables constructed above,
        # so that they can be used to populate the database that is ultimately used
        # to do the filtering. The check for the loading already having been completed
        # is required due to the really nasty workaround used to compensate for bokeh
        # not having a nice way to trigger a CustomJS callback programatically from
        # Python code (at least not on/as the page is loading).  
        self.load_tables_callback = CustomJS(args=dict(dummyDiv=self.dummydiv, tableInfo=self.table_info), code="""
            if (window.shared.hasOwnProperty('linker'))
            {
                console.log("Linker init already done!");
                return;
            }

            // create linker object and store in window.shared
            const linker = new Linker(tableInfo);
            window.shared['linker'] = linker;

            // update div to indicate loading completed (this is used to
            // cancel the periodic callbacks begun by the code in main.py)
            // yes it's a mess
            dummyDiv.text = 'LOADED';

            // hide the CSS 'loading...' overlay
            document.getElementById("overlay").style.display = "none";
        """)

        # finally attach the callback to the text property of the empty div
        # so that it will be triggered when the main.py callback updates 
        # that property
        self.dummydiv.js_on_change('text', self.load_tables_callback)
