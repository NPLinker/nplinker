import pandas as pd

from bokeh.models import Button, CustomJS, ColumnDataSource
from bokeh.models.widgets import DataTable, TableColumn

from tables_functions import create_links, get_table_info, NA

def setup_data_tables(nh):

    
    # 1. need to create a ColumnDataSource and DataTable per object type
    # and populate it with the relevant attributes (which are ??? TODO)
    #   BGCs: strain/name, class/product prediction
    #   GCFs: gcf_id, product_type
    #   Spectra: spectrum_id, family
    #   MolFams: family_id

    bgc_data = dict(bgc_pk=[NA], name=[NA], product_type=[NA], score=[NA])
    gcf_data = dict(gcf_pk=[NA], gcf_id=[NA], product_type=[NA], score=[NA])
    spec_data = dict(spec_pk=[NA], spectrum_id=[NA], family=[NA], score=[NA])
    molfam_data = dict(molfam_pk=[NA], family_id=[NA], score=[NA])

    for bgc in nh.nplinker.bgcs:
        bgc_data['bgc_pk'].append(bgc.id)
        bgc_data['name'].append(bgc.name)
        bgc_data['product_type'].append(bgc.product_prediction)
        bgc_data['score'].append(1)

    for gcf in nh.nplinker.gcfs:
        gcf_data['gcf_pk'].append(gcf.id)
        gcf_data['gcf_id'].append(gcf.gcf_id)
        gcf_data['product_type'].append(gcf.product_type)
        gcf_data['score'].append(1)

    for spec in nh.nplinker.spectra:
        spec_data['spec_pk'].append(spec.id)
        spec_data['spectrum_id'].append(spec.spectrum_id)
        spec_data['family'].append(spec.family.family_id)
        spec_data['score'].append(1)

    for molfam in nh.nplinker.molfams:
        molfam_data['molfam_pk'].append(molfam.id)
        molfam_data['family_id'].append(molfam.family_id)
        molfam_data['score'].append(1)

    bgc_df = pd.DataFrame(bgc_data)
    gcf_df = pd.DataFrame(gcf_data)
    spec_df = pd.DataFrame(spec_data)
    molfam_df = pd.DataFrame(molfam_data)

    bgc_cols = [TableColumn(field='bgc_pk', title='ID'), TableColumn(field='name', title='BGC name'), TableColumn(field='product_type', title='Product type'), TableColumn(field='score', title='Score')]
    gcf_cols = [TableColumn(field='gcf_pk', title='ID'), TableColumn(field='gcf_id', title='GCF ID'), TableColumn(field='product_type', title='Product type'), TableColumn(field='score', title='Score')]
    spec_cols = [TableColumn(field='spec_pk', title='ID'), TableColumn(field='spectrum_id', title='Spectrum ID'), TableColumn(field='family', title='Family ID'), TableColumn(field='score', title='Score')]
    molfam_cols = [TableColumn(field='molfam_pk', title='ID'), TableColumn(field='family_id', title='MolFam ID'), TableColumn(field='score', title='Score')]

    bgc_ds = ColumnDataSource(bgc_df)
    gcf_ds = ColumnDataSource(gcf_df)
    spec_ds = ColumnDataSource(spec_df)
    molfam_ds = ColumnDataSource(molfam_df)

    bgc_dt = DataTable(source=bgc_ds, columns=bgc_cols, width=300, name='table_bgcs')
    gcf_dt = DataTable(source=gcf_ds, columns=gcf_cols, width=300, name='table_gcfs')
    spec_dt = DataTable(source=spec_ds, columns=spec_cols, width=300, name='table_spectra')
    molfam_dt = DataTable(source=molfam_ds, columns=molfam_cols, width=300, name='table_molfams')

    widgets = []
    widgets.append(bgc_dt)
    widgets.append(gcf_dt)
    widgets.append(spec_dt)
    widgets.append(molfam_dt)

    data_sources = {
        'molfam_table': molfam_ds,
        'spec_table': spec_ds,
        'bgc_table': bgc_ds,
        'gcf_table': gcf_ds
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
    molfam_ds.selected.js_on_change('indices', CustomJS(args={'table_name': 'molfam_table', 'data_sources': data_sources}, code=code))
    spec_ds.selected.js_on_change('indices', CustomJS(args={'table_name': 'spec_table', 'data_sources': data_sources}, code=code))
    bgc_ds.selected.js_on_change('indices', CustomJS(args={'table_name': 'bgc_table', 'data_sources': data_sources}, code=code))
    gcf_ds.selected.js_on_change('indices', CustomJS(args={'table_name': 'gcf_table', 'data_sources': data_sources}, code=code))

    spec_links = nh.nplinker.get_links(nh.nplinker.spectra, scoring_method=nh.nplinker.scoring.metcalf)
    gcf_links = nh.nplinker.get_links(nh.nplinker.gcfs, scoring_method=nh.nplinker.scoring.metcalf)

    # create some example links
    # 0 is not used because it's for NA
    link_indices = [(1, 1), (1, 2), (1, 3), (2, 4), (2, 5), (3, 6), (4, 7), (5, 8), (5, 9), (5, 10)]
    mf_spectra = create_links(molfam_df, spec_df, 'molfam_pk', 'spec_pk', link_indices)

    link_indices = [(1, 1), (2, 2), (1, 3), (3, 4), (3, 5), (5, 7)]
    spectra_bgc = create_links(spec_df, bgc_df, 'spec_pk', 'bgc_pk', link_indices)

    # link_indices = [(1, 1), (2, 1), (3, 1), (4, 1), (5, 2), (6, 3), (7, 4), (8, 5), (9, 5), (10, 5)]
    # bgc_gcf = create_links(bgc_df, gcf_df, 'bgc_pk', 'gcf_pk', link_indices)

    # combine the data and link informations into a table info
    # this will be passed to the linker object when the Load button is clicked

    bgc_gcf = []
    print(spectra_bgc)
    bgc_gcf = [ {'bgc_pk': 1, 'gcf_pk': 10} ]

    table_info = get_table_info(molfam_df, mf_spectra, spec_df, spectra_bgc, bgc_df, bgc_gcf, gcf_df)

    load_tables = Button(label='Load table data', name='load_tables')
    load_tables.callback = CustomJS(args=dict(tableInfo=table_info), code="""
        // disable button
        cb_obj.disabled = true;

        // create linker object and store in window.shared 
        const linker = new Linker(tableInfo);
        window.shared['linker'] = linker;
    """)
    widgets.append(load_tables)

    return widgets
