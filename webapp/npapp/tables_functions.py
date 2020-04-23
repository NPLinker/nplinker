import numpy as np
import pandas as pd
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import DataTable, TableColumn
from natsort import natsorted, index_natsorted, order_by_index

NA = '-'


def get_data(label, max_items):
    # generate some random data
    items = ['%s_%d' % (label, i+1) for i in list(range(max_items))]
    x1 = np.random.randint(low=1, high=100, size=max_items).tolist()
    x2 = np.random.randint(low=1, high=100, size=max_items).tolist()
    x3 = np.random.randint(low=1, high=100, size=max_items).tolist()

    # insert NAs to the first row
    items.insert(0, NA)
    x1.insert(0, NA)
    x2.insert(0, NA)
    x3.insert(0, NA)

    # create pandas dataframe
    colname = '%s_pk' % label
    data = {
        colname: items,
        'x1': x1,
        'x2': x2,
        'x3': x3
    }
    df = pd.DataFrame(data)

    # https://stackoverflow.com/questions/29580978/naturally-sorting-pandas-dataframe
    df = df.reindex(index=order_by_index(df.index, index_natsorted(df[colname])))

    # create bokeh ColumnDataSource and DataTable
    columns = [
        TableColumn(field=colname, title='ID'),
        TableColumn(field='x1', title='x1'),
        TableColumn(field='x2', title='x2'),
        TableColumn(field='x3', title='x3'),
    ]
    ds = ColumnDataSource(df)
    dt = DataTable(source=ds, columns=columns, width=300, height=300)
    return df, ds, dt


def create_links(df1, df2, key1, key2, indices):
    # link entries in df1 and df2 to each other according to their indices
    # for idx1, idx2 in indices:
    #     val1 = df1.loc[idx1][key1]
    #     val2 = df2.loc[idx2][key2]
    #     link = {key1: val1, key2: val2}
    #     links.append(link)
    # this is *much* faster than the above and should do the same thing...
    links = [{key1: id1, key2: id2} for id1, id2 in indices]

    # make sure that there are no unlinked items
    # for items where no linking has been specified, we just link them to NA
    linked_1 = set([x[0] for x in indices])
    linked_2 = set([x[1] for x in indices])
    unlinked_1 = set(df1.index.tolist()) - linked_1
    unlinked_2 = set(df2.index.tolist()) - linked_2

    # add the linking for unlinked items in df1 to NA
    for idx in unlinked_1:
        val1 = df1.loc[idx][key1]
        val2 = NA
        if val1 != NA:
            link = {key1: val1, key2: val2}
            links.append(link)

    # add the linking for unlinked items in df2 to NA
    for idx in unlinked_2:
        val1 = NA
        val2 = df2.loc[idx][key2]
        if val2 != NA:
            link = {key1: val1, key2: val2}
            links.append(link)

    # finally link NA to NA
    link = {key1: NA, key2: NA}
    links.append(link)

    return links


def get_table_info(mf_df, mf_spectra, spectra_df, spectra_bgc, bgc_df, bgc_gcf, gcf_df):
    tables_info = [
        {
            'tableName': 'molfam_table',
            'tableData': mf_df.to_dict('records'),
            'options': {
                'visible': True,
                'pk': 'molfam_pk'
            },
            'relationship': {'with': 'mf_spectra', 'using': 'molfam_pk'}
        },
        {
            'tableName': 'mf_spectra',
            'tableData': mf_spectra,
            'options': {
                'visible': False
            },
            'relationship': {'with': 'spec_table', 'using': 'spec_pk'}
        },
        {
            'tableName': 'spec_table',
            'tableData': spectra_df.to_dict('records'),
            'options': {
                'visible': True,
                'pk': 'spec_pk'
            },
            'relationship': {'with': 'spectra_bgc', 'using': 'spec_pk'}
        },
        {
            'tableName': 'spectra_bgc',
            'tableData': spectra_bgc,
            'options': {
                'visible': False
            },
            'relationship': {'with': 'bgc_table', 'using': 'bgc_pk'}
        },
        {
            'tableName': 'bgc_table',
            'tableData': bgc_df.to_dict('records'),
            'options': {
                'visible': True,
                'pk': 'bgc_pk'
            },
            'relationship': {'with': 'bgc_gcf', 'using': 'bgc_pk'}
        },
        {
            'tableName': 'bgc_gcf',
            'tableData': bgc_gcf,
            'options': {
                'visible': False
            },
            'relationship': {'with': 'gcf_table', 'using': 'gcf_pk'}
        },
        {
            'tableName': 'gcf_table',
            'tableData': gcf_df.to_dict('records'),
            'options': {
                'visible': True,
                'pk': 'gcf_pk'
            }
        },

    ]
    return tables_info


def update_alert(alert_div, msg, alert_class='primary'):
    alert_div.text = '<div class="alert alert-{}" role="alert">{}</div>'.format(alert_class, msg)
