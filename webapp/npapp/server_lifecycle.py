import sys
import os

# add path to nplinker/prototype (for local dev)
sys.path.append(os.path.join(os.path.dirname(__file__), '../../prototype'))

from nplinker.nplinker import NPLinker
from nplinker.logconfig import LogConfig

from tables_init import TableData, TableSessionData

class NPLinkerHelper(object):
    """
    This is just a simple class to wrap up the various objects involved in 
    loading and plotting a dataset via nplinker/bokeh
    """
    def __init__(self):
        LogConfig.setLogLevelStr('DEBUG')
        self.data_dir = os.path.join(os.path.dirname(__file__), 'data')

    def load(self):
        self.load_nplinker()

        self.load_tables()

        self.load_rosetta()

    def load_rosetta(self):
        # trigger rosetta scoring setup during initial load as the first run 
        # can take a significant amount of time. the results are cached and
        # pickled for subsequent use.
        print('Preprocessing for rosetta scoring...')
        self.nplinker.scoring_method('rosetta')
        print('Finished preprocessing for rosetta scoring')

    def load_nplinker(self):
        # initialise nplinker and load the dataset, using a config file in the webapp dir
        datapath = os.getenv('NPLINKER_CONFIG', os.path.join(os.path.join(os.path.dirname(__file__), 'nplinker_webapp.toml')))
        print('DATAPATH: {}'.format(datapath))
        self.nplinker = NPLinker(datapath)
        if not self.nplinker.load_data():
            raise Exception('Failed to load data')

    def load_tables(self):
        print('Loading tables data')
        self.table_data = TableData(self)
        self.table_data.setup()

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
