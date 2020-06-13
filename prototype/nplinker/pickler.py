import pickle
import os

from .logconfig import LogConfig

logger = LogConfig.getLogger(__file__)

# basic wrapper for loading pickled data files. this is done in various places 
# in the webapp especially, and in most of those instances if the loading fails
# (e.g. because of a truncated file, or a change in data format) then the simple
# fix is to delete the file and recreate it. this just avoids having the same 
# snippet of code in many different places
def load_pickled_data(path, delete_on_error=True):
    if not os.path.exists(path):
        return None

    data = None
    try:
        data = pickle.load(open(path, 'rb'))
    except Exception as e:
        logger.warning('Failed to unpickle file "{}", deleting it'.format(path))
        try:
            os.unlink(path)
        except OSError as oe:
            pass

    return data

def save_pickled_data(data, path, protocol=4):
    dirpath = os.path.dirname(path)
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)

    pickle.dump(data, open(path, 'wb'), protocol=4)
