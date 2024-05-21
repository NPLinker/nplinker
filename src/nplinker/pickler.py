# Copyright 2021 The NPLinker Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import logging
import os
import pickle
from .genomics import BGC
from .genomics import GCF
from .metabolomics import MolecularFamily
from .metabolomics import Spectrum


logger = logging.getLogger(__name__)

# This is based on the example at https://docs.python.org/3/library/pickle.html#pickle-persistent
# Normally there are serious problems trying to pickle any object with a __hash__,
# so this is an attempt to workaround it and allow pickling of NPLinker API objects.
#
# The NPLinkerPickler class returns a "persistent ID" for selected types of object.
# This ID replaces the object itself during the pickle process. Other objects are
# treated normally. The persistent ID is just the object class name and the
# internal NPLinker ID in most cases.
#
# During unpickling, the NPLinkerUnpickler class performs the inverse operation,
# taking a persistent ID and replacing it with the object it references via an
# instance of the NPLinker class.


class NPLinkerPickler(pickle.Pickler):
    def persistent_id(self, obj):
        if isinstance(obj, BGC):
            return ("BGC", obj.bgc_id)
        elif isinstance(obj, GCF):
            return ("GCF", obj.gcf_id)
        elif isinstance(obj, Spectrum):
            return ("Spectrum", obj.spectrum_id)
        elif isinstance(obj, MolecularFamily):
            return ("MolecularFamily", obj.family_id)
        else:
            # TODO: ideally should use isinstance(obj, ScoringMethod) here
            # but it's currently a problem because it creates a circular
            # import situation
            name = type(obj).__name__
            if name == "RosettaScoring" or name == "MetcalfScoring":
                return ("ScoringMethod", obj.name)

        # pickle anything else as usual
        return None


class NPLinkerUnpickler(pickle.Unpickler):
    def __init__(self, file, nplinker, protocol=4):
        super().__init__(file)
        self.nplinker = nplinker

    def persistent_load(self, pid):
        # return the corresponding NPLinker API objects based on the
        # persistent ID that's been unpickled
        obj_type, obj_id = pid

        if obj_type == "BGC":
            return self.nplinker.bgcs[obj_id]
        elif obj_type == "GCF":
            return self.nplinker.gcfs[obj_id]
        elif obj_type == "Spectrum":
            return self.nplinker.spectra[obj_id]
        elif obj_type == "MolecularFamily":
            return self.nplinker.molfams[obj_id]
        elif obj_type == "ScoringMethod":
            return self.nplinker.scoring_method(obj_id)

        raise pickle.UnpicklingError(f"Unsupported persistent object: {pid}")


# basic wrapper for loading pickled data files. this is done in various places
# in the webapp especially, and in most of those instances if the loading fails
# (e.g. because of a truncated file, or a change in data format) then the simple
# fix is to delete the file and recreate it. this just avoids having the same
# snippet of code in many different places
def load_pickled_data(nplinker, path, delete_on_error=True):
    if not os.path.exists(path):
        return None

    unp = NPLinkerUnpickler(open(path, "rb"), nplinker)
    data = None
    try:
        data = unp.load()
    except Exception as e:
        logger.warning(
            'Failed to unpickle file "{}", deleting it (exception={})'.format(path, str(e))
        )
        try:
            os.unlink(path)
        except OSError:
            pass

    return data


def save_pickled_data(data, path):
    dirpath = os.path.dirname(path)
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)

    pic = NPLinkerPickler(open(path, "wb"))
    pic.dump(data)
