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

import pickle
import numpy
import scipy.io
from cdk_pywrapper.cdk_pywrapper import Compound


def main(datapath):
    data = scipy.io.loadmat(datapath + "/data_GNPS.mat")

    fp_size = 1024

    new_fps = []

    for i in data["inchi"]:
        print("Processing sample %s" % i)
        inchi = i[0][0]
        print(inchi)
        c = Compound(compound_string=inchi, identifier_type="inchi")
        fp = c.get_fingerprint()
        fp_array = numpy.zeros(fp_size)
        for fp_bit in range(fp_size):
            fp_array[fp_bit] = fp.get(fp_bit)
        new_fps.append(fp_array)

    with open(datapath + "/cdk_fingerprints.bin", "wb") as f:
        pickle.dump(new_fps, f)


def default_fingerprint_from_inchi(inchi):
    c = Compound(compound_string=inchi, identifier_type="inchi")
    # this should accommodate different types!
    fp = c.get_fingerprint()
    fp_size = 1024
    fp_array = numpy.zeros(fp_size)
    for fp_bit in range(fp_size):
        fp_array[fp_bit] = fp.get(fp_bit)
    return fp_array


def fingerprint_from_smiles(smiles, fingerprint_type=None):
    if fingerprint_type is None:
        fingerprint = numpy.array([])
        for fp_type in ("cdk_default", "substructure", "klekota-roth"):
            fingerprint = numpy.hstack((fingerprint, fingerprint_from_smiles(smiles, fp_type)))
        return fingerprint

    c = Compound(compound_string=smiles, identifier_type="smiles")

    if fingerprint_type == "cdk_default":
        fingerprinter = c.cdk.fingerprint.Fingerprinter()
    elif fingerprint_type == "substructure":
        fingerprinter = c.cdk.fingerprint.SubstructureFingerprinter()
    elif fingerprint_type == "klekota-roth":
        fingerprinter = c.cdk.fingerprint.KlekotaRothFingerprinter()
    else:
        raise SystemExit(f"Unknown fingerprint type: {fingerprint_type}")

    fp = fingerprinter.getBitFingerprint(c.mol_container)
    fp_size = fp.size()
    fp_array = numpy.zeros(fp_size)
    for fp_bit in range(fp_size):
        fp_array[fp_bit] = fp.get(fp_bit)
    return fp_array


def fingerprint_from_inchi(inchi, fingerprint_type=None):
    if fingerprint_type is None:
        fingerprint = numpy.array([])
        for fp_type in ("cdk_default", "substructure", "klekota-roth"):
            fingerprint = numpy.hstack((fingerprint, fingerprint_from_inchi(inchi, fp_type)))
        return fingerprint

    c = Compound(compound_string=inchi, identifier_type="inchi")

    if fingerprint_type == "cdk_default":
        fingerprinter = c.cdk.fingerprint.Fingerprinter()
    elif fingerprint_type == "substructure":
        fingerprinter = c.cdk.fingerprint.SubstructureFingerprinter()
    elif fingerprint_type == "klekota-roth":
        fingerprinter = c.cdk.fingerprint.KlekotaRothFingerprinter()
    else:
        raise SystemExit(f"Unknown fingerprint type: {fingerprint_type}")

    fp = fingerprinter.getBitFingerprint(c.mol_container)
    fp_size = fp.size()
    fp_array = numpy.zeros(fp_size)
    for fp_bit in range(fp_size):
        fp_array[fp_bit] = fp.get(fp_bit)
    return fp_array
