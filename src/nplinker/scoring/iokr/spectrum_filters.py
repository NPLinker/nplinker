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

import json
import os
import pickle
import numpy

# import sys
# sys.path.append('/home/grimur/git/lda')
# from lda.code.formula import Formula
from .formula import Formula


global _ALL_DAG_MASSES
global _ALL_TREE_MASSES
_ALL_DAG_MASSES = None
_ALL_TREE_MASSES = None

global datapath
datapath = ""


def load_formula_dag(filename):
    formula = set()

    with open(filename) as f:
        for line in f.readlines():
            if not line.startswith("v"):
                continue

            fields = line.strip().split()
            if fields[1] == "->":
                # network edge
                f = fields[3].split('="')[1][:-3]
            else:
                # network vertex
                f = line.split("\\n")[0].split('"')[-1]

            formula.add(f)

    return list(formula)


def load_formula_tree(filename):
    with open(filename, "rb") as f:
        data = json.load(f)

    formula = [data["molecularFormula"], data["root"]]
    formula.extend([x["molecularFormula"] for x in data["fragments"]])
    formula.extend([x["source"] for x in data["losses"]])
    formula.extend([x["target"] for x in data["losses"]])
    formula.extend([x["molecularFormula"] for x in data["losses"]])

    return list(set(formula))


def load_peaks_from_tree(filename):
    with open(filename, "rb") as f:
        data = json.load(f)

    spectrum = []
    for f in data["fragments"]:
        if len(f["peaks"]) > 0:
            for p in f["peaks"]:
                mass = p["mz"]
                intensity = p["int"]
                spectrum.append((mass, intensity))
        else:
            mass = f["mz"]
            intensity = f["intensity"]
            spectrum.append((mass, intensity))

    return numpy.array(spectrum)


# These functions should be hanged onto the MS object
def filter_by_tree(self):
    treepath = "/home/grimur/iokr/data/trees"
    tol = 0.005

    try:
        formula = load_formula_tree(treepath + os.path.sep + self.id + ".json")
    except:
        formula = load_formula_dag(treepath + os.path.sep + self.id + ".dot")
    formula_objects = [Formula(x) for x in formula]
    formula_masses = sorted(x.compute_exact_mass() for x in formula_objects)

    return filter_by_mass(self.shifted_spectrum, formula_masses, tol)


def filter_by_tree_unshifted(self):
    treepath = "/home/grimur/iokr/data/trees"

    try:
        spectrum = load_peaks_from_tree(treepath + os.path.sep + self.id + ".json")
    except FileNotFoundError:
        spectrum = self.raw_spectrum

    return spectrum


def filter_by_dag(self):
    treepath = "/home/grimur/iokr/data/trees"
    tol = 0.005

    formula = load_formula_dag(treepath + os.sep + self.id + ".dot")
    formula_objects = [Formula(x) for x in formula]
    formula_masses = sorted(x.compute_exact_mass() for x in formula_objects)

    return filter_by_mass(self.shifted_spectrum, formula_masses, tol)


def load_all_dag_masses(path):
    formula_masses_collected = []
    for filename in os.listdir(path):
        if filename.endswith(".dot"):
            formula = load_formula_dag(path + os.sep + filename)
            formula_objects = [Formula(x) for x in formula]
            formula_masses_collected.extend([x.compute_exact_mass() for x in formula_objects])

    return sorted(list(set(formula_masses_collected)))


def filter_by_collected_dag(self):
    treepath = "/home/grimur/iokr/data/trees"
    tol = 0.005
    global _ALL_DAG_MASSES

    if _ALL_DAG_MASSES is None:
        _ALL_DAG_MASSES = load_all_dag_masses(treepath)

    return filter_by_mass(self.shifted_spectrum, _ALL_DAG_MASSES, tol)


def filter_by_frozen_dag(self):
    with open(os.path.join(datapath, "dag_masses.bin"), "rb") as f:
        dag_masses = pickle.load(f)
    tol = 0.005
    return filter_by_mass(self.shifted_spectrum, dag_masses, tol)


def load_all_tree_masses(path):
    formula_masses_collected = []
    for filename in os.listdir(path):
        if filename.endswith(".json"):
            formula = load_formula_tree(path + os.sep + filename)
            formula_objects = [Formula(x) for x in formula]
            formula_masses_collected.extend([x.compute_exact_mass() for x in formula_objects])

    return sorted(list(set(formula_masses_collected)))


def filter_by_collected_tree(self):
    treepath = "/home/grimur/iokr/data/trees"
    tol = 0.005
    global _ALL_TREE_MASSES

    if _ALL_TREE_MASSES is None:
        _ALL_TREE_MASSES = load_all_tree_masses(treepath)

    return filter_by_mass(self.shifted_spectrum, _ALL_TREE_MASSES, tol)


def filter_by_mass(raw_spectrum, formula_masses, tol):
    filtered_spectrum = []
    spectrum_index = 0
    mass_index = 0
    while mass_index < len(formula_masses):
        while (
            spectrum_index < len(raw_spectrum)
            and raw_spectrum[spectrum_index, 0] < formula_masses[mass_index] + tol
        ):
            while (
                spectrum_index < len(raw_spectrum)
                and raw_spectrum[spectrum_index, 0] < formula_masses[mass_index] - tol
            ):
                spectrum_index += 1
            if (
                spectrum_index < len(raw_spectrum)
                and raw_spectrum[spectrum_index, 0] < formula_masses[mass_index] + tol
            ):
                filtered_spectrum.append(raw_spectrum[spectrum_index, :])
            spectrum_index += 1
        mass_index += 1

    return numpy.array(filtered_spectrum)
