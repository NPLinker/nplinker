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

import numpy
from numba import jit


class MSSpectrum:
    def __init__(self, mgf_dict=None, spec=None):
        self.compound = None
        self.formula = None
        self.ionisation = None
        self.raw_parentmass = None
        self.filename = None
        self.id = None

        self.raw_spectrum = None
        self.output_spectrum = None

        # self.filter is a function that can be used
        # to filter the raw spectrum, e.g. denoising,
        # removal of adduct, etc.
        self.filter = None
        self.normalise = False

        self.correct_for_ionisation = False

        if mgf_dict is not None:
            self.init_from_mgf(mgf_dict)

        if spec is not None:
            self.init_from_spec(spec)

    def init_from_spec(self, spec):
        self.id = spec.id
        self.raw_parentmass = spec.precursor_mz
        self.raw_spectrum = spec.peaks
        # TODO this is a temporary default for the Crusemann data
        # should check for it in the mgf in metabolomics.py and store
        # in the Spectrum object if found
        # (TODO what is the MGF field name for this??)
        self.ionisation = "[M+H]+"

    def init_from_mgf(self, mgf_dict):
        self.output_spectrum = None

        self.filename = mgf_dict["params"]["filename"]
        self.compound = None
        self.formula = None
        self.ionisation = mgf_dict["params"]["charge"]
        self.raw_parentmass = mgf_dict["params"]["pepmass"][0]
        self.inchi = None
        self.id = None
        if "smiles" in mgf_dict["params"]:
            self.smiles = mgf_dict["params"]["smiles"]

        spec = []
        for a in zip(mgf_dict["m/z array"], mgf_dict["intensity array"]):
            spec.append(a)

        self.raw_spectrum = numpy.array(spec)

    def load(self, filename):
        self.output_spectrum = None
        self.filename = filename
        spectrum = []
        with open(filename) as f:
            for line in f.readlines():
                line = line.strip()
                if len(line) == 0:
                    pass
                elif line.startswith(">compound"):
                    self.compound = strip_leading(line)
                elif line.startswith(">formula"):
                    self.formula = strip_leading(line)
                elif line.startswith(">ionization"):
                    self.ionisation = strip_leading(line)
                elif line.startswith(">parentmass"):
                    self.raw_parentmass = float(strip_leading(line))
                elif line.startswith(">"):
                    pass
                elif line.startswith("#inchi"):
                    self.inchi = strip_leading(line)
                elif line.startswith("#SpectrumID"):
                    self.id = strip_leading(line)
                elif line.startswith("#"):
                    pass
                else:
                    mass, charge = line.split()
                    mass = float(mass)
                    charge = float(charge)
                    spectrum.append((mass, charge))
        self.raw_spectrum = numpy.array(spectrum)

    @property
    def parentmass(self):
        if self.correct_for_ionisation:
            return self.raw_parentmass - self.ionisation_mass
        else:
            return self.raw_parentmass

    @property
    def spectrum(self):
        if self.normalise:
            return _normalise_spectrum(self.unnormalised_spectrum)
        else:
            return self.unnormalised_spectrum

    @property
    def unnormalised_spectrum(self):
        if self.filter is None:
            if self.correct_for_ionisation:
                return self.shifted_spectrum
            else:
                return self.raw_spectrum
        else:
            if self.output_spectrum is None:
                filtered_spectrum = self.filter(self)
                if len(filtered_spectrum) != 0:
                    self.output_spectrum = self.filter(self)
                else:
                    self.output_spectrum = self.shifted_spectrum
            return self.output_spectrum

    @property
    def shifted_spectrum(self):
        return self.raw_spectrum - [self.ionisation_mass, 0]

    @property
    def ionisation_mass(self):
        return IONISATION_MASSES[self.ionisation]


PROTON_MASS = 1.00727645199076
IONISATION_MASSES = {
    "[M+H]+": PROTON_MASS,
    "[M+H-H2O]+": PROTON_MASS - 18.01056468638,
    "[M+K]+": 38.963158,
    "[M+Na]+": 22.989218,
}


def _normalise_spectrum(spectrum, peak=100.0):
    _, max_peak = numpy.max(spectrum, axis=0)
    return spectrum * [1, peak / max_peak]


def strip_leading(line):
    return " ".join(line.split()[1:])


def _ppk(i_peaks, j_peaks, sm, si):
    X1 = i_peaks
    X2 = j_peaks
    # N1 = numpy.size(X1, 0); N2 = numpy.size(X2, 0)
    N1 = X1.shape[0]
    N2 = X2.shape[0]
    if N1 == 0 or N2 == 0:
        raise Exception("[ERROR]:No peaks when computing the kernel.(try not clean the peaks)")
    constant = 1.0 / (N1 * N2) * 0.25 / (numpy.pi * numpy.sqrt(sm * si))
    mass_term = (
        1.0
        / sm
        * numpy.power(
            numpy.kron(X1[:, 0].flatten(), numpy.ones(N2))
            - numpy.kron(numpy.ones(N1), X2[:, 0].flatten()),
            2,
        )
    )
    inte_term = (
        1.0
        / si
        * numpy.power(
            numpy.kron(X1[:, 1].flatten(), numpy.ones(N2))
            - numpy.kron(numpy.ones(N1), X2[:, 1].flatten()),
            2,
        )
    )
    return constant * numpy.sum(numpy.exp(-0.25 * (mass_term + inte_term)))


def ppk(*args):
    # t0 = time.time()
    # a = ppk_loop(*args)
    # = _ppk(*args)
    # t1 = time.time()
    b = ppk_limit(*args)
    # t2 = time.time()
    # print(t1-t0, t2-t1)
    return b


@jit(nopython=True)
def ppk_loop(spectrum_1, spectrum_2, sigma_mass, sigma_int):
    # the inputs are really sigma^2, though
    # sigma_mass = 0.00001
    # sigma_int = 100000
    sigma_array = numpy.array([[sigma_mass, 0], [0, sigma_int]])
    sigma_inv = numpy.linalg.inv(sigma_array)
    len_1 = spectrum_1.shape[0]
    len_2 = spectrum_2.shape[0]
    constant_term = 1.0 / (len_1 * len_2 * 4 * numpy.pi * numpy.sqrt(sigma_mass * sigma_int))
    sum_term = 0
    # for p_1, p_2 in itertools.product(spectrum_1, spectrum_2):
    for p_1_idx in range(len_1):
        p_1 = spectrum_1[p_1_idx, :]
        for p_2_idx in range(len_2):
            p_2 = spectrum_2[p_2_idx, :]
            d = p_1 - p_2
            sum_term += numpy.exp(-0.25 * numpy.sum(d * sigma_inv * d))
    # print(sum_term)
    # print(numpy.sum(sum_term))
    return constant_term * sum_term


@jit(nopython=True)
def ppk_limit(spectrum_1, spectrum_2, sigma_mass, sigma_int):
    # the inputs are really sigma^2, though
    # sigma_mass = 0.00001
    # sigma_int = 100000
    sigma_array = numpy.array([[sigma_mass, 0], [0, sigma_int]])
    sigma_inv = numpy.linalg.inv(sigma_array)
    len_1 = spectrum_1.shape[0]
    len_2 = spectrum_2.shape[0]
    constant_term = 1.0 / (len_1 * len_2 * 4 * numpy.pi * numpy.sqrt(sigma_mass * sigma_int))
    sum_term = 0

    tol = 5 * numpy.sqrt(sigma_mass)

    for p_1_idx, p_2_idx in find_pairs(spectrum_1, spectrum_2, tol):
        p_1 = spectrum_1[p_1_idx, :]
        p_2 = spectrum_2[p_2_idx, :]
        d = p_1 - p_2
        sum_term += numpy.exp(-0.25 * numpy.sum(d * sigma_inv * d))
    # print(sum_term)
    # print(numpy.sum(sum_term))
    # print(constant_term, sum_term)
    return constant_term * sum_term


@jit(nopython=True)
def find_pairs(spec1, spec2, tol, shift=0):
    matching_pairs = []
    spec2_lowpos = 0
    spec2_length = len(spec2)

    for idx in range(len(spec1)):
        mz, intensity = spec1[idx, :]
        while spec2_lowpos < spec2_length and spec2[spec2_lowpos][0] + shift < mz - tol:
            spec2_lowpos += 1
        if spec2_lowpos == spec2_length:
            break
        spec2_pos = spec2_lowpos
        while spec2_pos < spec2_length and spec2[spec2_pos][0] + shift < mz + tol:
            matching_pairs.append((idx, spec2_pos))
            spec2_pos += 1

    return matching_pairs


def ppk_nloss(spec1, spec2, prec1, prec2, sigma_mass, sigma_int):
    spec1_loss = ([prec1, 0] - spec1) * [1, -1]
    spec2_loss = ([prec2, 0] - spec2) * [1, -1]
    k_nloss = ppk(spec1_loss[::-1], spec2_loss[::-1], sigma_mass, sigma_int)
    return k_nloss
