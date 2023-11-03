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

import itertools
import os
import time
import iokrdata as data
import numpy
import scipy
import spectrum_filters
from pyteomics import mgf
from spectrum import MSSpectrum
from spectrum import ppk
from spectrum import ppk_nloss


def ppk_r(spec1, spec2, prec1, prec2, sigma_mass, sigma_int):
    k_peaks = ppk(spec1, spec2, sigma_mass, sigma_int)
    k_nloss = ppk_nloss(spec1, spec2, prec1, prec2, sigma_mass, sigma_int)
    k_diff = ppk_diff(spec1, spec2, sigma_mass, sigma_int)
    return k_peaks + k_nloss + k_diff


def ppk_diff(spec1, spec2, sigma_mass, sigma_int):
    spec1_diff = numpy.array([y - x for x, y in itertools.combinations(spec1, 2)])
    spec2_diff = numpy.array([y - x for x, y in itertools.combinations(spec2, 2)])
    k_diff = ppk(spec1_diff, spec2_diff, sigma_mass, sigma_int)
    return k_diff


def strip_leading(line):
    return " ".join(line.split()[1:])


def create_ppk_matrix():
    iokr_data_path = "/home/grimur/iokr/data"
    data_gnps = scipy.io.loadmat("/home/grimur/iokr/data/data_GNPS.mat")
    ms_path = "/home/grimur/iokr/data/SPEC"

    sigma_mass = 0.00001
    sigma_int = 100000.0

    iokrdata = data.IOKRDataServer(iokr_data_path, kernel="PPKr.txt")
    ker_size = len(iokrdata.spectra)

    kernel_matrix_peaks = numpy.zeros((ker_size, ker_size))
    kernel_matrix_nloss = numpy.zeros((ker_size, ker_size))
    kernel_matrix_ppkr = numpy.zeros((ker_size, ker_size))

    j_ms = MSSpectrum()
    j_ms.filter = spectrum_filters.filter_by_collected_dag
    i_ms = MSSpectrum()
    i_ms.filter = spectrum_filters.filter_by_collected_dag

    import time

    t0 = time.time()
    cnt = 0
    for i in range(len(iokrdata.spectra)):
        i_name = iokrdata.spectra[i][0]
        i_ms.load(ms_path + os.sep + i_name + ".ms")
        # for j in range(i + 1):
        if True:
            j = i
            cnt += 1

            j_name = iokrdata.spectra[j][0]
            j_ms.load(ms_path + os.sep + j_name + ".ms")

            # print('%s vs %s' % (len(i_ms.spectrum), len(j_ms.spectrum)))

            if len(i_ms.spectrum) == 0 or len(j_ms.spectrum) == 0:
                print("empty")
                ij_peaks = 0
                ij_nloss = 0
            else:
                ij_peaks = ppk(i_ms.spectrum, j_ms.spectrum, sigma_mass, sigma_int)
                ij_nloss = ppk_nloss(
                    i_ms.spectrum,
                    j_ms.spectrum,
                    i_ms.parentmass,
                    j_ms.parentmass,
                    sigma_mass,
                    sigma_int,
                )

            kernel_matrix_peaks[i, j] = ij_peaks
            kernel_matrix_peaks[j, i] = ij_peaks

            kernel_matrix_nloss[i, j] = ij_nloss
            kernel_matrix_nloss[j, i] = ij_nloss

            kernel_matrix_ppkr[i, j] = ij_peaks + ij_nloss
            kernel_matrix_ppkr[j, i] = ij_peaks + ij_nloss

            if cnt % 100 == 0:
                print("done {}/{}, {}".format(cnt, (ker_size**2) / 2, time.time() - t0))
                t0 = time.time()

    # numpy.savetxt('ppk_peaks.csv', kernel_matrix_peaks, delimiter=',')
    # numpy.savetxt('ppk_nloss.csv', kernel_matrix_nloss, delimiter=',')
    # numpy.savetxt('ppk_r.csv', kernel_matrix_ppkr, delimiter=',')
    numpy.save("ppk_dag_peaks.npy", kernel_matrix_peaks)
    numpy.save("ppk_dag_nloss.npy", kernel_matrix_nloss)


def create_ppk_matrix_stripe_serial(filter_func, shift, normalise, output_name):
    iokr_data_path = "/home/grimur/iokr/data"
    data_gnps = scipy.io.loadmat("/home/grimur/iokr/data/data_GNPS.mat")
    ms_path = "/home/grimur/iokr/data/SPEC"
    candidate_set = "/home/grimur/iokr/data/mibig/matched_mibig_gnps_2.0.mgf"
    candidate_set_size = 257

    iokrdata = data.IOKRDataServer(iokr_data_path)
    ker_size = len(iokrdata.spectra)

    kernel_matrix_peaks = numpy.zeros((candidate_set_size, ker_size))
    kernel_matrix_nloss = numpy.zeros_like(kernel_matrix_peaks)

    t0 = time.time()
    names = [x[0] for x in iokrdata.spectra]
    cnt = 0
    for i in mgf.read(candidate_set):
        i_ms = MSSpectrum(i)
        # active_jobs.append((i, p.apply_async(do_stripe, (i, names))))
        res = do_stripe(i_ms, names, filter_func, shift, normalise)

        for j_idx, values in enumerate(res):
            ij_peaks, ij_nloss = values

            kernel_matrix_peaks[cnt, j_idx] = ij_peaks
            kernel_matrix_nloss[cnt, j_idx] = ij_nloss

        cnt += 1
        print(f"done {cnt} / {candidate_set_size}, {time.time() - t0}")

    numpy.save(output_name + "_test_peaks.npy", kernel_matrix_peaks)
    numpy.save(output_name + "_test_nloss.npy", kernel_matrix_nloss)


def gather_results(active_jobs):
    while len(active_jobs) > 500:
        done_jobs = []
        remaining_jobs = []
        for i, j, job in active_jobs:
            if job.ready():
                res = job.get()
                done_jobs.append((i, j, res))
            else:
                remaining_jobs.append((i, j, job))
        active_jobs = remaining_jobs
        # time.sleep(1)
    return active_jobs, done_jobs


def gather_results_2(active_jobs, queue_length):
    while len(active_jobs) > queue_length:
        done_jobs = []
        remaining_jobs = []
        for i, job in active_jobs:
            if job.ready():
                res = job.get()
                done_jobs.append((i, res))
            else:
                remaining_jobs.append((i, job))
        active_jobs = remaining_jobs
        time.sleep(0.5)
    return active_jobs, done_jobs


def do_stripe(i_ms, names, filter_func, shift, normalise):
    iokr_data_path = "/home/grimur/iokr/data"
    data_gnps = scipy.io.loadmat("/home/grimur/iokr/data/data_GNPS.mat")
    ms_path = "/home/grimur/iokr/data/SPEC"

    sigma_mass = 0.00001
    sigma_int = 100000.0

    # i_ms = MSSpectrum()
    i_ms.correct_for_ionisation = shift
    i_ms.filter = filter_func
    i_ms.normalise = normalise

    j_ms = MSSpectrum()
    j_ms.correct_for_ionisation = shift
    j_ms.filter = filter_func
    j_ms.normalise = normalise

    ii_peaks = ppk(i_ms.spectrum, i_ms.spectrum, sigma_mass, sigma_int)
    ii_nloss = ppk_nloss(
        i_ms.spectrum, i_ms.spectrum, i_ms.parentmass, i_ms.parentmass, sigma_mass, sigma_int
    )

    # i_ms.load(ms_path + os.sep + names[i] + '.ms')
    results = []
    for name in names:
        j_ms.load(ms_path + os.sep + name + ".ms")

        jj_peaks = ppk(j_ms.spectrum, j_ms.spectrum, sigma_mass, sigma_int)
        jj_nloss = ppk_nloss(
            j_ms.spectrum, j_ms.spectrum, j_ms.parentmass, j_ms.parentmass, sigma_mass, sigma_int
        )

        ij_peaks = ppk(i_ms.spectrum, j_ms.spectrum, sigma_mass, sigma_int)
        ij_nloss = ppk_nloss(
            i_ms.spectrum, j_ms.spectrum, i_ms.parentmass, j_ms.parentmass, sigma_mass, sigma_int
        )

        results.append(
            (ij_peaks / numpy.sqrt(ii_peaks * jj_peaks), ij_nloss / numpy.sqrt(ii_nloss * jj_nloss))
        )

    return results


def do_pair(i_spectrum, j_spectrum, i_parentmass, j_parentmass, sigma_mass, sigma_int):
    if len(i_spectrum) == 0 or len(j_spectrum) == 0:
        print("empty spectrum!")
        return 0, 0
    ij_peaks = ppk(i_spectrum, j_spectrum, sigma_mass, sigma_int)
    ij_nloss = ppk_nloss(i_spectrum, j_spectrum, i_parentmass, j_parentmass, sigma_mass, sigma_int)
    # ij_diff = ppk_diff(i_spectrum, j_spectrum, sigma_mass, sigma_int)
    return ij_peaks, ij_nloss  # , ij_diff


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", dest="filter_type", help="filter type (dag, tree)", default=None)
    parser.add_argument(
        "-c",
        dest="collected",
        help="build filter from all input files",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-s", dest="shift", help="correct for ionisation", action="store_true", default=False
    )
    parser.add_argument("-o", dest="output", help="output name", required=True)
    parser.add_argument(
        "-n", dest="normalise", help="normalise", action="store_true", default=False
    )
    args = parser.parse_args()

    if args.filter_type == "dag":
        if args.collected:
            filter_func = spectrum_filters.filter_by_collected_dag
        else:
            filter_func = spectrum_filters.filter_by_dag
    elif args.filter_type == "tree":
        if args.collected:
            filter_func = spectrum_filters.filter_by_collected_tree
        else:
            if args.shift == False:
                filter_func = spectrum_filters.filter_by_tree_unshifted
            else:
                filter_func = spectrum_filters.filter_by_tree
    elif args.filter_type is None:
        filter_func = None
    else:
        raise SystemExit("Unknown filter: %s" % args.filter_type)

    # create_ppk_matrix_parallell()
    create_ppk_matrix_stripe_serial(filter_func, args.shift, args.normalise, args.output)
    # create_ppk_matrix()
