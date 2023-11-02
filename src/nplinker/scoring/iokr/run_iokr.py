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

import argparse
import multiprocessing
import os
import time
import iokr_opt
import iokrdata as iokrdataserver
import numpy


def normalise_kernel(matrix):
    return matrix / numpy.sqrt(numpy.outer(matrix.diagonal(), matrix.diagonal()))


def load_kernel_file(filename):
    kernel = numpy.load(filename)
    return normalise_kernel(kernel)


def load_kernels(kernel_files):
    kernel_matrices = [load_kernel_file(x) for x in kernel_files]
    kernel_sum = numpy.sum(kernel_matrices, axis=0)
    return normalise_kernel(kernel_sum)


def gather_results(active_jobs, limit=25):
    done_jobs = []
    while len(active_jobs) > limit:
        remaining_jobs = []
        for line in active_jobs:
            job = line[-1]
            if job.ready():
                res = job.get()
                new_line = list(line[:-1])
                new_line.append(res)
                done_jobs.append(new_line)
            else:
                remaining_jobs.append(line)
        active_jobs = remaining_jobs
        time.sleep(1)
    return active_jobs, done_jobs


def run_iokr(data):
    p = multiprocessing.Pool()
    active_jobs = []
    job_limit = 25

    excluded = []
    missing_candidate = []
    collected_rankings = []

    for label in sorted(list(set(data.folds))):
        print("label %s" % label)
        label_indices = data.get_indices(label, complement=True)

        iokr = iokr_opt.InputOutputKernelRegression(data)
        iokr.set_training_indices(label_indices, _lambda=0.001)
        iokr.fit()

        test_indices = data.get_indices(label)
        for i in test_indices:
            if i not in data.test_sample_indices:
                excluded.append(i)
                continue

            sample = data.get_sample(i)
            formula = sample["formula"]
            sample_inchi = sample["inchi"]

            candidates = data.get_candidates(formula)
            candidate_inchi = [x[0] for x in candidates]
            correct_index = candidate_inchi.index(sample_inchi)

            # TODO: Recalculate fingerprints
            candidates = data.get_candidates(formula)
            candidate_fingerprints = [numpy.array(x[1]) for x in candidates]

            total_count = len(candidate_inchi)

            print(f"iokr job idx {i}, cand.set size {total_count}")
            # ranking = iokr.rank_candidates_opt(i, candidate_fingerprints)
            # # print(ranking)
            # ranking = list(ranking)
            # correct_ranking = ranking.index(correct_index)
            # print('ranked {} / {}'.format(correct_ranking, total_count))

            latent, x_kernel_vector, latent_basis, gamma = iokr.get_data_for_candidate_ranking(i)
            args = (i, candidate_fingerprints, latent, x_kernel_vector, latent_basis, gamma)
            job = p.apply_async(iokr_opt.rank_candidates_opt, args)
            active_jobs.append((i, formula, correct_index, label, total_count, job))

            if len(active_jobs) > job_limit:
                active_jobs, results = gather_results(active_jobs, job_limit)
                for (
                    res_i,
                    res_formula,
                    res_correct_index,
                    res_label,
                    res_total_count,
                    res_output,
                ) in results:
                    res_ranking = list(res_output[0])
                    # print(res_ranking)
                    correct_ranking = res_ranking.index(res_correct_index)

                    collected_rankings.append((res_i, correct_ranking, res_total_count))

                    total = len(collected_rankings)
                    print(float([x[1] for x in collected_rankings].count(0)) / total, total)

                    # print(cr_b[res_i], cr_a[res_i], cr_a[res_i] == cr_b[res_i])

    print("Clean up remaining jobs")

    # clean up the last remaining jobs
    active_jobs, results = gather_results(active_jobs, 0)
    for res_i, res_formula, res_correct_index, res_label, res_total_count, res_output in results:
        res_ranking = list(res_output[0])
        correct_ranking = res_ranking.index(res_correct_index)
        collected_rankings.append((res_i, correct_ranking, res_total_count))
        total = len(collected_rankings)
        print(float([x[1] for x in collected_rankings].count(0)) / total, total)

    print("")
    print("IOKR test run done!")
    print(f"#samples: {len(collected_rankings)}")
    print(
        "top-1 acc: {}".format(
            float([x[1] for x in collected_rankings].count(0)) / total,
        )
    )

    return collected_rankings


def main():
    parser = argparse.ArgumentParser("Run IOKR test on a set")
    parser.add_argument("--kernel", dest="kernel", help="Kernel files", nargs="+")
    parser.add_argument(
        "--fp",
        dest="fingerprint",
        help="fingerprint type (substructure, cdk (default), klekota-roth",
        default="cdk_default",
    )
    parser.add_argument("--data", dest="datapath", help="data path", required=True)
    parser.add_argument("--output", dest="output", help="output label", required=True)
    args = parser.parse_args()

    # read from args
    # datapath = '/home/grimur/iokr/data'
    datapath = args.datapath
    # kernel_files = [datapath + os.sep + 'input_kernels_gh/ppk_dag_all_normalised_shifted_nloss.npy',
    #                 datapath + os.sep + 'input_kernels_gh/ppk_dag_all_normalised_shifted_peaks.npy']
    kernel_files = args.kernel
    # fingerprint = None

    fingerprint = args.fingerprint

    output_file = "IOKRranking_%s.bin" % args.output

    iokrdata = iokrdataserver.IOKRDataServer(datapath, kernel=None)
    kernel_matrix = load_kernels(kernel_files)
    iokrdata.kernel = kernel_matrix
    with open(datapath + os.sep + "ind_eval.txt") as f:
        raw_data = f.read()
        test_sample_indices = [int(x) - 1 for x in raw_data.strip().split()]
    iokrdata.test_sample_indices = test_sample_indices

    if fingerprint is not None:
        iokrdata.set_fingerprint(fingerprint)

    print("run iokr")
    rankings = run_iokr(iokrdata)

    numpy.save(output_file, rankings)


if __name__ == "__main__":
    main()
