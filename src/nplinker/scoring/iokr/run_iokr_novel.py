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
import os
import random
import time
import iokr_opt
import iokrdata as iokrdataserver
import numpy
from iokrdata import load_kernels
from mk_fprints import fingerprint_from_smiles


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
    collected_rankings = []
    collected_baseline_rankings = []

    all_indices = data.get_all_indices()
    filtered_indices = all_indices
    # filtered_indices = numpy.load('/home/grimur/iokr/data/mibig/mibig_idx_to_skip.bin.npy')
    # filtered_indices = list(filtered_indices)

    iokr = iokr_opt.InputOutputKernelRegression(data)
    iokr.set_training_indices(filtered_indices, _lambda=0.001)
    iokr.fit()

    latent, latent_basis, gamma = iokr.get_data_for_novel_candidate_ranking()

    # init the common candidate set
    candidate_set = [x for x in data.candidate_set]
    # init the sample set (same as candidate set?)
    spectra_kernels = data.spectra_kernels

    for i, sample_kernel_vector in enumerate(spectra_kernels):
        sample_kernel_vector = sample_kernel_vector[filtered_indices]
        # match the sample to the candidate set
        correct_indices = data.spectra_to_bgc_indices[i]
        candidate_fingerprints = numpy.array([x[1] for x in candidate_set]).copy()
        [x[0] for x in candidate_set]

        args = (0, candidate_fingerprints, latent, sample_kernel_vector, latent_basis, gamma)
        with open("secondary.bin", "wb") as f:
            import pickle

            pickle.dump(args, f)
        output = iokr_opt.rank_candidates_opt(*args)

        res_i = i
        res_correct_indices = correct_indices
        res_output = output

        res_ranking = list(res_output[0])

        # actual rankings
        correct_ranking = min(res_ranking.index(x) for x in res_correct_indices)
        collected_rankings.append((res_i, correct_ranking, len(res_ranking)))

        # correct_random_rankings = []
        # for random_iteration in range(100):
        #     random.shuffle(res_ranking)
        #     correct_random_ranking = min([res_ranking.index(x) for x in res_correct_indices])
        #     correct_random_rankings.append(correct_random_ranking)
        # collected_baseline_rankings.append((res_i, numpy.mean(correct_random_rankings), len(res_ranking)))

        random.shuffle(res_ranking)
        correct_random_ranking = min(res_ranking.index(x) for x in res_correct_indices)
        collected_baseline_rankings.append((res_i, correct_random_ranking, len(res_ranking)))

        total = len(collected_rankings)
        print(
            "%.06f, %.06f, %s"
            % (
                float([x[1] for x in collected_rankings].count(0)) / total,
                float([x[1] for x in collected_baseline_rankings].count(0)) / total,
                total,
            )
        )

        # print(cr_b[res_i], cr_a[res_i], cr_a[res_i] == cr_b[res_i])

    print("")
    print("IOKR test run done!")
    print(f"#samples: {len(collected_rankings)}")
    print(
        "top-1 acc: {}".format(
            float([x[1] for x in collected_rankings].count(0)) / total,
        )
    )
    print(
        "top-1 base: {}".format(
            float([x[1] for x in collected_baseline_rankings].count(0)) / total,
        )
    )

    return collected_rankings, collected_baseline_rankings


def get_test_set(fingerprint):
    import csv
    from pyteomics import mgf

    mapping = []
    # TODO: fix the path
    with open("/home/grimur/iokr/data/mibig/matched_mibig_gnps_2.0.csv") as f:
        for line in csv.reader(f):
            if not line[0].startswith("#"):
                mapping.append(line)

    mapping_dict = {}
    bgc_id_list = []
    bgc_fps = {}
    for line in mapping:
        line_id = line[0] + line[1]
        if line[0] in mapping_dict:
            mapping_dict[line_id].append(line)
        else:
            mapping_dict[line_id] = [line]

        bgc_id = line[4]
        bgc_smiles = line[7]
        bgc_fp = fingerprint_from_smiles(bgc_smiles, fingerprint)
        if bgc_id not in bgc_id_list:
            bgc_id_list.append(bgc_id)
            bgc_fps[bgc_id] = bgc_fp

    bgc_fp_list = [bgc_fps[x] for x in bgc_id_list]
    candidate_set = zip(bgc_id_list, bgc_fp_list)

    mgf_file = "/home/grimur/iokr/data/mibig/matched_mibig_gnps_2.0.mgf"
    spectra_to_bgc_indices = []
    for i, spectrum in enumerate(mgf.read(mgf_file)):
        spec_file = spectrum["params"]["filename"]
        spec_id = spectrum["params"]["feature_id"]
        spec_filename = spec_file + spec_id
        bgc_matches = mapping_dict[spec_filename]
        bgc_indices = []
        for line in bgc_matches:
            bgc_id = line[4]
            bgc_idx = bgc_id_list.index(bgc_id)
            bgc_indices.append(bgc_idx)
        spectra_to_bgc_indices.append(bgc_indices)

    kernel_files = [
        #'/home/grimur/iokr/sml_test_test_nloss.npy',
        #'/home/grimur/iokr/sml_test_test_peaks.npy'
        "/home/grimur/iokr/mibig_test_nloss.npy",
        "/home/grimur/iokr/mibig_test_peaks.npy",
    ]
    spectra_kernels = load_kernels(kernel_files, normalise=False)

    return spectra_kernels, candidate_set, spectra_to_bgc_indices


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
    output_baseline = "IOKRranking_baseline_%s.bin" % args.output

    iokrdata = iokrdataserver.IOKRDataServer(datapath, kernel=None)
    kernel_matrix = load_kernels(kernel_files)
    iokrdata.kernel = kernel_matrix
    with open(datapath + os.sep + "ind_eval.txt") as f:
        raw_data = f.read()
        test_sample_indices = [int(x) - 1 for x in raw_data.strip().split()]
    iokrdata.test_sample_indices = test_sample_indices

    if fingerprint is not None:
        iokrdata.set_fingerprint(fingerprint)

    spectra_kernels, candidate_set, spectra_to_bgc_indices = get_test_set(fingerprint)
    with open("inputtest.bin", "wb") as f:
        import pickle

        pickle.dump((spectra_kernels, candidate_set, spectra_to_bgc_indices), f)
    # with open('inputtest.bin', 'rb') as f:
    #     import pickle
    #     spectra_kernels, candidate_set, spectra_to_bgc_indices = pickle.load(f)

    # print(spectra_kernels)
    # print(candidate_set)
    # print(spectra_to_bgc_indices)

    iokrdata.spectra_kernels = spectra_kernels
    iokrdata.candidate_set = candidate_set
    iokrdata.spectra_to_bgc_indices = spectra_to_bgc_indices

    print("run iokr")
    rankings, baseline = run_iokr(iokrdata)

    numpy.save(output_file, rankings)
    numpy.save(output_baseline, baseline)


if __name__ == "__main__":
    main()
