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
import numpy as np
import pylab as plt


# Method to perform the projection from the input object that
# Grimur provided
def bgc_projection(input_pickle, bgc_list, weights=[0, 1, 0]):
    # javcard, sequence, adjacency
    with open(input_pickle, "rb") as f:
        data = pickle.load(f, encoding="latin1")

    bgc_name_list = data[0]
    families = data[2]
    similarity_matrices = np.array(data[1])
    weights = np.array(weights)[:, None, None]
    d = (weights * similarity_matrices).sum(axis=0)

    L = np.diag(d) - d
    w, v = np.linalg.eig(L)

    bgc_dict = {}
    family_dict = {}
    for i, b in enumerate(bgc_name_list):
        bgc = [ba for ba in bgc_list if ba.name == b][0]
        bgc_dict[bgc] = (v[i, 0], v[i, 1])
        for family, bl in families.items():
            if b in bl:
                family_dict[bgc] = family
                break

    return bgc_dict, family_dict


def plot_and_highlight(bgc_dict, strains):
    if not isinstance(strains, list):
        strains = [strains]
    plt.figure(figsize=(10, 10))
    for b, pos in bgc_dict.items():
        plt.plot(pos[0], pos[1], "ko")
        for s in strains:
            if b.strain == s:
                plt.plot(pos[0], pos[1], "ro", markersize=10)


def plot_families(bgc_dict, family_dict, min_size=0):
    plt.figure(figsize=(10, 10))
    fam_data = {}
    for b, fam in family_dict.items():
        if fam not in fam_data:
            fam_data[fam] = []
        fam_data[fam].append(bgc_dict[b])
    for f, dat in fam_data.items():
        if len(dat) < min_size:
            plt.plot(dat[0][0], dat[0][1], "ko", markersize=5)
        else:
            x, y = zip(*dat)
            plt.plot(x, y, "o", markersize=10)


def k_means(bgc_dict, K=5, n_starts=10):
    N = len(bgc_dict)
    M = len(bgc_dict[list(bgc_dict.keys())[0]])

    bgc_list = list(bgc_dict.keys())
    X = []
    for b in bgc_list:
        X.append(bgc_dict[b])
    X = np.array(X)
    best = None
    for iteration in range(n_starts):
        mu = X[np.random.choice(N, K, replace=False), :]
        finished = False
        old_z = np.zeros((N, K))
        while not finished:
            # compute z
            some_empty = True
            while some_empty:
                di = (((X[None, :, :] - mu[:, None, :]) ** 2).sum(axis=2)).T
                z = np.zeros((N, K))
                for i, d in enumerate(di):
                    z[i, d.argmin()] = 1
                pos = np.where(z.sum(axis=0) == 0)[0]
                if len(pos) == 0:
                    some_empty = False
                for p in pos:
                    mu[p, :] = X[np.random.choice(N), :]

            # compute mu
            for k in range(K):
                mu[k, :] = (X * z[:, k][:, None]).mean(axis=0)

            if (old_z - z).sum() == 0:
                break

            old_z = z.copy()
        if iteration == 0:
            best_z = z.copy()
            best_mu = mu.copy()
        else:
            dis = 0.0
            dis_current = 0.0
            for n in range(N):
                dis += np.sqrt(((X[n, :] - mu[z[n, :].argmax(), :]) ** 2).sum())
                dis_current += np.sqrt(((X[n, :] - best_mu[best_z[n, :].argmax(), :]) ** 2).sum())
            if dis < dis_current:
                best_z = z.copy()
                best_mu = mu.copy()
                print(dis)

    families = {}
    for i, b in enumerate(bgc_list):
        families[b] = best_z[i, :].argmax()

    return families


def plot_families_strain_highlight(bgc_dict, gcf_list, strain_list):
    plt.figure(figsize=(20, 20))
    for b, pos in bgc_dict.items():
        plt.plot(pos[0], pos[1], "ko")
        plt.text(pos[0] + 0.001, pos[1] + 0.001, b.strain, fontsize=10)
    fam_success = {s: set() for s in strain_list}
    for s in strain_list:
        for f in gcf_list:
            if f.has_strain(s):
                fam_success[s].add(f)

    all_success = fam_success[strain_list[0]]
    for s in strain_list[1:]:
        all_success = all_success.intersection(fam_success[s])
    for f in all_success:
        x = []
        y = []
        for b in f.bgc_list:
            x.append(bgc_dict[b][0])
            y.append(bgc_dict[b][1])
        plt.plot(x, y, "o", markersize=10)


def hierarchical_clustering(bgc_dict, K=10):
    from genomics import GCF
    from sklearn.cluster import AgglomerativeClustering

    X = []
    bgc_list = []
    for b, pos in bgc_dict.items():
        bgc_list.append(b)
        X.append(pos)
    X = np.array(X)
    clustering = AgglomerativeClustering(n_clusters=K).fit(X)

    families = {}
    for i, b in enumerate(bgc_list):
        families[b] = clustering.labels_[i]

    gcf_dict = {i: GCF(i) for i in set(clustering.labels_)}
    for i, b in enumerate(bgc_list):
        gcf_dict[clustering.labels_[i]].add_bgc(b)
    return families, list(gcf_dict.values())
