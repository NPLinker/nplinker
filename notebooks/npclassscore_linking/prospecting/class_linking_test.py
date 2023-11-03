#!/usr/bin/env python3
"""Initial testing for NPClassScore integration into NPLinker."""

from class_linking import NPLinker_classes


if __name__ == "__main__":
    # load local crusemann data for testing
    test_data_path = (
        "/mnt/scratch/louwe015/NPLinker/own/nplinker_shared/crusemann_3ids_AS6-AS3_30-11/"
    )
    npl = NPLinker_classes({"dataset": {"root": test_data_path}})
    npl.load_data()
    npl.read_class_info()

    # 1. print general config info about dataset
    # - a copy of the configuration as parsed from the .toml file (dict)
    print(npl.config)
    # - the path to the directory where various nplinker data files are located (e.g. the
    #   default configuration file template) (str)
    print(npl.data_dir)
    # - a dataset ID, derived from the path for local datasets or the paired platform ID
    #   for datasets loaded from that source (str)
    print(npl.dataset_id)
    # - the root directory for the current dataset (str)
    print(npl.root_dir)

    # objects
    # - you can directly access lists of each of the 4 object types:
    print("BGCs:", len(npl.bgcs))
    print("GCFs:", len(npl.gcfs))  # contains GCF objects
    print("Spectra:", len(npl.spectra))  # contains Spectrum objects
    print("Molecular Families:", len(npl.molfams))  # contains MolecularFamily objects

    # 2. check chemical compound predictions from canopus and molnetenhancer
    test_spec = list(npl.spectra)[500]

    print(npl.canopus.spectra_classes.get(str(test_spec.spectrum_id)))
    print(npl.molnetenhancer.spectra_classes(str(test_spec.spectrum_id)))

    # 3. example of a good score, (predicted) NRP linking to a (predicted) peptide like spectrum
    print(npl.class_linking_score(list(npl.gcfs)[0], test_spec))
    for class_method in ("main", "mix", "canopus", "molnetenhancer"):
        print(class_method)
        print(npl.npclass_score(list(npl.gcfs)[0], test_spec, method=class_method))

    # 4. do linking
    mc = npl.scoring_method("metcalf")
    mc.cutoff = 2.5
    mc.standardised = True

    results = npl.get_links(npl.gcfs, mc, and_mode=True)
    print(f"Number of results: {len(results)}")
    print(results.methods)
