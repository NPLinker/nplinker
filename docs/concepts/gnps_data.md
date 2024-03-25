NPLinker requires GNPS molecular networking data as input. It currently accepts data from the following 
GNPS workflows:

- `METABOLOMICS-SNETS` (data should be downloaded from the option `Download Clustered Spectra as MGF`)
- `METABOLOMICS-SNETS-V2` (`Download Clustered Spectra as MGF`)
- `FEATURE-BASED-MOLECULAR-NETWORKING` (`Download Cytoscape Data`)


## Mappings from GNPS data to NPLinker input

=== "`METABOLOMICS-SNETS` workflow"

    | NPLinker input         | GNPS file in the archive of `Download Clustered Spectra as MGF`  |
    | ---------------------- | ---------------------------------------------------------------- |
    | spectra.mgf            | METABOLOMICS-SNETS*.mgf                                          |
    | molecular_families.tsv | networkedges_selfloop/*.pairsinfo                                |
    | annotations.tsv        | result_specnets_DB/*.tsv                                         |
    | file_mappings.tsv      | clusterinfosummarygroup_attributes_withIDs_withcomponentID/*.tsv |

    For example, the file `METABOLOMICS-SNETS*.mgf` from the downloaded zip archive is used as 
    the `spectra.mgf` input file of NPLinker. 
    
    When manually preparing GNPS data for NPLinker, the `METABOLOMICS-SNETS*.mgf` must be renamed to
    `spectra.mgf` and placed in the `gnps` sub-directory of the NPLinker working directory.
    

=== "`METABOLOMICS-SNETS-V2`"

    | NPLinker input         | GNPS file in the archive of `Download Clustered Spectra as MGF`             |
    | ---------------------- | --------------------------------------------------------------------------- |
    | spectra.mgf            | METABOLOMICS-SNETS-V2*.mgf                                                  |
    | molecular_families.tsv | networkedges_selfloop/*.selfloop                                            |
    | annotations.tsv        | result_specnets_DB/*.tsv                                                    |
    | file_mappings.tsv      | clusterinfosummarygroup_attributes_withIDs_withcomponentID/*.clustersummary |


=== "`FEATURE-BASED-MOLECULAR-NETWORKING`"

    | NPLinker input         | GNPS file in the archive of `Download Cytoscape Data` |
    | ---------------------- | ----------------------------------------------------- |
    | spectra.mgf            | spectra/*.mgf                                         |
    | molecular_families.tsv | networkedges_selfloop/*.selfloop                      |
    | annotations.tsv        | DB_result/*.tsv                                       |
    | file_mappings.csv      | quantification_table/*.csv                            |

    Note that `file_mappings.csv` is a CSV file, not a TSV file, different from the other workflows.