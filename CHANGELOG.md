# Changelog

This file describes changes made to the NPLinker codebase in each tagged release. 

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

Starting from v1.1, each release listed here will also have a corresponding tagged [Docker image version](https://hub.docker.com/r/andrewramsay/nplinker/tags?page=1&ordering=last_updated). If you want to use the latest version of NPLinker, use `docker pull andrewramsay/nplinker:latest`. If you want a specific version, you can instead use `docker pull andrewramsay/nplinker:v1.1` etc. 

## [v1.2.0](https://github.com/sdrogers/nplinker/compare/v1.1.3...v1.2)

To use this release, run `docker pull andrewramsay:nplinker:v1.2.0`

### Changed

 - NPClassScore is added as a scoring method. Demo notebook can be found [here](notebooks/npclassscore_linking/NPClassScore_demo.ipynb).
 - CANOPUS is included in the workflow if `run_canopus = true` is set in the toml file/config.
 - Chemical class predictions from CANOPUS and MolNetEnhancer are read from directories 'canopus' and 'molnetenhancer'.
 - To perform class-based linking in NPClassScore, [MIBiG classes](prototype/nplinker/data/MIBiG2.0_compounds_with_AS_BGC_CF_NPC_classes.txt) are added from which [scoring tables](notebooks/npclassscore_linking/class_matching_tables.ipynb) are calculated depicting the matching class ontologies.

## [[v1.1.3]](https://github.com/sdrogers/nplinker/compare/v1.1.2...v1.1.3) - 2021-11-26

To use this release, run `docker pull andrewramsay:nplinker:v1.1.3`

### Changed

 - Reverted a change to inclusion of singleton MolFams in Metcalf scoring from v1.1.2 until [issue #57](https://github.com/sdrogers/nplinker/issues/57) can be resolved
 - Updated RefSeq accession ID process via NCBI website when downloading datasets from the paired data platform (website updates had broken the previous version)
 - Download process better handles genome records with no usable genome ID
 - Added more logging for information/debugging

## [[v1.1.2]](https://github.com/sdrogers/nplinker/compare/v1.1.1...v1.1.2) - 2021-09-23

To use this release, run `docker pull andrewramsay:nplinker:v1.1.2`

### Changed

 - BiG-SCAPE parameter handling behaviour has been updated. See the [default configuration file](https://github.com/sdrogers/nplinker/blob/666669e32724139bcc27d6869a986f891c1dc0cf/prototype/nplinker/data/nplinker.toml#L200) for details. The default settings are the same as the last release
 
## [[v1.1.1]](https://github.com/sdrogers/nplinker/compare/v1.1...v1.1.1) - 2021-09-13

To use this release, run `docker pull andrewramsay:nplinker:v1.1.1`

### Changed

 - BiG-SCAPE runs using the NPLinker Docker image now use the ``--mibig`` parameter (previously this was not a default option)
 - bugfix for NPLinker.molfams property (was returning a set instead of a list)
 - bugfix for parameter handling in NPLinker.get_links method
 

## [[v1.1]](https://github.com/sdrogers/nplinker/compare/v1.0...v1.1) - 2021-09-04

To use this release, run `docker pull andrewramsay:nplinker:v1.1`

### Added

 - Apache 2.0 license applied to project files and GitHub repo
 - AUTHORS file added
 - Support for NPLinker to better handle hybrid BGC objects. BGCs may now have multiple GCF parents corresponding to different BiG-SCAPE classes. The tables in the web application now have extra columns to allow you to check if a particular BGC is hybrid or not, and if a particular GCF is "pure" (only containing non-hybrid BGCs)
 - New configuration file option "extended_metadata_table_parsing". This applies only to datasets with FBMN metabolomics data from GNPS, and if enabled allows NPLinker to try and extract growth media labels for each strain automatically
 - .gitattributes 

### Changed

 - improvements to parsing of strain labels from GNPS FBMN workflows (previous implementation was not flexible enough)
 - bugfixes in Rosetta scoring
 - bugfix for parsing GNPS metadata

### Removed

 - old `strain_id_mapping.csv` file is now removed to avoid name clashes with uesr-provided strain mappings (this was an early development feature anyway)
 - various unused/older files cluttering the repository

## [[v1.0]](https://github.com/sdrogers/nplinker/compare/v.0.1...v1.0) - 2021-04-12

### Added
 
 - "ignore_spaces" setting in configuration file. If enabled, this allows NPLinker to automatically rename BiG-SCAPE input files if they contain spaces
 - more optional diagnostic output
 - automatically handle empty columns that may appear in `strain_mappings.csv`
 - new parser for antiSMASH .json files, to provide an alternative to relying on the deprecated text files which contain the same information
 - caching of Metcalf scoring data structures (provides a speedup on larger datasets)
 - ability to provide an "include_strains" file, listing strain IDs which should be explicitly included (all others are discarded)

### Changed

 - major update for handling downloads of [Paired Omics Platform](https://pairedomicsdata.bioinformatics.nl/) datasets: fix several bugs; increase ability of NPLinker to resolve accession IDs; store JSON data using UTF-8 encoding; improve detection of previously downloaded files to reduce network traffic/time required to retrieve a dataset
 - bugfix for a crash that occurred when trying to access a disabled scoring method
 - major update to Rosetta scoring: improve caching of data structures; save version number with cached data to detect changes; reorder preprocessing code so it will fail faster if required files are not all available; bugfixes
 - web application will now use latest available version of the [bokeh](https://github.com/bokeh/bokeh) module instead of an older customised version (previous issues that affected NPLinker are fixed)
 - expanded .gbk file parsing
 - reduce number of files extracted from downloaded .zip files when using Paired Platform datasets to save disk space
 - updated comments/content of default configuration file

### Removed

 - plotting interface in the web application (removed in favour of the newer tables-based interface)
